//This program is used to simulate the superradiant laser using the cNumber Langevin method.
#include "superradiantLaser.hpp"
#include "config.hpp"

//Changes required subject to the definition of Param 
void getParam(const char* filename, Param *param) 
{
  std::ifstream configInput(filename);
  std::string dummy;

  while (!configInput.eof()) {
    configInput >> dummy;
    if (configInput.eof()) break;
    if (!configInput.good()) {
      std::cout << "Bad read in input file" << std::endl;
      exit(-1);
    }
    if (dummy.compare("dt") == 0)
      configInput >> param->dt;
    else if (dummy.compare("tmax") == 0)
      configInput >> param->tmax;
    else if (dummy.compare("nstore") == 0)
      configInput >> param->nstore;
    else if (dummy.compare("nTrajectory") == 0)
      configInput >> param->nTrajectory;
    else if (dummy.compare("nAtom") == 0)
      configInput >> param->nAtom;
    else if (dummy.compare("gammac") == 0)
      configInput >> param->gammac;
    else if (dummy.compare("repumping") == 0)
      configInput >> param->repumping;
    else if (dummy.compare("name") == 0)
      configInput >> param->name;
    else {
      std::cout << "Error: invalid label " << dummy << " in "
          << filename << std::endl;
      exit(-1);
    }
  }
}

void generateInitialAtoms(Ensemble& ensemble, const Param& param)
{ 
  //For convenience
  const int nTrajectory = param.nTrajectory;
  const int nAtom = param.nAtom;
  
  //Generate nAtom atoms
  for (int i = 0; i < nAtom; i++) {
    Atom newAtom;
    
    //sx and sy
    newAtom.sx = VectorXd::Zero(nTrajectory);
     for (int j = 0; j < nTrajectory; j++) 
       newAtom.sx[j] = double(rng.get_binomial_int(0.5, 1))*2-1; //50percent giving 1 or -1
    newAtom.sy = VectorXd::Zero(nTrajectory);
     for (int j = 0; j < nTrajectory; j++)   
       newAtom.sy[j] = double(rng.get_binomial_int(0.5, 1))*2-1;     //50percent giving 1 or -1
   
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Ground state
    newAtom.sz = -VectorXd::Ones(nTrajectory);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Excited state
    //newAtom.sz = VectorXd::Ones(nTrajectory);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    ensemble.atoms.push_back(newAtom);
  }
}


void getDiffusionMatrix(const Ensemble& ensemble, MatrixXd& diffusion, const Param& param)
{ 
  //For convenience
  const double gc = param.gammac;
  const double w = param.repumping;
  const int nAtom = param.nAtom;
  const int nTrajectory = param.nTrajectory;

  //Diagonal block matrices
  for (int i = 0; i < nAtom; i++) {
    diffusion(NVAR*i, NVAR*i) = gc+w;
    diffusion(NVAR*i, NVAR*i+1) = 0;
    diffusion(NVAR*i, NVAR*i+2) = (gc-w)*ensemble.atoms[i].sx.sum()/nTrajectory;
    diffusion(NVAR*i+1, NVAR*i+1) = gc+w;
    diffusion(NVAR*i+1, NVAR*i+2) = gc*ensemble.atoms[i].sy.sum()/nTrajectory;
    diffusion(NVAR*i+2, NVAR*i+2) = 2*gc*(1+ensemble.atoms[i].sz.sum()/nTrajectory)
                                   +2*w*(1-ensemble.atoms[i].sz.sum()/nTrajectory);
  }

  //Off-diagonal block matrices
  for (int i = 0; i < nAtom; i++) {
    for (int j = i+1; j < nAtom; j++) {
      diffusion(NVAR*i, NVAR*j) = gc/nTrajectory
        *(ensemble.atoms[i].sz.cwiseProduct(ensemble.atoms[j].sz)).sum();
      diffusion(NVAR*i, NVAR*j+1) = 0;
      diffusion(NVAR*i, NVAR*j+2) = -gc/nTrajectory
        *(ensemble.atoms[i].sz.cwiseProduct(ensemble.atoms[j].sx)).sum();
      diffusion(NVAR*i+1, NVAR*j) = 0;
      diffusion(NVAR*i+1, NVAR*j+1) = gc/nTrajectory
        *(ensemble.atoms[i].sz.cwiseProduct(ensemble.atoms[j].sz)).sum();
      diffusion(NVAR*i+1, NVAR*j+2) = -gc/nTrajectory
        *(ensemble.atoms[i].sz.cwiseProduct(ensemble.atoms[j].sy)).sum();
      diffusion(NVAR*i+2, NVAR*j) = -gc/nTrajectory
        *(ensemble.atoms[j].sz.cwiseProduct(ensemble.atoms[i].sx)).sum();
      diffusion(NVAR*i+2, NVAR*j+1) = -gc/nTrajectory
        *(ensemble.atoms[j].sz.cwiseProduct(ensemble.atoms[i].sy)).sum();
      diffusion(NVAR*i+2, NVAR*j+2) = gc/nTrajectory
        *((ensemble.atoms[i].sx.cwiseProduct(ensemble.atoms[j].sx)).sum()
         +(ensemble.atoms[i].sy.cwiseProduct(ensemble.atoms[j].sy)).sum());
    }
  }

  //The lower half of the symmetric matrix
  for (int i = 1; i < NVAR*nAtom; i++) {
    for (int j = 0; j < i; j++) {
      diffusion(i,j) = diffusion(j,i);
    }
  }
}

void getSqrtMatrix(const MatrixXd& diffusion, MatrixXd& bMatrix)
{
  SelfAdjointEigenSolver<MatrixXd> es(diffusion);
  bMatrix = es.operatorSqrt();
  // MatrixXd v = es.eigenvectors();
  // MatrixXd vt = v.transpose();
  // VectorXd eigen = es.eigenvalues();
  // for (int i=0; i<eigen.size(); i++)
  //   if(abs(eigen[i]) < 1.0E-10)
  //     eigen[i] = 0;
  // MatrixXd l = eigen.asDiagonal();
  // bMatrix = v*l.cwiseSqrt()*vt;
  // //Quit if NaN
  // std::cout << "bMatrix is " << std::endl << bMatrix << std::endl << std::endl;
  // std::cout << "Eigenvalues of mMatrix is " << std::endl << es.eigenvalues() << std::endl << std::endl;
  if (bMatrix(1,1) != bMatrix(1,1)) {
    std::cout << "Not a number\n" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void getDriftVector(const VectorXd& sAtoms, VectorXd& drift, const Param& param) 
{
  //For convenience
  const double gc = param.gammac;
  const double w = param.repumping;
  const int nAtom = param.nAtom;
  const int size = NVAR*nAtom;
  
  //Definition of Jx and Jy
  double jx = 0, jy = 0;
  for (int j = 0; j < nAtom; j++) {
    jx += sAtoms[NVAR*j];
    jy += sAtoms[NVAR*j+1];
  }

  //Drift vector terms. Dimension 3*nAtom, structure {D1+,D1-,D1z, D2+, D2-, D2z,....}
  for (int j = 0; j < nAtom; j++) {
    drift[NVAR*j] = gc/2*(jx*sAtoms[NVAR*j+2]-sAtoms[NVAR*j])-w/2*sAtoms[NVAR*j];
    drift[NVAR*j+1] = gc/2*(jy*sAtoms[NVAR*j+2]-sAtoms[NVAR*j+1])-w/2*sAtoms[NVAR*j+1];
    drift[NVAR*j+2] = -gc/2*(jx*sAtoms[NVAR*j]+jy*sAtoms[NVAR*j+1]-pow(sAtoms[NVAR*j],2)-pow(sAtoms[NVAR*j+1],2)
                     +2+2*sAtoms[NVAR*j+2])+w*(1-sAtoms[NVAR*j+2]);
  }
  
}

void advanceInterval(Ensemble& ensemble, const Param& param)
{
  //For convenience
  const double dt = param.dt;
  const double gc = param.gammac;
  const int nTrajectory = param.nTrajectory;
  const int nAtom = ensemble.atoms.size();
  const int size = NVAR*nAtom;

  //Diffusion matrix 2*m_ij. The data structure is {1x,1y,1z,2x,2y,2z,...}^2;
  MatrixXd diffusion = MatrixXd::Zero(size, size);
  getDiffusionMatrix(ensemble, diffusion, param);
  //std::cout << diffusion <<std::endl << std::endl;
  //Get the sqrt matrix.
  MatrixXd bMatrix = MatrixXd::Zero(size, size);
  getSqrtMatrix(diffusion, bMatrix);
  

  //Loop over all trajectories. "n" stands for the number of the current trajectory.
  for (int n = 0; n < nTrajectory; n++) {

    //Y_n. Notations see "Beam laser/11.1" on ipad pro.
    VectorXd sAtoms = VectorXd::Zero(size);;//a vector of spins for all the atoms
    for (int i = 0; i < nAtom; i++) {
      sAtoms[NVAR*i] = ensemble.atoms[i].sx[n];
      sAtoms[NVAR*i+1] = ensemble.atoms[i].sy[n];
      sAtoms[NVAR*i+2] = ensemble.atoms[i].sz[n];
    }

    //drift as a function of sAtoms
    VectorXd drift = VectorXd::Zero(size);
    getDriftVector(sAtoms, drift, param);


    //dW
    // double dw = rng.get_gaussian_rn(sqrt(dt));
    // VectorXd dW = VectorXd::Ones(size)*dw;
    
    VectorXd dW = VectorXd::Zero(size);
    for (int i = 0; i < size; i++) {
      dW[i] = rng.get_gaussian_rn(sqrt(dt));
    }

    //Y'. Notations see "Beam laser/11.1" on ipad pro.
    VectorXd sAtomsTemp = sAtoms+drift*dt+bMatrix*dW;//Temporary value for sAtoms

    //driftTemp as a function of sAtomsTemp
    VectorXd driftTemp = VectorXd::Zero(size);
    getDriftVector(sAtomsTemp, driftTemp, param);

    //Y_{n+1}
    sAtoms += (driftTemp+drift)*dt/2+bMatrix*dW;

    //Put back
    for (int i = 0; i < nAtom; i++) {
      ensemble.atoms[i].sx[n] = sAtoms[NVAR*i];
      ensemble.atoms[i].sy[n] = sAtoms[NVAR*i+1];
      ensemble.atoms[i].sz[n] = sAtoms[NVAR*i+2];
    }
  }
}

void storeObservables(Observables& observables, int s, Ensemble& ensemble, 
    const Param& param)
{
  //For convenience
  const double gc = param.gammac;
  const int nTrajectory = param.nTrajectory;
  const int nAtom = ensemble.atoms.size();
  
  //The order of the definitions of observables matters. 
  
  //inversion
  double inversion = 0;
  for (int i = 0; i < nAtom; i++) {
    inversion += ensemble.atoms[i].sz.sum();
  }
  observables.inversion(s) = inversion/nAtom/nTrajectory;

  //intensityUnCor
  double intensityUnCor = 0;
  observables.intensityUnCor(s) = gc/2*(nAtom+inversion/nTrajectory);

  //spinSpinCor
  double spinSpinCor =0;
  for (int n = 0; n < nTrajectory; n++) {
    for (int i = 0; i < nAtom; i++) {
      for (int j = 0; j < nAtom; j++) {
        spinSpinCor += (ensemble.atoms[i].sx[n]*ensemble.atoms[j].sx[n]
                       +ensemble.atoms[i].sy[n]*ensemble.atoms[j].sy[n])/4;
      }
    }
  }
  observables.spinSpinCor(s) = spinSpinCor/nAtom/(nAtom-1)/nTrajectory;

  //intensity
  observables.intensity(s) = gc*spinSpinCor/nTrajectory+observables.intensityUnCor(s);
}
  

void evolve(Ensemble& ensemble, const Param& param, Observables& observables)
{
  //Integration conditions
  double dt = param.dt;
  double tmax = param.tmax;
  
  //evolve
  int nTimeStep = tmax/dt+0.5;
  double tstep = dt, t=0;

  for (int n = 0, s = 0; n <= nTimeStep; n++, t += tstep) {
    if ((long)(n+1)*param.nstore/(nTimeStep+1) > s) {
      storeObservables(observables, s++, ensemble, param);
      //debug
      std::cout << "This is timestep " << n << "/" << nTimeStep << std::endl << std::endl;
      //debug
    }
    if (n != nTimeStep)
      advanceInterval(ensemble, param);
  }
}

void writeObservables(ObservableFiles& observableFiles, 
    Observables& observables)
{
  observableFiles.intensity << observables.intensity << std::endl;
  observableFiles.intensityUnCor << observables.intensityUnCor << std::endl;
  observableFiles.inversion << observables.inversion << std::endl;
  observableFiles.spinSpinCor << observables.spinSpinCor << std::endl;

}

int main(int argc, char *argv[])
{
  //Count time
  clock_t t1,t2;
  t1=clock();

  //Configuration
  CmdLineArgs config;
  getOptions(argc, argv, &config);

  //Set up parameters
  Param param;
  getParam (config.configFile, &param);
	
  //Set up initial conditions
  Ensemble ensemble;
  generateInitialAtoms(ensemble, param);
  Observables observables(param.nstore);

  //Start simulation
  evolve(ensemble, param, observables);

  //Write Observables
  ObservableFiles observableFiles;
  writeObservables(observableFiles, observables);

  //Move .dat files into the directory named "name"
  mkdir(param);
  
  //Count time
  t2=clock();
  float diff = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
  std::cout << "\nThis program takes " << diff << " seconds." << std::endl << std::endl;

  return 0;
}