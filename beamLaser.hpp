#ifndef __BEAMLASER__HPP__
#define __BEAMLASER__HPP__
//This program is used to simulate the beam laser using the cNumber Langevin method.

//Include Eigen package
//Work in Eigen namespace
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Eigenvalues>
using namespace Eigen;

//Include standard packages
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <cmath>
#include <getopt.h>
#include <time.h>
#include <stdlib.h> 

//Include and define RNG
#include "RNG.hpp"
RNG rng(time(NULL));

//Define the complex I and ONE
static const std::complex<double> I = std::complex<double>(0.0,1.0);
static const std::complex<double> ONE = std::complex<double>(1.0,0.0);

#define NVAR 3 // 3 variables for each atom for this code;

//Define data structures

//Atom external states
typedef struct {
  Vector3d X;     //position
  Vector3d P;     //momentum. We suppose mass is one, so momentum is velocity.
} External;

//Atom internal states
typedef struct {
  VectorXd sx;       //sigma_x. Dim: nTrajectory
  VectorXd sy;       //sigma_y. Dim: nTrajectory
  VectorXd sz;       //sigma_z. Dim: nTrajectory
} Internal;

//Atom total states
typedef struct {
  External external;  //The position and velocity of an atom.
  Internal internal;  //The internal states of an atom. We keep track of sx, sy, and sz
                        //of a single atom at all times for all trajectories.
} Atom;

//Ensemble of atoms
typedef struct {
  std::vector<Atom> atoms;
} Ensemble;

//Configuration setup; copied from Murray's old codes
typedef struct {
  const char* configFile;
} CmdLineArgs;
const char* usageHeader = "\nBeam Laser Simulation.\n";
const char* usageMessage =
  "\n"
  "Usage:         "
  "beamLaser "
  "--file"
  "\n"
  "--file, -f     : Configuration file to setup the simulation\n"
  "--help, -h     : Print this usage message\n"
  "\n\n";

//Simulation parameters
typedef struct Param {
  //simulation specification
  double meanAtomNumber; //the desired average number of intracavity atoms N0
  int nTrajectory; //number of trajectories
  int nstore; // number of times to store observables
  int steadyTime;  // tmax/tau
  //beam parameters
  double yWall; //position of the wall where atoms are destroyed
                //The coordinated are chosen s.t. atoms are created at -yWall.
                //The walls are assumed to be in xz plane.
  Vector2d sigmaX;  //standard deviation of position in xz. y deviation is 
                    //taken care of by the Poisson distribution.
  double transitTime;     //the transit time tau1>0, unit 1/gammaC.
  Vector3d sigmaP;  //standard deviation of momentum
  int meanAtomGeneratingNumber;   //the mean number of atoms dN generated in time dt, dN << N0;
  double gammac;    //collective decay rate

  //Other parameters
  std::string name; //name of the directory to store results

  //Set up initial values of the parameters
  Param() : meanAtomNumber(100), nTrajectory(1), nstore(10), steadyTime(10), yWall(5.0e0), 
    sigmaX(0.0e0,0.0e0), transitTime(1.0e0),
    sigmaP(0.0e0,0.0e0,0.0e0), meanAtomGeneratingNumber(1), gammac(0.1e0), name("abracadabra")  {}
} Param;

std::ostream& operator<< (std::ostream& o, const Param& s)
{
  o << s.meanAtomNumber << std::endl;
  o << s.steadyTime << std::endl;
  o << s.nTrajectory << std::endl;
  o << s.nstore << std::endl;
  o << s.yWall << std::endl;
  o << s.sigmaX << std::endl;
  o << s.transitTime << std::endl;
  o << s.sigmaP << std::endl;
  o << s.meanAtomGeneratingNumber << std::endl;
  o << s.gammac << std::endl;
  return o;
}

//Observables; n is the nTimeStep
typedef struct Observables {
  Observables(const int n/*,const int m*/) : nAtom(n), 
                                          intensity(n), intensityUnCor(n),
                                          inversion(n), 
                                          spinSpinCor(n)
                                          //,g1(m)
  {}
  Matrix <unsigned long int, 1, Dynamic> nAtom; 
  VectorXd intensity;
  VectorXd intensityUnCor;
  VectorXd inversion;
  VectorXd spinSpinCor;
  //VectorXd g1;
} Observables;

typedef struct ObservableFiles {
  ObservableFiles() : nAtom("nAtom.dat"), 
                      intensity("intensity.dat"), intensityUnCor("intensityUnCor.dat"), 
                      inversion("inversion.dat"), 
                      spinSpinCor("spinSpinCor.dat")
                      //,g1("g1.dat")
  {}
  ~ObservableFiles() {
    nAtom.close();
    intensity.close();
    intensityUnCor.close();
    inversion.close();
    spinSpinCor.close();
    //g1.close();
  }
  std::ofstream nAtom, intensity, intensityUnCor, inversion, spinSpinCor;//, g1;
} ObservableFiles;

#endif
