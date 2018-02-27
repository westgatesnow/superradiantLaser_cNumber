#ifndef __SUPERRADIANTLASER__HPP__
#define __SUPERRADIANTLASER__HPP__
//This program is used to simulate the superradiant laser using the cNumber Langevin method.

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


//Define data structures

//Atom Internal states
typedef struct {
  VectorXd sx;       //sigma_x. Dim: nTrajectory
  VectorXd sy;       //sigma_y. Dim: nTrajectory
  VectorXd sz;       //sigma_z. Dim: nTrajectory
} Atom;

#define NVAR 3 // 3 variables for each atom for this code;

//Ensemble of atoms
typedef struct {
  std::vector<Atom> atoms;
} Ensemble;

//Simulation parameters
typedef struct Param {
  //simulation specification
  double dt; 
  double tmax;
  int nstore; //number of times to store observables
  int nTrajectory; //number of trajectories
  //beam parameters
  int nAtom; //number of intracavity atoms
  double gammac;    //collective decay rate
  double repumping; //repumping rate w
  //Other parameters
  std::string name; //name of the directory to store results

  //Set up initial values of the parameters
  Param() : dt(1.0e-3), tmax(1), nstore(100), nTrajectory(1), nAtom(10), gammac(0.1e0), repumping(10), 
            name("abracadabra")  {}
} Param;

std::ostream& operator<< (std::ostream& o,const Param& s)
{
  o << s.dt << std::endl;
  o << s.tmax << std::endl;
  o << s.nstore << std::endl;
  o << s.nTrajectory << std::endl; 
  o << s.nAtom << std::endl;
  o << s.gammac << std::endl;
  o << s.repumping << std::endl;

  return o;
}

typedef struct Observables {
  Observables(const int n/*,const int m*/) : intensity(n), intensityUnCor(n),
                                          inversion(n), spinSpinCor(n)
                                          //,g1(m)
  {}
  VectorXd intensity;
  VectorXd intensityUnCor;
  VectorXd inversion;
  VectorXd spinSpinCor;
  //VectorXd g1;
} Observables;

typedef struct ObservableFiles {
  ObservableFiles() : intensity("intensity.dat"), intensityUnCor("intensityUnCor.dat"), 
                      inversion("inversion.dat"), spinSpinCor("spinSpinCor.dat")
                      //,g1("g1.dat")
  {}
  ~ObservableFiles() {
    intensity.close();
    intensityUnCor.close();
    inversion.close();
    spinSpinCor.close();
    //g1.close();
  }
  std::ofstream intensity, intensityUnCor, inversion, spinSpinCor;//, g1;
} ObservableFiles;

#endif