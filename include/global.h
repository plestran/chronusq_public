/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#ifndef INCLUDED_GLOBAL
#define INCLUDED_GLOBAL
// Hello

#define MAXNAMELEN 50
#define MAXANGULARMOMENTUM 6
#define MAXCONTRACTION 10
#define MAXATOMS 1000

// CMake Compilation Configuration
#include <config_chronusq.h>

// IO
#include <iostream>
#include <fstream>
#include <iomanip>

// Math
#include <cmath>
#include <complex>
#define EIGEN_MATRIXBASE_PLUGIN "eigenplugin.h" // ChronusQ plugin for Eigen
#include <Eigen/Core> // Eigen Linear Algebra
#include <unsupported/Eigen/MatrixFunctions>
#ifdef USE_LIBINT
#  include <libint2.hpp> // Libint Gaussian Integrals library
#endif
#include <eiginterface.h> // ChronusQ interface (TODO consolidate into plugin)
#include <btas/btas.h> // BTAS Tensor Algebra library (header only)

// Parallelization
#ifdef USE_OMP
#include <omp.h>
#endif
//#include "oompi.h"
//#include <pthread.h>

// Misc
#include <stdlib.h>
//#include <sys/stat.h>
#include <cstring>
#include <vector>
#include <time.h>
#include <chrono>
#include <memory>
#include <exception>
#include <stdexcept>

//using namespace std;
/* Things from STD that we need always */
using std::cout;
using std::endl;
using std::cin;
using std::fstream;
using std::ostream;
using std::ios;
using std::nothrow;
using std::ifstream;

/* Things from Eigen that we always need */
using Eigen::Infinity;
using Eigen::Dynamic;
using Eigen::ColMajor;
using Eigen::RowMajor;
using Eigen::Upper;
using Eigen::Lower;

/* Things from BTAS that we always need */
using btas::Tensor;

/* Things from Libint that we always need */
using libint2::OneBodyEngine;
using libint2::TwoBodyEngine;

// Useful typedefs
typedef std::complex<double> dcomplex;
typedef Eigen::Matrix<double,Dynamic,Dynamic,RowMajor>     RealMatrix;    // RowMajor BC Libint
typedef Eigen::Matrix<dcomplex,Dynamic,Dynamic,RowMajor>   ComplexMatrix; // RowMajor BC Libint
typedef Eigen::MatrixExponentialReturnValue<RealMatrix>    RealMatExp;
typedef Eigen::MatrixExponentialReturnValue<ComplexMatrix> ComplexMatExp;
typedef TwoBodyEngine<libint2::Coulomb> coulombEngine;

//----------------//
//number constants//
//----------------//
struct Math {
  double zero, one, two, three, four, five, six, seven, eight, nine, ten, half, quarter;
  double sqrt2;
  double pi,pi32,sqrt2pi54; //pi, pi^{3/2} sqrt(2)*pi^{5/4}
  double small;
  dcomplex ii;
};
const Math math = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 0.5, 0.25,
		   1.4142135623731,
		   3.14159265358979,5.56832799683171,5.91496717279561,
		   1.0e-10, dcomplex(0,1.0)};
//factorials n!
static double Factorial[8] ={
1,1,2, 6, 24, 120, 720, 5040 };

//double factorials n!!
static double dFactorial[21] = {
  1.00000000,
  1.00000000,
  3.00000000,
  15.00000000,
  105.00000000,
  945.00000000,
  10395.00000000,
  135135.00000000,
  2027025.00000000,
  34459425.00000000,
  654729075.00000000,
  13749310575.00000000,
  316234143225.00000000,
  7905853580625.00000000,
  213458046676875.00000000,
  6190283353629375.00000000,
  191898783962510624.00000000,
  6332659870762850304.00000000,
  221643095476699758592.00000000,
  8200794532637890838528.00000000,
  319830986772877752139776.00000000,
};

//------------------//
//Physical constants//
//------------------//
struct Phys {
  //Bohr radius per Angstrom
  double bohr;
  //number of cartesian AO's in a shell
};
const Phys phys = {0.5291772083000001};


//------------------//
// IO block numbers //
//------------------//
enum {blockControlFlags,blockMolecule,blockBasisSet, // 1,2,3
      blockSingleSlater,blockIntegrals};  // 4,5,6

//----------//
//MPI global//
//----------//

struct GlobalMPI {
  int  myid;
  int  size;
//  char nodeName[MPI_MAX_PROCESSOR_NAME];;
  int  nodeNameLen;
};
//---------//
//MPI Tags //
//---------//
enum {tagMolecule,tagMatrix,tagBasisSet,tagSingleSlater,tagIntegrals,tagSDResponse};


#endif
