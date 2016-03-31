/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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

#define MAXNAMELEN 50 ///< Define maximum length for char arrays
#define MAXANGULARMOMENTUM 6 ///< Define maximum allowed total angular momentum for basis functions
#define MAXCONTRACTION 10 ///< Define maximum contraction depth for gaussian basis functions
#define MAXATOMS 1000 ///< Define maximum number of allowed nuclei

// Maximum number of FileIO Sratch Partitions (can override)
#ifndef CQ_MAX_SCRATCH_PARTITIONS
  #define CQ_MAX_SCRATCH_PARTITIONS 1000
#endif

// CMake Compilation Configuration
#include <config_chronusq.h>

// Boost Headers
#include <boost/python.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/algorithm/string.hpp>


// IO
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

// Math
#include <cmath>
#include <complex>
#define EIGEN_MATRIXBASE_PLUGIN "eigenplugin.h" ///< ChronusQ plugin for Eigen
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO // Eröffne die Matrizen als Null.
#include <Eigen/Core> // Eigen Linear Algebra
#include <Eigen/SparseCore> // Eigen Sparse Linear Algebra
#include <unsupported/Eigen/MatrixFunctions>
#ifdef CQ_ENABLE_MPI
#  define EIGEN_DONT_PARALLELIZE
#endif
#ifdef USE_LIBINT
#  include <libint2.hpp> // Libint Gaussian Integrals library
#endif
#include <eiginterface.h> // ChronusQ interface (TODO consolidate into plugin)
#include <btas/btas.h> // BTAS Tensor Algebra library (header only)

// Parallelization
#ifdef _OPENMP
#  include <omp.h>
#endif
#ifdef CQ_ENABLE_MPI
#  include <mpi.h>
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
#include <algorithm>
#include <gauinterface.h>
#include <iterator>
#include <typeinfo>
#include <H5Cpp.h>

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
using Eigen::Upper;
using Eigen::Lower;
using Eigen::VectorXd;
using Eigen::VectorXcd;

/* Things from BTAS that we always need */
using btas::Tensor;

// Alias for Boost::Geometry
namespace bg = boost::geometry;

// Useful typedefs
typedef std::complex<double> dcomplex; ///< Support for complex numbers (double precision)
typedef Eigen::Matrix<double,Dynamic,Dynamic,ColMajor>     RealMatrix;    ///< Dynamically allocated Real (double) matrix. Row major for integration with Libint
typedef Eigen::Matrix<dcomplex,Dynamic,Dynamic,ColMajor>   ComplexMatrix; ///< Dynamically allocated Complex (dcomplex) matrix. Row major for integration with Libint
typedef Eigen::SparseMatrix<double> RealSparseMatrix;
typedef Eigen::SparseMatrix<dcomplex> ComplexSparseMatrix;
typedef Eigen::Map<VectorXd> RealVecMap;
//AP
typedef Eigen::SelfAdjointView<RealMatrix,Lower> RealSelfAdjLow;
//APE
typedef Eigen::Map<VectorXcd> ComplexVecMap;
typedef Eigen::Map<RealMatrix> RealMap; ///< Map double precision real array onto RealMatrix object
typedef Eigen::Map<const RealMatrix> ConstRealMap; ///< Map double precision real array onto const RealMatrix object
typedef Eigen::Map<ComplexMatrix> ComplexMap; ///< Map double precision complex array onto ComplexMatrix object
typedef Eigen::Map<const ComplexMatrix> ConstComplexMap; ///< Map double precision complex array onto const ComplexMatrix object
typedef Eigen::MatrixExponentialReturnValue<RealMatrix>    RealMatExp; ///< Driver for matrix exponentaial (RealMatrix)
typedef Eigen::MatrixExponentialReturnValue<ComplexMatrix> ComplexMatExp; ///< Driver for matrix exponential (ComplexMatrix)
typedef btas::RangeNd<CblasColMajor,std::array<long,4>> Range4d; ///< BTAS range specification for rank-4 tensors
typedef btas::RangeNd<CblasColMajor,std::array<long,3>> Range3d; ///< BTAS range specification for rank-3 tensors
typedef btas::RangeNd<CblasColMajor,std::array<long,2>> Range2d; ///< BTAS range specification for rank-2 tensors (isomorphic with matrix)
typedef btas::RangeNd<CblasColMajor,std::array<long,1>> Range1d; ///< BTAS range specification for rank-1 tensors (isomorphic with vector)
typedef Tensor<double,Range4d> RealTensor4d; ///< Support for real-valued rank-4 tensors using BTAS
typedef Tensor<double,Range3d> RealTensor3d; ///< Support for real-valued rank-3 tensors using BTAS
typedef Tensor<double,Range2d> RealTensor2d; ///< Support for real-values rank-2 tensors (aka Matricies) using BTAS
typedef Tensor<double,Range1d> RealTensor1d; ///< Support for real-values rank-1 tensors (aka Vectors) using BTAS
typedef Tensor<dcomplex,Range4d> ComplexTensor4d; ///< Support for complex-valued rank-4 tensors using BTAS
typedef Tensor<dcomplex,Range3d> ComplexTensor3d; ///< Support for complex-valued rank-3 tensors using BTAS
typedef Tensor<dcomplex,Range2d> ComplexTensor2d; ///< Support for complex-values rank-2 tensors (aka Matricies) using BTAS
typedef Tensor<dcomplex,Range1d> ComplexTensor1d; ///< Support for complex-values rank-1 tensors (aka Vectors) using BTAS


typedef bg::model::point< double, 3, bg::cs::spherical<bg::radian> > sph3GP; ///< 3 Coordinate Spherical (w varying radius) (phi,theta,radius)
typedef bg::model::point< double, 2, bg::cs::spherical<bg::radian> > sph2GP; ///< 2 Coordinate Spherical (unit sphere)      (phi,theta) (Azimut[0,2pi],Elevation[0,pi]
typedef bg::model::point< double, 3, bg::cs::cartesian > cartGP;             ///< 3 Coordinate Carteaisn                    (x,y,z)

//----------------//
//number constants//
//----------------//
/**
 * Various numerical constants.
 */
struct Math {
  double zero, one, two, three, four, five, six, seven, eight, nine, ten, half, quarter;
  double sqrt2;
  double pi; ///< Mathematical constant \f$\pi\f$
  double pi32; ///< Mathematical constant \f$\pi^{3/2}\f$
  double sqrt2pi54; ///< Mathematical constant \f$\sqrt{2}\pi^{5/4}\f$
  double small;
  dcomplex ii; ///< Imaginary unit (0 + 1i)
};
const Math math = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 0.5, 0.25,
		   1.4142135623731,
		   boost::math::constants::pi<double>(),5.56832799683171,5.91496717279561,
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

/**
 *  Physical constants (considering phasing out to use Boost phyiscal constants)
 */
struct Phys {
  double bohr; ///< Bohr radii per Angstrom
  double debye; ///< e*bohr in 1 Debye
  double eVPerHartree;
  double nmPerHartree;
  double AuToFs;
  double SPEED_OF_LIGHT;
  //number of cartesian AO's in a shell
};
const Phys phys = {0.5291772083000001,0.393430307,27.211396132,45.56335,0.02418884326505,
                   137.035999139};

#include <clapack.h> // Extern "C" defs for LAPACK routines (require "_" extension)

#endif
