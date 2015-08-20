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
#ifndef INCLUDED_SDRESPONSE
#define INCLUDED_SDRESPONSE
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <controls.h>
#include <mointegrals.h>
#include <singleslater.h>
#include <basisset.h>

/****************************/
/* Error Messages 5000-5999 */
/****************************/

namespace ChronusQ {
template<typename T>
class SDResponse {
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TCMMatrix;
  typedef Eigen::Matrix<T,Dynamic,Dynamic,RowMajor> TMatrix;
  typedef Tensor<T,Range3d> TTensor3d;
  typedef Tensor<T,Range4d> TTensor4d;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TVec> TVecMap;
  typedef Eigen::Map<TCMMatrix> TCMMap;

  int       nBasis_;
  int       nTCS_;
  int       Ref_;
  int       nSek_;
  int       nGuess_;
  int       iMeth_;
  int       iPPRPA_;
  bool      haveDag_;
  int       nOA_;
  int       nVA_;
  int       nOB_;
  int       nVB_;
  int       nOAVA_;      // NOA * NVA
  int       nOBVB_;      // NOB * NVB
  int       nOAVB_;      // NOA * NVB
  int       nOBVA_;      // NOB * NVA
  int       nVAVA_SLT_;  // NVA * (NVA - 1) / 2
  int       nVBVB_SLT_;  // NVB * (NVB - 1) / 2
  int       nVAVA_LT_;   // NVA * (NVA + 1) / 2
  int       nVAVA_;      // NVA * NVA
  int       nVBVB_;      // NVB * NVB
  int       nOAOA_SLT_;  // NOA * (NOA - 1) / 2
  int       nOBOB_SLT_;  // NOB * (NOB - 1) / 2
  int       nOAOA_LT_;   // NOA * (NOA + 1) / 2
  int       nOAOA_;      // NOA * NOA
  int       nOBOB_;      // NOB * NOB
  int       nVAVB_;      // NVA * NVB
  int       nOAOB_;      // NOA * NOB
  int       nO_;         // NOA + NOB
  int       nV_;         // NVA + NVB
  int       nOV_;        // NO * NV
  int       nVV_SLT_;    // NV * (NV - 1) / 2
  int       nVV_LT_;     // NV * (NV + 1) / 2
  int       nVV_;        // NV * NV
  int       nOO_SLT_;    // NO * (NO - 1) / 2
  int       nOO_LT_;     // NO * (NO + 1) / 2
  int       nOO_;        // NO * NO

  int       nSingleDim_; // Single dimension of response matrix
  double    rMu_;

  std::unique_ptr<TCMMatrix>  transDen_;
  std::unique_ptr<RealMatrix> oscStrength_;
  std::unique_ptr<VectorXd>   omega_;
  std::unique_ptr<TTensor3d>  transDipole_;
  BasisSet        * basisSet_;
  Molecule        * molecule_;
  FileIO          * fileio_;
  Controls        * controls_;
  MOIntegrals<T>  * mointegrals_;
  SingleSlater<T> * singleSlater_;
  TTensor4d       * aoERI_;
  TTensor3d       * elecDipole_;

  std::unique_ptr<RealCMMatrix> rmDiag_;
  std::unique_ptr<TMatrix>  davGuess_;

public:
  enum{
    __invalid,
    CIS,
    RPA,
    PPRPA,
    PPATDA,
    PPCTDA,
    STAB,
    CCSD
  };
 
  // constructor & destructor
  SDResponse(){;};
  ~SDResponse() {;};
  // pseudo-constructor
  void iniSDResponse(Molecule *,BasisSet *,MOIntegrals<T> *,FileIO *,
                     Controls *, SingleSlater<T> *);

  #include <sdresponse_getset.h>
  #include <sdresponse_qnrelated.h>
  #include <sdresponse_io.h>
  #include <sdresponse_misc.h>
  #include <sdresponse_prop.h>

};
} // namespace ChronusQ
#endif
