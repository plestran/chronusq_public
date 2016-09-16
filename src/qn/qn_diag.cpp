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
#include <qn.h>

namespace ChronusQ {

template<>
void QuasiNewton2<double>::checkImaginary(const int N, double *EI){
  double xSmall = 1e-08;
  for(auto i = 0; i < N; i++)
    if(std::abs(EI[i]) < xSmall)
      CErr("Imaginary Eigenroot has been found in QuasiNewton",(*this->out_));
}; // QuasiNewton2<double>::checkImaginary

template<>
void QuasiNewton2<double>::reducedDimDiag(const int NTrial){
  char JOBVR = 'V';
  char JOBVL = 'N';
  char UPLO  = 'L';
  bool INFO;

  (*this->out_) << "Diagonalizing Projected Matrix in QuasiNewton" << endl;
  int N = NTrial;
  if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL) N *= 2;

  double *A  = this->XTSigmaRMem_;
  double *VR = this->XTSigmaRMem_;
  double *VL = this->XTSigmaRMem_;

  if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL) {
    CErr();
    this->formNHrProd(NTrial);
    A = this->NHrProdMem_;
    VR = this->SSuperMem_;
    VL = this->SSuperMem_;
    CErr();
  }


  // Non Hermetian, real diagonalization requires 2x space for eigenvalues
  // for real and imaginary parts. This is separated in stdNonHermetianDiag
  double *EMem;
  if(this->qnObj_->matrixType_ == HERMETIAN)
    EMem = this->memManager_->malloc<double>(N);
  else
    EMem = this->memManager_->malloc<double>(2*N);
  
  if(this->qnObj_->matrixType_ == HERMETIAN)
  //dsyev_(&JOBVR,&UPLO,&N,A,&N,this->ERMem_,this->WORK,&this->LWORK,&INFO);
    INFO = this->stdHermetianDiag(JOBVR,UPLO,N,A,EMem);
  else {
  //dgeev_(&JOBVR,&JOBVL,&N,A,&N,this->ERMem_,this->EIMem_,VL,&N,VR,&N,
  //  this->WORK,&this->LWORK,&INFO);
    INFO = this->stdNonHermetianDiag(JOBVR,JOBVL,N,A,EMem,VR,VL);
    this->checkImaginary(N,EMem+N);
    RealVecMap EMap(EMem,N);
    RealMap    VRMap(VR,N,N);
    this->eigSrt(VRMap,EMap);
  };

  std::copy_n(EMem,N,this->EPersist_);

  if(this->qnObj_->matrixType_ == HERMETIAN)
    this->memManager_->free(EMem,N);
  else
    this->memManager_->free(EMem,2*N);

  if(!INFO) CErr("Diagonalization in Reduced Dimension Failed!",
                  (*this->out_));
}; //QuasiNewton2<double>::reducedDimDiag
}; // namespace ChronusQ
