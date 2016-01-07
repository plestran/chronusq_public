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
#include <qn.h>

namespace ChronusQ {

  template<>
  void QuasiNewton2<double>::checkLinearDependence(int &NTrial){
    auto N = this->qnObj_->nSingleDim();
    char JOBU  = 'N';
    char JOBVT = 'N';
    int  INFO;

    (*this->out_) << "Performing SVD on Trial Vectors in QuasiNewton" << endl;
    this->writeTrialVectors(NTrial);
  //std::memcpy(this->URMem_,this->TRMem_,NTrial * N * sizeof(double));

    // URMem is also passed into U and VT as dummy variables, shouldn't
    // be touched, this assumes that JOBU = JOBVT = 'N'

  //dgesvd_(&JOBU,&JOBVT,&N,&NTrial,this->URMem_,&N,this->ERMem_,this->URMem_,
  //  &N,this->URMem_,&N,this->WORK,&this->LWORK,&INFO);
    dgesvd_(&JOBU,&JOBVT,&N,&NTrial,this->TRMem_,&N,this->ERMem_,this->URMem_,
      &N,this->URMem_,&N,this->WORK,&this->LWORK,&INFO);
    this->readTrialVectors(NTrial);

    if(INFO != 0) 
      CErr("DGESVD Failed in QuasiNewton2<double>::checkLinearDependence",
        (*this->out_));

    for(auto i = 0; i < NTrial; i++){
      if(std::abs(this->ERMem_[i]) < 1e-08) 
        CErr("QuasiNewton2 FATAL: Linear Dependency found in Trial Vectors");
    }
  }; // QuasiNewton2<double>::checkLinearDependence

  template<>
  void QuasiNewton2<double>::orthogonalize(int NTrial){
    auto N = this->qnObj_->nSingleDim();
    int INFO;

    (*this->out_) 
      << "Performing QR Decomposition on Trial Vectors in QuasiNewton" << endl;

    dgeqrf_(&N,&NTrial,this->TRMem_,&N,this->ERMem_,this->WORK,&this->LWORK,
      &INFO);

    if(INFO != 0) 
      CErr("DGEQRF Failed in QuasiNewton2<double>::orthogonalize",
        (*this->out_));

    dorgqr_(&N,&NTrial,&NTrial,this->TRMem_,&N,this->ERMem_,this->WORK,
      &this->LWORK,&INFO);

    if(INFO != 0) 
      CErr("DORGQR Failed in QuasiNewton2<double>::orthogonalize",
        (*this->out_));
    
  }; // QuasiNewton2<double>::orthogonalize

}; // namespace ChronusQ
