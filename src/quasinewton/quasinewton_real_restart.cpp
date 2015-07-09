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
#include <quasinewton.h>

namespace ChronusQ {
  template<>
  void QuasiNewton<double>::setupRestart(){
    this->allocGuess();
    RealCMMap UR(this->URMem,this->N_,this->nGuess_);
    (*this->guessR_) = UR;
    if(this->symmetrizedTrial_){
      RealCMMap UL(this->ULMem,this->N_,this->nGuess_);
      (*this->guessL_) = UL;
    } 
    // Zero out scratch space
    std::memset(this->SCR,0.0,this->LenScr*sizeof(double));

    // Ensure that the the new guess vectors are orthonormal
    int N = this->guessR_->cols();
    int M = this->guessR_->rows();
    int LDA = this->guessR_->rows();
    int INFO;
    double *AMATR = this->guessR_->data();
    double *AMATL;
    if(this->symmetrizedTrial_ || !this->isHermetian_) AMATL = this->guessL_->data();
    double *TAU = this->LAPACK_SCR;
    this->WORK = TAU + N;
    
    dgeqrf_(&M,&N,AMATR,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    dorgqr_(&M,&N,&N,AMATR,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    if(this->symmetrizedTrial_ || !this->isHermetian_){
      dgeqrf_(&M,&N,AMATL,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
      dorgqr_(&M,&N,&N,AMATL,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    }

    /** DO NOT RESET doRestart_ here! next iteration needs to know that we
        restarted **/
  } // setupRestart

}; // namespace ChronusQ
