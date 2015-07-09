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
  /** Form Residual Based Guess **/
  template<>
  void QuasiNewton<double>::formResidualGuess(double Omega, 
                                     const RealCMMap & ResR, RealCMMap & QR, 
                                     const RealCMMap & ResL, RealCMMap & QL){
    if(this->symmetrizedTrial_) {
      for(auto i = 0; i < this->N_; i++){
        QR(i) = ResR(i) * (*this->diag_)(i,0);
        QL(i) = ResL(i) * (*this->diag_)(i,0);
      }
      for(auto i = 0; i < this->N_/2; i++){
        QR(i) += Omega*ResL(i);
        QR(this->N_/2 + i) -= Omega*ResL(this->N_/2 + i);
        QL(i) += Omega*ResR(i);
        QL(this->N_/2 + i) -= Omega*ResR(this->N_/2 + i);
      }
      for(auto i = 0; i < this->N_; i++){
        QR(i) = QR(i) / (std::pow((*this->diag_)(i,0),2.0)-std::pow(Omega,2.0));
        QL(i) = QL(i) / (std::pow((*this->diag_)(i,0),2.0)-std::pow(Omega,2.0));
      }
    } else {
      for(auto i = 0; i < this->N_; i++)
        QR(i) = -ResR(i) / ((*this->diag_)(i,0) - Omega);
    }
  }
    /** Form New Perturbed Guess **/
  template<>
  void QuasiNewton<double>::formNewGuess(std::vector<bool> &resConv,int &NTrial, 
                                    int NNotConv, int &NOld, int &NNew){

    RealCMMap TrialVecR(this->TVecRMem,0,0);
    RealCMMap TrialVecL(this->TVecLMem,0,0);
    RealCMMap ResR     (this->ResRMem, 0,0);
    RealCMMap ResL     (this->ResLMem, 0,0);
    RealCMMap QR       (this->TVecRMem,0,0);
    RealCMMap RR       (this->ResRMem ,0,0);
    RealCMMap QL       (this->TVecLMem,0,0);
    RealCMMap RL       (this->ResLMem ,0,0);

    new (&TrialVecR) RealCMMap(this->TVecRMem,this->N_,NTrial+NNotConv);
    new (&ResR) RealCMMap(this->ResRMem,this->N_,NTrial);
    if(this->symmetrizedTrial_ || !this->isHermetian_){
      new (&TrialVecL) RealCMMap(this->TVecLMem,this->N_,NTrial+NNotConv);
      new (&ResL) RealCMMap(this->ResLMem,this->N_,NTrial);
    }

    RealVecMap ER(this->ERMem,NTrial);
    int INDX = 0;
    for(auto k = 0; k < this->nSek_; k++) {
      // If the residual for root "k" is not converged, construct
      // a perturbed guess vector according to Davidson's diagonal
      // preconditioning scheme.
      //
      // **WARNING** This is only valid for strongly diagonal dominant
      //             matricies. Convergence will be slow (maybe infinitely)
      //             if this criteria is not met.
      if(!resConv[k]) {
        new (&RR) RealCMMap(this->ResRMem + k*this->N_,this->N_,1);
        new (&QR) RealCMMap(this->TVecRMem+(NTrial+INDX)*this->N_,this->N_,1);
        if(this->symmetrizedTrial_ || !this->isHermetian_){
          new (&RL) RealCMMap(this->ResLMem + k*this->N_,this->N_,1);
          new (&QL) RealCMMap(this->TVecLMem+(NTrial+INDX)*this->N_,this->N_,1);
        }
         
        this->formResidualGuess(ER(k),RR,QR,RL,QL);
        INDX++;
      }
    }
    // Normalize and orthogonalize the new guess vectors to the
    // existing set using QR factorization
    int N,M,LDA,INFO;
    double * AMATR, * AMATL;
    N = TrialVecR.cols();
    M = TrialVecR.rows();
    LDA = TrialVecR.rows();
    AMATR = TrialVecR.data();
    if(this->symmetrizedTrial_ || !this->isHermetian_) AMATL = TrialVecL.data();
    double *TAU = this->LAPACK_SCR;
    this->WORK = TAU + N;
  
    dgeqrf_(&M,&N,AMATR,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    dorgqr_(&M,&N,&N,AMATR,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    if(this->symmetrizedTrial_ || !this->isHermetian_){
      dgeqrf_(&M,&N,AMATL,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
      dorgqr_(&M,&N,&N,AMATL,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    }
    // Update number of vectors
    NOld = NTrial;
    NNew = NNotConv;
    NTrial += NNew;
  }

}; // namespace ChronusQ
