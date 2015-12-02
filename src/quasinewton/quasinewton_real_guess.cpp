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
  void QuasiNewton<double>::genStdHerResGuess(double Omega, const RealMap &Res, RealMap &Q){
    for(auto i = 0; i < this->N_; i++)
      Q(i) = -Res(i) / ((*this->diag_)(i,0) - Omega);
  }

  template<>
  void QuasiNewton<double>::genSymmResGuess(double Omega, 
                                     const RealMap & ResR, RealMap & QR, 
                                     const RealMap & ResL, RealMap & QL){
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
  }

  /** Form Residual Based Guess **/
  template<>
  void QuasiNewton<double>::formResidualGuess(double Omega, 
                                     const RealMap & ResR, RealMap & QR, 
                                     const RealMap & ResL, RealMap & QL){
    if(this->symmetrizedTrial_) this->genSymmResGuess(  Omega,ResR,QR,ResL,QL);
    else                        this->genStdHerResGuess(Omega,ResR,QR        );
  }
    /** Form New Perturbed Guess **/
  template<>
  void QuasiNewton<double>::formNewGuess(std::vector<bool> &resConv,int &NTrial, 
                                    int NNotConv, int &NOld, int &NNew){

    RealMap TrialVecR(this->TVecRMem,0,0);
    RealMap TrialVecL(this->TVecLMem,0,0);
    RealMap ResR     (this->ResRMem, 0,0);
    RealMap ResL     (this->ResLMem, 0,0);
    RealMap QR       (this->TVecRMem,0,0);
    RealMap RR       (this->ResRMem ,0,0);
    RealMap QL       (this->TVecLMem,0,0);
    RealMap RL       (this->ResLMem ,0,0);

    new (&TrialVecR) RealMap(this->TVecRMem,this->N_,NTrial+NNotConv);
    new (&ResR) RealMap(this->ResRMem,this->N_,NTrial);
    if(this->symmetrizedTrial_ || !this->isHermitian_){
      new (&TrialVecL) RealMap(this->TVecLMem,this->N_,NTrial+NNotConv);
      new (&ResL) RealMap(this->ResLMem,this->N_,NTrial);
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
        new (&RR) RealMap(this->ResRMem + k*this->N_,this->N_,1);
        new (&QR) RealMap(this->TVecRMem+(NTrial+INDX)*this->N_,this->N_,1);
        if(this->symmetrizedTrial_ || !this->isHermitian_){
          new (&RL) RealMap(this->ResLMem + k*this->N_,this->N_,1);
          new (&QL) RealMap(this->TVecLMem+(NTrial+INDX)*this->N_,this->N_,1);
        }
         
        this->formResidualGuess(ER(k),RR,QR,RL,QL);
        INDX++;
      }
    }
    // Normalize and orthogonalize the new guess vectors to the
    // existing set using QR factorization
    this->Orth(TrialVecR);
    if(this->symmetrizedTrial_) this->Orth(TrialVecL);
    // Update number of vectors
    NOld = NTrial;
    NNew = NNotConv;
    NTrial += NNew;
  }

}; // namespace ChronusQ
