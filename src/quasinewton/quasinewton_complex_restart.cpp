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
  void QuasiNewton<dcomplex>::setupRestart(){
    this->allocGuess();
    ComplexMap UR(this->URMem,this->N_,this->nGuess_);
    (*this->guessR_) = UR;
    if(this->symmetrizedTrial_){
      ComplexMap UL(this->ULMem,this->N_,this->nGuess_);
      (*this->guessL_) = UL;
    } 
    // Zero out scratch space
    std::memset(this->SCR,0.0,this->LenScr*sizeof(dcomplex));

    // Ensure that the the new guess vectors are orthonormal
    this->Orth(*this->guessR_);
    if(this->symmetrizedTrial_) this->Orth(*this->guessL_);  

    /** DO NOT RESET doRestart_ here! next iteration needs to know that we
        restarted **/
  } // setupRestart

}; // namespace ChronusQ
