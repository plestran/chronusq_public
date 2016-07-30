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
  void QuasiNewton2<double>::readGuess(){
    (*this->out_) << "Reading the Guess in QuasiNetwon" << endl;
    auto N    = this->qnObj_->nSingleDim();
    auto NVec = this->qnObj_->nGuess();

    H5::DataSpace dataspace = this->qnObj_->guessFile()->getSpace();
    this->qnObj_->guessFile()->read(this->TRMem_,H5PredType<double>(),
      dataspace,dataspace);

    if(this->qnObj_->needsLeft())
      this->qnObj_->guessFile()->read(this->TLMem_,H5PredType<double>(),
        dataspace,dataspace);
  /*
    RealMap GUESS(this->TRMem_,N,NVec);
    prettyPrint(cout,GUESS,"GUESS");
  */
  };
}; // namespace ChronusQ
