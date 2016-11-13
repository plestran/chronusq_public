/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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
#include <singleslater.h>

namespace ChronusQ{
/*
  template<>
  void SingleSlater<double>::printDensity(){
    if(this->nTCS_ == 1) {
      prettyPrint(this->fileio_->out,(*this->onePDMA_),"\u03B1-Density");
      if(!this->isClosedShell) 
        prettyPrint(this->fileio_->out,(*this->onePDMB_),"\u03B2-Density");
    } else
      prettyPrintTCS(this->fileio_->out,(*this->onePDMA_),"Density");
  };
  template<>
  void SingleSlater<double>::printFock(){
    if(this->nTCS_ == 1) {
      prettyPrint(this->fileio_->out,(*this->fockA_),"\u03B1-Fock");
      if(!this->isClosedShell) 
        prettyPrint(this->fileio_->out,(*this->fockB_),"\u03B2-Fock");
    } else
      prettyPrintTCS(this->fileio_->out,(*this->fockA_),"Fock");
  };
  template<>
  void SingleSlater<double>::printPT(){
    if(this->nTCS_ == 1) {
      prettyPrint(this->fileio_->out,(*this->PTA_),"\u03B1-Perturbation Tensor");
      if(!this->isClosedShell) 
        prettyPrint(this->fileio_->out,(*this->PTB_),"\u03B2-Perturbation Tensor");
    } else
      prettyPrintTCS(this->fileio_->out,(*this->PTA_),"Perturbation Tensor");
  };
*/

}; // namespace ChronusQ
