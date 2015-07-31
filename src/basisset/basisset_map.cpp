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
#include <basisset.h>
namespace ChronusQ{
/**
 *  Generate shell index -> starting basis function map
 */
void BasisSet::makeMapSh2Bf(){
  auto n = 0;
  for(auto shell : this->shells_){
     this->mapSh2Bf_.push_back(n);
     n += shell.size();
  } // loop shell
  this->haveMapSh2Bf = true;
} // BasisSet::makeMapSh2Bf

/**
 *  Generate shell index -> atomic center map
 */
void BasisSet::makeMapSh2Cen(Molecule *mol){
  for(auto shell : this->shells_){
    for(auto iAtom = 0; iAtom < mol->nAtoms(); iAtom++){
      std::array<double,3> center = {{ (*mol->cart())(0,iAtom),
                                       (*mol->cart())(1,iAtom),
                                       (*mol->cart())(2,iAtom) }};
      if(shell.O == center){
        this->mapSh2Cen_.push_back(iAtom+1);
        break;
      }
    } // loop iAtom 
  } // loop shell
  this->haveMapSh2Cen = true;
} // BasisSet::makeMapSh2Cen

/**
 *  Generate atomic center index -> starting basis function map
 */
void BasisSet::makeMapCen2Bf(Molecule *mol){
  if(!this->haveMapSh2Bf ) this->makeMapSh2Bf();
  if(!this->haveMapSh2Cen) this->makeMapSh2Cen(mol);

  for(auto iAtm = 0; iAtm < mol->nAtoms(); iAtm++){
    auto nSize = 0;
    for(auto iShell = 0; iShell < this->nShell_; iShell++){
      if((iAtm+1) == this->mapSh2Cen_[iShell]) nSize += this->shells_[iShell].size();
    } // loop iShell
    auto iSt = -1;
    for(auto iShell = 0; iShell < this->nShell_; iShell++){
      if((iAtm+1) == this->mapSh2Cen_[iShell]){
       iSt = this->mapSh2Bf_[iShell];
       break;
      }
    } // loop iShell
    if(iSt == -1) CErr("Could not find Center in Basis definition",this->fileio_->out);
    this->mapCen2Bf_.push_back({{iSt,nSize}});
  } // loop iAtm
  this->haveMapCen2Bf = true;
} // BasisSet::makeMapCen2Bf

}; // namespace ChronusQ
