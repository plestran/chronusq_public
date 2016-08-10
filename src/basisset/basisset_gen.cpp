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
#include <basisset.h>
namespace ChronusQ{
/**
 *  Construct a local basis defintion using libint2::Shell struct
 *  from the reference shells
 */
void BasisSet::constructLocal(Molecule * mol){
  //cout << "Reference Shells" << endl;
  for(auto iAtom = 0; iAtom < mol->nAtoms(); iAtom++){
    bool found = false;
    for(auto iRef = this->refShells_.begin(); iRef != this->refShells_.end(); ++iRef){
      if(mol->index(iAtom) == (*iRef).index){
        for(auto iShell = (*iRef).shells.begin(); iShell != (*iRef).shells.end(); ++iShell){
          //cout << (*iShell) << endl;
/*
          this->shells_.push_back(
            libint2::Shell{ 
              iShell->alpha, 
              {iShell->contr[0]}, 
              {{ (*mol->cart())(0,iAtom),
                 (*mol->cart())(1,iAtom),
                 (*mol->cart())(2,iAtom)}}
            }
          );
*/
          this->shells_.push_back(libint2::Shell::unit());
          this->shells_.back() = (*iShell);
          this->shells_.back().O = 
            { (*mol->cart())(0,iAtom),
              (*mol->cart())(1,iAtom),
              (*mol->cart())(2,iAtom)};

          this->shellsCQ.push_back(this->shells_.back());
        };

        for(auto iCons = (*iRef).unNormalizedCons.begin(); 
            iCons != (*iRef).unNormalizedCons.end(); ++iCons){
          this->unNormCons_.push_back(*iCons);
        }
        found = true;
      }
    } // loop iRef
    if(!found)  
      CErr("Atomic Number " + 
             std::to_string(elements[mol->index(iAtom)].atomicNumber) +
             " not found in current Basis Set",
           this->fileio_->out);
  } // loop iAtom

  this->computeMeta();
//cout << "Construct Local Shells" << endl;
//for(auto i = 0; i < this->shells_.size(); i++) cout << this->shells_[i] << endl;
} // BasisSet::constructLocal

/**
 *  Construct an external basis defintion using the local reference shells
 */
void BasisSet::constructExtrn(Molecule * mol, BasisSet *genBasis){
  int nShell = 0;
  genBasis->fileio_ = this->fileio_;
  for(auto iAtom = 0; iAtom < mol->nAtoms(); iAtom++){
    bool found = false;
    for(auto iRef = this->refShells_.begin(); iRef != this->refShells_.end(); ++iRef){
      if(mol->index(iAtom) == (*iRef).index){
        for(auto iShell = (*iRef).shells.begin(); iShell != (*iRef).shells.end(); ++iShell){
/*
          genBasis->shells_.push_back(
            libint2::Shell{ 
              iShell->alpha, 
              {iShell->contr[0]}, 
              {{ (*mol->cart())(0,iAtom),
                 (*mol->cart())(1,iAtom),
                 (*mol->cart())(2,iAtom)}}
            }
          );
*/
          genBasis->shells_.push_back(libint2::Shell::unit());
          genBasis->shells_.back() = (*iShell);
          genBasis->shells_.back().O = 
            { (*mol->cart())(0,iAtom),
              (*mol->cart())(1,iAtom),
              (*mol->cart())(2,iAtom)};
          genBasis->shellsCQ.push_back(ChronusQ::ShellCQ{genBasis->shells_[nShell]});
          nShell++;
        }
        found = true;
      }
    } // Loop iRef
    if(!found)  
      CErr("Atomic Number " + 
             std::to_string(elements[mol->index(iAtom)].atomicNumber) +
             " not found in current Basis Set",
           this->fileio_->out);
  } // loop IAtom
  genBasis->computeMeta();
  genBasis->forceCart_ = this->forceCart_;

} // BasisSet::constructExtrn

/**
 *  Generate uncontracted basis definition from local shells
 */

void BasisSet::genUCvomLocal(BasisSet *genBasis){
  genBasis->fileio_ = this->fileio_;
  genBasis->shells_ = this->uncontractBasis();
  genBasis->computeMeta();
} // BasisSet::genUCausLocal

/**
 *  Compute BasisSet metadata
 */
void BasisSet::computeMeta(){
  this->nBasis_ = 0;
  this->nShell_     = this->shells_.size();
  this->nShellPair_ = this->nShell_ * (this->nShell_ + 1) / 2;
  
  for(auto iShell = this->shells_.begin(); iShell != this->shells_.end(); ++iShell){
    this->nBasis_ += (*iShell).size();
    auto L = (*iShell).contr[0].l;
    auto shPrim = (*iShell).alpha.size();  
    if( L      > this->maxL_   ) this->maxL_    = L     ;
    if( shPrim > this->maxPrim_) this->maxPrim_ = shPrim;
    this->nPrimitive_ += shPrim * iShell->size();
  } // loop iShell

  this->nLShell_ = std::vector<int>(this->maxL_+1,0);
  for(auto shell : this->shells_){
    this->nLShell_[shell.contr[0].l]++;
  } // loop shell

  // Always allocate scratch for basis eval and first derivatives
  this->basisEvalScr_.resize(4*(2*this->maxL_+2));

} // BasisSet::computeMeta
}; // namespace ChronusQ

