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
 *  Compute the shell block norm matrix for some AO quantity to
 *  be contracted to AO ERIs w/ screening
 *
 *  Each element is the infinity norm of the shell pair block
 */
template<>
void BasisSet::computeShBlkNorm(bool doBeta, int nTCS, const RealMatrix *DAlpha, 
                                   const RealMatrix *DBeta){
  // If map doesnt exist, make it
  if(!this->haveMapSh2Bf) this->makeMapSh2Bf(nTCS);

  // Allocate Matricies
  this->shBlkNormAlpha = 
    std::unique_ptr<RealMatrix>(new RealMatrix(this->nShell_,this->nShell_));
  if(doBeta)
    this->shBlkNormBeta = 
      std::unique_ptr<RealMatrix>(new RealMatrix(this->nShell_,this->nShell_));

  for(int s1 = 0; s1 < this->nShell_; s1++) {
    int bf1 = this->mapSh2Bf_[s1];
    int n1  = this->shells_[s1].size();
    for(int s2 = 0; s2 < this->nShell_; s2++) {
      int bf2 = this->mapSh2Bf_[s2];
      int n2  = this->shells_[s2].size();
     
      (*this->shBlkNormAlpha)(s1,s2) = DAlpha->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
      if(doBeta)
        (*this->shBlkNormBeta)(s1,s2) = DBeta->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
    }
  }
} // BasisSet::computeShBlkNorm (TMat = RealMatrix)

/**
 *  Compute the shell block norm matrix for some AO quantity to
 *  be contracted to AO ERIs w/ screening
 *
 *  Each element is the infinity norm of the shell pair block
 */
template<>
void BasisSet::computeShBlkNorm(bool doBeta,int nTCS,const ComplexMatrix *DAlpha, 
                                   const ComplexMatrix *DBeta){
  // If map doesnt exist, make it
  if(!this->haveMapSh2Bf) this->makeMapSh2Bf(nTCS);

  // Allocate Matricies
  this->shBlkNormAlpha = 
    std::unique_ptr<RealMatrix>(new RealMatrix(this->nShell_,this->nShell_));
  if(doBeta)
    this->shBlkNormBeta = 
      std::unique_ptr<RealMatrix>(new RealMatrix(this->nShell_,this->nShell_));

  for(int s1 = 0; s1 < this->nShell_; s1++) {
    int bf1 = this->mapSh2Bf_[s1];
    int n1  = this->shells_[s1].size();
    for(int s2 = 0; s2 < this->nShell_; s2++) {
      int bf2 = this->mapSh2Bf_[s2];
      int n2  = this->shells_[s2].size();
     
      (*this->shBlkNormAlpha)(s1,s2) = DAlpha->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
      if(doBeta)
        (*this->shBlkNormBeta)(s1,s2) = DBeta->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
    }
  }
} // BasisSet::computeShBlkNorm (TMat = ComplexMatrix)

/**
 *  Renormalize the libint2::Shell vector (this is important)
 */
void BasisSet::renormShells(){
  libint2::Engine engine(libint2::Operator::overlap,this->maxPrim_,this->maxL_,0);

  for(auto iShell = this->shells_.begin(); iShell != this->shells_.end(); ++iShell){
  //iShell->renorm();
    auto buff = engine.compute(*iShell,*iShell);
    for(auto k = 0; k < iShell->alpha.size(); k++)
      iShell->contr[0].coeff[k] /= std::sqrt(buff[0]);
  }
} // BasisSet::renormShells

/**
 * Uncontract the basis
 */
std::vector<libint2::Shell> BasisSet::uncontractBasis(){
  std::vector<libint2::Shell> newShells;

  for(auto iShell = this->shells_.begin(); iShell != this->shells_.end(); ++iShell) {
//  cout << " New Shell " << endl;
    for(auto i = 0; i < iShell->alpha.size(); ++i){
      newShells.push_back(libint2::Shell{
        { iShell->alpha[i] },
        { {iShell->contr[0].l,iShell->contr[0].pure,{1.0}  }},
        { {iShell->O[0],iShell->O[1],iShell->O[2]}}
      } );
//    cout << iShell->alpha[i] << " " << iShell->contr[0].l << endl;
    }
  }

  return newShells;
} // BasisSet::uncontractBasis


}; //namespace ChronusQ
