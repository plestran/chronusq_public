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

#ifdef USE_LIBINT

using ChronusQ::BasisSet;
using ChronusQ::Molecule;
using ChronusQ::HashL;
typedef ChronusQ::Shell CShell;
typedef libint2::Shell LIShell;

void BasisSet::convShell(Molecule* mol) {
  std::vector<double> coeff;
  std::vector<double> exp;
  std::array<double,3> center;
  for(auto i = 0; i < this->nShell_; ++i) {
    for(auto iPrim = 0; iPrim < this->shells[i].nPGTOs; ++iPrim){
      coeff.push_back(this->shells[i].coef[iPrim]);
      exp.push_back(this->shells[i].expo[iPrim]);
    }
    center = {{(*mol->cart())(0,this->shells[i].center),
              (*mol->cart())(1,this->shells[i].center),
              (*mol->cart())(2,this->shells[i].center)}};

    int L = HashL(&this->shells[i].name[0]);

    this->shells_libint.push_back(
       LIShell{
         exp,
	 {
           {L, false, coeff}
	 },
	 center
       }
    );
    coeff.resize(0);
    exp.resize(0);
    shells_libint[i].renorm();
    if(i==0) {
      this->maxPrim = this->shells[i].nPGTOs;
      this->maxL = L;
    } else {
      if(this->shells[i].nPGTOs > this->maxPrim)this->maxPrim = this->shells[i].nPGTOs;
      if(L > this->maxL) this->maxL = L;
    }
  }
  this->convToLI = true;
}

void BasisSet::makeMap(Molecule * mol) {
  if(!this->convToLI) this->convShell(mol);
  int n = 0;
  for( auto shell: this->shells_libint) {
    this->mapSh2Bf.push_back(n);
    n += shell.size();
  }
  this->haveMap = true;
}

void BasisSet::computeShBlkNorm(Molecule *mol, RealMatrix *D){
  // This will be much easier in Eigen
  if(!this->convToLI) this->convShell(mol);
  if(!this->haveMap)  this->makeMap(mol);

  this->shBlkNorm = new RealMatrix(this->nShell(),this->nShell());
  for(int s1 = 0; s1 < this->nShell(); s1++) {
    int bf1 = this->mapSh2Bf[s1];
    int n1  = this->shells_libint[s1].size();
    for(int s2 = 0; s2 < this->nShell(); s2++) {
      int bf2 = this->mapSh2Bf[s2];
      int n2  = this->shells_libint[s2].size();
     
      (*this->shBlkNorm)(s1,s2) = D->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
    }
  }
}
#endif
