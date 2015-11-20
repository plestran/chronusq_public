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
#include <aointegrals.h>
using ChronusQ::AOIntegrals;
//-----------------------------------------------//
// compute atomic orbital two-electron integrals //
//-----------------------------------------------//
void AOIntegrals::computeAOTwoE(){
  int ijShell,klShell, ij, kl, i, j, k, l, lA[3], lB[3], lC[3], lD[3], jBasis, kBasis, lBasis, iMax, jMax, kMax, lMax;
  int nPGTOs[4],totalL,nBasis;
  int nn,mm,nFmT;
  double X,C,divConst,twoeijkl,temp;
//  ConstQuartetNew *constQuartetNew = new ConstQuartetNew;
  nBasis = this->basisSet_->nBasis();
  ShellPair *ijShellPair, *klShellPair, *tmpShellPair;
  this->twoEX_->setZero();
  // integration over shell quartet on same shell pairs
  for(ijShell=0;ijShell<this->basisSet_->nShellPair();ijShell++) {
    ijShellPair = &(this->basisSet_->shellPairs[ijShell]);
    if((*ijShellPair).LTotal>=1) this->iniQuartetConstants(ijShellPair, ijShellPair);
    // first pass to integrate (ij|ij)
    for(ij=0;ij<(*ijShellPair).nAOPair;ij++) {
      i = (*ijShellPair).aoPairIndex[ij][0];
      j = (*ijShellPair).aoPairIndex[ij][1];
      k = (*ijShellPair).aoPairIndex[ij][0];
      l = (*ijShellPair).aoPairIndex[ij][1];
      nPGTOs[0] = ijShellPair->nPGTOs[0];
      nPGTOs[1] = ijShellPair->nPGTOs[1];
      nPGTOs[2] = ijShellPair->nPGTOs[0];
      nPGTOs[3] = ijShellPair->nPGTOs[1];
      twoeijkl = this->twoehRRabcd(nPGTOs,ijShellPair,ijShellPair,
  		 (*ijShellPair).L[0],(basisSet_->ao[i]).l,(*ijShellPair).L[1],(basisSet_->ao[j]).l,
	 	 (*ijShellPair).L[0],(basisSet_->ao[k]).l,(*ijShellPair).L[1],(basisSet_->ao[l]).l);
//		 /((*ijShellPair).divConst[ij]*(*ijShellPair).divConst[ij]);
      this->twoEC(i,j,i,j) = twoeijkl;

      if(i==k|j==l){
        this->twoEX(i,k,j,l)+= twoeijkl*math.two;
        this->twoEX(j,l,i,k)+= twoeijkl*math.two;
      }else{
        this->twoEX(i,k,j,l)+= twoeijkl;
        this->twoEX(j,l,i,k)+= twoeijkl;
      };
      if(i!=j&&k!=l) {
        if(i==l|j==k){
          this->twoEX(i,l,j,k)+= twoeijkl*math.two;
          this->twoEX(j,k,i,l)+= twoeijkl*math.two;
	}else{
          this->twoEX(i,l,j,k)+= twoeijkl;
          this->twoEX(j,k,i,l)+= twoeijkl;
	};
      };

    };
    // second pass to integrate (ij|kl) in a same shell
    for(ij=0;ij<(*ijShellPair).nAOPair;ij++) for(kl=ij+1;kl<(*ijShellPair).nAOPair;kl++) {
      i = (*ijShellPair).aoPairIndex[ij][0];
      j = (*ijShellPair).aoPairIndex[ij][1];
      k = (*ijShellPair).aoPairIndex[kl][0];
      l = (*ijShellPair).aoPairIndex[kl][1];
      if(std::abs(this->twoEC(i,j,i,j)*this->twoEC(k,l,k,l))>this->controls_->thresholdSchawrtz) {
        nPGTOs[0] = ijShellPair->nPGTOs[0];
        nPGTOs[1] = ijShellPair->nPGTOs[1];
        nPGTOs[2] = ijShellPair->nPGTOs[0];
	nPGTOs[3] = ijShellPair->nPGTOs[1];
        twoeijkl = this->twoehRRabcd(nPGTOs,ijShellPair,ijShellPair,
	  	   (*ijShellPair).L[0],(basisSet_->ao[i]).l,(*ijShellPair).L[1],(basisSet_->ao[j]).l,
		   (*ijShellPair).L[0],(basisSet_->ao[k]).l,(*ijShellPair).L[1],(basisSet_->ao[l]).l);
//		   /((*ijShellPair).divConst[ij]*(*ijShellPair).divConst[kl]);
        this->twoEC(i,j,k,l) = twoeijkl;
        this->twoEC(k,l,i,j) = twoeijkl;

	if(i==k|j==l){
  	  this->twoEX(i,k,j,l)+= twoeijkl*math.two;
  	  this->twoEX(j,l,i,k)+= twoeijkl*math.two;
	}else{
  	  this->twoEX(i,k,j,l)+= twoeijkl;
  	  this->twoEX(j,l,i,k)+= twoeijkl;
	};
	if(i!=j&&k!=l){
	  if(i==l|j==k){
            this->twoEX(i,l,j,k)+= twoeijkl*math.two;
            this->twoEX(j,k,i,l)+= twoeijkl*math.two;
	  }else{
            this->twoEX(i,l,j,k)+= twoeijkl;
            this->twoEX(j,k,i,l)+= twoeijkl;
	  };
	};

      };
    };
  };

  int iniFlag;
  // integration over shell quartet on different shell pairs
  for(ijShell=0;ijShell<this->basisSet_->nShellPair();ijShell++) for(klShell=ijShell+1;klShell<this->basisSet_->nShellPair();klShell++) {
    ijShellPair = &(this->basisSet_->shellPairs[ijShell]);
    klShellPair = &(this->basisSet_->shellPairs[klShell]);
    if((*klShellPair).LTotal>(*ijShellPair).LTotal){
      tmpShellPair = ijShellPair;
      ijShellPair  = klShellPair;
      klShellPair  = tmpShellPair;
    };
    iniFlag = 0;
    totalL = (*ijShellPair).LTotal+(*klShellPair).LTotal;
    for(ij=0;ij<(*ijShellPair).nAOPair;ij++) for(kl=0;kl<(*klShellPair).nAOPair;kl++) {
      i = (*ijShellPair).aoPairIndex[ij][0];
      j = (*ijShellPair).aoPairIndex[ij][1];
      k = (*klShellPair).aoPairIndex[kl][0];
      l = (*klShellPair).aoPairIndex[kl][1];
      if(std::abs(this->twoEC(i,j,i,j)*this->twoEC(k,l,k,l))>this->controls_->thresholdSchawrtz) {
        if(totalL>=1&&iniFlag==0) {
          this->iniQuartetConstants(ijShellPair, klShellPair);
          iniFlag = 1;
        };
	nPGTOs[0] = ijShellPair->nPGTOs[0];
	nPGTOs[1] = ijShellPair->nPGTOs[1];
	nPGTOs[2] = klShellPair->nPGTOs[0];
	nPGTOs[3] = klShellPair->nPGTOs[1];
	twoeijkl = this->twoehRRabcd(nPGTOs,ijShellPair,klShellPair,
		   (*ijShellPair).L[0],(basisSet_->ao[i]).l,(*ijShellPair).L[1],(basisSet_->ao[j]).l,
		   (*klShellPair).L[0],(basisSet_->ao[k]).l,(*klShellPair).L[1],(basisSet_->ao[l]).l);
//		   /((*ijShellPair).divConst[ij]*(*klShellPair).divConst[kl]);
        this->twoEC(i,j,k,l) = twoeijkl;
        this->twoEC(k,l,i,j) = twoeijkl;

	if(i==k|j==l){
  	  this->twoEX(i,k,j,l)+= twoeijkl*math.two;
  	  this->twoEX(j,l,i,k)+= twoeijkl*math.two;
	}else{
  	  this->twoEX(i,k,j,l)+= twoeijkl;
  	  this->twoEX(j,l,i,k)+= twoeijkl;
	};
	if(i!=j&&k!=l) {
  	  if(i==l|j==k){
            this->twoEX(i,l,j,k)+= twoeijkl*math.two;
            this->twoEX(j,k,i,l)+= twoeijkl*math.two;
	  }else{
            this->twoEX(i,l,j,k)+= twoeijkl;
            this->twoEX(j,k,i,l)+= twoeijkl;
	  };
	};
      };
    };
  };

  // Taking care of the diagonal elements of the exchange integral matrix
  for(ij=0;ij<nTT_;ij++) this->twoEX(ij,ij)=this->twoEX(ij,ij)*math.half;

  this->haveAOTwoE = true;
};

//-----------------------------------//
// compute one-electron integrals    //
//   overlap matrix                  //
//   kinetic energy matrix           //
//   nuclear-electron potential      //
//-----------------------------------//
void AOIntegrals::computeAOOneE(){
  /* overlap, electron-nuclear potential and kinetic energy integrals */
  double divSTV,S,T,V;
  int i,j,m,n,ij,ijShell,nPGTOs[2];
  ShellPair *ijShellPair;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(controls_->printLevel>=1) start = std::chrono::high_resolution_clock::now();
  this->iniMolecularConstants();
  for(ijShell=0;ijShell<basisSet_->nShellPair();ijShell++) {
    ijShellPair = &(basisSet_->shellPairs[ijShell]);
    this->iniPairConstants(ijShellPair);
    for(ij=0;ij<(*ijShellPair).nAOPair;ij++){
      i = (*ijShellPair).aoPairIndex[ij][0];
      j = (*ijShellPair).aoPairIndex[ij][1];
      divSTV = (*ijShellPair).divConst[ij];
      T = math.zero;
      if(std::abs(this->pairConstants_->ssPairTotal)<this->pairConstants_->intSmall) {
	S = math.zero;
	V = math.zero;
	T = math.zero;
      } else {
	if(i==j) S = math.one*divSTV;
	else S = this->oneehRRSab((*ijShellPair).L[0],(basisSet_->ao[i]).l,(*ijShellPair).L[1],(basisSet_->ao[j]).l);
	V = this->oneehRRVab((*ijShellPair).L[0],(basisSet_->ao[i]).l,(*ijShellPair).L[1],(basisSet_->ao[j]).l);
	if(std::abs(S)>this->pairConstants_->intSmall) for(m=0;m<this->pairConstants_->nPGTOs[0];m++) for(n=0;n<this->pairConstants_->nPGTOs[1];n++) 
	  if(this->pairConstants_->ssNonzero[m][n]) T += this->pairConstants_->ssPair[m][n]*this->oneevRRTab((*ijShellPair).L[0],(basisSet_->ao[i]).l,(*ijShellPair).L[1],(basisSet_->ao[j]).l,&m,&n);
	S = S/divSTV;
	V = V/divSTV;
	T = T/divSTV;
      };
      (*this->overlap_)(i,j) = S;
      (*this->potential_)(i,j) = V;
      (*this->kinetic_)(i,j) = T;
    };
  };
  (*this->oneE_) = (*this->kinetic_)-(*this->potential);
  finish = std::chrono::high_resolution_clock::now();
  this->OneED = finish - start;
  if(controls_->printLevel>=2) {
    prettyPrint(this->fileio_->out,(*this->overlap_),"Overlap");
    prettyPrint(this->fileio_->out,(*this->kinetic_),"Kinetic");
    prettyPrint(this->fileio_->out,(*this->potential_),"Potential");
    prettyPrint(this->fileio_->out,(*this->oneE_),"Core Hamiltonian");
  };
  if(controls_->printLevel>=1) {
//    finish = clock();
    this->fileio_->out<<"\nCPU time for one-electron integral:  "<< this->OneED.count() <<" seconds."<<endl;
  };
  this->haveAOOneE = true;
};


void AOIntegrals::computeOverlapS(){
  double S;
  int i,j,ijShell,iAOP;
  ShellPair_New *ijS;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(controls_->printLevel>=1) start = std::chrono::high_resolution_clock::now();
  this->overlap_->setZero();
  this->createShellPairs();
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);
    for(iAOP=0;iAOP<ijS->nAOPair;iAOP++) {
      i = ijS->aoPairIndex[0][iAOP];
      j = ijS->aoPairIndex[1][iAOP];
      if(i==j) S = math.one;
      else S = this->hRRSab(ijS,ijS->iShell->L,(basisSet_->ao[i]).l,ijS->jShell->L,(basisSet_->ao[j]).l);
      (*this->overlap_)(i,j) = S;
    };
  };
  if(controls_->printLevel>=1) {
    prettyPrint(this->fileio_->out,(*this->overlap_),"Overlap");
    this->fileio_->out<<"\nCPU time for overlap integral:  "<< this->OneED.count() <<" seconds."<<endl;
  };
};

void AOIntegrals::computeKineticT(){
  double T;
  int i,j,ijShell,iAOP;
  ShellPair_New *ijS;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(controls_->printLevel>=1) start = std::chrono::high_resolution_clock::now();
  this->kinetic_->setZero();
  this->createShellPairs();
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);
    for(iAOP=0;iAOP<ijS->nAOPair;iAOP++) {
      i = ijS->aoPairIndex[0][iAOP];
      j = ijS->aoPairIndex[1][iAOP];
      T = this->hRRSab(ijS,ijS->iShell->L,(basisSet_->ao[i]).l,ijS->jShell->L,(basisSet_->ao[j]).l);
      (*this->kinetic_)(i,j) = T;
    };
  };
  if(controls_->printLevel>=1) {
    prettyPrint(this->fileio_->out,(*this->kinetic_),"Kinetic");
    this->fileio_->out<<"\nCPU time for kinetic energy integral:  "<< this->OneED.count() <<" seconds."<<endl;
  };
};
