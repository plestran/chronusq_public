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
#include "aointegrals.h"
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
  this->twoEX_->clearAll();
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

#ifndef USE_LIBINT // Only use Libint OneE for now
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
//clock_t start,finish;
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
  oneE_->sub(kinetic_,potential_);
  if(controls_->printLevel>=2) {
    this->overlap_->printAll(5,fileio_->out);
    this->kinetic_->printAll(5,fileio_->out);
    this->potential_->printAll(5,fileio_->out);
    this->oneE_->printAll(5,fileio_->out);
  };
  if(controls_->printLevel>=1) {
//    finish = clock();
    finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> OneED = finish - start;
    this->fileio_->out<<"\nCPU time for one-electron integral:  "<< OneED.count() <<" seconds."<<endl;
  };
  this->haveAOOneE = true;
};
#else // USE_LIBINT
using libint2::OneBodyEngine;

void AOIntegrals::OneEDriver(OneBodyEngine::integral_type iType) {
  // Not parallelizing yet
  //


  Matrix<double> *mat;
  if(iType == OneBodyEngine::overlap){
    mat = this->overlap_;
  } else if(iType == OneBodyEngine::kinetic) {
    mat = this->kinetic_;
  } else if(iType == OneBodyEngine::nuclear) {
    mat = this->potential_;
  } else {
    cout << "OneBodyEngine type not recognized" << endl;
    exit(EXIT_FAILURE);
  }
 
  // Check to see if the basisset had been converted
  if(!this->basisSet_->convToLI) this->basisSet_->convShell(this->molecule_);

  // Define integral Engine
  OneBodyEngine engine = OneBodyEngine(iType,this->basisSet_->maxPrim,
                                       this->basisSet_->maxL,0);
  // If engine is V, define nuclear charges
  if(iType == OneBodyEngine::nuclear){
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(int i = 0; i < this->molecule_->nAtoms(); i++) {
      q.push_back(
        {
          static_cast<double>((*this->molecularConstants_).atomZ[i]), 
          {
            {
	      (*this->molecularConstants_).cart[0][i],
	      (*this->molecularConstants_).cart[1][i],
	      (*this->molecularConstants_).cart[2][i]
	    }
	  }
	}
      );
    }
    engine.set_q(q);
  }
 
 std::vector<int> mapSh2Bf;
 int n = 0;
 for( auto shell: this->basisSet_->shells_libint) {
   mapSh2Bf.push_back(n);
   n += shell.size();
 }

  for(int s1=0; s1 < this->basisSet_->nShell(); s1++){
    int bf1 = mapSh2Bf[s1];
    int n1  = this->basisSet_->shells_libint[s1].size();
    for(int s2=0; s2 <= s1; s2++){
      int bf2 = mapSh2Bf[s2];
      int n2  = this->basisSet_->shells_libint[s2].size();
 
      const double* buff = engine.compute(
        this->basisSet_->shells_libint[s1],
        this->basisSet_->shells_libint[s2]
      );
      int ij = 0;
      for(int i = 0; i < n1; i++) {
        for(int j = 0; j < n2; j++) {
          (*mat)(bf1+i,bf2+j) = buff[ij];
 	 ij++;
        }
      }
    }
  }
  if(this->controls_->printLevel>=2)  mat->printAll(5,fileio_->out);

}

void AOIntegrals::computeAOOneE(){
  // Collect Relevant data into a struct (odd, but convienient) 
  this->iniMolecularConstants();

  // Start timer for one-electron integral evaluation
  auto oneEStart = std::chrono::high_resolution_clock::now();

  // Compute and time overlap integrals
  auto OStart = std::chrono::high_resolution_clock::now();
  OneEDriver(OneBodyEngine::overlap);
  auto OEnd = std::chrono::high_resolution_clock::now();

  // Compute and time kinetic integrals
  auto TStart = std::chrono::high_resolution_clock::now();
  OneEDriver(OneBodyEngine::kinetic);
  auto TEnd = std::chrono::high_resolution_clock::now();

  // Compute and time nuclear attraction integrals (negative sign is factored in)
  auto VStart = std::chrono::high_resolution_clock::now();
  OneEDriver(OneBodyEngine::nuclear);
  auto VEnd = std::chrono::high_resolution_clock::now();
  this->oneE_->add(this->kinetic_,this->potential_);

  // Get end time of one-electron integral evaluation
  auto oneEEnd = std::chrono::high_resolution_clock::now();
  if(this->controls_->printLevel>=2) this->oneE_->printAll(5,fileio_->out);

  // Compute time differenes
  std::chrono::duration<double> OneED = oneEEnd - oneEStart;
  std::chrono::duration<double> SED = OEnd - OStart;
  std::chrono::duration<double> TED = TEnd - TStart;
  std::chrono::duration<double> VED = VEnd - VStart;
  if(this->controls_->printLevel >= 1) {
    this->fileio_->out << endl;
    this->fileio_->out << std::left << std::setw(60) << "CPU time for Overlap evaluation:" 
                       << std::left << std::setw(15) << SED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << "CPU time for Kinetic evaluation:" 
                       << std::left << std::setw(15) << TED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << "CPU time for Nuclear Attraction Potential evaluation:" 
                       << std::left << std::setw(15) << VED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << " "
                       << std::left << std::setw(15) << "---------------" << "----" << endl;
    this->fileio_->out << std::left << std::setw(60) << "Total CPU time for one-electron integral evaluation:" 
                       << std::left << std::setw(15) << OneED.count() << " sec" << endl;
  this->haveAOOneE = true;
  }
}

using libint2::TwoBodyEngine;
void AOIntegrals::computeSchwartz(){
  ChronusQ::Matrix<double> *ShBlk; 
  this->schwartz_->clearAll();
  // Check to see if the basisset had been converted
  if(!this->basisSet_->convToLI) this->basisSet_->convShell(this->molecule_);

  // Define Integral Engine
  TwoBodyEngine<libint2::Coulomb> engine = 
    TwoBodyEngine<libint2::Coulomb>(this->basisSet_->maxPrim,
                                    this->basisSet_->maxL,0);
  engine.set_precision(0.); // Don't screen primitives during schwartz

  this->fileio_->out << "Computing Schwartz Bound Tensor ... " << endl;

  for(int s1=0; s1 < this->basisSet_->nShell(); s1++){
    int n1  = this->basisSet_->shells_libint[s1].size();
    for(int s2=0; s2 <= s1; s2++){
      int n2  = this->basisSet_->shells_libint[s2].size();
 
      const double* buff = engine.compute(
        this->basisSet_->shells_libint[s1],
        this->basisSet_->shells_libint[s2],
        this->basisSet_->shells_libint[s1],
        this->basisSet_->shells_libint[s2]
      );
      
      try{
        ShBlk = new ChronusQ::Matrix<double>(n1,n2);
      } catch (int msg) {cout << "couldn't allocate shblk" << endl;}
   
      int ij = 0;
      for(int i = 0; i < n1; i++) {
        for(int j = 0; j < n2; j++) {
          (*ShBlk)(i,j) = buff[ij*n1*n2 + ij];
 	 ij++;
        }
      }

      (*this->schwartz_)(s1,s2) = std::sqrt(ShBlk->infNorm());
      
      delete ShBlk;
    }
  }

  this->schwartz_->printAll();
  this->haveSchwartz = true;
  exit(EXIT_FAILURE);
}

#endif


