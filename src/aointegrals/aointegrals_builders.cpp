#include "aointegrals.h"
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
  clock_t start,finish;
  if(controls_->printLevel>=1) start = clock();
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
    finish = clock();
    this->fileio_->out<<"\nCPU time for one-electron integral:"<<(finish-start)/CLOCKS_PER_SEC<<" seconds."<<endl;
  };
  this->haveAOOneE = true;
};


