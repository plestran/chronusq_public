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
#ifndef USE_LIBINT // Libint is default integral driver
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
//oneE_->sub(kinetic_,potential_);
  (*this->oneE_) = (*this->kinetic_)-(*this->potential);
  finish = std::chrono::high_resolution_clock::now();
  this->OneED = finish - start;
  if(controls_->printLevel>=2) {
/*
    this->overlap_->printAll(5,fileio_->out);
    this->kinetic_->printAll(5,fileio_->out);
    this->potential_->printAll(5,fileio_->out);
    this->oneE_->printAll(5,fileio_->out);
*/
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
#else // USE_LIBINT

void AOIntegrals::OneEDriver(OneBodyEngine::integral_type iType) {

  std::vector<RealMap> mat;
  int NB = this->nBasis_;
  int NBSq = NB*NB;
  if(iType == OneBodyEngine::overlap){
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
  } else if(iType == OneBodyEngine::kinetic) {
    mat.push_back(RealMap(this->kinetic_->data(),NB,NB));
  } else if(iType == OneBodyEngine::nuclear) {
    mat.push_back(RealMap(this->potential_->data(),NB,NB));
  } else if(iType == OneBodyEngine::emultipole1) {
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
    mat.push_back(RealMap(&this->elecDipole_->storage()[0],NB,NB));
    mat.push_back(RealMap(&this->elecDipole_->storage()[NBSq],NB,NB));
    mat.push_back(RealMap(&this->elecDipole_->storage()[2*NBSq],NB,NB));
  } else if(iType == OneBodyEngine::emultipole2) {
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
    mat.push_back(RealMap(&this->elecDipole_->storage()[0],NB,NB));
    mat.push_back(RealMap(&this->elecDipole_->storage()[NBSq],NB,NB));
    mat.push_back(RealMap(&this->elecDipole_->storage()[2*NBSq],NB,NB));
    mat.push_back(RealMap(&this->elecQuadpole_->storage()[0],NB,NB));
    mat.push_back(RealMap(&this->elecQuadpole_->storage()[NBSq],NB,NB));
    mat.push_back(RealMap(&this->elecQuadpole_->storage()[2*NBSq],NB,NB));
    mat.push_back(RealMap(&this->elecQuadpole_->storage()[3*NBSq],NB,NB));
    mat.push_back(RealMap(&this->elecQuadpole_->storage()[4*NBSq],NB,NB));
    mat.push_back(RealMap(&this->elecQuadpole_->storage()[5*NBSq],NB,NB));
  } else {
    cout << "OneBodyEngine type not recognized" << endl;
    exit(EXIT_FAILURE);
  }
 
  // Check to see if the basisset had been converted
  if(!this->basisSet_->convToLI) this->basisSet_->convShell(this->molecule_);
  // Define integral Engine
  std::vector<OneBodyEngine> engines(this->controls_->nthreads);
  engines[0] = OneBodyEngine(iType,this->basisSet_->maxPrim,this->basisSet_->maxL,0);
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
    engines[0].set_q(q);
  }
  for(size_t i = 1; i < this->controls_->nthreads; i++) engines[i] = engines[0];

  if(!this->basisSet_->haveMap) this->basisSet_->makeMap(this->molecule_); 
#ifdef USE_OMP
  #pragma omp parallel
#endif
  {
#ifdef USE_OMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
    for(auto s1=0l, s12=0l; s1 < this->basisSet_->nShell(); s1++){
      int bf1 = this->basisSet_->mapSh2Bf[s1];
      int n1  = this->basisSet_->shells_libint[s1].size();
      for(int s2=0; s2 <= s1; s2++, s12++){
        if(s12 % this->controls_->nthreads != thread_id) continue;
        int bf2 = this->basisSet_->mapSh2Bf[s2];
        int n2  = this->basisSet_->shells_libint[s2].size();
  
        const double* buff = engines[thread_id].compute(
          this->basisSet_->shells_libint[s1],
          this->basisSet_->shells_libint[s2]
        );
/*
        for(int i = 0, ij=0; i < n1; i++) {
          for(int j = 0; j < n2; j++, ij++) {
            (*mat)(bf1+i,bf2+j) = buff[ij];
          }
        }
*/
        int IOff = 0;
        for(auto nMat = 0; nMat < mat.size(); nMat++) {
          ConstRealMap bufMat(&buff[IOff],n1,n2); // Read only map
          mat[nMat].block(bf1,bf2,n1,n2) = bufMat;
          IOff += n1*n2;
        }
      }
    }
  } // end openmp parallel
  for(auto nMat = 0; nMat < mat.size(); nMat++) 
    mat[nMat] = mat[nMat].selfadjointView<Lower>();

//if(this->controls_->printLevel>=2)  mat->printAll(5,fileio_->out);
  if(this->controls_->printLevel>=2){
    if(iType == OneBodyEngine::overlap){
      prettyPrint(this->fileio_->out,(mat[0]),"Overlap");
    } else if(iType == OneBodyEngine::kinetic) {
      prettyPrint(this->fileio_->out,(mat[0]),"Kinetic");
    } else if(iType == OneBodyEngine::nuclear) {
      prettyPrint(this->fileio_->out,(mat[0]),"Potential");
    } else if(iType == OneBodyEngine::emultipole1) {
      prettyPrint(this->fileio_->out,(mat[0]),"Overlap");
      prettyPrint(this->fileio_->out,(mat[1]),"Electric Dipole (x)");
      prettyPrint(this->fileio_->out,(mat[2]),"Electric Dipole (y)");
      prettyPrint(this->fileio_->out,(mat[3]),"Electric Dipole (z)");
    } else if(iType == OneBodyEngine::emultipole2) {
      prettyPrint(this->fileio_->out,(mat[0]),"Overlap");
      prettyPrint(this->fileio_->out,(mat[1]),"Electric Dipole (x)");
      prettyPrint(this->fileio_->out,(mat[2]),"Electric Dipole (y)");
      prettyPrint(this->fileio_->out,(mat[3]),"Electric Dipole (z)");
      prettyPrint(this->fileio_->out,(mat[4]),"Electric Quadrupole (xx)");
      prettyPrint(this->fileio_->out,(mat[5]),"Electric Quadrupole (xy)");
      prettyPrint(this->fileio_->out,(mat[6]),"Electric Quadrupole (xz)");
      prettyPrint(this->fileio_->out,(mat[7]),"Electric Quadrupole (yy)");
      prettyPrint(this->fileio_->out,(mat[8]),"Electric Quadrupole (yz)");
      prettyPrint(this->fileio_->out,(mat[9]),"Electric Quadrupole (zz)");
    } else {
      cout << "OneBodyEngine type not recognized" << endl;
      exit(EXIT_FAILURE);
    }
  }

}

void AOIntegrals::computeAOOneE(){
  // Collect Relevant data into a struct (odd, but convienient) 
  this->iniMolecularConstants();

  // Start timer for one-electron integral evaluation
  auto oneEStart = std::chrono::high_resolution_clock::now();

  // Compute and time overlap integrals
  auto OStart = std::chrono::high_resolution_clock::now();
  if(this->controls_->doQuadpole) OneEDriver(OneBodyEngine::emultipole2);
  else if(this->controls_->doDipole) OneEDriver(OneBodyEngine::emultipole1);
  else OneEDriver(OneBodyEngine::overlap);
  auto OEnd = std::chrono::high_resolution_clock::now();

  // Compute and time kinetic integrals
  auto TStart = std::chrono::high_resolution_clock::now();
  OneEDriver(OneBodyEngine::kinetic);
  auto TEnd = std::chrono::high_resolution_clock::now();

  // Compute and time nuclear attraction integrals (negative sign is factored in)
  auto VStart = std::chrono::high_resolution_clock::now();
  OneEDriver(OneBodyEngine::nuclear);
  auto VEnd = std::chrono::high_resolution_clock::now();
//this->oneE_->add(this->kinetic_,this->potential_);
  (*this->oneE_) = (*this->kinetic_) + (*this->potential_);

  // Get end time of one-electron integral evaluation
  auto oneEEnd = std::chrono::high_resolution_clock::now();
//if(this->controls_->printLevel>=2) this->oneE_->printAll(5,fileio_->out);
  if(this->controls_->printLevel>=2) 
    prettyPrint(this->fileio_->out,(*this->oneE_),"Core Hamiltonian");

  // Compute time differenes
  this->OneED = oneEEnd - oneEStart;
  this->SED = OEnd - OStart;
  this->TED = TEnd - TStart;
  this->VED = VEnd - VStart;
/*
  if(this->controls_->printLevel >= 1) {
    this->fileio_->out << endl;
    this->fileio_->out << std::left << std::setw(60) << "CPU time for Overlap evaluation:" 
                       << std::left << std::setw(15) << this->SED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << "CPU time for Kinetic evaluation:" 
                       << std::left << std::setw(15) << this->TED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << "CPU time for Nuclear Attraction Potential evaluation:" 
                       << std::left << std::setw(15) << this->VED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << " "
                       << std::left << std::setw(15) << "---------------" << "----" << endl;
    this->fileio_->out << std::left << std::setw(60) << "Total CPU time for one-electron integral evaluation:" 
                       << std::left << std::setw(15) << this->OneED.count() << " sec" << endl;
  }
*/
  this->haveAOOneE = true;
}

using libint2::TwoBodyEngine;
void AOIntegrals::computeSchwartz(){
  RealMatrix *ShBlk; 
  this->schwartz_->setZero();
  // Check to see if the basisset had been converted
  if(!this->basisSet_->convToLI) this->basisSet_->convShell(this->molecule_);

  // Define Integral Engine
  TwoBodyEngine<libint2::Coulomb> engine = 
    TwoBodyEngine<libint2::Coulomb>(this->basisSet_->maxPrim,
                                    this->basisSet_->maxL,0);
  engine.set_precision(0.); // Don't screen primitives during schwartz

  this->fileio_->out << "Computing Schwartz Bound Tensor ... ";
  auto start =  std::chrono::high_resolution_clock::now();
  for(int s1=0; s1 < this->basisSet_->nShell(); s1++){
    int n1  = this->basisSet_->shells_libint[s1].size();
    for(int s2=0; s2 <= s1; s2++){
      int n2  = this->basisSet_->shells_libint[s2].size();
 
      const auto* buff = engine.compute(
        this->basisSet_->shells_libint[s1],
        this->basisSet_->shells_libint[s2],
        this->basisSet_->shells_libint[s1],
        this->basisSet_->shells_libint[s2]
      );

      
      ShBlk = new RealMatrix(n1,n2);
      ShBlk->setZero();
   
      int ij = 0;
      for(int i = 0; i < n1; i++) {
        for(int j = 0; j < n2; j++) {
          (*ShBlk)(i,j) = buff[ij*n1*n2 + ij];
 	 ij++;
        }
      }

      (*this->schwartz_)(s1,s2) = std::sqrt(ShBlk->lpNorm<Infinity>());
      
      delete ShBlk;
    }
  }
  auto finish =  std::chrono::high_resolution_clock::now();
  this->SchwartzD = finish - start;
  (*this->schwartz_) = this->schwartz_->selfadjointView<Lower>();

  this->fileio_->out << "done (" << this->SchwartzD.count() << ")" << endl;
  this->haveSchwartz = true;
}
void AOIntegrals::computeAOTwoE(){
  if(!this->haveSchwartz) this->computeSchwartz();


  std::vector<coulombEngine> engines(this->controls_->nthreads);
  engines[0] = coulombEngine(this->basisSet_->maxPrim,this->basisSet_->maxL,0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<this->controls_->nthreads; i++) engines[i] = engines[0];
  if(!this->basisSet_->haveMap) this->basisSet_->makeMap(this->molecule_); 

#ifdef USE_OMP
  #pragma omp parallel
#endif
  {
#ifdef USE_OMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
  for(int s1 = 0, s1234=0; s1 < this->basisSet_->nShell(); s1++) {
    int bf1_s = this->basisSet_->mapSh2Bf[s1];
    int n1    = this->basisSet_->shells_libint[s1].size();
    for(int s2 = 0; s2 <= s1; s2++) {
      int bf2_s = this->basisSet_->mapSh2Bf[s2];
      int n2    = this->basisSet_->shells_libint[s2].size();
      for(int s3 = 0; s3 <= s1; s3++) {
        int bf3_s = this->basisSet_->mapSh2Bf[s3];
        int n3    = this->basisSet_->shells_libint[s3].size();
        int s4_max = (s1 == s3) ? s2 : s3;
        for(int s4 = 0; s4 <= s4_max; s4++, s1234++) {

          if(s1234 % this->controls_->nthreads != thread_id) continue;

          int bf4_s = this->basisSet_->mapSh2Bf[s4];
          int n4    = this->basisSet_->shells_libint[s4].size();
    
          // Schwartz and Density screening
          if((*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4)
              < this->controls_->thresholdSchawrtz ) continue;
 
          const double* buff = engines[thread_id].compute(
            this->basisSet_->shells_libint[s1],
            this->basisSet_->shells_libint[s2],
            this->basisSet_->shells_libint[s3],
            this->basisSet_->shells_libint[s4]);


          std::vector<std::array<int,4>> lower;
          std::vector<std::array<int,4>> upper;
/*
          auto lower_p1 = {bf1,bf2,bf3,bf4};
          auto lower_p2 = {bf1,bf2,bf4,bf3};
          auto lower_p3 = {bf2,bf1,bf3,bf4};
          auto lower_p4 = {bf2,bf1,bf4,bf3};
          auto lower_p5 = {bf3,bf4,bf1,bf2};
          auto lower_p6 = {bf4,bf3,bf1,bf2};
          auto lower_p7 = {bf3,bf4,bf2,bf1};
          auto lower_p8 = {bf4,bf3,bf2,bf1};


          auto upper_p1 = {bf1+n1,bf2+n2,bf3+n3,bf4+n4};
          auto upper_p2 = {bf1+n1,bf2+n2,bf4+n4,bf3+n3};
          auto upper_p3 = {bf2+n2,bf1+n1,bf3+n3,bf4+n4};
          auto upper_p4 = {bf2+n2,bf1+n1,bf4+n4,bf3+n3};
          auto upper_p5 = {bf3+n3,bf4+n4,bf1+n1,bf2+n2};
          auto upper_p6 = {bf4+n4,bf3+n3,bf1+n1,bf2+n2};
          auto upper_p7 = {bf3+n3,bf4+n4,bf2+n2,bf1+n1};
          auto upper_p8 = {bf4+n4,bf3+n3,bf2+n2,bf1+n1};
*/
/*        TODO Why doesn't this work?
          lower.push_back({{bf1_s,bf2_s,bf3_s,bf4_s}});
          lower.push_back({{bf1_s,bf2_s,bf4_s,bf3_s}});
          lower.push_back({{bf2_s,bf1_s,bf3_s,bf4_s}});
          lower.push_back({{bf2_s,bf1_s,bf4_s,bf3_s}});
          lower.push_back({{bf3_s,bf4_s,bf1_s,bf2_s}});
          lower.push_back({{bf4_s,bf3_s,bf1_s,bf2_s}});
          lower.push_back({{bf3_s,bf4_s,bf2_s,bf1_s}});
          lower.push_back({{bf4_s,bf3_s,bf2_s,bf1_s}});


          upper.push_back({{bf1_s+n1,bf2_s+n2,bf3_s+n3,bf4_s+n4}});
          upper.push_back({{bf1_s+n1,bf2_s+n2,bf4_s+n4,bf3_s+n3}});
          upper.push_back({{bf2_s+n2,bf1_s+n1,bf3_s+n3,bf4_s+n4}});
          upper.push_back({{bf2_s+n2,bf1_s+n1,bf4_s+n4,bf3_s+n3}});
          upper.push_back({{bf3_s+n3,bf4_s+n4,bf1_s+n1,bf2_s+n2}});
          upper.push_back({{bf4_s+n4,bf3_s+n3,bf1_s+n1,bf2_s+n2}});
          upper.push_back({{bf3_s+n3,bf4_s+n4,bf2_s+n2,bf1_s+n1}});
          upper.push_back({{bf4_s+n4,bf3_s+n3,bf2_s+n2,bf1_s+n1}});

          for(auto p = 0; p < lower.size(); ++p) {
            auto view = btas::make_view(
              this->aoERI_->range().slice(lower[p],upper[p]),
              this->aoERI_->storage());
            std::copy(buff,buff+n1*n2*n3*n4,view.begin());
          }
*/
//        lower.clear();
//        upper.clear();
          for(int i = 0, ijkl = 0 ; i < n1; ++i) {
            int bf1 = bf1_s + i;
            for(int j = 0; j < n2; ++j) {
              int bf2 = bf2_s + j;
              for(int k = 0; k < n3; ++k) {
                int bf3 = bf3_s + k;
                for(int l = 0; l < n4; ++l, ++ijkl) {
                  int bf4 = bf4_s + l;
                  (*this->aoERI_)(bf1,bf2,bf3,bf4) = buff[ijkl];
                  (*this->aoERI_)(bf1,bf2,bf4,bf3) = buff[ijkl];
                  (*this->aoERI_)(bf2,bf1,bf3,bf4) = buff[ijkl];
                  (*this->aoERI_)(bf2,bf1,bf4,bf3) = buff[ijkl];
                  (*this->aoERI_)(bf3,bf4,bf1,bf2) = buff[ijkl];
                  (*this->aoERI_)(bf4,bf3,bf1,bf2) = buff[ijkl];
                  (*this->aoERI_)(bf3,bf4,bf2,bf1) = buff[ijkl];
                  (*this->aoERI_)(bf4,bf3,bf2,bf1) = buff[ijkl];
		}
	      }
	    }
	  }

	}
      }
    }
  }
  } // OMP Parallel

}

void AOIntegrals::computeAORII(){
  if(!this->haveSchwartz) this->computeSchwartz();
  if(!this->haveRIS)      this->computeAORIS(); 


  std::vector<coulombEngine> engines(this->controls_->nthreads);
  engines[0] = coulombEngine(
    std::max(this->basisSet_->maxPrim,this->DFbasisSet_->maxPrim),
    std::max(this->basisSet_->maxL,this->DFbasisSet_->maxL),0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<this->controls_->nthreads; i++) engines[i] = engines[0];
  if(!this->basisSet_->haveMap) this->basisSet_->makeMap(this->molecule_); 
  if(!this->DFbasisSet_->haveMap) this->DFbasisSet_->makeMap(this->molecule_); 

#ifdef USE_OMP
  #pragma omp parallel
#endif
  {
#ifdef USE_OMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
  for(int s1 = 0, s123=0; s1 < this->basisSet_->nShell(); s1++) {
    int bf1_s = this->basisSet_->mapSh2Bf[s1];
    int n1    = this->basisSet_->shells_libint[s1].size();
    for(int s2 = 0; s2 < this->basisSet_->nShell(); s2++) {
      int bf2_s = this->basisSet_->mapSh2Bf[s2];
      int n2    = this->basisSet_->shells_libint[s2].size();
      for(int dfs = 0; dfs < this->DFbasisSet_->nShell(); dfs++,s123++) {
        if(s123 % this->controls_->nthreads != thread_id) continue;
        int dfbf3_s = this->DFbasisSet_->mapSh2Bf[dfs];
        int dfn3    = this->DFbasisSet_->shells_libint[dfs].size();

        // Schwartz and Density screening
        if((*this->schwartz_)(s1,s2) * (*this->aoRIS_)(dfs,dfs)
            < this->controls_->thresholdSchawrtz ) continue;
 
        const double* buff = engines[thread_id].compute(
          this->basisSet_->shells_libint[s1],
          this->basisSet_->shells_libint[s2],
          this->DFbasisSet_->shells_libint[dfs],
          libint2::Shell::unit());

        auto lower = {bf1_s,bf2_s,dfbf3_s};
        auto upper = {bf1_s+n1,bf2_s+n2,dfbf3_s+dfn3};
        auto view  = btas::make_view(
          this->aoRII_->range().slice(lower,upper),
          this->aoRII_->storage());
        std::copy(buff,buff+n1*n2*dfn3,view.begin());
      }
    }
  }
  } // omp parallel scope
}

void AOIntegrals::computeAORIS(){
  this->haveRIS= true;
  std::vector<coulombEngine> engines(this->controls_->nthreads);
  engines[0] = coulombEngine( this->DFbasisSet_->maxPrim, 
                              this->DFbasisSet_->maxL,0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<this->controls_->nthreads; i++) engines[i] = engines[0];
  if(!this->DFbasisSet_->haveMap) this->DFbasisSet_->makeMap(this->molecule_); 

  RealMap aoRISMap(&this->aoRIS_->storage()[0],
    this->DFbasisSet_->nBasis(),this->DFbasisSet_->nBasis());

#ifdef USE_OMP
  #pragma omp parallel
#endif
  {
#ifdef USE_OMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
  for(int s1 = 0, s12=0; s1 < this->DFbasisSet_->nShell(); s1++) {
    int bf1_s = this->DFbasisSet_->mapSh2Bf[s1];
    int n1    = this->DFbasisSet_->shells_libint[s1].size();
    for(int s2 = 0; s2 < this->DFbasisSet_->nShell(); s2++,s12++) {
      int bf2_s = this->DFbasisSet_->mapSh2Bf[s2];
      int n2    = this->DFbasisSet_->shells_libint[s2].size();
 
      if(s12 % this->controls_->nthreads != thread_id) continue;

      const double* buff = engines[thread_id].compute(
        this->DFbasisSet_->shells_libint[s1],
        libint2::Shell::unit(),
        this->DFbasisSet_->shells_libint[s2],
        libint2::Shell::unit());

      ConstRealMap buffMat(buff,n1,n2);
      aoRISMap.block(bf1_s,bf2_s,n1,n2) = buffMat; 
    }
  }
  } // omp parallel scope

  aoRISMap = aoRISMap.selfadjointView<Lower>(); // Symmetrize
}
#endif


