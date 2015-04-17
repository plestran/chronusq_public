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
#include "singleslater.h"
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;
//------------------------------//
// allocate memory for matrices //
//------------------------------//
void SingleSlater::iniSingleSlater(Molecule *molecule, BasisSet *basisset, AOIntegrals *aointegrals, FileIO *fileio, Controls *controls) {
  int nBasis  = basisset->nBasis();
  int nTotalE = molecule->nTotalE();
  int spin = molecule->spin();
  if(spin!=1) this->RHF_ = 0;
  else this->RHF_ = 1;
  try {
    this->densityA_  = new RealMatrix(nBasis,nBasis); // Alpha Density
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::densityA_! E#:"<<msg<<endl;
    exit(1);
  };
  try{
    this->fockA_     = new RealMatrix(nBasis,nBasis); // Alpha Fock
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::fockA_! E#:"<<msg<<endl;
    exit(1);
  };
#ifndef USE_LIBINT
  try{
    this->coulombA_  = new RealMatrix(nBasis,nBasis); // Alpha Coulomb Integral
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::coulombA_! E#:"<<msg<<endl;
    exit(1);
  };
  try{
    this->exchangeA_ = new RealMatrix(nBasis,nBasis); // Alpha Exchange Integral
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::exchangeA_! E#:"<<msg<<endl;
    exit(1);
  };
#else // USE_LIBINT
  try{
    this->PTA_  = new RealMatrix(nBasis,nBasis); // Alpha Perturbation Tensor
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::PTA_! E#:"<<msg<<endl;
    exit(1);
  };
#endif
  try{
    this->moA_       = new RealMatrix(nBasis,nBasis); // Alpha Molecular Orbital Coefficients
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::moA_! E#:"<<msg<<endl;
    exit(1);
  };

  if(!this->RHF_) {
    try{
      this->densityB_  = new RealMatrix(nBasis,nBasis); // Beta Density
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::densityB_! E#:"<<msg<<endl;
      exit(1);
    };
    try{
      this->fockB_     = new RealMatrix(nBasis,nBasis); // Beta Fock
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::fockB_! E#:"<<msg<<endl;
      exit(1);
    };
#ifndef USE_LIBINT
    try{
      this->coulombB_  = new RealMatrix(nBasis,nBasis); // Beta Coulomb Integral
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::coulombB_! E#:"<<msg<<endl;
      exit(1);
    };
    try{
      this->exchangeB_ = new RealMatrix(nBasis,nBasis); // Beta Exchange Integral
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::exchangeB_! E#:"<<msg<<endl;
      exit(1);
    };
#else // USE_LIBINT
  try{
    this->PTB_  = new RealMatrix(nBasis,nBasis); // Beta Perturbation Tensor
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::PTB_! E#:"<<msg<<endl;
    exit(1);
  };
#endif
    try{
      this->moB_       = new RealMatrix(nBasis,nBasis); // Beta Molecular Orbital Coefficients
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::moB_! E#:"<<msg<<endl;
      exit(1);
    };
  };

  this->nBasis_= nBasis;
  this->nTT_   = this->nBasis_*(this->nBasis_+1)/2;
  this->spin_  = spin;
  int nSingleE = spin - 1;
  this->nOccB_ = (nTotalE - nSingleE)/2;
  this->nVirB_ = this->nBasis_ - this->nOccB_;
  this->nOccA_ = this->nOccB_ + nSingleE;
  this->nVirA_ = this->nBasis_ - this->nOccA_;
  this->energyNuclei = molecule->energyNuclei();

  this->molecule_ = molecule;
  this->basisset_ = basisset;
  this->fileio_   = fileio;
  this->controls_ = controls;
  this->aointegrals_= aointegrals;

  int i,j,ij;
  this->R2Index_ = new int*[nBasis];
  for(i=0;i<nBasis;i++) this->R2Index_[i] = new int[nBasis];
  for(i=0;i<nBasis;i++) for(j=0;j<nBasis;j++) {
    if(i>=j) ij=j*(nBasis)-j*(j-1)/2+i-j;
    else ij=i*(nBasis)-i*(i-1)/2+j-i;
    this->R2Index_[i][j] = ij;
  };

  this->haveCoulomb = false;
  this->haveExchange= false;
  this->haveDensity = false;
  this->haveMO	    = false;
#ifdef USE_LIBINT
  this->havePT = false;
#endif
};
//-----------------------------------//
// print a wave function information //
//-----------------------------------//
void SingleSlater::printInfo() {
  this->fileio_->out<<"\nSingle Slater Determinent Wave Function Information:"<<endl;
  this->fileio_->out<<std::setw(15)<<"nOccA ="<<std::setw(8)<<this->nOccA_<<std::setw(5)<<" "<<std::setw(20)<<"nVirA ="<<std::setw(8)<<this->nVirA_<<endl;
  this->fileio_->out<<std::setw(15)<<"nOccB ="<<std::setw(8)<<this->nOccB_<<std::setw(5)<<" "<<std::setw(20)<<"nVirB ="<<std::setw(8)<<this->nVirB_<<endl;
  this->fileio_->out<<std::setw(15)<<"Multiplicity ="<<std::setw(8)<<this->spin_<<endl;
};
//----------------------//
// compute energies     //
//----------------------//
void SingleSlater::computeEnergy(){
  this->energyOneE = (this->aointegrals_->oneE_)->cwiseProduct(*this->densityA_).sum();
#ifndef USE_LIBINT
  this->energyTwoE = (this->coulombA_->cwiseProduct(*this->densityA_).sum() - this->exchangeA_->cwiseProduct(*this->densityA_).sum());
#else
  this->energyTwoE = (this->PTA_)->cwiseProduct(*this->densityA_).sum();
#endif
  this->totalEnergy= this->energyOneE + this->energyTwoE + this->energyNuclei;
  this->printEnergy();
};
//--------------------//
// print energies     //
//--------------------//
void SingleSlater::printEnergy(){
  this->fileio_->out<<"\nEnergy Information:"<<endl;
  this->fileio_->out<<std::setw(30)<<"E(one electron) = "<<std::setw(15)<<this->energyOneE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::setw(30)<<"E(two electron) = "<<std::setw(15)<<this->energyTwoE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::setw(30)<<"E(nuclear repulsion) = "<<std::setw(15)<<this->energyNuclei<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::setw(30)<<"E(total) = "<<std::setw(15)<<this->totalEnergy<<std::setw(5)<<" Eh "<<endl;
};
//-------------------------//
// form the density matrix //
//-------------------------//
void SingleSlater::formDensity(){
  if(!this->haveMO) {
    this->fileio_->out<<"No MO coefficients available to form one-particle density matrix!"<<endl;
    exit(1);
  };

  this->densityA_->setZero();
  int i,j,k,nE=0;
  for(i=0;i<this->nBasis_;i++){
    for(j=i;j<this->nBasis_;j++){
      for(k=0;k<this->nOccA_;k++) (*(this->densityA_))(i,j)+=(*(this->moA_))(i,k)*(*(this->moA_))(j,k);
      if(this->RHF_) (*(this->densityA_))(i,j) *= math.two;
    };
  };
  (*this->densityA_) = this->densityA_->selfadjointView<Upper>();

  if(!this->RHF_) {
    this->densityB_->setZero();
    for(i=0;i<this->nBasis_;i++) for(j=i;j<this->nBasis_;j++)
      for(k=0;k<this->nOccB_;k++) (*(this->densityB_))(i,j)+=(*(this->moB_))(i,k)*(*(this->moB_))(j,k);
  };

  if(this->controls_->printLevel>=2) {
/*
    this->densityA_->printAll(5,this->fileio_->out);
    if(!this->RHF_) this->densityB_->printAll(5,this->fileio_->out);
*/
    prettyPrint(this->fileio_->out,(*this->densityA_),"Alpha Density");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->densityB_),"Beta Density");
  };
  this->haveDensity = true;
};
//-------------------------//
// form the Fock matrix    //
//-------------------------//
void SingleSlater::formFock(){
  if(!this->haveDensity) this->formDensity();
#ifndef USE_LIBINT
  if(!this->haveCoulomb) this->formCoulomb();
  if(!this->haveExchange) this->formExchange();
#else
  if(!this->havePT) this->formPT();
#endif
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

  this->fockA_->setZero();
  *(fockA_)+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
  *(fockA_)+=(*this->coulombA_);
  *(fockA_)-=(*this->exchangeA_);
#else
  *(fockA_)+=(*this->PTA_);
#endif
  if(!this->RHF_){
    this->fockB_->setZero();
    *(fockB_)+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
    *(fockB_)+=(*this->coulombB_);
    *(fockB_)-=(*this->exchangeB_);
#else
    *(fockB_)+=(*this->PTB_);
#endif
  };

#ifndef USE_LIBINT
  RealMatrix *tmp = new RealMatrix(this->nBasis_,this->nBasis_); // NO
  tmp->setZero();
  (*tmp) = (*this->coulombA_) - (*this->exchangeA_);
  prettyPrint(this->fileio_->out,(*tmp),"Coul - Exch");
#endif
  if(this->controls_->printLevel>=2) {
/*
    this->fockA_->printAll(5,this->fileio_->out);
    if(!this->RHF_) this->fockB_->printAll(5,this->fileio_->out);
*/
    prettyPrint(this->fileio_->out,(*this->fockA_),"Alpha Fock");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->fockB_),"Beta Fock");
  };
};
#ifndef USE_LIBINT
//----------------------------//
// form the Coulomb matrix    //
//----------------------------//
void SingleSlater::formCoulomb(){
//clock_t start,finish;
//start = clock();
  std::chrono::high_resolution_clock::time_point start, finish;
  start = std::chrono::high_resolution_clock::now();
  int i;
  
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOTwoE) this->aointegrals_->computeAOTwoE();
  this->coulombA_->setZero();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.half;

  this->densityA_->vectorize();
  this->coulombA_->vectorize();
  (*this->coulombA_) = (*this->aointegrals_->twoEC_)*(*this->densityA_);
  this->coulombA_->unvectorize();
  this->densityA_->unvectorize();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.two;

  finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> CoulD = finish - start; 
  this->fileio_->out<<"\nCPU time for building the Coulomb matrix:  "<< CoulD.count() <<" seconds."<<endl;

//if(this->controls_->printLevel>=2) this->coulombA_->printAll(5,this->fileio_->out);
  if(this->controls_->printLevel>=2) prettyPrint(this->fileio_->out,(*this->coulombA_),"Alpha Coulomb");
  this->haveCoulomb = true;
};
//----------------------------//
// form the exchange matrix    //
//----------------------------//
void SingleSlater::formExchange(){
//clock_t start,finish;
//start = clock();
  std::chrono::high_resolution_clock::time_point start, finish;
  start = std::chrono::high_resolution_clock::now();
  int i;
  
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOTwoE) this->aointegrals_->computeAOTwoE();
  this->exchangeA_->setZero();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.half;

  this->densityA_->vectorize();
  this->exchangeA_->vectorize();
  (*this->exchangeA_)= (*this->aointegrals_->twoEX_)*(*this->densityA_);
  this->exchangeA_->scale(math.quarter);
  this->exchangeA_->unvectorize();
  this->densityA_->unvectorize();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.two;

  finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> ExchD = finish - start; 
  this->fileio_->out<<"\nCPU time for building the Exchange matrix:  "<< ExchD.count() <<" seconds."<<endl;

//if(this->controls_->printLevel>=2) this->exchangeA_->printAll(5,this->fileio_->out);
  if(this->controls_->printLevel>=2) prettyPrint(this->fileio_->out,(*this->exchangeA_),"Alpha Exchange");

  this->haveExchange = true;
};
#endif
//dbwys
#ifdef USE_LIBINT
using libint2::TwoBodyEngine;
typedef TwoBodyEngine<libint2::Coulomb> coulombEngine;
// Form perturbation tensor (G)
void SingleSlater::formPT() {
  if(!this->aointegrals_->haveSchwartz) this->aointegrals_->computeSchwartz();
  if(!this->haveDensity) this->formDensity();
  this->PTA_->setZero();

  coulombEngine engine = coulombEngine(this->basisset_->maxPrim,
                                       this->basisset_->maxL,0);
  engine.set_precision(std::numeric_limits<double>::epsilon());
  this->fileio_->out << "Computing Two Electron Integrals with " <<
    std::scientific << engine.precision() << " precision" << endl;

  if(!this->basisset_->haveMap) this->basisset_->makeMap(this->molecule_); 
  this->basisset_->computeShBlkNorm(this->molecule_,this->densityA_);
  int ijkl = 0;
  for(int s1 = 0; s1 < this->basisset_->nShell(); s1++) {
    int bf1_s = this->basisset_->mapSh2Bf[s1];
    int n1    = this->basisset_->shells_libint[s1].size();
    for(int s2 = 0; s2 <= s1; s2++) {
      int bf2_s = this->basisset_->mapSh2Bf[s2];
      int n2    = this->basisset_->shells_libint[s2].size();
      for(int s3 = 0; s3 <= s1; s3++) {
        int bf3_s = this->basisset_->mapSh2Bf[s3];
        int n3    = this->basisset_->shells_libint[s3].size();
        int s4_max = (s1 == s3) ? s2 : s3;
        for(int s4 = 0; s4 <= s4_max; s4++) {
          int bf4_s = this->basisset_->mapSh2Bf[s4];
          int n4    = this->basisset_->shells_libint[s4].size();
    
          if( std::max((*this->basisset_->shBlkNorm)(s1,s4),
                 std::max((*this->basisset_->shBlkNorm)(s2,s4),
                    std::max((*this->basisset_->shBlkNorm)(s3,s4),
                       std::max((*this->basisset_->shBlkNorm)(s1,s3),
                          std::max((*this->basisset_->shBlkNorm)(s2,s3),
                                   (*this->basisset_->shBlkNorm)(s1,s2))
                          )
                       )      
                    )
                 ) * (*this->aointegrals_->schwartz_)(s1,s2)
                   * (*this->aointegrals_->schwartz_)(s3,s4)
                 < this->controls_->thresholdSchawrtz ) continue;

          const double* buff = engine.compute(
            this->basisset_->shells_libint[s1],
            this->basisset_->shells_libint[s2],
            this->basisset_->shells_libint[s3],
            this->basisset_->shells_libint[s4]);
       // cout << "SHELL :" << s1 << " " << s2 << " " << s3 << " " <<s4<<endl; 
       // for(int i = 0; i < n1*n2*n3*n4; i++) cout << buff[i] << endl;
    
          double s12_deg = (s1 == s2) ? 1.0 : 2.0;
          double s34_deg = (s3 == s4) ? 1.0 : 2.0;
          double s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
          double s1234_deg = s12_deg * s34_deg * s12_34_deg;
          for(int i = 0, ijkl = 0 ; i < n1; ++i) {
            int bf1 = bf1_s + i;
            for(int j = 0; j < n2; ++j) {
              int bf2 = bf2_s + j;
              for(int k = 0; k < n3; ++k) {
                int bf3 = bf3_s + k;
                for(int l = 0; l < n4; ++l, ++ijkl) {
                  int bf4 = bf4_s + l;
                  double v = buff[ijkl]*s1234_deg;

                  // Coulomb
                  (*this->PTA_)(bf1,bf2) += (*this->densityA_)(bf3,bf4)*v;
                  (*this->PTA_)(bf3,bf4) += (*this->densityA_)(bf1,bf2)*v;

                  // Exchange
                  (*this->PTA_)(bf1,bf3) -= 0.25*(*this->densityA_)(bf2,bf4)*v;
                  (*this->PTA_)(bf2,bf4) -= 0.25*(*this->densityA_)(bf1,bf3)*v;
                  (*this->PTA_)(bf1,bf4) -= 0.25*(*this->densityA_)(bf2,bf3)*v;
                  (*this->PTA_)(bf2,bf3) -= 0.25*(*this->densityA_)(bf1,bf4)*v;
                }
              }
            }
          }
        }
      }
    }
  }
//exit(EXIT_FAILURE);
  RealMatrix Tmp = 0.5*((*this->PTA_) + (*this->PTA_).transpose());
  (*this->PTA_) = 0.25*Tmp; // Can't consolidate where this comes from?
//this->PTA_->printAll(5,this->fileio_->out);
  prettyPrint(this->fileio_->out,(*this->PTA_),"Alpha Perturbation Tensor");
  
}
#endif
//dbwye
//--------------------------------//
// form the initial guess of MO's //
//--------------------------------//
void SingleSlater::formGuess() {
  if(this->controls_->printLevel>=3) {
/*
    this->moA_->printAll(5,this->fileio_->out);
    if(!this->RHF_) this->moB_->printAll(5,this->fileio_->out);
*/
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
};
//------------------------------------------//
// form the initial guess of MOs from input //
//------------------------------------------//
void SingleSlater::readGuessIO() {
  int i,j;
  this->fileio_->in.clear();
  char readString[MAXNAMELEN];
  this->fileio_->in.seekg(0,ios::beg);
  this->fileio_->in>>readString;
  strupr(readString);
  while(!this->fileio_->in.eof()&&strcmp(readString,"$GUESS")) {
    this->fileio_->in>>readString;
    strupr(readString);
  };
  this->fileio_->in >> readString;
  for(i=0;i<this->nBasis_;i++) for(j=0;j<this->nBasis_;j++){
    this->fileio_->in>>(*(this->moA_))(j,i);
  };
  this->fileio_->in.close();
  if(this->controls_->printLevel>=3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
};
//-----------------------------------------------------------------------//
// form the initial guess of MOs from Gaussian formatted checkpoint file //
//-----------------------------------------------------------------------//
void SingleSlater::readGuessGauFChk(char *filename) {
  this->fileio_->out<<"reading formatted checkpoint file "<<filename<<endl;
  char readString[MAXNAMELEN];
  int i,j,nBasis;
  double data;
  ifstream *fchk = new ifstream(filename);

  *fchk>>readString;
  while((!(fchk->eof()))&&(strcmp(readString,"basis"))) *fchk>>readString;
  if(!strcmp(readString,"basis")) {
    *fchk>>readString;
    *fchk>>readString;
    *fchk>>nBasis;
    if(nBasis!=this->nBasis_) {
      this->fileio_->out<<"Basis set mismatch! "<<nBasis<<"!="<<this->nBasis_<<endl;
      exit(1);
    };
  } else {
    this->fileio_->out<<"No basis set found in the formatted checkpoint file! "<<endl;
  };

  while(!(fchk->eof())&&strcmp(readString,"MO")) *fchk>>readString;
  if(!strcmp(readString,"MO")) {
    *fchk>>readString;
    *fchk>>readString;
    *fchk>>readString;
    *fchk>>readString;
    for(i=0;i<nBasis;i++) 
      for(j=0;j<nBasis;j++){
        *fchk>>data;
        (*(this->moA_))(j,i) = data;
    };
  } else {
    this->fileio_->out<<"No basis set found in the formatted checkpoint file! "<<endl;
  };

  fchk->close();
  if(this->controls_->printLevel>=3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
};

/*************************/
/* MPI Related Routines  */
/*************************/
void SingleSlater::mpiSend(int toID,int tag) {
  //OOMPI_COMM_WORLD[toID].Send(this->nAtoms_,tag);
  //OOMPI_COMM_WORLD[toID].Send(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiSend(toID,tag);
};
void SingleSlater::mpiRecv(int fromID,int tag) {
  //OOMPI_COMM_WORLD[fromID].Recv(this->nAtoms_,tag);
  //this->index_=new int[this->nAtoms_];
  //this->cart_ =new RealMatrix(3, this->nAtoms_, "Molecule");
  //OOMPI_COMM_WORLD[fromID].Recv(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiRecv(fromID,tag);
};

