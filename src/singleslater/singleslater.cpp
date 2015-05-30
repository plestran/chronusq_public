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
#include <singleslater.h>
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;
//------------------------------//
// allocate memory for matrices //
//------------------------------//
void SingleSlater::iniSingleSlater(Molecule * molecule, BasisSet * basisset, 
                                   AOIntegrals * aointegrals, FileIO * fileio, 
                                   Controls * controls) {
  int nTotalE = molecule->nTotalE();
  this->nBasis_  = basisset->nBasis();
  this->nTT_   = this->nBasis_*(this->nBasis_+1)/2;
  this->spin_  = molecule->spin();
  int nSingleE = this->spin_ - 1;
  this->nOccB_ = (nTotalE - nSingleE)/2;
  this->nVirB_ = this->nBasis_ - this->nOccB_;
  this->nOccA_ = this->nOccB_ + nSingleE;
  this->nVirA_ = this->nBasis_ - this->nOccA_;
  this->energyNuclei = molecule->energyNuclei();
  if(this->spin_!=1) this->RHF_ = 0;
  else this->RHF_ = 1;

  // FIXME Nedd try statements for allocation
  try { this->densityA_  = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Alpha Density
  catch (...) { CErr(std::current_exception(),"Alpha Density Matrix Allocation"); }
  try { this->fockA_     = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Alpha Fock
  catch (...) { CErr(std::current_exception(),"Alpha Fock Matrix Allocation"); }
#ifndef USE_LIBINT
  try { this->coulombA_  = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Alpha Coulomb Integral
  catch (...) { CErr(std::current_exception(),"Alpha Coulomb Tensor (R2) Allocation"); }
  try { this->exchangeA_ = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); }// Alpha Exchange Integral
  catch (...) { CErr(std::current_exception(),"Alpha Exchange Tensor (R2) Allocation"); }
#else // USE_LIBINT
  try { this->PTA_  = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Alpha Perturbation Tensor
  catch (...) { CErr(std::current_exception(),"Alpha Perturbation Tensor (G[P]) Allocation"); }
#endif
  try { this->moA_       = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Alpha Molecular Orbital Coefficients
  catch (...) { CErr(std::current_exception(),"Alpha MO Coefficients Allocation"); }
  

  if(!this->RHF_) {
    try { this->densityB_  = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Beta Density
    catch (...) { CErr(std::current_exception(),"Beta Density Matrix Allocation"); }
    try { this->fockB_     = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Beta Fock
    catch (...) { CErr(std::current_exception(),"Beta Fock Matrix Allocation"); }
#ifndef USE_LIBINT
    try { this->coulombB_  = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Beta Coulomb Integral
    catch (...) { CErr(std::current_exception(),"Beta Coulomb Tensor (R2) Allocation"); }
    try { this->exchangeB_ = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Beta Exchange Integral
    catch (...) { CErr(std::current_exception(),"Beta Exchange Tensor (R2) Allocation"); }
#else // USE_LIBINT
    try { this->PTB_  = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Beta Perturbation Tensor
    catch (...) { CErr(std::current_exception(),"Beta Perturbation Tensor (G[P]) Allocation"); }
#endif
    try { this->moB_       = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_)); } // Beta Molecular Orbital Coefficients
    catch (...) { CErr(std::current_exception(),"Beta MO Coefficients Allocation"); }
  };

  this->dipole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,1));
  this->quadpole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,3));
  this->tracelessQuadpole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,3));
  this->octpole_  = std::unique_ptr<RealTensor3d>(new RealTensor3d(3,3,3));
  this->molecule_ = molecule;
  this->basisset_ = basisset;
  this->fileio_   = fileio;
  this->controls_ = controls;
  this->aointegrals_= aointegrals;
/* Leaks memory
  int i,j,ij;
  this->R2Index_ = new int*[nBasis];
  for(i=0;i<nBasis;i++) this->R2Index_[i] = new int[nBasis];
  for(i=0;i<nBasis;i++) for(j=0;j<nBasis;j++) {
    if(i>=j) ij=j*(nBasis)-j*(j-1)/2+i-j;
    else ij=i*(nBasis)-i*(i-1)/2+j-i;
    this->R2Index_[i][j] = ij;
  };
*/

  this->haveCoulomb = false;
  this->haveExchange= false;
  this->haveDensity = false;
  this->haveMO	    = false;
  this->havePT = false;
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
/*
  this->energyOneE = (this->aointegrals_->oneE_)->cwiseProduct(*this->densityA_).sum();
#ifndef USE_LIBINT
  this->energyTwoE = (this->coulombA_->cwiseProduct(*this->densityA_).sum() - this->exchangeA_->cwiseProduct(*this->densityA_).sum());
#else
  this->energyTwoE = (this->PTA_)->cwiseProduct(*this->densityA_).sum();
#endif
*/
  this->energyOneE = (*this->aointegrals_->oneE_).frobInner(*this->densityA_);
#ifndef USE_LIBINT
  this->energyTwoE = ((*this->coulombA_)-(*this->exchangeA_)).frobInner(*this->densityA_);
#else
  this->energyTwoE = (*this->PTA_).frobInner(*this->densityA_);
#endif
  this->totalEnergy= this->energyOneE + this->energyTwoE + this->energyNuclei;
  this->printEnergy();
};
//--------------------//
// print energies     //
//--------------------//
void SingleSlater::printEnergy(){
  this->fileio_->out<<"\nEnergy Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"E(one electron) = "<<std::setw(15)<<this->energyOneE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"E(two electron) = "<<std::setw(15)<<this->energyTwoE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"E(nuclear repulsion) = "<<std::setw(15)<<this->energyNuclei<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"E(total) = "<<std::setw(15)<<this->totalEnergy<<std::setw(5)<<" Eh "<<endl;
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
  this->formPT();
#endif
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

  this->fockA_->setZero();
  *(fockA_)+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
  *(fockA_)+=2*(*this->coulombA_);
  *(fockA_)-=2*(*this->exchangeA_);
#else
  *(fockA_)+=2*(*this->PTA_);
#endif
  if(!this->RHF_){
    this->fockB_->setZero();
    *(fockB_)+=(*this->aointegrals_->oneE_);
#ifndef USE_LIBINT
    *(fockB_)+=2*(*this->coulombB_);
    *(fockB_)-=2*(*this->exchangeB_);
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
  this->aointegrals_->CoulD = finish - start; 
  this->fileio_->out<<"\nCPU time for building the Coulomb matrix:  "<< this->aointegrals_->CoulD.count() <<" seconds."<<endl;

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
  this->aointegrals_->ExchD = finish - start; 
  this->fileio_->out<<"\nCPU time for building the Exchange matrix:  "<< this->aointegrals_->ExchD.count() <<" seconds."<<endl;

//if(this->controls_->printLevel>=2) this->exchangeA_->printAll(5,this->fileio_->out);
  if(this->controls_->printLevel>=2) prettyPrint(this->fileio_->out,(*this->exchangeA_),"Alpha Exchange");

  this->haveExchange = true;
};
#endif
//dbwys
#ifdef USE_LIBINT
// Form perturbation tensor (G)
void SingleSlater::formPT() {
  if(!this->haveDensity) this->formDensity();
  this->aointegrals_->twoEContract(true,*this->densityA_,*this->PTA_);
  if(this->controls_->printLevel >= 3) prettyPrint(this->fileio_->out,(*this->PTA_),"Alpha Perturbation Tensor");
}
/*
void SingleSlater::formPT() {
  if(!this->aointegrals_->haveSchwartz) this->aointegrals_->computeSchwartz();
  if(!this->haveDensity) this->formDensity();
  this->PTA_->setZero();
  std::vector<RealMatrix> 
    G(this->controls_->nthreads,RealMatrix::Zero(this->nBasis_,this->nBasis_));

  std::vector<coulombEngine> engines(this->controls_->nthreads);
  engines[0] = coulombEngine(this->basisset_->maxPrim,this->basisset_->maxL,0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<this->controls_->nthreads; i++) engines[i] = engines[0];

  if(!this->basisset_->haveMap) this->basisset_->makeMap(this->molecule_); 
  auto start = std::chrono::high_resolution_clock::now();
  this->basisset_->computeShBlkNorm(this->molecule_,this->densityA_.get());
  auto finish = std::chrono::high_resolution_clock::now();
  this->aointegrals_->DenShBlkD = finish - start;
  int ijkl = 0;
  start = std::chrono::high_resolution_clock::now();

  auto lambda = [&] (int thread_id) {
    coulombEngine &engine = engines[thread_id];
    RealMatrix &g = G[thread_id];
    for(int s1 = 0, s1234=0; s1 < this->basisset_->nShell(); s1++) {
      int bf1_s = this->basisset_->mapSh2Bf[s1];
      int n1    = this->basisset_->shells_libint[s1].size();
      for(int s2 = 0; s2 <= s1; s2++) {
        int bf2_s = this->basisset_->mapSh2Bf[s2];
        int n2    = this->basisset_->shells_libint[s2].size();
        for(int s3 = 0; s3 <= s1; s3++) {
          int bf3_s = this->basisset_->mapSh2Bf[s3];
          int n3    = this->basisset_->shells_libint[s3].size();
          int s4_max = (s1 == s3) ? s2 : s3;
          for(int s4 = 0; s4 <= s4_max; s4++, s1234++) {
            if(s1234 % this->controls_->nthreads != thread_id) continue;
            int bf4_s = this->basisset_->mapSh2Bf[s4];
            int n4    = this->basisset_->shells_libint[s4].size();
      
            // Schwartz and Density screening
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
                    g(bf1,bf2) += (*this->densityA_)(bf3,bf4)*v;
                    g(bf3,bf4) += (*this->densityA_)(bf1,bf2)*v;
 
                    // Exchange
                    g(bf1,bf3) -= 0.25*(*this->densityA_)(bf2,bf4)*v;
                    g(bf2,bf4) -= 0.25*(*this->densityA_)(bf1,bf3)*v;
                    g(bf1,bf4) -= 0.25*(*this->densityA_)(bf2,bf3)*v;
                    g(bf2,bf3) -= 0.25*(*this->densityA_)(bf1,bf4)*v;
                  }
                }
              }
            }
          }
        }
      }
    }
  };

#ifdef USE_OMP
  #pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#else
  lambda(0);
#endif
  for(int i = 0; i < this->controls_->nthreads; i++) (*this->PTA_) += G[i];
//exit(EXIT_FAILURE);
  RealMatrix Tmp = 0.5*((*this->PTA_) + (*this->PTA_).transpose());
  (*this->PTA_) = 0.25*Tmp; // Can't consolidate where this comes from?
  finish = std::chrono::high_resolution_clock::now();
  this->aointegrals_->PTD = finish - start;
//this->PTA_->printAll(5,this->fileio_->out);
  if(this->controls_->printLevel >= 3) prettyPrint(this->fileio_->out,(*this->PTA_),"Alpha Perturbation Tensor");
  
}
*/
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
  std::string readString;
  this->fileio_->in.seekg(0,ios::beg);
  this->fileio_->in>>readString;
  readString=stringupper(readString);
  while(!this->fileio_->in.eof()&&readString.compare("$GUESS")) {
    this->fileio_->in>>readString;
    readString=stringupper(readString);
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
// form the initial guess of MOs from Gaussian raw matrix element file   //
//-----------------------------------------------------------------------//
void SingleSlater::readGuessGauMatEl(GauMatEl& matEl){
  this->fileio_->out << "Reading MO coefficients from " <<matEl.fname()<< endl;
  if(matEl.nBasis()!=this->nBasis_) CErr("Basis Set mismatch",this->fileio_->out);
  double *scr = matEl.readRec(GauMatEl::moa); 
  for(auto i = 0; i < this->nBasis_*this->nBasis_; i++)
    this->moA_->data()[i] = scr[i];
  this->moA_->transposeInPlace();
  if(this->controls_->printLevel>=3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  delete [] scr;
  this->haveMO = true;
}
//-----------------------------------------------------------------------//
// form the initial guess of MOs from Gaussian formatted checkpoint file //
//-----------------------------------------------------------------------//
void SingleSlater::readGuessGauFChk(std::string &filename) {
  this->fileio_->out<<"Reading formatted checkpoint file "<<filename<<endl;
  std::string readString;
  int i,j,nBasis;
  double data;
  std::unique_ptr<ifstream> fchk(new ifstream(filename));
  if(fchk->fail()) 
    CErr("Could not find "+filename,this->fileio_->out);

  *fchk>>readString;
  while((!(fchk->eof()))&&(readString.compare("basis"))) *fchk>>readString;
  if(!readString.compare("basis")) {
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

  while(!(fchk->eof())&&readString.compare("MO")) *fchk>>readString;
  if(!readString.compare("MO")) {
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

void SingleSlater::computeMultipole(){
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  if(!this->controls_->doDipole && !this->controls_->doQuadpole) return;

  int NB = this->nBasis_;
  int NBSq = NB*NB;
  int iBuf = 0;
  for(auto ixyz = 0; ixyz < 3; ixyz++){
    ConstRealMap mu(&this->aointegrals_->elecDipole_->storage()[iBuf],NB,NB);
    (*dipole_)(ixyz,0) = -this->densityA_->frobInner(mu);
    iBuf += NBSq;
  }
  for(int iA = 0; iA < this->molecule_->nAtoms(); iA++)
    *this->dipole_ += elements[this->molecule_->index(iA)].atomicNumber *
          this->molecule_->cart()->col(iA);
  
  if(this->controls_->doQuadpole){
    iBuf = 0;
    for(auto jxyz = 0; jxyz < 3; jxyz++)
    for(auto ixyz = jxyz; ixyz < 3; ixyz++){
      ConstRealMap 
        mu(&this->aointegrals_->elecQuadpole_->storage()[iBuf],NB,NB);
      (*quadpole_)(jxyz,ixyz) = -this->densityA_->frobInner(mu);
      iBuf += NBSq;
    }
    *this->quadpole_ = this->quadpole_->selfadjointView<Upper>();
//  this->quadpole_->setZero();
    for(int iA = 0; iA < this->molecule_->nAtoms(); iA++)
      *this->quadpole_ += elements[this->molecule_->index(iA)].atomicNumber *
            this->molecule_->cart()->col(iA) * 
            this->molecule_->cart()->col(iA).transpose();

    (*this->tracelessQuadpole_) = (*this->quadpole_) - RealMatrix::Identity(3,3)*this->quadpole_->trace()/3;
  }

  if(this->controls_->doOctpole){
    iBuf = 0;
    for(auto kxyz = 0;    kxyz < 3; kxyz++)
    for(auto jxyz = kxyz; jxyz < 3; jxyz++)
    for(auto ixyz = jxyz; ixyz < 3; ixyz++){
      ConstRealMap 
        mu(&this->aointegrals_->elecOctpole_->storage()[iBuf],NB,NB);
      (*octpole_)(kxyz,jxyz,ixyz) = -this->densityA_->frobInner(mu);
      iBuf += NBSq;
    }
    for(auto kxyz = 0;    kxyz < 3; kxyz++)
    for(auto jxyz = kxyz; jxyz < 3; jxyz++)
    for(auto ixyz = jxyz; ixyz < 3; ixyz++){
      (*octpole_)(ixyz,jxyz,kxyz) = (*octpole_)(kxyz,jxyz,ixyz);
      (*octpole_)(ixyz,kxyz,jxyz) = (*octpole_)(kxyz,jxyz,ixyz);
      (*octpole_)(jxyz,ixyz,kxyz) = (*octpole_)(kxyz,jxyz,ixyz);
      (*octpole_)(jxyz,kxyz,ixyz) = (*octpole_)(kxyz,jxyz,ixyz);
      (*octpole_)(kxyz,ixyz,jxyz) = (*octpole_)(kxyz,jxyz,ixyz);
    }

    for(auto iA = 0; iA < this->molecule_->nAtoms(); iA++)
    for(auto ixyz = 0; ixyz < 3; ixyz++)
    for(auto jxyz = 0; jxyz < 3; jxyz++)
    for(auto kxyz = 0; kxyz < 3; kxyz++)
      (*this->octpole_)(ixyz,jxyz,kxyz) += elements[this->molecule_->index(iA)].atomicNumber *
            (*this->molecule_->cart())(ixyz,iA)*
            (*this->molecule_->cart())(jxyz,iA)*
            (*this->molecule_->cart())(kxyz,iA);
    
  }
  this->printMultipole();

}
void SingleSlater::printMultipole(){
  this->fileio_->out << bannerTop << endl;
  this->fileio_->out << std::setw(50) << std::left <<"Electric Dipole Moment"
                        << "(Debye)" << endl;
  this->fileio_->out << std::left << std::setw(5) <<"X=" 
                     << std::fixed << std::right << std::setw(20) 
                     << (*this->dipole_)(0,0)/phys.debye << endl;
  this->fileio_->out << std::left << std::setw(5) <<"Y=" 
                     << std::fixed << std::right << std::setw(20) 
                     << (*this->dipole_)(1,0)/phys.debye << endl;
  this->fileio_->out << std::left << std::setw(5) <<"Z=" 
                     << std::fixed << std::right << std::setw(20) 
                     << (*this->dipole_)(2,0)/phys.debye << endl;
  if(this->controls_->doQuadpole) {
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left << "Electric Quadrupole Moment" 
                       <<  "(Debye-\u212B)" << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(0,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(0,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(0,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(1,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(1,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(1,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(2,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(2,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(2,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left << "Electric Quadrupole Moment (Traceless)" 
                       <<  "(Debye-\u212B)" << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(0,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(0,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(0,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(1,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(1,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(1,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(2,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(2,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(2,2)*phys.bohr/phys.debye << endl;
  }
  if(this->controls_->doOctpole) {
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left << "Electric Octupole Moment" 
                       << "(Debye-\u212B\u00B2)" << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,0,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,0,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,0,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,1,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,1,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,1,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,2,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,2,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,2,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,0,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,0,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,0,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,1,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,1,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,1,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,2,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,2,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,2,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,0,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,0,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,0,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,1,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,1,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,1,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,2,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,2,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,2,2)*phys.bohr*phys.bohr/phys.debye << endl;
  }
  this->fileio_->out << endl << bannerEnd << endl;
}
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

