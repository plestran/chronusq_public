#include "singleslater.h"
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
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
    this->densityA_  = new Matrix<double>(nBasis,nBasis,"Alpha Density","LT");
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::densityA_! E#:"<<msg<<endl;
    exit(1);
  };
  try{
    this->fockA_     = new Matrix<double>(nBasis,nBasis,"Alpha Fock","LT");
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::fockA_! E#:"<<msg<<endl;
    exit(1);
  };
  try{
    this->coulombA_  = new Matrix<double>(nBasis,nBasis,"Alpha Coulomb Integral","LT");
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::coulombA_! E#:"<<msg<<endl;
    exit(1);
  };
  try{
    this->exchangeA_ = new Matrix<double>(nBasis,nBasis,"Alpha Exchange Integral","LT");
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::exchangeA_! E#:"<<msg<<endl;
    exit(1);
  };
  try{
    this->moA_       = new Matrix<double>(nBasis,nBasis,"Alpha Molecular Orbital Coefficients");
  } catch (int msg) {
    fileio->out<<"Unable to allocate memory for SingleSlater::moA_! E#:"<<msg<<endl;
    exit(1);
  };

  if(!this->RHF_) {
    try{
      this->densityB_  = new Matrix<double>(nBasis,nBasis,"Beta Density","LT");
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::densityB_! E#:"<<msg<<endl;
      exit(1);
    };
    try{
      this->fockB_     = new Matrix<double>(nBasis,nBasis,"Beta Fock","LT");
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::fockB_! E#:"<<msg<<endl;
      exit(1);
    };
    try{
      this->coulombB_  = new Matrix<double>(nBasis,nBasis,"Beta Coulomb Integral","LT");
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::coulombB_! E#:"<<msg<<endl;
      exit(1);
    };
    try{
      this->exchangeB_ = new Matrix<double>(nBasis,nBasis,"Beta Exchange Integral","LT");
    } catch (int msg) {
      fileio->out<<"Unable to allocate memory for SingleSlater::exchangeB_! E#:"<<msg<<endl;
      exit(1);
    };
    try{
      this->moB_       = new Matrix<double>(nBasis,nBasis,"Alpha Molecular Orbital Coefficients","STD");
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
  this->energyOneE = (this->aointegrals_->oneE_)->scalarProd(this->densityA_);
  this->energyTwoE = (this->coulombA_->scalarProd(this->densityA_) - this->exchangeA_->scalarProd(this->densityA_));
  this->totalEnergy= this->energyOneE + this->energyTwoE + this->energyNuclei;
  this->printEnergy();
};
//--------------------//
// print energies     //
//--------------------//
void SingleSlater::printEnergy(){
  this->fileio_->out<<"\nEnergy Information:"<<endl;
  this->fileio_->out<<std::setw(30)<<"E(one electron) = "<<std::setw(15)<<this->energyOneE<<std::setw(5)<<" a.u. "<<endl;
  this->fileio_->out<<std::setw(30)<<"E(two electron) = "<<std::setw(15)<<this->energyTwoE<<std::setw(5)<<" a.u. "<<endl;
  this->fileio_->out<<std::setw(30)<<"E(nuclear repulsion) = "<<std::setw(15)<<this->energyNuclei<<std::setw(5)<<" a.u. "<<endl;
  this->fileio_->out<<std::setw(30)<<"E(total) = "<<std::setw(15)<<this->totalEnergy<<std::setw(5)<<" a.u. "<<endl;
};
//-------------------------//
// form the density matrix //
//-------------------------//
void SingleSlater::formDensity(){
  if(!this->haveMO) {
    this->fileio_->out<<"No MO coefficients available to form one-particle density matrix!"<<endl;
    exit(1);
  };

  this->densityA_->clearAll();
  int i,j,k,nE=0;
  for(i=0;i<this->nBasis_;i++){
    for(j=i;j<this->nBasis_;j++){
      for(k=0;k<this->nOccA_;k++) (*(this->densityA_))(i,j)+=(*(this->moA_))(i,k)*(*(this->moA_))(j,k);
      if(this->RHF_) (*(this->densityA_))(i,j) *= math.two;
    };
  };

  if(!this->RHF_) {
    this->densityB_->clearAll();
    for(i=0;i<this->nBasis_;i++) for(j=i;j<this->nBasis_;j++)
      for(k=0;k<this->nOccB_;k++) (*(this->densityB_))(i,j)+=(*(this->moB_))(i,k)*(*(this->moB_))(j,k);
  };

  if(this->controls_->printLevel>=2) {
    this->densityA_->printAll(5,this->fileio_->out);
    if(!this->RHF_) this->densityB_->printAll(5,this->fileio_->out);
  };
  this->haveDensity = true;
};
//-------------------------//
// form the Fock matrix    //
//-------------------------//
void SingleSlater::formFock(){
  if(!this->haveDensity) this->formDensity();
  if(!this->haveCoulomb) this->formCoulomb();
  if(!this->haveExchange) this->formExchange();
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

  this->fockA_->clearAll();
  *(fockA_)+=this->aointegrals_->oneE_;
  *(fockA_)+=this->coulombA_;
  *(fockA_)-=this->exchangeA_;
  if(!this->RHF_){
    this->fockB_->clearAll();
    *(fockB_)+=this->aointegrals_->oneE_;
    *(fockB_)+=this->coulombB_;
    *(fockB_)-=this->exchangeB_;
  };
  if(this->controls_->printLevel>=2) {
    this->fockA_->printAll(5,this->fileio_->out);
    if(!this->RHF_) this->fockB_->printAll(5,this->fileio_->out);
  };
};
//----------------------------//
// form the Coulomb matrix    //
//----------------------------//
void SingleSlater::formCoulomb(){
  clock_t start,finish;
  start = clock();
  int i;
  
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOTwoE) this->aointegrals_->computeAOTwoE();
  this->coulombA_->clearAll();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.half;

  this->densityA_->vectorize();
  this->coulombA_->vectorize();
  (*this->coulombA_) = (*this->aointegrals_->twoEC_)*(*this->densityA_);
  this->coulombA_->unvectorize();
  this->densityA_->unvectorize();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.two;

  finish = clock();
  this->fileio_->out<<"\nCPU time for building the Coulomb matrix:"<<(finish-start)/CLOCKS_PER_SEC<<" seconds."<<endl;

  if(this->controls_->printLevel>=2) this->coulombA_->printAll(5,this->fileio_->out);
  this->haveCoulomb = true;
};
//----------------------------//
// form the exchange matrix    //
//----------------------------//
void SingleSlater::formExchange(){
  clock_t start,finish;
  start = clock();
  int i;
  
  if(!this->haveDensity) this->formDensity();
  if(!this->aointegrals_->haveAOTwoE) this->aointegrals_->computeAOTwoE();
  this->exchangeA_->clearAll();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.half;

  this->densityA_->vectorize();
  this->exchangeA_->vectorize();
  (*this->exchangeA_)= (*this->aointegrals_->twoEX_)*(*this->densityA_);
  this->exchangeA_->scale(math.quarter);
  this->exchangeA_->unvectorize();
  this->densityA_->unvectorize();

  for(i=0;i<this->nBasis_;i++) (*densityA_)(i,i)*=math.two;

  finish = clock();
  this->fileio_->out<<"\nCPU time for building the Exchange matrix:"<<(finish-start)/CLOCKS_PER_SEC<<" seconds."<<endl;

  if(this->controls_->printLevel>=2) this->exchangeA_->printAll(5,this->fileio_->out);

  this->haveExchange = true;
};
//--------------------------------//
// form the initial guess of MO's //
//--------------------------------//
void SingleSlater::formGuess() {
  if(this->controls_->printLevel>=3) {
    this->moA_->printAll(5,this->fileio_->out);
    if(!this->RHF_) this->moB_->printAll(5,this->fileio_->out);
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
    this->moA_->printAll(5,this->fileio_->out);
    if(!this->RHF_) this->moB_->printAll(5,this->fileio_->out);
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
    this->moA_->printAll(5,this->fileio_->out);
    if(!this->RHF_) this->moB_->printAll(5,this->fileio_->out);
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
  //this->cart_ =new Matrix<double>(3, this->nAtoms_, "Molecule");
  //OOMPI_COMM_WORLD[fromID].Recv(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiRecv(fromID,tag);
};

