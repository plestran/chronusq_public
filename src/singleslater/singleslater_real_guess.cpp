/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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
#include <singleslater.h>
#include <aointegrals.h>
#include <basisset.h>
#include <workers.h>
#ifdef USE_LIBINT
using ChronusQ::BasisSet;
using ChronusQ::Molecule;
using ChronusQ::HashNAOs;
namespace ChronusQ {
/*
template<>
void SingleSlater<double>::placeAtmDen(std::vector<int> atomIndex, 
  SingleSlater<double> &hfA){
  // Place atomic SCF densities in the right place of the total density
  // ** Note: ALWAYS spin average, even for UHF **
  for(auto iAtm : atomIndex){
    auto iBfSt = this->basisset_->mapCen2Bf(iAtm)[0];
    auto iSize = this->basisset_->mapCen2Bf(iAtm)[1]; 
    if(this->nTCS_ == 1){
      this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize)= (*hfA.onePDMA_);
      if(this->isClosedShell){
        if(hfA.isClosedShell) 
          this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.onePDMA_);
        else
          this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.onePDMB_);
      } else {
        this->onePDMB_->block(iBfSt,iBfSt,iSize,iSize)= (*hfA.onePDMA_);
        if(hfA.isClosedShell){
          this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.onePDMA_);
          this->onePDMB_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.onePDMA_);
        } else {
          this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.onePDMB_);
          this->onePDMB_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.onePDMB_);
        }
      }
    } else {
      for(auto I = iBfSt, i = 0; I < (iBfSt +iSize); I += 2, i++)
      for(auto J = iBfSt, j = 0; J < (iBfSt +iSize); J += 2, j++){
        (*this->onePDMA_)(I,J)     = 
          (*hfA.onePDMA_)(i,j) + (*hfA.onePDMB_)(i,j);
        (*this->onePDMA_)(I+1,J+1) = 
          (*hfA.onePDMA_)(i,j) + (*hfA.onePDMB_)(i,j);
      }
    }
  } // loop iAtm
}
*/
/*
template<>
void SingleSlater<double>::scaleDen(){
  // Scale UHF densities according to desired multiplicity
  if(!this->isClosedShell && this->nTCS_ == 1){
    int nE = this->molecule_->nTotalE();
    (*this->onePDMA_) *= (double)this->nAE_/(double)nE ;
    (*this->onePDMB_) *= (double)this->nBE_/(double)nE ;
  } else if(this->nTCS_ == 2) {
    int nE = this->molecule_->nTotalE();
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i += 2)
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j += 2){
      (*this->onePDMA_)(i,j)      *= (double)this->nAE_/(double)nE ;
      (*this->onePDMA_)(i+1,j+1)  *= (double)this->nBE_/(double)nE ;
    }
    
//  (*this->onePDMA_) *= (double)(this->nAE_+this->nBE_)/(double)nE ;
  }
//CErr();
}; // SingleSlater::scaleDen [T=double]
*/

/*
template<>
void SingleSlater<double>::RandomGuess() {
//JJG make random guess
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  this->haveMO = true;
  if(this->molecule_->nAtoms() > 1) this->haveDensity = true;
  // Need to init random number otherwise Eigen is not truly random
  srand((unsigned int) time(0));
  *this->onePDMA_ = RealMatrix::Random(NTCSxNBASIS,NTCSxNBASIS);
  *this->onePDMA_ = this->onePDMA_->selfadjointView<Lower>();
  if(!this->isClosedShell && this->nTCS_ == 1){
    *this->onePDMB_ = RealMatrix::Random(NTCSxNBASIS,NTCSxNBASIS);
    *this->onePDMB_ = this->onePDMB_->selfadjointView<Lower>();
  }  
};
*/


/*
//------------------------------------------//
// form the initial guess of MOs from input //
//------------------------------------------//
template<>
void SingleSlater<double>::readGuessIO() {
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
  if(this->printLevel_ >= 3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(this->Ref_ != RHF) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
};
*/

/*
//-----------------------------------------------------------------------//
// form the initial guess of MOs from Gaussian raw matrix element file   //
//-----------------------------------------------------------------------//
template<>
void SingleSlater<double>::readGuessGauMatEl(GauMatEl& matEl){
  this->fileio_->out << "Reading MO coefficients from " <<matEl.fname()<< endl;
  if(matEl.nBasis()!=this->nBasis_) CErr("Basis Set mismatch",this->fileio_->out);
  double *scr = matEl.readRec(GauMatEl::moa); 
  for(auto i = 0; i < this->nBasis_*this->nBasis_; i++)
    this->moA_->data()[i] = scr[i];
//this->moA_->transposeInPlace(); // Row Major
  delete [] scr;
  if(this->Ref_ != RHF) {
    scr = matEl.readRec(GauMatEl::mob); 
    for(auto i = 0; i < this->nBasis_*this->nBasis_; i++)
      this->moB_->data()[i] = scr[i];
//  this->moB_->transposeInPlace(); // Row Major
    delete [] scr;
  }
  this->matchord();
  if(this->printLevel_ >= 3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(this->Ref_ != RHF) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
}
*/
//-----------------------------------------------------------------------//
// form the initial guess of MOs from Gaussian formatted checkpoint file //
//-----------------------------------------------------------------------//
/*
template<>
void SingleSlater<double>::readGuessGauFChk(std::string &filename) {
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
  if(this->printLevel_ >= 3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(this->Ref_ != RHF) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
};
*/

/*
template<>
void SingleSlater<double>::READGuess(){
  if(getRank() == 0) {
    this->fileio_->out << "Reading SCF Density from disk" << endl;
    H5::DataSpace dataspace = this->fileio_->alphaSCFDen->getSpace();
    this->fileio_->alphaSCFDen->read(this->onePDMA_->data(),H5::PredType::NATIVE_DOUBLE,dataspace,dataspace);
    this->fileio_->alphaMO->read(this->moA_->data(),H5::PredType::NATIVE_DOUBLE,dataspace,dataspace);
    if(!this->isClosedShell && this->nTCS_ == 1){
      this->fileio_->betaSCFDen->read(this->onePDMB_->data(),H5::PredType::NATIVE_DOUBLE,dataspace,dataspace);
      this->fileio_->betaMO->read(this->moB_->data(),H5::PredType::NATIVE_DOUBLE,dataspace,dataspace);
    }
  }
  this->haveMO = true;
  if(this->molecule_->nAtoms() > 1) this->haveDensity = true;
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(this->onePDMA_->data(),
    this->nTCS_*this->nTCS_*this->nBasis_*this->nBasis_,MPI_DOUBLE,0,
    MPI_COMM_WORLD);
  if(!this->isClosedShell && this->nTCS_ == 1)
    MPI_Bcast(this->onePDMB_->data(),
      this->nTCS_*this->nTCS_*this->nBasis_*this->nBasis_,MPI_DOUBLE,0,
      MPI_COMM_WORLD);
#endif
  
}
*/
}; //namespace ChronusQ
#endif
