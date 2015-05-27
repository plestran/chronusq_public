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
namespace ChronusQ {
//--------------------------------//
// form the initial guess of MO's //
//--------------------------------//
template<>
void SingleSlater<double>::formGuess() {
  this->moA_->setZero();
  this->haveMO = true;
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
  if(this->controls_->printLevel>=3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
};
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
  if(this->controls_->printLevel>=3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(!this->RHF_) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
};
}; //namespace ChronusQ
