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
#include "molecule.h"
using ChronusQ::Molecule;
//-------------//
// initializer //
//-------------//
void Molecule::iniMolecule(int nAtoms, FileIO *fileio) {
  if(nAtoms<1) throw 2000;
  if(fileio==NULL) throw 2001;
  this->nAtoms_= nAtoms;
  this->cart_  = new Matrix<double>(3, this->nAtoms_,"Molecular Cartesian");
  this->index_ = new int[this->nAtoms_];
  this->size_  = fileio->sizeInt()*5 + fileio->sizeInt()*nAtoms_ + cart_->size();
};
//--------------------------------------------//
// Read molecular information from input file //
//--------------------------------------------//
void Molecule::readMolecule(FileIO *fileio){
  int i, j, n, readInt;
  char readString[MAXNAMELEN];
  fileio->in >> readInt;
  try{ iniMolecule(readInt,fileio);}
  catch(int msg) {
    fileio->out<<"Molecule initialization failed! E#: "<<msg<<endl;
    exit(1);
  };
  nTotalE_ = 0;
  for(i=0;i<nAtoms_;i++) {
    fileio->in >> readString;
    fileio->in >> readInt;
    if((n=HashAtom(readString,readInt))!=-1) index_[i] = n;
    else {
      fileio->out <<"Error: invalid atomic symbol or mass number!"<<endl;
      throw 2002;
    };
    nTotalE_ += atom[n].atomicNumber;
    fileio->in >> (*cart_)(0,i);
    (*cart_)(0,i) = (*cart_)(0,i)/phys.bohr;
    fileio->in >> (*cart_)(1,i);
    (*cart_)(1,i) = (*cart_)(1,i)/phys.bohr;
    fileio->in >> (*cart_)(2,i);
    (*cart_)(2,i) = (*cart_)(2,i)/phys.bohr;
  };
  double sqrAB;
  energyNuclei_ = 0.0;
  for(i=0;i<nAtoms_;i++) 
    for(j=i+1;j<nAtoms_;j++) {
      sqrAB = 0.0;
      for(n=0;n<3;n++) sqrAB += ((*cart_)(n,i)-(*cart_)(n,j))*((*cart_)(n,i)-(*cart_)(n,j));
      energyNuclei_ += atom[index_[i]].atomicNumber*atom[index_[j]].atomicNumber/sqrt(sqrAB);
  };

};
//---------------------------------------------------//
// Print out molecular carteisan coordinates in bohr //
//---------------------------------------------------//
void Molecule::printInfo(FileIO *fileio,Controls *controls) {
  fileio->out.precision(8);
  fileio->out.fill(' ');
  fileio->out.setf(ios::right,ios::adjustfield);
  fileio->out.setf(ios::fixed,ios::floatfield);
  fileio->out<<"\nCartesian coordinates (bohr):"<<endl;
  fileio->out<<bannerTop<<endl;
  fileio->out<<std::setw(24)<<" "<<std::setw(15)<<"X"<<std::setw(15)<<"Y"<<std::setw(15)<<"Z"<<endl;
  fileio->out<<bannerMid<<endl;
  for(int i=0;i<nAtoms_;i++)
    fileio->out<<std::setw(8)<<i+1<<std::setw(8)<<atom[index_[i]].symbol<<std::setw(8)<<atom[index_[i]].atomicNumber
	       <<std::setw(15)<<(*cart_)(0,i)<<std::setw(15)<<(*cart_)(1,i)<<std::setw(15)<<(*cart_)(2,i)<<endl;
  fileio->out<<bannerEnd<<endl;
  fileio->out<<"\nMolecular Information:"<<endl;
  fileio->out<<std::setw(15)<<"nAtoms ="<<std::setw(8)<<nAtoms_<<std::setw(5)<<" "
	     <<std::setw(20)<<"Charge ="<<std::setw(8)<<charge_<<endl;
  fileio->out<<std::setw(15)<<"nElectrons ="<<std::setw(8)<<nTotalE_<<endl;
};
//--------------------------------//
// read from binary files //
//--------------------------------//
void Molecule::ioRead(FileIO *fileio) {
  const int nInteger = 5;
  int storage[nInteger];
  try { fileio->io("R",blockMolecule,"BIN",storage,nInteger,0);}
  catch (int msg) {
    fileio->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);    
  };
  nAtoms_ = storage[0];
  charge_ = storage[1];
  spin_   = storage[2];
  nTotalE_= storage[3];
  size_   = storage[4];
  if(cart_!=NULL) delete[] cart_;
  if(index_!=NULL) delete[] index_;
  cart_  = new Matrix<double>(3, nAtoms_,"Molecular Cartesian");
  index_ = new int[nAtoms_];
  try { fileio->io("R",blockMolecule,"BIN",index_,nAtoms_);}
  catch (int msg) {
    fileio->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);    
  };
  cart_->ioRead(fileio,blockMolecule,"BIN");
};
//-----------------------//
// write to binary files //
//----------------------//
void Molecule::ioWrite(FileIO *fileio) {
  if(!fileio->isOpen(blockMolecule)) fileio->iniBlock(blockMolecule);
  const int nInteger = 5;
  int storage[nInteger];
  storage[0] = nAtoms_;
  storage[1] = charge_;
  storage[2] = spin_;
  storage[3] = nTotalE_;
  storage[4] = size_;
  try { fileio->io("W",blockMolecule,"BIN",storage,nInteger,0);}
  catch (int msg) {
    fileio->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);    
  };
  try { fileio->io("W",blockMolecule,"BIN",index_,nAtoms_);}
  catch (int msg) {
    fileio->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);    
  };
  cart_->ioWrite(fileio,blockMolecule,"BIN");
};
/*************************/
/* MPI Related Routines  */
/*************************/
/*
void Molecule::mpiSend(int toID,int tag) {
  OOMPI_COMM_WORLD[toID].Send(this->nAtoms_,tag);
  OOMPI_COMM_WORLD[toID].Send(this->index_,this->nAtoms_,tag);
  this->cart_->mpiSend(toID,tag);
};

void Molecule::mpiRecv(int fromID,int tag) {
  OOMPI_COMM_WORLD[fromID].Recv(this->nAtoms_,tag);
  this->index_=new int[this->nAtoms_];
  this->cart_ =new Matrix<double>(3, this->nAtoms_, "Molecule");
  OOMPI_COMM_WORLD[fromID].Recv(this->index_,this->nAtoms_,tag);
  this->cart_->mpiRecv(fromID,tag);
};
*/
