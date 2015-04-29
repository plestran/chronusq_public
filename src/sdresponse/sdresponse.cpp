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
#include <sdresponse.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
using ChronusQ::SDResponse;
using ChronusQ::SingleSlater;
//------------------------------//
// allocate memory for matrices //
//------------------------------//
void SDResponse::iniSDResponse( Molecule *molecule, BasisSet *basisSet, MOIntegrals *mointegrals, 
                                FileIO *fileio, Controls *controls, SingleSlater *singleSlater) {
  this->nBasis_  = basisSet->nBasis();

  this->molecule_       = molecule;
  this->basisSet_       = basisSet;
  this->fileio_         = fileio;
  this->controls_       = controls;
  this->mointegrals_    = mointegrals;
  this->singleSlater_   = singleSlater;

  this->aoERI_ = singleSlater->aointegrals()->aoERI_;
};
//-----------------------------------//
// print a wave function information //
//-----------------------------------//
void SDResponse::printInfo() {
};
//----------------------//
// compute energies     //
//----------------------//
void SDResponse::computeExcitedStates(){
};
//--------------------//
// print energies     //
//--------------------//
void SDResponse::printExcitedStateEnergies(){
};
//--------------------//
//        form        //
//    < i j | a b >   //
//--------------------//
void SDResponse::formRM(){
  int nOA = this->singleSlater_->nOccA();
  int nO = nOA;
  int nVA = this->singleSlater_->nVirA();
  int nV  = nVA;

  // Copy this->singleSlater_->moA_ (Eigen) to local BTAS tensor
  // Split MO_O & MO_V
  Tensor<double> LocMoAO(this->nBasis_,nO);
  Tensor<double> LocMoAV(this->nBasis_,nV);
  for(auto ii = 0; ii < this->nBasis_; ii++) {
  for(auto jj = 0; jj < nO; jj++) {
    LocMoAO(ii,jj) = (*this->singleSlater_->moA())(ii,jj);
  }
  for(auto kk = nO; kk< this->nBasis_; kk++) {
    LocMoAV(ii,kk-nO) = (*this->singleSlater_->moA())(ii,kk);
  }
  }

  // Create 4 Tensor<double> objects as intermetiates
  Tensor<double> Inter1(nO,this->nBasis_,this->nBasis_,this->nBasis_); // <i nu | lam sig>
  Tensor<double> Inter2(nO,nO,this->nBasis_,this->nBasis_); // <i j  | lam sig>
  Tensor<double> Inter3(nO,nO,nV,this->nBasis_); // <i j  | a   sig>
  Tensor<double> Inter4(nO,nO,nV,nV); // <i j  | a   b  >
  Tensor<double> Inter5(nO,nO,nV,nV);
  Tensor<double> dbbar(nO,nO,nV,nV);

  enum{i,j,a,b,mu,nu,lam,sig};
  // <i nu | lam sig>
  contract(1.0,LocMoAO,{mu,i},(*this->aoERI_),{mu,nu,lam,sig},0.0,Inter1,{i,nu,lam,sig});
  // <i j  | lam sig>
  contract(1.0,LocMoAO,{nu,j},Inter1,{i,nu,lam,sig},0.0,Inter2,{i,j,lam,sig});
  // <i j  | a   sig>
  contract(1.0,Inter2,{i,j,lam,sig},LocMoAV,{lam,a},0.0,Inter3,{i,j,a,sig});
  // <i j  | a   b  >
  contract(1.0,Inter3,{i,j,a,sig},LocMoAV,{sig,b},0.0,Inter4,{i,j,a,b});

    for(auto b=0;b<nV;b++) 
    for(auto a=0;a<nV;a++) 
    for(auto j=0;j<nO;j++) 
    for(auto i=0;i<nO;i++) {
      cout << "( "<< (i+1) << " " 
           << (j+1) << " " 
           << (a+1) << " " 
           << (b+1) << " ) "
           << Inter4(i,j,a,b)<<"\n"; 
    }
//
//  enum{i,j,a,b,mu,nu,lam,sig};
//  // <i nu | lam sig>
//  contract(1.0,LocMoAO(mu,i),this->aoERI_(mu,nu,lam,sig),0.0,Inter1);
//  // <i j  | lam sig>
//  contract(1.0,LocMoAO(nu,j),Inter1(i,nu,lam,sig),0.0,Inter2);
//  // <i j  | a   sig>
//  contract(1.0,Inter2(i,j,lam,sig),LocMoAV(lam,b),0.0,Inter3);
//  // <i j  | a   b  >
//  contract(1.0,Inter3(i,j,a,sig),LocMoAV(sig,a),0.0,Inter5);
//  
//  for(auto i=0;i<nO;i++) 
//  for(auto j=0;j<nO;j++) 
//  for(auto a=0;a<nV;a++) 
//  for(auto b=0;b<nV;b++) {
//    dbbar(i,j,a,b)=Inter4(i,j,a,b)-Inter5(i,j,a,b);
//  }
//  
}

/*************************/
/* MPI Related Routines  */
/*************************/
void SDResponse::mpiSend(int toID,int tag) {
  //OOMPI_COMM_WORLD[toID].Send(this->nAtoms_,tag);
  //OOMPI_COMM_WORLD[toID].Send(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiSend(toID,tag);
};
void SDResponse::mpiRecv(int fromID,int tag) {
  //OOMPI_COMM_WORLD[fromID].Recv(this->nAtoms_,tag);
  //this->index_=new int[this->nAtoms_];
  //this->cart_ =new Matrix(3, this->nAtoms_, "Molecule");
  //OOMPI_COMM_WORLD[fromID].Recv(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiRecv(fromID,tag);
};

