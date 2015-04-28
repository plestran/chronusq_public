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
//------------------------------//
// allocate memory for matrices //
//------------------------------//
void SDResponse::iniSDResponse( Molecule *molecule, BasisSet *basisSet, MOIntegrals *mointegrals, 
                                FileIO *fileio, Controls *controls, SingleSlater *singleSlater) {
  int nBasis  = basisSet->nBasis();

  this->molecule_       = molecule;
  this->basisSet_       = basisSet;
  this->fileio_         = fileio;
  this->controls_       = controls;
  this->mointegrals_    = mointegrals;
  this->singleSlater_   = singleSlater;
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
//    < i j | a b >   //
//--------------------//
void SDResponse::formRM(){
  int i,j,a,b,mu,nu,lam,sig;
  int nBasis = this->basisSet->nBasis();
  // Copy this->singleslater->moA_ (Eigen) to local BTAS tensor
  Tensor<double> this->LocMoA(nBasis,nBasis);
  for(auto i = 0; i < this->basisSet_->nBasis(); i++)
  for(auto j = 0; j < this->basisSet_->nBasis(); j++) {
    //    this->LocalMoA(i,j) = oldMOA(i,j);
    LocMoA(i,j) = this->singleslater->moA_(i,j);
  }
  // Create 4 Tensor<double> objects as intermetiates
  Tensor<double> this->Inter1(nBasis,nBasis,nBasis,nBasis); // <i nu | lam sig>
  Tensor<double> this->Inter2(nBasis,nBasis,nBasis,nBasis); // <i j  | lam sig>
  Tensor<double> this->Inter3(nBasis,nBasis,nBasis,nBasis); // <i j  | a   sig>
  Tensor<double> this->Inter4(nBasis,nBasis,nBasis,nBasis); // <i j  | a   b  >
  // <i nu | lam sig>
  contract(1.0,LocMoA(i,mu),AO,0.0,Inter1);
  // <i j  | lam sig>
  contract(1.0,LocMoA(j,nu),Inter1(i,nu,lam,sig),0.0,Inter2);
  // <i j  | a   sig>
  contract(1.0,LocMoA(lam,a),Inter2(i,j,lam,sig),0.0,Inter3);
  // <i j  | a   b  >
  contract(1.0,LocMoA(sig,b),Inter3(i,j,a,sig),0.0,Inter4);
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

