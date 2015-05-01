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
void SDResponse::iniSDResponse( std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisSet, std::shared_ptr<MOIntegrals> mointegrals, 
                                std::shared_ptr<FileIO> fileio, std::shared_ptr<Controls> controls, std::shared_ptr<SingleSlater> singleSlater) {
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
  cout << "Number of Occupied: " << nO << ", Number of Virtual " << nV << " Number of Basis: "<<this->nBasis_<< ".\n";
  // Copy this->singleSlater_->moA_ (Eigen) to local BTAS tensor
  // Split MO_O & MO_V
  Tensor<double> LocMoAO(this->nBasis_,nO);
  Tensor<double> LocMoAV(this->nBasis_,nV);
  for(auto ii = 0; ii < this->nBasis_; ii++) {
  for(auto jj = 0; jj < nO; jj++) {
    LocMoAO(ii,jj) = (*this->singleSlater_->moA())(ii,jj);
    cout << "The MO_NO: "<< ii<<jj<< " " <<LocMoAO(ii,jj) << "\n";
  }
  for(auto kk = nO; kk< this->nBasis_; kk++) {
    LocMoAV(ii,kk-nO) = (*this->singleSlater_->moA())(ii,kk);
    cout << "The MO_NV: "<< ii<<(kk-nO)<< LocMoAO(ii,kk-nO) << "\n";

  }
  }
  // Create the A matrix <aj||ib>
  // Create 4 Tensor<double> objects as intermetiates
  Tensor<double> InterA1(nV,this->nBasis_,this->nBasis_,this->nBasis_); // <a nu | lam sig>
  Tensor<double> InterA2(nV,nO,this->nBasis_,this->nBasis_); // <a j  | lam sig>
  Tensor<double> InterA3(nV,nO,nO,this->nBasis_); // <a j  | i   sig>
  Tensor<double> InterA4(nV,nO,nO,nV); // < a j | i b >
  Tensor<double> InterA5(nV,nO,nV,this->nBasis_); // < a j | b sig >
  Tensor<double> InterA6(nV,nO,nV,nO); //< a j | b i >
  Tensor<double> dbbarA(nV,nO,nO,nV); // <aj||ib>

  enum{a,j,i,b,mu,nu,lam,sig};

  // <aj|ib>
  // <a nu | lam sig>
  contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,InterA1,{a,nu,lam,sig});
  cout <<"finish 1\n";
  // <a j  | lam sig>
  contract(1.0,LocMoAO,{nu,j},InterA1,{a,nu,lam,sig},0.0,InterA2,{a,j,lam,sig});
  cout <<"finish 2\n";
  // <i j  | a   sig>
  contract(1.0,InterA2,{a,j,lam,sig},LocMoAO,{lam,i},0.0,InterA3,{a,j,i,sig});
  cout <<"finish 3\n";
  // <i j  | a   b  >
  contract(1.0,InterA3,{a,j,i,sig},LocMoAV,{sig,b},0.0,InterA4,{a,j,i,b});
  cout <<"finish 4\n";

  for(auto b=0;b<nV;b++) 
  for(auto i=0;i<nO;i++) 
  for(auto j=0;j<nO;j++) 
  for(auto a=0;a<nV;a++) {
    cout << "( "<< (a+1) << " " 
         << (j+1) << " " 
         << (i+1) << " " 
         << (b+1) << " ) "
         << InterA4(a,j,i,b)<<"\n"; 
  }

  // <aj|bi>
  // <a nu | lam sig>
  contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,InterA1,{a,nu,lam,sig});
  // <a j  | lam sig>
  contract(1.0,LocMoAO,{nu,j},InterA1,{a,nu,lam,sig},0.0,InterA2,{a,j,lam,sig});
  // <a j  | b   sig>
  contract(1.0,InterA2,{a,j,lam,sig},LocMoAV,{lam,b},0.0,InterA5,{a,j,b,sig});
  // <a j  | b   i  >
  contract(1.0,InterA5,{a,j,b,sig},LocMoAO,{sig,i},0.0,InterA6,{a,j,b,i});

  for(auto a=0;a<nV;a++) 
  for(auto j=0;j<nO;j++) 
  for(auto i=0;i<nO;i++) 
  for(auto b=0;b<nV;b++) {
    dbbarA(a,j,i,b)=InterA4(a,j,i,b)-InterA6(a,j,b,i);
  }
  //Print A matrix
  cout << "Here is Matrix A: <aj||ib>\n";
  for(auto b=0;b<nV;b++) 
  for(auto i=0;i<nO;i++) 
  for(auto j=0;j<nO;j++) 
  for(auto a=0;a<nV;a++) {
    cout << "( "<< (a+1) << " " 
         << (j+1) << " " 
         << (i+1) << " " 
         << (b+1) << " ) "
         << dbbarA(a,j,i,b)<<"\n"; 
  }

  // Create the B matrix <ij||ab>
  // Create 4 Tensor<double> objects as intermetiates
  Tensor<double> InterB1(nO,this->nBasis_,this->nBasis_,this->nBasis_); // <i nu | lam sig>
  Tensor<double> InterB2(nO,nO,this->nBasis_,this->nBasis_); // <i j  | lam sig>
  Tensor<double> InterB3(nO,nO,nV,this->nBasis_); // <i j  | a   sig>
  Tensor<double> InterB4(nO,nO,nV,nV); // < i j | a b >
  Tensor<double> InterB5(nO,nO,nV,nV); // < i j | b a >
  Tensor<double> dbbarB(nO,nO,nV,nV); // <ij||ab>
  
  // <ij|ab>
  // <i nu | lam sig>
  contract(1.0,LocMoAO,{mu,i},(*this->aoERI_),{mu,nu,lam,sig},0.0,InterB1,{i,nu,lam,sig});
  cout << "finish 1\n";
  // <i j  | lam sig>
  contract(1.0,LocMoAO,{nu,j},InterB1,{i,nu,lam,sig},0.0,InterB2,{i,j,lam,sig});
  cout << "finish 2\n";
  // <i j  | a   sig>
  contract(1.0,InterB2,{i,j,lam,sig},LocMoAV,{lam,a},0.0,InterB3,{i,j,a,sig});
  cout << "finish 3\n";
  // <i j  | a   b  >
  contract(1.0,InterB3,{i,j,a,sig},LocMoAV,{sig,b},0.0,InterB4,{i,j,a,b});
  cout << "finish 4\n";

//  for(auto b=0;b<nV;b++) 
//  for(auto a=0;a<nV;a++) 
//  for(auto j=0;j<nO;j++) 
//  for(auto i=0;i<nO;i++) {
//    cout << "( "<< (i+1) << " " 
//         << (j+1) << " " 
//         << (a+1) << " " 
//         << (b+1) << " ) "
//         << Inter4(i,j,a,b)<<"\n"; 
//  }
  
  // <ij|ba>
  // <i nu | lam sig>
  contract(1.0,LocMoAO,{mu,i},(*this->aoERI_),{mu,nu,lam,sig},0.0,InterB1,{i,nu,lam,sig});
  // <i j  | lam sig>
  contract(1.0,LocMoAO,{nu,j},InterB1,{i,nu,lam,sig},0.0,InterB2,{i,j,lam,sig});
  // <i j  | b   sig>
  contract(1.0,InterB2,{i,j,lam,sig},LocMoAV,{lam,b},0.0,InterB3,{i,j,b,sig});
  // <i j  | b   a  >
  contract(1.0,InterB3,{i,j,b,sig},LocMoAV,{sig,a},0.0,InterB4,{i,j,b,a});

  for(auto i=0;i<nO;i++) 
  for(auto j=0;j<nO;j++) 
  for(auto a=0;a<nV;a++) 
  for(auto b=0;b<nV;b++) {
    dbbarB(i,j,a,b)=InterB4(i,j,a,b)-InterB5(i,j,b,a);
  }
  //Print B matrix
  cout << "Here is Matrix B: <ij||ab>\n";
  for(auto b=0;b<nV;b++) 
  for(auto a=0;a<nV;a++) 
  for(auto j=0;j<nO;j++) 
  for(auto i=0;i<nO;i++) {
    cout << "( "<< (i+1) << " " 
         << (j+1) << " " 
         << (a+1) << " " 
         << (b+1) << " ) "
         << dbbarB(i,j,a,b)<<"\n"; 
  }
 
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

