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
using std::cout;
using std::setw;
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
  this->aoERI_          = singleSlater->aointegrals()->aoERI_;
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
  }
  for(auto kk = nO; kk< this->nBasis_; kk++) {
    LocMoAV(ii,kk-nO) = (*this->singleSlater_->moA())(ii,kk);
  }
  }
  // Create the A matrix
  // Ixxxx for intermediate Sxxxx for Mulliken Notation single bar
  // Dxxxx for Dirac Notation double bar, dxxxx for Dirac Notation single bar
  Tensor<double> Ianls(nV,this->nBasis_,this->nBasis_,this->nBasis_); // (a nu | lam sig)
  Tensor<double> Iails(nV,nO,this->nBasis_,this->nBasis_); // (a j  | lam sig)
  Tensor<double> Iaijs(nV,nO,nO,this->nBasis_); // (a j  | i   sig)
  Tensor<double> Saijb(nV,nO,nO,nV); // ( a j | i b )
  Tensor<double> Iabls(nV,nV,this->nBasis_,this->nBasis_); // ( a j | b sig )
  Tensor<double> Iabjs(nV,nV,nO,this->nBasis_); // ( a j | b i )
  Tensor<double> Sabji(nV,nV,nO,nO); // ( a b | j i )
  Tensor<double> Dajib(nV,nO,nO,nV); // <aj||ib>
  Tensor<double> dajib(nV,nO,nO,nV); // <aj|ib>
  Tensor<double> Dabij(nV,nV,nO,nO); // <ab||ij>
  Tensor<double> dabij(nV,nV,nO,nO); // <ab|ij>

  enum{a,j,i,b,mu,nu,lam,sig};

  // (ai|jb)
  // (a nu | lam sig)
  contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,Ianls,{a,nu,lam,sig});
  // (a i  | lam sig)
  contract(1.0,LocMoAO,{nu,i},Ianls,{a,nu,lam,sig},0.0,Iails,{a,i,lam,sig});
  // (a i  | j   sig)
  contract(1.0,Iails,{a,i,lam,sig},LocMoAO,{lam,j},0.0,Iaijs,{a,i,j,sig});
  // (a i  | j   b  )
  contract(1.0,Iaijs,{a,i,j,sig},LocMoAV,{sig,b},0.0,Saijb,{a,i,j,b});

  // (ab|ji)
  // (a nu | lam sig)
  contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,Ianls,{a,nu,lam,sig});
  // (a b  | lam sig)
  contract(1.0,LocMoAV,{nu,b},Ianls,{a,nu,lam,sig},0.0,Iabls,{a,b,lam,sig});
  // (a b  | j   sig)
  contract(1.0,Iabls,{a,b,lam,sig},LocMoAO,{lam,j},0.0,Iabjs,{a,b,j,sig});
  // (a b  | j   i  )
  contract(1.0,Iabjs,{a,b,j,sig},LocMoAO,{sig,i},0.0,Sabji,{a,b,j,i});

  // Build <aj||ib>
  for(auto a=0;a<nV;a++) 
  for(auto j=0;j<nO;j++) 
  for(auto i=0;i<nO;i++) 
  for(auto b=0;b<nV;b++) {
    Dajib(a,j,i,b) = Saijb(a,i,j,b)-Sabji(a,b,j,i);
    dajib(a,j,i,b) = Saijb(a,i,j,b);
  }
  //Build <ab||ij>
  for (auto a=0;a<nV;a++)
  for (auto b=0;b<nV;b++)
  for (auto i=0;i<nO;i++)
  for (auto j=0;j<nO;j++){
    Dabij(a,b,i,j) = Saijb(a,i,j,b) - Saijb(a,j,i,b);
    dabij(a,b,i,j) = Saijb(a,i,j,b);
  }

  // Build A & B matrix
  int nOV = nO*nV;
  int ia,jb;
  RealMatrix ABBA(4*nOV,4*nOV);
  RealMatrix A(2*nOV,2*nOV);
  RealMatrix B(2*nOV,2*nOV);
  RealMatrix Ad(nOV,nOV);
  RealMatrix Bd(nOV,nOV);
  RealMatrix Aod(nOV,nOV);
  RealMatrix Bod(nOV,nOV);
  RealMatrix EigO(nO,1);
  RealMatrix EigV(nV,1);

  for (auto i=0;i<nO;i++){
    EigO(i,0) = (*this->singleSlater_->epsA())(i,0);
    cout << "The " << (i+1) << " eigenvalue in Occupied is: " << EigO(i,0) << endl;
  }
  for (auto j=0;j<nV;j++){
    EigV(j,0) = (*this->singleSlater_->epsA())((j+nO),0);
    cout << "The " << (j+1) << " eigenvalue in Virtual is: " << EigV(j,0) << endl;
  }

  ia = 0;
  for (auto i=0;i<nO;i++)
  for (auto a=0;a<nV;a++){
    jb =0;
    for (auto j=0;j<nO;j++)
    for (auto b=0;b<nV;b++){
      Ad(ia,jb) = 0.0;
      if ((a==b)&&(i==j)){
	Ad(ia,jb) = EigV(a,0)-EigO(i,0);
      }
      Ad(ia,jb) = Ad(ia,jb) + Dajib(a,j,i,b);
      Bd(ia,jb) = Dabij(a,b,i,j);
      Aod(ia,jb) = dajib(a,j,i,b);
      Bod(ia,jb) = dabij(a,b,i,j);
      jb = jb+1;
    }
    ia = ia+1;
  }
  A.block(0,0,nOV,nOV) = Ad;
  A.block(nOV,nOV,nOV,nOV) = Ad;
  A.block(0,nOV,nOV,nOV) = Aod;
  A.block(nOV,0,nOV,nOV) = Aod;
  B.block(0,0,nOV,nOV) = Bd;
  B.block(nOV,nOV,nOV,nOV) = Bd;
  B.block(0,nOV,nOV,nOV) = Aod;
  B.block(nOV,0,nOV,nOV) = Aod;

  for (auto i=0;i<2*nOV;i++)
  for (auto j=0;j<2*nOV;j++){
    cout << "A symmetry " << A(i,j)-A(j,i) <<endl;
    cout << "B symmetry " << B(i,j)-B(j,i) <<endl;
  }

  // Build the ABBA matrix

  ABBA.block(0,0,2*nOV,2*nOV) = A;
  ABBA.block(2*nOV,2*nOV,2*nOV,2*nOV) = -A;
  ABBA.block(0,2*nOV,2*nOV,2*nOV) = B;
  ABBA.block(2*nOV,0,2*nOV,2*nOV) = -B;

  // Print the A & B matrix
  cout << "Print the A matrix" << endl;
  for (auto i=0;i<2*nOV;i++)
  for (auto j=0;j<2*nOV;j++){
    cout << "(" << (i+1) << ", "
	 << (j+1) << ")" << "  "
	 << A(i,j) << endl;
  }
  cout << "Print the B matrix" << endl;
  for (auto i=0;i<2*nOV;i++)
  for (auto j=0;j<2*nOV;j++){
    cout << "(" << (i+1) << ", "
	 << (j+1) << ")" << "  "
	 << B(i,j) << endl;
  }
  // CIS routine
  // Diagonalize the A matrix
  Eigen::SelfAdjointEigenSolver<RealMatrix> CIS;
  CIS.compute(A);
  CIS.eigenvalues();
  CIS.eigenvectors();
  
  // Print the CIS Excitation Energies
  for (auto i=0;i<2*nOV;i++){
    cout << "The " << (i+1) << " CIS Exicitation Energy is: "
         << (CIS.eigenvalues())(i) << endl;
  }
  
  RealMatrix AXs(4*nOV,4*nOV);
  AXs = A * CIS.eigenvectors();
  // Print the AX matrix
  for (auto i=0;i<4*nOV;i++)
  for (auto j=0;j<4*nOV;j++){
    cout << "( "<< i+1 << ", "<< j+1 << " )      "
         << AXs(i,j)<<endl;
  }

  // LR TDHF routine
  Eigen::EigenSolver<RealMatrix> TD;
  TD.compute(ABBA);
  TD.eigenvalues();
  //TD.eigenvectors();

  // Print the LR-TDHF Excitation Energies
  for (auto i=0;i<4*nOV;i++){
    cout << "The " << (i+1) << " LR-TDHF Exicitation Energy is: "
         << (TD.eigenvalues())(i) << endl;
  }

}

RealMatrix SDResponse::formRM2(int NTrial){

  // Convert X from MO to AO X(AO) = C(occ) * X(MO) * C(vir)^T
  // Contract X(AO) with (mu nu || lambda sigma)  IX_{ij} = (ij||kl) X(AO)_{kl}
  // Transform IX from AO to MO  IX(MO) = C(occ)^T * IX(AO) * C(vir)
  // Contract IX(MO) with eigenvalue differences
  //   AX_{ia} = IX(ia) / (e(a) - e(i))

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

