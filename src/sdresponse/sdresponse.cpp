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
#include <davidson.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
using ChronusQ::SDResponse;
using ChronusQ::SingleSlater;
using ChronusQ::Davidson; 
using std::cout;
using std::setw;
//------------------------------//
// allocate memory for matrices //
//------------------------------//
void SDResponse::iniSDResponse( Molecule * molecule, BasisSet * basisSet, MOIntegrals * mointegrals, 
                                FileIO * fileio, Controls * controls, SingleSlater<double> * singleSlater) {
  this->nBasis_  = basisSet->nBasis();

  this->molecule_       = molecule;
  this->basisSet_       = basisSet;
  this->fileio_         = fileio;
  this->controls_       = controls;
  this->mointegrals_    = mointegrals;
  this->singleSlater_   = singleSlater;
  this->aoERI_          = singleSlater->aointegrals()->aoERI_.get();
  this->elecDipole_     = singleSlater->aointegrals()->elecDipole_.get();
  cout << "Allocate Memory for CISEnergy_ and CISTransDen_" << endl;
  int nOVA = singleSlater_->nOVA();
  int nOVB = singleSlater_->nOVB();
  this->CISEnergy_ = std::unique_ptr<RealMatrix>(new RealMatrix(nOVA+nOVB,1));
  this->CISTransDen_ = std::unique_ptr<RealMatrix>(new RealMatrix(nOVA+nOVB,nOVA+nOVB));
  this->TransDipole_ = std::unique_ptr<RealMatrix>(new RealMatrix(1,3));
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
  // Get info from SCF and make Local copy of MO
  int nOA = this->singleSlater_->nOccA();
  int nOB = this->singleSlater_->nOccB();
  int nVA = this->singleSlater_->nVirA();
  int nVB = this->singleSlater_->nVirB();
  int nOVA = this->singleSlater_->nOVA();
  int nOVB = this->singleSlater_->nOVB();
  cout << "Number of Occupied Alpha: " << nOA << ", Number of Virtual Alpha" << nVA << " Number of Basis: "<<this->nBasis_<< ".\n";
  cout << "Number of Occupied Beta: " << nOB << ", Number of Virtual Beta" << nVB << " Number of Basis: "<<this->nBasis_<< ".\n";
  Tensor<double> LocMoAO(this->nBasis_,nOA);
  Tensor<double> LocMoAV(this->nBasis_,nVA);
  Tensor<double> LocMoBO(this->nBasis_,nOB);
  Tensor<double> LocMoBV(this->nBasis_,nVB);

  for(auto ii = 0; ii < this->nBasis_; ii++) {
  for(auto jj = 0; jj < nOA; jj++) {
    LocMoAO(ii,jj) = (*this->singleSlater_->moA())(ii,jj);
  }
  for(auto kk = nOA; kk< this->nBasis_; kk++) {
    LocMoAV(ii,kk-nOA) = (*this->singleSlater_->moA())(ii,kk);
  }
  }
  cout << "MoA" << endl;
  cout << (*this->singleSlater_->moA())<<endl;
  cout << "MoB" << endl;
  cout << (*this->singleSlater_->moB())<<endl;

  for(auto ii = 0; ii < this->nBasis_; ii++) {
  for(auto jj = 0; jj < nOB; jj++) {
    LocMoBO(ii,jj) = (*this->singleSlater_->moB())(ii,jj);
  }
  for(auto kk = nOB; kk< this->nBasis_; kk++) {
    LocMoBV(ii,kk-nOB) = (*this->singleSlater_->moB())(ii,kk);
  }
  }
  cout << "Storge local copy of MO" << endl;

  // Prepare the 2e integrals
  // Ixxxx for intermediate Sxxxx for Mulliken Notation single bar
  // Dxxxx for Dirac Notation double bar, dxxxx for Dirac Notation single bar
  Tensor<double> IanlsA(nVA,this->nBasis_,this->nBasis_,this->nBasis_); // (a nu | lam sig)
  Tensor<double> IanlsB(nVB,this->nBasis_,this->nBasis_,this->nBasis_);
  Tensor<double> IailsA(nVA,nOA,this->nBasis_,this->nBasis_); // (a j  | lam sig)
  Tensor<double> IailsB(nVB,nOB,this->nBasis_,this->nBasis_);
  Tensor<double> IaijsAA(nVA,nOA,nOA,this->nBasis_); // (a j  | i   sig)
  Tensor<double> IaijsAB(nVA,nOA,nOB,this->nBasis_);
  Tensor<double> IaijsBA(nVB,nOB,nOA,this->nBasis_);
  Tensor<double> IaijsBB(nVB,nOB,nOB,this->nBasis_);
  Tensor<double> SaijbAA(nVA,nOA,nOA,nVA); // ( a j | i b )
  Tensor<double> SaijbAB(nVA,nOA,nOB,nVB);
  Tensor<double> SaijbBA(nVB,nOB,nOA,nVA);
  Tensor<double> SaijbBB(nVB,nOB,nOB,nVB);
  Tensor<double> IablsA(nVA,nVA,this->nBasis_,this->nBasis_); // ( a j | b sig )
  Tensor<double> IablsB(nVB,nVB,this->nBasis_,this->nBasis_);
  Tensor<double> IabjsA(nVA,nVA,nOA,this->nBasis_); // ( a j | b i )
  Tensor<double> IabjsB(nVB,nVB,nOB,this->nBasis_);
  Tensor<double> SabjiAA(nVA,nVA,nOA,nOA); // ( a b | j i )
  Tensor<double> SabjiBB(nVB,nVB,nOB,nOB);
  Tensor<double> DajibAA(nVA,nOA,nOA,nVA); // <aj||ib>
  Tensor<double> DajibBB(nVB,nOB,nOB,nVB);
  Tensor<double> dajibAB(nVA,nOB,nOA,nVB); // <aj|ib>
  Tensor<double> dajibBA(nVB,nOA,nOB,nVA);
  Tensor<double> DabijAA(nVA,nVA,nOA,nOA); // <ab||ij>
  Tensor<double> DabijBB(nVB,nVB,nOB,nOB);
  Tensor<double> dabijAB(nVA,nVB,nOA,nOB); // <ab|ij>
  Tensor<double> dabijBA(nVB,nVA,nOB,nOA);

  enum{a,j,i,b,mu,nu,lam,sig};

  // (ai|jb)AAAA
  // (a nu | lam sig)
  contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsA,{a,nu,lam,sig});
  // (a i  | lam sig)
  contract(1.0,LocMoAO,{nu,i},IanlsA,{a,nu,lam,sig},0.0,IailsA,{a,i,lam,sig});
  // (a i  | j   sig)
  contract(1.0,IailsA,{a,i,lam,sig},LocMoAO,{lam,j},0.0,IaijsAA,{a,i,j,sig});
  // (a i  | j   b  )
  contract(1.0,IaijsAA,{a,i,j,sig},LocMoAV,{sig,b},0.0,SaijbAA,{a,i,j,b});
  cout << "1"<<endl;
  // (ai|jb)AABB
  // (a nu | lam sig)
  contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsA,{a,nu,lam,sig});
  // (a i  | lam sig)
  contract(1.0,LocMoAO,{nu,i},IanlsA,{a,nu,lam,sig},0.0,IailsA,{a,i,lam,sig});
  // (a i  | j   sig)
  contract(1.0,IailsA,{a,i,lam,sig},LocMoBO,{lam,j},0.0,IaijsAB,{a,i,j,sig});
  // (a i  | j   b  )
  contract(1.0,IaijsAB,{a,i,j,sig},LocMoBV,{sig,b},0.0,SaijbAB,{a,i,j,b});
  cout << "2" <<endl;
  // (ai|jb)BBAA
  // (a nu | lam sig)
  contract(1.0,LocMoBV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsB,{a,nu,lam,sig});
  // (a i  | lam sig)
  contract(1.0,LocMoBO,{nu,i},IanlsB,{a,nu,lam,sig},0.0,IailsB,{a,i,lam,sig});
  // (a i  | j   sig)
  contract(1.0,IailsB,{a,i,lam,sig},LocMoAO,{lam,j},0.0,IaijsBA,{a,i,j,sig});
  // (a i  | j   b  )
  contract(1.0,IaijsBA,{a,i,j,sig},LocMoAV,{sig,b},0.0,SaijbBA,{a,i,j,b});
  cout << "3" << endl;
  // (ai|jb)BBBB
  // (a nu | lam sig)
  contract(1.0,LocMoBV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsB,{a,nu,lam,sig});
  // (a i  | lam sig)
  contract(1.0,LocMoBO,{nu,i},IanlsB,{a,nu,lam,sig},0.0,IailsB,{a,i,lam,sig});
  // (a i  | j   sig)
  contract(1.0,IailsB,{a,i,lam,sig},LocMoBO,{lam,j},0.0,IaijsBB,{a,i,j,sig});
  // (a i  | j   b  )
  contract(1.0,IaijsBB,{a,i,j,sig},LocMoBV,{sig,b},0.0,SaijbBB,{a,i,j,b});
  cout << "4" <<endl;
  // (ab|ji)AAAA
  // (a nu | lam sig)
  contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsA,{a,nu,lam,sig});
  // (a b  | lam sig)
  contract(1.0,LocMoAV,{nu,b},IanlsA,{a,nu,lam,sig},0.0,IablsA,{a,b,lam,sig});
  // (a b  | j   sig)
  contract(1.0,IablsA,{a,b,lam,sig},LocMoAO,{lam,j},0.0,IabjsA,{a,b,j,sig});
  // (a b  | j   i  )
  contract(1.0,IabjsA,{a,b,j,sig},LocMoAO,{sig,i},0.0,SabjiAA,{a,b,j,i});
  cout << "5" << endl;
  // (ab|ji)BBBB
  // (a nu | lam sig)
  contract(1.0,LocMoBV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsB,{a,nu,lam,sig});
  // (a b  | lam sig)
  contract(1.0,LocMoBV,{nu,b},IanlsB,{a,nu,lam,sig},0.0,IablsB,{a,b,lam,sig});
  // (a b  | j   sig)
  contract(1.0,IablsB,{a,b,lam,sig},LocMoBO,{lam,j},0.0,IabjsB,{a,b,j,sig});
  // (a b  | j   i  )
  contract(1.0,IabjsB,{a,b,j,sig},LocMoBO,{sig,i},0.0,SabjiBB,{a,b,j,i});
  cout << "6" << endl;

  // Build <aj||ib>AAAA
  for(auto a=0;a<nVA;a++) 
  for(auto j=0;j<nOA;j++) 
  for(auto i=0;i<nOA;i++) 
  for(auto b=0;b<nVA;b++) {
    DajibAA(a,j,i,b) = SaijbAA(a,i,j,b)-SabjiAA(a,b,j,i);
  }
  cout << "7" << endl;
  // Build <aj||ib>ABAB
  for(auto a=0;a<nVA;a++) 
  for(auto j=0;j<nOB;j++) 
  for(auto i=0;i<nOA;i++) 
  for(auto b=0;b<nVB;b++) {
    dajibAB(a,j,i,b) = SaijbAB(a,i,j,b);
  }
  cout << "8" << endl;
  // Build <aj||ib>BABA
  for(auto a=0;a<nVB;a++)
  for(auto j=0;j<nOA;j++)
  for(auto i=0;i<nOB;i++)
  for(auto b=0;b<nVA;b++) {
    dajibBA(a,j,i,b) = SaijbBA(a,i,j,b);
  }
  cout << "9" << endl;
  // Build <aj||ib>BBBB
  for(auto a=0;a<nVB;a++)
  for(auto j=0;j<nOB;j++)
  for(auto i=0;i<nOB;i++)
  for(auto b=0;b<nVB;b++) {
    DajibBB(a,j,i,b) = SaijbBB(a,i,j,b)-SabjiBB(a,b,j,i);
  }
  cout << "10" << endl;
  //Build <ab||ij>AAAA
  for (auto a=0;a<nVA;a++)
  for (auto b=0;b<nVA;b++)
  for (auto i=0;i<nOA;i++)
  for (auto j=0;j<nOA;j++){
    DabijAA(a,b,i,j) = SaijbAA(a,i,j,b) - SaijbAA(a,j,i,b);
  }
  cout << "11" << endl;
  //Build <ab||ij> ABAB
  for (auto a=0;a<nVA;a++)
  for (auto b=0;b<nVB;b++)
  for (auto i=0;i<nOA;i++)
  for (auto j=0;j<nOB;j++){
    dabijAB(a,b,i,j) = SaijbAB(a,i,j,b);
  }
  cout << "12" << endl;
  //Build <ab||ij> BABA
  for (auto a=0;a<nVB;a++)
  for (auto b=0;b<nVA;b++)
  for (auto i=0;i<nOB;i++)
  for (auto j=0;j<nOA;j++){
    dabijBA(a,b,i,j) = SaijbBA(a,i,j,b);
  }
  cout << "13"<< endl;
  //Build <ab||ij> BBBB
  for (auto a=0;a<nVB;a++)
  for (auto b=0;b<nVB;b++)
  for (auto i=0;i<nOB;i++)
  for (auto j=0;j<nOB;j++){
    DabijBB(a,b,i,j) = SaijbBB(a,i,j,b) - SaijbBB(a,j,i,b);
  }
  cout << "14" << endl;
  // Build A & B matrix
  int ia,jb;
  RealMatrix ABBA(2*(nOVA+nOVB),2*(nOVA+nOVB));
  RealMatrix A(nOVA+nOVB,nOVA+nOVB);
  RealMatrix B(nOVA+nOVB,nOVA+nOVB);
  RealMatrix Aud(nOVA,nOVA);
  RealMatrix Add(nOVB,nOVB);
  RealMatrix Bud(nOVA,nOVA);
  RealMatrix Bdd(nOVB,nOVB);
  RealMatrix Auod(nOVA,nOVB);
  RealMatrix Adod(nOVB,nOVA);
  RealMatrix Buod(nOVA,nOVB);
  RealMatrix Bdod(nOVB,nOVA);
  RealMatrix EigAO(nOA,1);
  RealMatrix EigAV(nVA,1);
  RealMatrix EigBO(nOB,1);
  RealMatrix EigBV(nVB,1);

  for (auto i=0;i<nOA;i++){
    EigAO(i) = (*this->singleSlater_->epsA())(i);
    cout << "The " << (i+1) << " eigenvalue in Occupied Alpha is: " << EigAO(i) << endl;
  }
  for (auto j=0;j<nVA;j++){
    EigAV(j) = (*this->singleSlater_->epsA())((j+nOA));
    cout << "The " << (j+1) << " eigenvalue in Virtual Alpha is: " << EigAV(j) << endl;
  }
  for (auto i=0;i<nOB;i++){
    EigBO(i) = (*this->singleSlater_->epsB())(i);
    cout << "The " << (i+1) << " eigenvalue in Occupied Beta is: " << EigBO(i) << endl;
  }
  for (auto j=0;j<nVB;j++){
    EigBV(j) = (*this->singleSlater_->epsB())((j+nOB));
    cout << "The " << (j+1) << " eigenvalue in Virtual Beta is: " << EigBV(j) << endl;
  }


  ia = 0;
  for (auto a=0;a<nVA;a++)
  for (auto i=0;i<nOA;i++)
  {
    jb =0;
    for (auto b=0;b<nVA;b++)
    for (auto j=0;j<nOA;j++)
    {
      Aud(ia,jb) = 0.0;
      if ((a==b)&&(i==j)){
	Aud(ia,jb) = EigAV(a)-EigAO(i);
      }
      Aud(ia,jb) = Aud(ia,jb) + DajibAA(a,j,i,b);
      Bud(ia,jb) = DabijAA(a,b,i,j);
      //Aod(ia,jb) = dajib(a,j,i,b);
      //Bod(ia,jb) = dabij(a,b,i,j);
      jb = jb+1;
    }
    ia = ia+1;
  }

  ia = 0;
  for (auto a=0;a<nVB;a++)
  for (auto i=0;i<nOB;i++)
  {
    jb =0;
    for (auto b=0;b<nVB;b++)
    for (auto j=0;j<nOB;j++)
    {
      Add(ia,jb) = 0.0;
      if ((a==b)&&(i==j)){
	Add(ia,jb) = EigBV(a)-EigBO(i);
      }
      Add(ia,jb) = Add(ia,jb) + DajibBB(a,j,i,b);
      Bdd(ia,jb) = DabijBB(a,b,i,j);
      //Aod(ia,jb) = dajib(a,j,i,b);
      //Bod(ia,jb) = dabij(a,b,i,j);
      jb = jb+1;
    }
    ia = ia+1;
  }

  ia = 0;
  for (auto a=0;a<nVA;a++)
  for (auto i=0;i<nOA;i++)
  {
    jb =0;
    for (auto b=0;b<nVB;b++)
    for (auto j=0;j<nOB;j++)
    {
      Auod(ia,jb) = dajibAB(a,j,i,b);
      Buod(ia,jb) = dabijAB(a,b,i,j);
      //Aod(ia,jb) = dajib(a,j,i,b);
      //Bod(ia,jb) = dabij(a,b,i,j);
      jb = jb+1;
    }
    ia = ia+1;
  }

  ia = 0;
  for (auto a=0;a<nVB;a++)
  for (auto i=0;i<nOB;i++)
  {
    jb =0;
    for (auto b=0;b<nVA;b++)
    for (auto j=0;j<nOA;j++)
    {
      Adod(ia,jb) = dajibBA(a,j,i,b);
      Bdod(ia,jb) = dabijBA(a,b,i,j);
      //Aod(ia,jb) = dajib(a,j,i,b);
      //Bod(ia,jb) = dabij(a,b,i,j);
      jb = jb+1;
    }
    ia = ia+1;
  }


  A.block(0,0,nOVA,nOVA) = Aud;
  A.block(nOVA,nOVA,nOVB,nOVB) = Add;
  A.block(0,nOVA,nOVA,nOVB) = Auod;
  A.block(nOVA,0,nOVB,nOVA) = Adod;
  B.block(0,0,nOVA,nOVA) = Bud;
  B.block(nOVA,nOVA,nOVB,nOVB) = Bdd;
  B.block(0,nOVA,nOVA,nOVB) = Auod;
  B.block(nOVA,0,nOVB,nOVA) = Adod;
  
  // Build the ABBA matrix

  ABBA.block(0,0,nOVA+nOVB,nOVA+nOVB) = A;
  ABBA.block(nOVA+nOVB,nOVA+nOVB,nOVA+nOVB,nOVA+nOVB) = -A;
  ABBA.block(0,nOVA+nOVB,nOVA+nOVB,nOVA+nOVB) = B;
  ABBA.block(nOVA+nOVB,0,nOVA+nOVB,nOVA+nOVB) = -B;

  // CIS routine
  // Diagonalize the A matrix
  Eigen::SelfAdjointEigenSolver<RealMatrix> CIS;
  CIS.compute(A);
  CIS.eigenvalues();
  CIS.eigenvectors();
  // Print the CIS Excitation Energies
//  for (auto i=0;i<2*nOV;i++){
//    cout << "The " << (i+1) << " CIS Exicitation Energy is: "
//         << (CIS.eigenvalues())(i) << endl;
//  }
  RealMatrix XMO = CIS.eigenvectors().col(1);
  formRM2(XMO);
  cout << "True AX" << endl;
  cout << A*CIS.eigenvectors() << endl;
  cout << "***********" << endl;
  cout << "Here " << endl;
  cout << "Output Energy" << endl;
  *this->CISEnergy_   = CIS.eigenvalues();
  cout << *this->CISEnergy_ << endl;
  cout << "Output Transition Density" << endl;
  *this->CISTransDen_ = CIS.eigenvectors();
  cout << *this->CISTransDen_ << endl;
  cout << "Here " << endl;
  TransDipole();
  double Oscstr = OscStrength();
  cout << "f = " << Oscstr << endl;

  // LR TDHF routine
  //Eigen::EigenSolver<RealMatrix> TD;
  //TD.compute(ABBA);
  //TD.eigenvalues();
  //TD.eigenvectors();

  // Print the LR-TDHF Excitation Energies
  //for (auto i=0;i<4*nOV;i++){
  //  cout << "The " << (i+1) << " LR-TDHF Exicitation Energy is: "
  //       << (TD.eigenvalues())(i) << endl;
  //}

}

void SDResponse::DavidsonCIS(){
  int nOVA = this->singleSlater_->nOVA();
  int nOVB = this->singleSlater_->nOVB();
  RealMatrix PDiag(nOVA+nOVB,1);
  PDiag = ReturnDiag();
  RealMatrix GVec(nOVA+nOVB,nOVA+nOVB);
  GVec = Guess(PDiag);
  RealMatrix Gpass = GVec.block(0,0,(nOVA+nOVB),3);
  cout << Gpass << endl;
  Davidson<double> davA(this,Davidson<double>::CIS,3,&Gpass,3,&PDiag);
  davA.run(this->fileio_->out);
  cout << "The lowest 3 eigenvalue solved by Davidson Algorithm:" <<endl;
  cout << *davA.eigenvalues() << endl;

}

RealMatrix SDResponse::formRM2(RealMatrix &XMO){
  int nOA = this->singleSlater_->nOccA();
  int nVA = this->singleSlater_->nVirA();
  int nOB = this->singleSlater_->nOccB();
  int nVB = this->singleSlater_->nVirB();

  Tensor<double> LocMoAO(this->nBasis_,nOA);
  Tensor<double> LocMoAV(this->nBasis_,nVA);
  Tensor<double> LocMoBO(this->nBasis_,nOB);
  Tensor<double> LocMoBV(this->nBasis_,nVB);

  for(auto ii = 0; ii < this->nBasis_; ii++) {
    for(auto jj = 0; jj < nOA; jj++) {
      LocMoAO(ii,jj) = (*this->singleSlater_->moA())(ii,jj);
    }
    for(auto kk = nOA; kk< this->nBasis_; kk++) {
      LocMoAV(ii,kk-nOA) = (*this->singleSlater_->moA())(ii,kk);
    }
  }

  for(auto ii = 0; ii < this->nBasis_; ii++) {
    for(auto jj = 0; jj < nOB; jj++) {
      LocMoBO(ii,jj) = (*this->singleSlater_->moB())(ii,jj);
    }
    for(auto kk = nOB; kk< this->nBasis_; kk++) {
      LocMoBV(ii,kk-nOB) = (*this->singleSlater_->moB())(ii,kk);
    }
  }

  int nOVA = this->singleSlater_->nOVA();
  int nOVB = this->singleSlater_->nOVB();
  RealMatrix EigAO(nOA,1);
  RealMatrix EigAV(nVA,1);
  RealMatrix EigBO(nOB,1);
  RealMatrix EigBV(nVB,1);

  for (auto i=0;i<nOA;i++){
    EigAO(i) = (*this->singleSlater_->epsA())(i);
    cout << "The " << (i+1) << " eigenvalue in Occupied is: " << EigAO(i) << endl;
  }
  for (auto j=0;j<nVA;j++){
    EigAV(j) = (*this->singleSlater_->epsA())((j+nOA));
    cout << "The " << (j+1) << " eigenvalue in Virtual is: " << EigAV(j) << endl;
  }

  for (auto i=0;i<nOB;i++){
    EigBO(i) = (*this->singleSlater_->epsB())(i);
    cout << "The " << (i+1) << " eigenvalue in Occupied is: " << EigBO(i) << endl;
  }
  for (auto j=0;j<nVB;j++){
    EigBV(j) = (*this->singleSlater_->epsB())((j+nOB));
    cout << "The " << (j+1) << " eigenvalue in Virtual is: " << EigBV(j) << endl;
  }

  enum{a,j,i,b,mu,nu,lam,sig};

  int nCol = XMO.cols();
  cout << "The dimension of input Matrix XMO:  " << (nOVA+nOVB) << " * "<< nCol << endl;
  RealMatrix AX(nOVA+nOVB,nCol);
  Tensor<double> XAAOTsr(this->nBasis_,this->nBasis_);
  Tensor<double> XBAOTsr(this->nBasis_,this->nBasis_);
  RealMatrix X(nOVA+nOVB,1);
  RealMatrix XAAO(this->nBasis_,this->nBasis_);
  RealMatrix XBAO(this->nBasis_,this->nBasis_);
  Tensor<double> IXAO1(this->nBasis_,this->nBasis_);
  Tensor<double> IIXMO1(nVA,this->nBasis_);
  Tensor<double> IXMOTsr1(nVA,nOA);
  Tensor<double> IXAO2(this->nBasis_,this->nBasis_);
  Tensor<double> IIXMO2(nVA,this->nBasis_);
  Tensor<double> IXMOTsr2(nVA,nOA);
  Tensor<double> IXAO3(this->nBasis_,this->nBasis_);
  Tensor<double> IIXMO3(nVB,this->nBasis_);
  Tensor<double> IXMOTsr3(nVB,nOB);
  Tensor<double> IXAO4(this->nBasis_,this->nBasis_);
  Tensor<double> IIXMO4(nVB,this->nBasis_);
  Tensor<double> IXMOTsr4(nVB,nOB);
  RealMatrix IXMO1(nVA,nOA);
  RealMatrix IXMO2(nVA,nOA);
  RealMatrix IXMO3(nVB,nOB);
  RealMatrix IXMO4(nVB,nOB);
  RealMatrix IXMOA(nVA,nOA);
  RealMatrix IXMOB(nVB,nOB);
  // Build <mn||ls> and <mn|ls>
  Tensor<double> Dmnls(this->nBasis_,this->nBasis_,this->nBasis_,this->nBasis_);
  Tensor<double> dmnls(this->nBasis_,this->nBasis_,this->nBasis_,this->nBasis_);
  for (auto m=0;m<this->nBasis_;m++)
  for (auto n=0;n<this->nBasis_;n++)
  for (auto l=0;l<this->nBasis_;l++)
  for (auto s=0;s<this->nBasis_;s++)
  {
    Dmnls(m,n,l,s) = (*this->aoERI_)(m,l,n,s)-(*this->aoERI_)(m,s,n,l);
    dmnls(m,n,l,s) = (*this->aoERI_)(m,l,n,s);
  }

  for (auto idx=0;idx<nCol;idx++)
  {
    // Build AX by column 
    X = XMO.col(idx);
    //cout << "Contract the " << idx+1 << " column vector" <<endl;
    //cout << X << endl;
    RealMap XA(X.data(),nVA,nOA);
    //cout << "Print the XA" <<endl;
    //cout << XA << endl;
    RealMap XB(X.data()+nOVA,nVB,nOB);
    //cout << "Print the XB" <<endl;
    //cout << XB << endl;
    XAAO = this->singleSlater_->moA()->block(0,nOA,this->nBasis_,nVA)*XA*this->singleSlater_->moA()->block(0,0,this->nBasis_,nOA).transpose();
    XBAO = this->singleSlater_->moB()->block(0,nOB,this->nBasis_,nVB)*XB*this->singleSlater_->moB()->block(0,0,this->nBasis_,nOB).transpose();
    // XAAO,XBAO in tensor form
    for (auto i=0;i<this->nBasis_;i++)
    for (auto j=0;j<this->nBasis_;j++)
    {
      XAAOTsr(i,j) = XAAO(i,j);
      XBAOTsr(i,j) = XBAO(i,j);
    }
    RealMatrix IXAO1t(this->nBasis_,this->nBasis_);
    // Contract A_AAAA( <mn||ls> ) with XA
    cout << "****Test twoEContractN4" << endl;
    RealMatrix AXA(this->nBasis_,this->nBasis_);
    RealMatrix AXB(this->nBasis_,this->nBasis_);
    RealMatrix MYAXA(this->nBasis_,this->nBasis_);
    RealMatrix MYAXB(this->nBasis_,this->nBasis_);
    this->singleSlater_->aointegrals()->twoEContractN4(false, true, XAAO,AXA,XBAO,AXB);
    cout << "AXA" << endl;
    cout << AXA << endl;
    cout << "AXB" << endl;
    cout << AXB   << endl;
    contract(1.0,XAAOTsr,{sig,nu},Dmnls,{mu,nu,lam,sig},0.0,IXAO1,{mu,lam});
    for (auto a=0;a<this->nBasis_;a++)
    for (auto i=0;i<this->nBasis_;i++)
    {
      MYAXA(a,i)=IXAO1(a,i);
    }
    cout << "Show my AXAA" << endl;
    cout << MYAXA << endl;
    contract(1.0,LocMoAV,{mu,a},IXAO1,{mu,lam},0.0,IIXMO1,{a,lam});
    contract(1.0,LocMoAO,{lam,i},IIXMO1,{a,lam},0.0,IXMOTsr1,{a,i});
    for (auto a=0;a<nVA;a++)
    for (auto i=0;i<nOA;i++)
    {
      IXMO1(a,i)=IXMOTsr1(a,i);
      cout << IXMO1(a,i) << endl;
      IXMO1(a,i)= IXMO1(a,i) + XA(a,i)*(EigAV(a)-EigAO(i));
    }
    RealMatrix MYAXAB(this->nBasis_,this->nBasis_);
    // Contract A_AABB( <mn|ls> ) with XB
    contract(1.0,XBAOTsr,{sig,nu},dmnls,{mu,nu,lam,sig},0.0,IXAO2,{mu,lam});
    for (auto a=0;a<this->nBasis_;a++)
    for (auto i=0;i<this->nBasis_;i++)
    {
      MYAXAB(a,i) = IXAO2(a,i);
    }
    cout << "Show my AXAB" << endl;
    cout << MYAXAB << endl;
    contract(1.0,LocMoAV,{mu,a},IXAO2,{mu,lam},0.0,IIXMO2,{a,lam});
    contract(1.0,LocMoAO,{lam,i},IIXMO2,{a,lam},0.0,IXMOTsr2,{a,i});
    for (auto a=0;a<nVA;a++)
    for (auto i=0;i<nOA;i++)
    {
      IXMO2(a,i)=IXMOTsr2(a,i);
    }
    RealMatrix MYAXBA(this->nBasis_,this->nBasis_);
    // Contract A_BBAA( <mn|ls> ) with XA
    contract(1.0,XAAOTsr,{sig,nu},dmnls,{mu,nu,lam,sig},0.0,IXAO3,{mu,lam});
    for (auto a=0;a<this->nBasis_;a++)
    for (auto i=0;i<this->nBasis_;i++)
    {
      MYAXBA(a,i) = IXAO3(a,i);
    }
    cout << "Show my AXBA" << endl;
    cout << MYAXBA << endl;
    contract(1.0,LocMoBV,{mu,a},IXAO3,{mu,lam},0.0,IIXMO3,{a,lam});
    contract(1.0,LocMoBO,{lam,i},IIXMO3,{a,lam},0.0,IXMOTsr3,{a,i});
    for (auto a=0;a<nVB;a++)
    for (auto i=0;i<nOB;i++)
    {
      IXMO3(a,i)=IXMOTsr3(a,i);
    }
    // Contract A_BBBB( <mn||ls> ) with XB
    contract(1.0,XBAOTsr,{sig,nu},Dmnls,{mu,nu,lam,sig},0.0,IXAO4,{mu,lam});
    for (auto a=0;a<this->nBasis_;a++)
    for (auto i=0;i<this->nBasis_;i++)
    {
      MYAXB(a,i)=IXAO4(a,i);
    }
    cout << "Show my AXBB" << endl;
    cout << MYAXB << endl;
    contract(1.0,LocMoBV,{mu,a},IXAO4,{mu,lam},0.0,IIXMO4,{a,lam});
    contract(1.0,LocMoBO,{lam,i},IIXMO4,{a,lam},0.0,IXMOTsr4,{a,i});
    for (auto a=0;a<nVB;a++)
    for (auto i=0;i<nOB;i++)
    {
      IXMO4(a,i) = IXMOTsr4(a,i);
      cout << IXMO4(a,i) << endl;
      IXMO4(a,i) = IXMO4(a,i) + XB(a,i)*(EigBV(a)-EigBO(i));
    }
    cout << "MyAXA" << endl;
    cout << MYAXA+MYAXAB << endl;
    cout << "MyAXB" << endl;
    cout << MYAXB+MYAXBA << endl;
  
    // Get the Final AX matrix
    IXMOA = IXMO1+IXMO2;
    IXMOB = IXMO3+IXMO4;
 
    // Print AX(i)
    //cout << "Print AX(i) (direct build)" <<endl;
    //cout << IXMOA <<endl;
    //cout << IXMOB << endl;
    RealMap IXMOAV(IXMOA.data(),nOVA,1);
    RealMap IXMOBV(IXMOB.data(),nOVB,1);
    AX.block(0,idx,nOVA,1) = IXMOAV;
    AX.block(nOVA,idx,nOVB,1) = IXMOBV; 
  }
  cout << "AX"  <<endl;
  cout << AX    <<endl;
  return AX;
}

RealMatrix SDResponse::ReturnDiag(){
  int nOA = this->singleSlater_->nOccA();
  int nVA = this->singleSlater_->nVirA();
  int nOB = this->singleSlater_->nOccB();
  int nVB = this->singleSlater_->nVirB();
  int nOVA = this->singleSlater_->nOVA();
  int nOVB = this->singleSlater_->nOVB();
  RealMatrix EigAO(nOA,1);
  RealMatrix EigAV(nVA,1);
  RealMatrix EigBO(nOB,1);
  RealMatrix EigBV(nVB,1);
  for (auto i=0;i<nOA;i++){
    EigAO(i) = (*this->singleSlater_->epsA())(i);
  }
  for (auto j=0;j<nVA;j++){
    EigAV(j) = (*this->singleSlater_->epsA())((j+nOA));
  }
  for (auto i=0;i<nOB;i++){
    EigBO(i) = (*this->singleSlater_->epsB())(i);
  }
  for (auto j=0;j<nVB;j++){
    EigBV(j) = (*this->singleSlater_->epsB())((j+nOB));
  }

  RealMatrix PDiag(nOVA+nOVB,1);
  for (auto a=0;a<nVA;a++)
  for (auto i=0;i<nOA;i++)
  {
    PDiag(a*nOA+i,0) = EigAV(a)-EigAO(i);
  }

  for (auto a=0;a<nVB;a++)
  for (auto i=0;i<nOB;i++)
  {
    PDiag(a*nOB+i+nOVA,0) = EigBV(a)-EigBO(i);
  }
  return PDiag;
}

RealMatrix SDResponse::Guess(RealMatrix &PDiag){
  int nOA = this->singleSlater_->nOccA();
  int nOB = this->singleSlater_->nOccB();
  int nVA = this->singleSlater_->nVirA();
  int nVB = this->singleSlater_->nVirB();
  int nOVA = this->singleSlater_->nOVA();
  int nOVB = this->singleSlater_->nOVB();
  double temp1;
  double temp2;                                                                                                    
  bool isSorted;
  RealMatrix Permt(nOVA+nOVB,1);                                                                                       
  for (auto i=0;i<(nOVA+nOVB);i++)                                                                                       
  { 
    Permt(i,0) = i;
  }   
  for (auto i=0;i<(nOVA+nOVB-1);i++)
  {   
    isSorted = true;
    for (auto j=0;j<(nOVA+nOVB-1);j++)                                                                                   
    { 
      if (PDiag(j,0)>PDiag(j+1,0))                                                                                 
      {
        isSorted=false;
        temp1=PDiag(j,0);
        temp2=Permt(j,0);                                                                                          
        PDiag(j,0)=PDiag(j+1,0);
        Permt(j,0)=Permt(j+1,0);                                                                                   
        PDiag(j+1,0)=temp1;
        Permt(j+1,0)=temp2;                                                                                        
      } 
    }                                                                                                              
  }
  RealMatrix GVec(nOVA+nOVB,nOVA+nOVB);                                                                                    
  GVec.Zero(nOVA+nOVB,nOVA+nOVB);
  for (auto i=0;i<(nOVA+nOVB);i++)                                                                                       
  {   
    GVec(Permt(i,0),i) = 1.0;                                                                                      
  }   
  return GVec;
}

void SDResponse::TransDipole(){
  int nOA   = this->singleSlater_->nOccA();
  int nVA   = this->singleSlater_->nVirA();
  int nOVA  = this->singleSlater_->nOVA();
  int nOB   = this->singleSlater_->nOccB();
  int nVB   = this->singleSlater_->nVirB();
  int nOVB  = this->singleSlater_->nOVB();

  int NBSq = this->nBasis_*this->nBasis_;
  double transdipole;
  double Tmax=0.0;
  int order=0;
  (*this->CISTransDen_).transposeInPlace();
  RealMap TransDen(this->CISTransDen_->data()+nOVA+nOVB,(nOVA+nOVB),1);
  cout << "The transition density matrix we use: " << endl;
  cout << TransDen << endl;
  Tmax = TransDen(0,0);
  cout << "Tmax = " << Tmax << endl;
  for (auto i=0;i<(nOVA+nOVB);i++)
  {
    if (fabs(TransDen(i,0))>fabs(Tmax))
    {
      Tmax = TransDen(i,0);
      order = i;
    }
  }
  cout << "Maximum contribution is " <<endl;
  for (auto i=0,IOff=0;i<3;i++,IOff+=NBSq)
  { 
    transdipole = 0.0;
    RealMap Dipole(&this->elecDipole_->storage()[IOff],this->nBasis_,this->nBasis_);
    RealMap TDenMO1(TransDen.data(),nVA,nOA);
    RealMatrix TDenAO = this->singleSlater_->moA()->block(0,nOA,this->nBasis_,nVA)*TDenMO1*this->singleSlater_->moA()->block(0,0,this->nBasis_,nOA).transpose();
    TDenAO = TDenAO.transpose()*Dipole;
    transdipole += TDenAO.trace();
    RealMap TDenMO2(TransDen.data()+nOVA,nVB,nOB);
    TDenAO = this->singleSlater_->moB()->block(0,nOB,this->nBasis_,nVB)*TDenMO2*this->singleSlater_->moB()->block(0,0,this->nBasis_,nOB).transpose();
    TDenAO = TDenAO.transpose()*Dipole;
    transdipole += TDenAO.trace();
    (*this->TransDipole_)(0,i) = transdipole;
  }
  
}

double SDResponse::OscStrength(){
  double Oscstr = 0.0;
  double Omega  = (*this->CISEnergy_)(6);
  for (auto i=0;i<3;i++)
  {
    Oscstr += (2.0/3.0)*Omega*pow((*this->TransDipole_)(0,i),2);
  }
  return Oscstr;
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

