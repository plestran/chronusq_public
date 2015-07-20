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
#include <quasinewton.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
using ChronusQ::SDResponse;
using ChronusQ::SingleSlater;
using ChronusQ::QuasiNewton; 
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
  this->RHF_            = singleSlater->RHF();
/*
  if (this->RHF_)
  {
    cout << "The calculation is based on RHF" << endl;
  }
  else
  {
    cout << "The calculation is based on UHF" << endl;
  }
  cout << "Allocate Memory for CISEnergy_ and CISTransDen_" << endl;
*/
  int nOVA = singleSlater_->nOVA();
  int nOVB = singleSlater_->nOVB();
/*
  this->CISEnergy_ = std::unique_ptr<RealMatrix>(new RealMatrix(nOVA+nOVB,1));
  this->CISTransDen_ = std::unique_ptr<RealMatrix>(new RealMatrix(nOVA+nOVB,nOVA+nOVB));
  this->TransDipole_ = std::unique_ptr<RealMatrix>(new RealMatrix(1,3));
*/
//dbwys
  this->haveDag_ = false;
  this->nOA_ = this->singleSlater_->nOccA();
  this->nOB_ = this->singleSlater_->nOccB();
  this->nVA_ = this->singleSlater_->nVirA();
  this->nVB_ = this->singleSlater_->nVirB();
  this->nOAVA_ = this->nOA_*this->nVA_;
  this->nOBVB_ = this->nOB_*this->nVB_;
  this->nOAVB_ = this->nOA_*this->nVB_;
  this->nOBVA_ = this->nOB_*this->nVA_;
  this->setNSek(this->controls_->SDNSek);
  this->setMeth(this->controls_->SDMethod);
  this->omega_ = std::unique_ptr<VectorXd>(new VectorXd(this->nSek_));
  this->transDen_ = std::unique_ptr<RealCMMatrix>(new RealCMMatrix(this->nSingleDim_,this->nSek_));
  this->oscStrength_ = std::unique_ptr<RealMatrix>(new RealMatrix(this->nSek_+1,this->nSek_+1));
  this->transDipole_ = std::unique_ptr<RealTensor3d>(new RealTensor3d(this->nSek_+1,this->nSek_+1,3));
//dbwye
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
  this->fileio_->out << bannerTop << endl;
  if(this->iMeth_ == CIS)
    this->fileio_->out << "CIS";
  else if(this->iMeth_ == RPA)
    this->fileio_->out << "RPA";
  this->fileio_->out << " Diagonalization for lowest " << this->nSek_ << " eigenstates" << endl;
  this->fileio_->out << bannerMid << endl << endl;
  
  for(auto iSt = 0; iSt < this->nSek_; iSt++){
    double Omega = (*this->omega_)(iSt);
    this->fileio_->out << "Excited State " << iSt+1 << ":" << endl;
    this->fileio_->out << "  \u03C9 = " << std::setw(10) << std::setprecision(7) << std::fixed << Omega                   << " Eh   ";
    this->fileio_->out << "  \u03C9 = " << std::setw(10) << std::setprecision(7) << std::fixed << Omega*phys.eVPerHartree << " eV   ";
    this->fileio_->out << "  \u03C9 = " << std::setw(10) << std::setprecision(7) << std::fixed << Omega*phys.nmPerHartree << " nm   " << endl;
    this->fileio_->out << "  f(" << 0 << "," << iSt+1 << ") = " << (*this->oscStrength_)(0,iSt+1) << endl;
    this->printPrinciple(iSt);
  }
};
//--------------------//
//        form        //
//    < i j | a b >   //
//--------------------//
void SDResponse::formRM(){
  // Get info from SCF and make Local copy of MO
  int nOA = this->singleSlater_->nOccA();
  int nVA = this->singleSlater_->nVirA();
  int nOVA = this->singleSlater_->nOVA();
  int nOB = this->singleSlater_->nOccB();
  int nVB = this->singleSlater_->nVirB();
  int nOVB = this->singleSlater_->nOVB();
//cout << "Number of Occupied Alpha: " << nOA << ", Number of Virtual Alpha" << nVA << " Number of Basis: "<<this->nBasis_<< ".\n";
//cout << "Number of Occupied Beta: " << nOB << ", Number of Virtual Beta" << nVB << " Number of Basis: "<<this->nBasis_<< ".\n";
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

  // Prepare the 2e integrals
  // Ixxxx for intermediate Sxxxx for Mulliken Notation single bar
  // Dxxxx for Dirac Notation double bar, dxxxx for Dirac Notation single bar
  Tensor<double> IanlsA(nVA,this->nBasis_,this->nBasis_,this->nBasis_); // (a nu | lam sig)
  Tensor<double> IailsA(nVA,nOA,this->nBasis_,this->nBasis_); // (a j  | lam sig)
  Tensor<double> IaijsAA(nVA,nOA,nOA,this->nBasis_); // (a j  | i   sig)
  Tensor<double> SaijbAA(nVA,nOA,nOA,nVA); // ( a j | i b )
  Tensor<double> IablsA(nVA,nVA,this->nBasis_,this->nBasis_); // ( a b | lam sig )
  Tensor<double> IabjsA(nVA,nVA,nOA,this->nBasis_); // ( a b | j sig )
  Tensor<double> SabjiAA(nVA,nVA,nOA,nOA); // ( a b | j i )
  Tensor<double> DajibAA(nVA,nOA,nOA,nVA); // <aj||ib>
  Tensor<double> dajibAB(nVA,nOB,nOA,nVB); // <aj|ib>
  Tensor<double> DabijAA(nVA,nVA,nOA,nOA); // <ab||ij>
  Tensor<double> dabijAB(nVA,nVB,nOA,nOB); // <ab|ij>
  Tensor<double> IanlsB(nVB,this->nBasis_,this->nBasis_,this->nBasis_);
  Tensor<double> IailsB(nVB,nOB,this->nBasis_,this->nBasis_);
  Tensor<double> IaijsAB(nVA,nOA,nOB,this->nBasis_);
  Tensor<double> IaijsBA(nVB,nOB,nOA,this->nBasis_);
  Tensor<double> IaijsBB(nVB,nOB,nOB,this->nBasis_);
  Tensor<double> SaijbAB(nVA,nOA,nOB,nVB);
  Tensor<double> SaijbBA(nVB,nOB,nOA,nVA);
  Tensor<double> SaijbBB(nVB,nOB,nOB,nVB);
  Tensor<double> IablsB(nVB,nVB,this->nBasis_,this->nBasis_);
  Tensor<double> IabjsB(nVB,nVB,nOB,this->nBasis_);
  Tensor<double> SabjiBB(nVB,nVB,nOB,nOB);
  Tensor<double> DajibBB(nVB,nOB,nOB,nVB);
  Tensor<double> dajibBA(nVB,nOA,nOB,nVA);
  Tensor<double> DabijBB(nVB,nVB,nOB,nOB);
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


  // (ab|ji)AAAA
  // (a nu | lam sig)
  contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsA,{a,nu,lam,sig});
  // (a b  | lam sig)
  contract(1.0,LocMoAV,{nu,b},IanlsA,{a,nu,lam,sig},0.0,IablsA,{a,b,lam,sig});
  // (a b  | j   sig)
  contract(1.0,IablsA,{a,b,lam,sig},LocMoAO,{lam,j},0.0,IabjsA,{a,b,j,sig});
  // (a b  | j   i  )
  contract(1.0,IabjsA,{a,b,j,sig},LocMoAO,{sig,i},0.0,SabjiAA,{a,b,j,i});

  // Build <aj||ib>AAAA
  for(auto a=0;a<nVA;a++) 
  for(auto j=0;j<nOA;j++) 
  for(auto i=0;i<nOA;i++) 
  for(auto b=0;b<nVA;b++) {
    DajibAA(a,j,i,b) = SaijbAA(a,i,j,b)-SabjiAA(a,b,j,i);
  }

  // Build <aj||ib>ABAB
  for(auto a=0;a<nVA;a++) 
  for(auto j=0;j<nOB;j++) 
  for(auto i=0;i<nOA;i++) 
  for(auto b=0;b<nVB;b++) {
    dajibAB(a,j,i,b) = SaijbAA(a,i,j,b);
  }
   
  //Build <ab||ij>AAAA
  for (auto a=0;a<nVA;a++)
  for (auto b=0;b<nVA;b++)
  for (auto i=0;i<nOA;i++)
  for (auto j=0;j<nOA;j++){
    DabijAA(a,b,i,j) = SaijbAA(a,i,j,b) - SaijbAA(a,j,i,b);
  }
  //Build <ab||ij> ABAB
  for (auto a=0;a<nVA;a++)
  for (auto b=0;b<nVB;b++)
  for (auto i=0;i<nOA;i++)
  for (auto j=0;j<nOB;j++){
    dabijAB(a,b,i,j) = SaijbAA(a,i,j,b);
  }

  // Build A & B matrix
  int ia,jb;
  RealMatrix ABBA(2*(nOVA+nOVB),2*(nOVA+nOVB));
  RealMatrix A(nOVA+nOVB,nOVA+nOVB);
  RealMatrix B(nOVA+nOVB,nOVA+nOVB);
  RealMatrix Aud(nOVA,nOVA);
  RealMatrix Bud(nOVA,nOVA);
  RealMatrix Auod(nOVA,nOVB);
  RealMatrix Buod(nOVA,nOVB);
  RealMatrix EigAO(nOA,1);
  RealMatrix EigAV(nVA,1);
  RealMatrix Add(nOVB,nOVB);
  RealMatrix Bdd(nOVB,nOVB);
  RealMatrix Adod(nOVB,nOVA);
  RealMatrix Bdod(nOVB,nOVA);
  RealMatrix EigBO(nOB,1);
  RealMatrix EigBV(nVB,1);
/*
  for (auto i=0;i<nOA;i++){
    EigAO(i) = (*this->singleSlater_->epsA())(i);
    cout << "The " << (i+1) << " eigenvalue in Occupied Alpha is: " << EigAO(i) << endl;
  }
  for (auto j=0;j<nVA;j++){
    EigAV(j) = (*this->singleSlater_->epsA())((j+nOA));
    cout << "The " << (j+1) << " eigenvalue in Virtual Alpha is: " << EigAV(j) << endl;
  }
*/

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
      jb = jb+1;
    }
    ia = ia+1;
  }
  if (this->RHF_)
  {
    A.block(0,0,nOVA,nOVA) = Aud;
    A.block(nOVA,nOVA,nOVB,nOVB) = Aud;
    A.block(0,nOVA,nOVA,nOVB) = Auod;
    A.block(nOVA,0,nOVB,nOVA) = Auod;
    B.block(0,0,nOVA,nOVA) = Bud;
    B.block(nOVA,nOVA,nOVB,nOVB) = Bud;
    B.block(0,nOVA,nOVA,nOVB) = Buod;
    B.block(nOVA,0,nOVB,nOVA) = Buod;
  }

  if (!this->RHF_)
  {
    for(auto ii = 0; ii < this->nBasis_; ii++) {
      for(auto jj = 0; jj < nOB; jj++) {
        LocMoBO(ii,jj) = (*this->singleSlater_->moB())(ii,jj);
      }
      for(auto kk = nOB; kk< this->nBasis_; kk++) {
        LocMoBV(ii,kk-nOB) = (*this->singleSlater_->moB())(ii,kk);
      }
    }
    // (ai|jb)AABB
    // (a nu | lam sig)
    contract(1.0,LocMoAV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsA,{a,nu,lam,sig});
    // (a i  | lam sig)
    contract(1.0,LocMoAO,{nu,i},IanlsA,{a,nu,lam,sig},0.0,IailsA,{a,i,lam,sig});
    // (a i  | j   sig)
    contract(1.0,IailsA,{a,i,lam,sig},LocMoBO,{lam,j},0.0,IaijsAB,{a,i,j,sig});
    // (a i  | j   b  )
    contract(1.0,IaijsAB,{a,i,j,sig},LocMoBV,{sig,b},0.0,SaijbAB,{a,i,j,b});
    // (ai|jb)BBAA
    // (a nu | lam sig)
    contract(1.0,LocMoBV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsB,{a,nu,lam,sig});
    // (a i  | lam sig)
    contract(1.0,LocMoBO,{nu,i},IanlsB,{a,nu,lam,sig},0.0,IailsB,{a,i,lam,sig});
    // (a i  | j   sig)
    contract(1.0,IailsB,{a,i,lam,sig},LocMoAO,{lam,j},0.0,IaijsBA,{a,i,j,sig});
    // (a i  | j   b  )
    contract(1.0,IaijsBA,{a,i,j,sig},LocMoAV,{sig,b},0.0,SaijbBA,{a,i,j,b});
    // (ai|jb)BBBB
    // (a nu | lam sig)
    contract(1.0,LocMoBV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsB,{a,nu,lam,sig});
    // (a i  | lam sig)
    contract(1.0,LocMoBO,{nu,i},IanlsB,{a,nu,lam,sig},0.0,IailsB,{a,i,lam,sig});
    // (a i  | j   sig)
    contract(1.0,IailsB,{a,i,lam,sig},LocMoBO,{lam,j},0.0,IaijsBB,{a,i,j,sig});
    // (a i  | j   b  )
    contract(1.0,IaijsBB,{a,i,j,sig},LocMoBV,{sig,b},0.0,SaijbBB,{a,i,j,b});
    // (ab|ji)BBBB
    // (a nu | lam sig)
    contract(1.0,LocMoBV,{mu,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsB,{a,nu,lam,sig});
    // (a b  | lam sig)
    contract(1.0,LocMoBV,{nu,b},IanlsB,{a,nu,lam,sig},0.0,IablsB,{a,b,lam,sig});
    // (a b  | j   sig)
    contract(1.0,IablsB,{a,b,lam,sig},LocMoBO,{lam,j},0.0,IabjsB,{a,b,j,sig});
    // (a b  | j   i  )
    contract(1.0,IabjsB,{a,b,j,sig},LocMoBO,{sig,i},0.0,SabjiBB,{a,b,j,i});

    // Build <aj||ib>ABAB
    for(auto a=0;a<nVA;a++)
    for(auto j=0;j<nOB;j++)
    for(auto i=0;i<nOA;i++)
    for(auto b=0;b<nVB;b++) {
      dajibAB(a,j,i,b) = SaijbAB(a,i,j,b);
    }
    // Build <aj||ib>BABA
    for(auto a=0;a<nVB;a++)
    for(auto j=0;j<nOA;j++)
    for(auto i=0;i<nOB;i++)
    for(auto b=0;b<nVA;b++) {
      dajibBA(a,j,i,b) = SaijbBA(a,i,j,b);
    }
    // Build <aj||ib>BBBB
    for(auto a=0;a<nVB;a++)
    for(auto j=0;j<nOB;j++)
    for(auto i=0;i<nOB;i++)
    for(auto b=0;b<nVB;b++) {
      DajibBB(a,j,i,b) = SaijbBB(a,i,j,b)-SabjiBB(a,b,j,i);
    }

    //Build <ab||ij> ABAB
    for (auto a=0;a<nVA;a++)
    for (auto b=0;b<nVB;b++)
    for (auto i=0;i<nOA;i++)
    for (auto j=0;j<nOB;j++){
      dabijAB(a,b,i,j) = SaijbAB(a,i,j,b);
    }

    //Build <ab||ij> BABA
    for (auto a=0;a<nVB;a++)
    for (auto b=0;b<nVA;b++)
    for (auto i=0;i<nOB;i++)
    for (auto j=0;j<nOA;j++){
      dabijBA(a,b,i,j) = SaijbBA(a,i,j,b);
    }
    //Build <ab||ij> BBBB
    for (auto a=0;a<nVB;a++)
    for (auto b=0;b<nVB;b++)
    for (auto i=0;i<nOB;i++)
    for (auto j=0;j<nOB;j++){
      DabijBB(a,b,i,j) = SaijbBB(a,i,j,b) - SaijbBB(a,j,i,b);
    }
/*
    for (auto i=0;i<nOB;i++){
      EigBO(i) = (*this->singleSlater_->epsB())(i);
      cout << "The " << (i+1) << " eigenvalue in Occupied Beta is: " << EigBO(i) << endl;
    }
    for (auto j=0;j<nVB;j++){
      EigBV(j) = (*this->singleSlater_->epsB())((j+nOB));
      cout << "The " << (j+1) << " eigenvalue in Virtual Beta is: " << EigBV(j) << endl;
    }
*/
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
  }

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
/*
  *this->CISEnergy_   = CIS.eigenvalues();
  *this->CISTransDen_ = CIS.eigenvectors(); 
  (*this->CISTransDen_).transposeInPlace(); 
  int nstate =5;
  cout << "Output the lowest "<< nstate << " states" << endl;
  for (auto st_rank=0;st_rank<nstate;st_rank++)
  {
    double Omega  = (*this->CISEnergy_)(st_rank);
    RealMap TransDen(this->CISTransDen_->data()+st_rank*(nOVA+nOVB),(nOVA+nOVB),1);
    TransDipole(st_rank,TransDen);
    double Oscstr = OscStrength(st_rank,Omega);
    cout << "Excitation energy is: " << " " << CIS.eigenvalues().row(st_rank) << " f = "<< Oscstr << endl << endl;
  }
*/

  // LR TDHF routine
  Eigen::EigenSolver<RealMatrix> TD;
  TD.compute(ABBA);
  TD.eigenvalues();
  TD.eigenvectors();
  cout <<"ABBA" << ABBA << endl;

  // Print the LR-TDHF Excitation Energies
  cout << "Linear response energy" << endl;
  cout << TD.eigenvalues() << endl;
  RealMatrix ReE(ABBA.rows(),1);
  ReE = TD.eigenvalues().real();
  std::sort(ReE.data(),ReE.data()+ReE.size());
  cout << ReE*phys.eVPerHartree << endl;
  cout << TD.eigenvectors().col(0) << endl;
  RealMatrix EVec = TD.eigenvectors().real();
  RealMatrix T(this->nSingleDim_,1);
  RealMatrix sigMOA(T);
  RealMatrix rhoMOA(T);
  sigMOA.setZero();
  rhoMOA.setZero();
  T = EVec.col(0);
  ABBA.block(nOVA+nOVB,nOVA+nOVB,nOVA+nOVB,nOVA+nOVB) = A;
  ABBA.block(nOVA+nOVB,0,nOVA+nOVB,nOVA+nOVB) = B;
  cout << "SIG" << endl;
  RealCMMap sMap(sigMOA.data(),this->nSingleDim_,1);
  RealCMMap tMap(T.data(),this->nSingleDim_,1);
  RealCMMap rMap(rhoMOA.data(),this->nSingleDim_,1);
  formRM3(tMap,sMap,rMap);
  cout << endl << ABBA*T-sigMOA << endl;
  T.block(this->nSingleDim_/2,0,this->nSingleDim_/2,1) = -T.block(this->nSingleDim_/2,0,this->nSingleDim_/2,1);
  cout << "RHO" << endl;
  cout << endl << T-rhoMOA << endl << endl;
  T.block(this->nSingleDim_/2,0,this->nSingleDim_/2,1) = -T.block(this->nSingleDim_/2,0,this->nSingleDim_/2,1);

  cout << sigMOA -  TD.eigenvalues()(0).real()*rhoMOA<< endl;
}

void SDResponse::IterativeRPA(){
/*
  int nOVA = this->singleSlater_->nOVA();
  int nOVB = this->singleSlater_->nOVB();
*/
//int nstate =3;
/*
  RealMatrix PDiag(nOVA+nOVB,1);
  PDiag = ReturnDiag();
  RealMatrix GVec(nOVA+nOVB,nOVA+nOVB);
  GVec = Guess(PDiag);
  RealMatrix Gpass = GVec.block(0,0,(nOVA+nOVB),3);
  cout << Gpass << endl;
//CErr();
*/
  this->formGuess();
  QuasiNewton<double> davA(this);
  davA.run(this->fileio_->out);
//davA.run(this->fileio_->out);
//this->fileio_->out << "The lowest " << this->nSek_ << " eigenstates solved by Davidson Algorithm:" <<endl;
  this->formTransDipole();
  this->formOscStrength();
/*
  this->fileio_->out << bannerTop << endl << "CIS Diagonalization for lowest " << this->nSek_ << " eigenstates" << endl;
  this->fileio_->out << bannerMid << endl;
  for (auto st_rank = 0;st_rank < this->nSek_; st_rank++)
  {
    double Omega = (*this->omega_)(st_rank);
    this->fileio_->out << "Excited State " << st_rank+1 << ":";
    this->fileio_->out << endl << "  \u03C9 = " << std::setw(15) << std::setprecision(7) << std::fixed << Omega                   << " Eh";
    this->fileio_->out << endl << "  \u03C9 = " << std::setw(15) << std::setprecision(7) << std::fixed << Omega*phys.eVPerHartree << " eV";
    this->fileio_->out << endl << "  \u03C9 = " << std::setw(15) << std::setprecision(7) << std::fixed << Omega*phys.nmPerHartree << " nm" << endl;
*/
/*
    RealMap TransDen(DavEvector.data()+st_rank*(this->nSingleDim_),(this->nSingleDim_),1);
    TransDipole(st_rank,TransDen);
    double Oscstr = OscStrength(st_rank,Omega);
    this->fileio_->out << "Excitation energy is: " << " " << Omega*phys.eVPerHartree << " eV   f = "<< Oscstr << endl << endl;
*/
//}
  this->printExcitedStateEnergies();

}

RealMatrix SDResponse::formRM2(RealMatrix &XMO){
  RealMatrix AX(this->nSingleDim_,XMO.cols());
  RealMatrix X(this->nSingleDim_,1);
  RealMatrix XAAO(this->nBasis_,this->nBasis_);
  RealMatrix XBAO(this->nBasis_,this->nBasis_);
  RealMatrix AXA(this->nVA_,this->nOA_);
  RealMatrix AXB(this->nVB_,this->nOB_);
  RealMatrix IXA(this->nBasis_,this->nBasis_);
  RealMatrix IXB(this->nBasis_,this->nBasis_);

  // Build AX by column 
  for (auto idx = 0; idx < XMO.cols(); idx++)
  {

    X = XMO.col(idx);
    RealMap XA(X.data(),this->nVA_,this->nOA_);
    RealMap XB(X.data()+this->nOAVA_,this->nVB_,this->nOB_);
    /*
     *  XAO(s) = Cv(s) * XMO(s) * Co(s)**H
     *
     *  XAO(s) - s-spin block of X in the AO basis
     *  XMO(s) - s-spin block of X in the MO basis
     *  Cv(s)  - s-spin block of the virtual block of the MO coefficients
     *  Co(s)  - s-spin block of the occupied block of the MO coefficients
     *  H      - Adjoint
     */ 
    XAAO = 
      this->singleSlater_->moA()->block(0,this->nOA_,this->nBasis_,this->nVA_)*
      XA*
      this->singleSlater_->moA()->block(0,0,this->nBasis_,this->nOA_).adjoint();
    if (this->RHF_)
    {
      XBAO = 
        this->singleSlater_->moA()->block(0,this->nOB_,this->nBasis_,this->nVB_)*
        XB*
        this->singleSlater_->moA()->block(0,0,this->nBasis_,this->nOB_).adjoint();
    }
    else
    {
      XBAO = 
        this->singleSlater_->moB()->block(0,this->nOB_,this->nBasis_,this->nVB_)*
        XB*
        this->singleSlater_->moB()->block(0,0,this->nBasis_,this->nOB_).adjoint();
    }


    AXA.setZero();  
    AXB.setZero();
    IXA.setZero();
    IXB.setZero();

    /*
     *  IXAO(s)_{i,j} = [ (ij,s|kl,s') + delta_{s,s'}*(il,s|kj,s) ] * XAO(s')_{l,k}
     */ 
    if(this->controls_->directTwoE && !this->controls_->doDF)
      this->singleSlater_->aointegrals()->twoEContractDirect(false,false,false,XAAO,IXA,XBAO,IXB);
    else
      this->singleSlater_->aointegrals()->twoEContractN4(false, true, XAAO,IXA,XBAO,IXB);
    

    /*
     *  AX(s)   += IXMO(s)
     *  IXMO(s) =  Cv(s)**H * IXAO(s) * Co(s)
     */  
    AXA = 
      this->singleSlater_->moA()->block(0,this->nOA_,this->nBasis_,this->nVA_).adjoint()*
      IXA*
      this->singleSlater_->moA()->block(0,0,this->nBasis_,this->nOA_);
    if (this->RHF_)
    {
      AXB = 
        this->singleSlater_->moA()->block(0,this->nOB_,this->nBasis_,this->nVB_).adjoint()*
        IXB*
        this->singleSlater_->moA()->block(0,0,this->nBasis_,this->nOB_);
    }
    else
    {
      AXB = 
        this->singleSlater_->moB()->block(0,this->nOB_,this->nBasis_,this->nVB_).adjoint()*
        IXB*
        this->singleSlater_->moB()->block(0,0,this->nBasis_,this->nOB_);
    }

/*  Old inefficient way of adding in the eigenenergie differences...
 *  Keep around just in case
    for (auto a = 0, ia = 0; a < this->nVA_; a++)
    for (auto i = 0; i  < this->nOA_;  i++, ia++)
    {
      AXA(a,i) += XA(a,i) * (*this->rmDiag_)(ia,0);
    }
    for (auto a = 0, ia = this->nOAVA_; a < this->nVB_; a++)
    for (auto i = 0; i  < this->nOB_;   i++,           ia++)
    {
      AXB(a,i) += XB(a,i) * (*this->rmDiag_)(ia,0);
    }
*/
    
    RealMap AXu(AXA.data(),this->nOAVA_,1);
    RealMap AXd(AXB.data(),this->nOBVB_,1);
    RealMap  Xu(XA.data(),this->nOAVA_,1);
    RealMap  Xd(XB.data(),this->nOBVB_,1);

   
    /*
     *  AX(s)_{a,i} += [Eps(s)_{a} - Eps(s)_{i}] * XMO(s)_{a,i}
     *
     *  Eps(s)_{a/i} - s-spin virtual / occupied eigenenergies (cannonical)
     */ 
    AXu += Xu.cwiseProduct(this->rmDiag_->block(0,0,this->nOAVA_,1));
    AXd += Xd.cwiseProduct(this->rmDiag_->block(this->nOAVA_,0,this->nOBVB_,1));
    
    AX.block(0,idx,this->nOAVA_,1) = AXu;
    AX.block(this->nOAVA_,idx,this->nOBVB_,1) = AXd; 
  }
  return AX;
}
//dbwys
void SDResponse::formRM3(RealCMMap &XMO, RealCMMap &Sigma, RealCMMap &Rho){
/*
 *  Forms sigma (and possibly rho) for the linear transfomation of E^(2)
 *  (and possibly S^(2) ) onto trial vectors (or any vector in general)
 *
 *  Adapted from Helgaker, et al. JCP 113, 8908 (2000)
 *
 *  DBWY (2015)
 */
  std::vector<RealMatrix> CommA(XMO.cols(),RealMatrix::Zero(this->nBasis_,this->nBasis_));
  std::vector<RealMatrix> CommB(XMO.cols(),RealMatrix::Zero(this->nBasis_,this->nBasis_));
  std::vector<RealMatrix> GCommA(XMO.cols(),RealMatrix::Zero(this->nBasis_,this->nBasis_));
  std::vector<RealMatrix> GCommB(XMO.cols(),RealMatrix::Zero(this->nBasis_,this->nBasis_));
  RealMatrix XAAO(this->nBasis_,this->nBasis_);
  RealMatrix XBAO(this->nBasis_,this->nBasis_);
  RealMatrix SigAOA(this->nBasis_,this->nBasis_);
  RealMatrix SigAOB(this->nBasis_,this->nBasis_);
  RealMatrix RhoAOA, RhoAOB;
  if(this->iMeth_ == RPA){
    RhoAOA = RealMatrix::Zero(this->nBasis_,this->nBasis_);
    RhoAOB = RealMatrix::Zero(this->nBasis_,this->nBasis_);
  }
  RealMatrix SDA,SDB;
  SDA = (*this->singleSlater_->aointegrals()->overlap_) * (*this->singleSlater_->densityA());
  if(!this->RHF_) SDB = (*this->singleSlater_->aointegrals()->overlap_) * (*this->singleSlater_->densityB());

  double fact = 1.0;
  if(this->RHF_) fact = 0.5;
  int iOff = this->nOAVA_ + this->nOBVB_;

  // Build Sigma by column 
  for (auto idx = 0; idx < XMO.cols(); idx++){
    /*
     *  XAO(s) = C(s) * XMO(s) * C(s)**H
     *
     *  XAO(s) - s-spin block of X in the AO basis
     *  XMO(s) - s-spin block of X in the MO basis
     *  C(s)   - s-spin block of the MO coefficients
     *  H      - Adjoint
     */ 
    RealVecMap X(XMO.data()+idx*this->nSingleDim_,this->nSingleDim_);
    this->formAOTDen(X,XAAO,XBAO);

    CommA[idx] =  fact * XAAO * SDA; 
    CommA[idx] += fact * SDA.adjoint() * XAAO;
    if(this->RHF_){ 
      CommB[idx] =  fact * XBAO * SDA;
      CommB[idx] += fact * SDA.adjoint() * XBAO;
    } else {
      CommB[idx] =  fact * XBAO * SDB;
      CommB[idx] += fact * SDB.adjoint() * XBAO;
    }
  }

  this->singleSlater_->aointegrals()->multTwoEContractDirect(XMO.cols(),false, false, CommA,GCommA,CommB,GCommB);


  for(auto idx = 0; idx < XMO.cols(); idx++){
    SigAOA =  (*this->singleSlater_->fockA()) * CommA[idx] * (*this->singleSlater_->aointegrals()->overlap_);
    SigAOA -= (*this->singleSlater_->aointegrals()->overlap_) * CommA[idx] * (*this->singleSlater_->fockA());
    SigAOA += fact * GCommA[idx] * SDA.adjoint();
    SigAOA -= fact * SDA * GCommA[idx];

    if(this->RHF_) {
      SigAOB =  (*this->singleSlater_->fockA()) * CommB[idx] * (*this->singleSlater_->aointegrals()->overlap_);
      SigAOB -= (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * (*this->singleSlater_->fockA());
      SigAOB += fact * GCommB[idx] * SDA.adjoint();
      SigAOB -= fact * SDA * GCommB[idx];
    } else {
      SigAOB =  (*this->singleSlater_->fockB()) * CommB[idx] * (*this->singleSlater_->aointegrals()->overlap_);
      SigAOB -= (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * (*this->singleSlater_->fockB());
      SigAOB += fact * GCommB[idx] * SDB.adjoint();
      SigAOB -= fact * SDB * GCommB[idx];
    }

    RealVecMap SVec(Sigma.data()+idx*this->nSingleDim_,this->nSingleDim_);
    this->formMOTDen(SVec,SigAOA,SigAOB);

    if(this->iMeth_ == RPA){
      RhoAOA =  (*this->singleSlater_->aointegrals()->overlap_) * CommA[idx] * (*this->singleSlater_->aointegrals()->overlap_);
      RhoAOB =  (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * (*this->singleSlater_->aointegrals()->overlap_);
      RealVecMap RVec(Rho.data()+idx*this->nSingleDim_,this->nSingleDim_);
      this->formMOTDen(RVec,RhoAOA,RhoAOB);
    }
  }
}
//dbwye

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

  if (!this->RHF_)
  {
    for (auto i=0;i<nOB;i++){
      EigBO(i) = (*this->singleSlater_->epsB())(i);
    }
    for (auto j=0;j<nVB;j++){
      EigBV(j) = (*this->singleSlater_->epsB())((j+nOB));
    }
  }

  RealMatrix PDiag(nOVA+nOVB,1);
  for (auto a=0;a<nVA;a++)
  for (auto i=0;i<nOA;i++)
  {
    PDiag(a*nOA+i,0) = EigAV(a)-EigAO(i);
    if (this->RHF_)
    {
      PDiag(a*nOB+i+nOVA,0) = EigBV(a)-EigBO(i);
    }
  }

  if (!this->RHF_)
  {
    for (auto a=0;a<nVB;a++)
    for (auto i=0;i<nOB;i++)
    {
      PDiag(a*nOB+i+nOVA,0) = EigBV(a)-EigBO(i);
    }
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

void SDResponse::TransDipole(int st_rank,RealMatrix TransDen){
  int nOA   = this->singleSlater_->nOccA();
  int nVA   = this->singleSlater_->nVirA();
  int nOVA  = this->singleSlater_->nOVA();
  int nOB   = this->singleSlater_->nOccB();
  int nVB   = this->singleSlater_->nVirB();
  int nOVB  = this->singleSlater_->nOVB();

  int NBSq = this->nBasis_*this->nBasis_;
  double transdipole;
  this->fileio_->out << "For the " << (st_rank+1) << " excited state transition, the contribution of MO are:"<< endl;
  for (auto i=0;i<nOVA;i++)
  {
    if (fabs(TransDen(i,0))>0.1)
    {
      this->fileio_->out <<(i%nOA+1) <<  "A -> " << (i/nOA+1+nOA) << "A" << "    " << std::fixed << TransDen(i,0) << endl;
    }
  }
  if (!this->RHF_)
  {
    for (auto i=0;i<nOVB;i++)
    {
      if (fabs(TransDen(i+nOVA,0))>0.1)
      {
        this->fileio_->out <<(i%nOB+1) <<  "B -> " << (i/nOB+1+nOB) << "B" << "    " << std::fixed << TransDen(i+nOVA,0) << endl;
      }
    }
  }

  for (auto i=0,IOff=0;i<3;i++,IOff+=NBSq)
  { 
    transdipole = 0.0;
    RealMap Dipole(&this->elecDipole_->storage()[IOff],this->nBasis_,this->nBasis_);
    RealMap TDenMO1(TransDen.data(),nVA,nOA);
    RealMatrix TDenAO = this->singleSlater_->moA()->block(0,nOA,this->nBasis_,nVA)*TDenMO1*this->singleSlater_->moA()->block(0,0,this->nBasis_,nOA).transpose();
    TDenAO = TDenAO.transpose()*Dipole;
    transdipole += TDenAO.trace();
    RealMap TDenMO2(TransDen.data()+nOVA,nVB,nOB);
    if (this->RHF_)
    {
      TDenAO = this->singleSlater_->moA()->block(0,nOA,this->nBasis_,nVA)*TDenMO2*this->singleSlater_->moA()->block(0,0,this->nBasis_,nOA).transpose();
    }
    else 
    {
      TDenAO = this->singleSlater_->moB()->block(0,nOB,this->nBasis_,nVB)*TDenMO2*this->singleSlater_->moB()->block(0,0,this->nBasis_,nOB).transpose();
    }
    TDenAO = TDenAO.transpose()*Dipole;
    transdipole += TDenAO.trace();
//  (*this->TransDipole_)(0,i) = transdipole;
  }
  
}

double SDResponse::OscStrength(int st_rank,double Omega){
  double Oscstr = 0.0;
  for (auto i=0;i<3;i++)
  {
//  Oscstr += (2.0/3.0)*Omega*pow((*this->TransDipole_)(0,i),2);
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


//dbwys
void SDResponse::formGuess(){
  this->checkValid();
  if(!this->haveDag_) this->getDiag();
  this->davGuess_ = 
    std::unique_ptr<RealMatrix>(
      new RealMatrix(this->nSingleDim_,this->nGuess_)
    ); 
  int nRHF = 1;
//if(this->RHF_) nRHF = 2;
  if(this->iMeth_==RPA) nRHF *= 2;
  RealMatrix dagCpy(this->nSingleDim_/nRHF,1);
  std::memcpy(dagCpy.data(),this->rmDiag_->data(),dagCpy.size()*sizeof(double));
  std::sort(dagCpy.data(),dagCpy.data()+dagCpy.size());
  std::vector<int> alreadyAdded; 
  for(auto i = 0; i < this->nGuess_; i++){
    int indx;
    for(auto k = 0; k < dagCpy.size(); k++){
      auto it = std::find(alreadyAdded.begin(),alreadyAdded.end(),k);
      if((dagCpy(i % (this->nSingleDim_/nRHF),0) == (*this->rmDiag_)(k,0)) && 
          it == alreadyAdded.end()){
        indx = k;
        alreadyAdded.push_back(indx);
        break;
      }
    }
    (*this->davGuess_)(indx,i) = 1.0;
  }
}

void SDResponse::formPerturbedGuess(double Omega, const RealCMMap & ResR, RealCMMap & QR, const RealCMMap & ResL, 
                        RealCMMap & QL){
  if(this->iMeth_ == CIS){
    for(auto i = 0; i < this->nSingleDim_; i++)
      QR(i) = -ResR(i) / ((*this->rmDiag())(i,0) - Omega);
  } else if(this->iMeth_ == RPA) {
    for(auto i = 0; i < this->nSingleDim_; i++){
      QR(i) = ResR(i) * (*this->rmDiag())(i,0);
      QL(i) = ResL(i) * (*this->rmDiag())(i,0);
    }
    for(auto i = 0; i < this->nSingleDim_/2; i++){
      QR(i) += Omega*ResL(i);
      QR(this->nSingleDim_/2 + i) -= Omega*ResL(this->nSingleDim_/2 + i);
      QL(i) += Omega*ResR(i);
      QL(this->nSingleDim_/2 + i) -= Omega*ResR(this->nSingleDim_/2 + i);
    }
    for(auto i = 0; i < this->nSingleDim_; i++){
      QR(i) = QR(i) / (std::pow((*this->rmDiag())(i,0),2.0) - std::pow(Omega,2.0));
      QL(i) = QL(i) / (std::pow((*this->rmDiag())(i,0),2.0) - std::pow(Omega,2.0));
    }
  } else {
    CErr("Perturbed guess vectors NYI for SDResponse::iMeth_ = " + std::to_string(this->iMeth_),this->fileio_->out);
  }
}

void SDResponse::checkValid(){
  if(this->nSek_ == 0)
    CErr("Specification of zero desired roots is not acceptable",
         this->fileio_->out);
  if(this->iMeth_ == 0)
    CErr("Invalid Method: SDResponse::iMeth_ = " + std::to_string(this->iMeth_),
         this->fileio_->out);
  if(this->nSingleDim_ == 0)
    CErr("Leading Dimenstion not defined for SDResponse::iMeth_ = " +
         std::to_string(this->iMeth_),this->fileio_->out);
}
void SDResponse::getDiag(){
  this->rmDiag_ = std::unique_ptr<RealCMMatrix>(new RealCMMatrix(nSingleDim_,1)); 

  for(auto aAlpha = 0; aAlpha < this->nVA_; aAlpha++)
  for(auto iAlpha = 0; iAlpha < this->nOA_; iAlpha++){
    auto iaAlpha = aAlpha*this->nOA_ + iAlpha; 
    auto eiAlpha = (*this->singleSlater_->epsA())(iAlpha);
    auto eaAlpha = (*this->singleSlater_->epsA())(aAlpha+this->nOA_);
    (*this->rmDiag_)(iaAlpha,0) = eaAlpha - eiAlpha;
    if(this->RHF_) 
      (*this->rmDiag_)(iaAlpha+this->nOAVA_,0) = eaAlpha - eiAlpha;
  }
  if(!this->RHF_){
    for(auto aBeta = 0; aBeta < this->nVB_; aBeta++)
    for(auto iBeta = 0; iBeta < this->nOB_; iBeta++){
      auto iaBeta = aBeta*this->nOB_ + iBeta + this->nOAVA_; 
      auto eiBeta = (*this->singleSlater_->epsB())(iBeta);
      auto eaBeta = (*this->singleSlater_->epsB())(aBeta+this->nOB_);
      (*this->rmDiag_)(iaBeta,0) = eaBeta - eiBeta;
    }
  }
  if(this->iMeth_ == RPA)
    this->rmDiag_->block(nSingleDim_/2,0,nSingleDim_/2,1)
      = this->rmDiag_->block(0,0,nSingleDim_/2,1);
  this->haveDag_ = true;
}
void SDResponse::initMeth(){
  if(this->nSek_ == 0) 
    CErr("Must set NSek before initializing a PSCF method",this->fileio_->out);
  if(this->iMeth_ == CIS){
    this->nSingleDim_ = this->nOAVA_ + this->nOBVB_;
  } else if(this->iMeth_ == RPA){
    this->nSingleDim_ = 2*(this->nOAVA_ + this->nOBVB_);
  } else {
    CErr("PSCF Method " + std::to_string(this->iMeth_) + " NYI",this->fileio_->out);
  }
}

void SDResponse::formAOTDen(const RealVecMap &TMOV, RealMatrix &TAOA, RealMatrix &TAOB){
  RealMatrix TMOA(this->nBasis_,this->nBasis_);
  RealMatrix TMOB(this->nBasis_,this->nBasis_);
  int iOff = this->nOAVA_ + this->nOBVB_;
  for(auto a = this->nOA_, ia = 0; a < this->nBasis_; a++)
  for(auto i = 0         ; i < this->nOA_; i++, ia++){
    TMOA(a,i) = TMOV(ia);
    if(this->iMeth_ == RPA) TMOA(i,a) = TMOV(ia+iOff);
  }
  for(auto a = this->nOB_, ia = this->nOAVA_; a < this->nBasis_; a++)
  for(auto i = 0         ; i < this->nOB_; i++, ia++){
    TMOB(a,i) = TMOV(ia);
    if(this->iMeth_ == RPA) TMOB(i,a) = TMOV(ia+iOff);
  }
  TAOA = (*this->singleSlater_->moA()) * TMOA * this->singleSlater_->moA()->adjoint();
  if(this->RHF_)
    TAOB = (*this->singleSlater_->moA()) * TMOB * this->singleSlater_->moA()->adjoint();
  else
    TAOB = (*this->singleSlater_->moB()) * TMOB * this->singleSlater_->moB()->adjoint();
}

void SDResponse::formMOTDen(RealVecMap &TMOV, const RealMatrix &TAOA, const RealMatrix &TAOB){
  RealMatrix TMOA(this->nBasis_,this->nBasis_);
  RealMatrix TMOB(this->nBasis_,this->nBasis_);
  TMOA = this->singleSlater_->moA()->adjoint() * TAOA * (*this->singleSlater_->moA());
  if(this->RHF_)
    TMOB = this->singleSlater_->moA()->adjoint() * TAOB * (*this->singleSlater_->moA());
  else 
    TMOB = this->singleSlater_->moB()->adjoint() * TAOB * (*this->singleSlater_->moB());
  int iOff = this->nOAVA_ + this->nOBVB_;
  for(auto a = this->nOA_, ia = 0; a < this->nBasis_; a++)
  for(auto i = 0         ; i < this->nOA_; i++, ia++){
    TMOV(ia) = TMOA(a,i);
    if(this->iMeth_ == RPA) TMOV(ia+iOff) = -TMOA(i,a);
  }
  for(auto a = this->nOB_, ia = this->nOAVA_; a < this->nBasis_; a++)
  for(auto i = 0         ; i < this->nOB_; i++, ia++){
    TMOV(ia) = TMOB(a,i);
    if(this->iMeth_ == RPA) TMOV(ia+iOff) = -TMOB(i,a);
  }
}

void SDResponse::formTransDipole(){
   RealMatrix TAOA(this->nBasis_,this->nBasis_);
   RealMatrix TAOB(this->nBasis_,this->nBasis_);
   auto NBSq = this->nBasis_*this->nBasis_;
   for(auto iSt = 0; iSt < this->nSek_; iSt++){
     RealVecMap TMOV(this->transDen_->data()+iSt*this->nSingleDim_,this->nSingleDim_);
     this->formAOTDen(TMOV,TAOA,TAOB);
     for(auto iXYZ = 0, iOff = 0; iXYZ < 3; iXYZ++, iOff += NBSq){
       RealMap dipole(&this->elecDipole_->storage()[iOff],this->nBasis_,this->nBasis_);
       (*this->transDipole_)(0,iSt+1,iXYZ) = (TAOA+TAOB).frobInner(dipole);
     }
   }
}

void SDResponse::formOscStrength(){
  this->oscStrength_->setZero();
  for(auto iSt  = 0; iSt  < this->nSek_; iSt++ )
  for(auto iXYZ = 0; iXYZ < 3;           iXYZ++){
    (*this->oscStrength_)(0,iSt+1) +=
      (2.0/3.0)*(*this->omega_)(iSt)*
      std::pow((*this->transDipole_)(0,iSt+1,iXYZ),2.0);
  }
}

void SDResponse::printPrinciple(int iSt){
  this->fileio_->out << "  Principle Transitions   ( tol = 0.1 )" << endl;
/*
  for(auto ia = 0; ia < this->nOAVA_; ia++){
    auto xIA = ia; auto yIA = ia + this->nSingleDim_/2;
    if(std::abs((*this->transDen_)(xIA,iSt)) > 0.1)
      this->fileio_->out << "    "
                       << (ia % this->nOA_) + 1              << "A -> "
                       << (ia / this->nOA_) + this->nOA_ + 1 << "A    "
                       << std::fixed << std::setw(10) << std::right << (*this->transDen_)(xIA,iSt)
                       << endl;
    if(this->iMeth_ == RPA)
      if(std::abs((*this->transDen_)(yIA,iSt)) > 0.1)
        this->fileio_->out << "    "
                         << (ia % this->nOA_) + 1              << "A <- "
                         << (ia / this->nOA_) + this->nOA_ + 1 << "A    "
                         << std::fixed << std::setw(10) << std::right << (*this->transDen_)(yIA,iSt)
                         << endl;
  }
  if(!this->RHF_){
    for(auto ia = this->nOAVA_; ia < this->nOAVA_ + this->nOBVB_; ia++){
      auto xIA = ia; auto yIA = ia + this->nSingleDim_/2;
 
      if(std::abs((*this->transDen_)(xIA,iSt)) > 0.1)
        this->fileio_->out << "    "
                         << (ia % this->nOA_) + 1              << "B -> "
                         << (ia / this->nOA_) + this->nOA_ + 1 << "B    "
                         << std::fixed << std::setw(10) << std::right << (*this->transDen_)(xIA,iSt)
                         << endl;
      if(this->iMeth_ == RPA)
        if(std::abs((*this->transDen_)(yIA,iSt)) > 0.1)
          this->fileio_->out << "    "
                           << (ia % this->nOA_) + 1              << "B <- "
                           << (ia / this->nOA_) + this->nOA_ + 1 << "B    "
                           << std::fixed << std::setw(10) << std::right << (*this->transDen_)(yIA,iSt)
                           << endl;
    }
  }
*/
  double printTol = 0.1;
  for(auto ia = 0; ia < this->nOAVA_; ia++){
    auto xIA_Alpha = ia; auto yIA_Alpha = ia + this->nSingleDim_/2;

    auto alphaOccOrb = (xIA_Alpha % this->nOA_) + 1;
    auto alphaVirOrb = (xIA_Alpha / this->nOA_) + this->nOA_ + 1;

    double absXIA_Alpha, absYIA_Alpha;

    absXIA_Alpha = std::abs((*this->transDen_)(xIA_Alpha,iSt));
    if(this->iMeth_ == RPA)
      absYIA_Alpha = std::abs((*this->transDen_)(yIA_Alpha,iSt));

    if(absXIA_Alpha > printTol)
        this->fileio_->out << "    "
                           << alphaOccOrb << "A -> " << alphaVirOrb << "A   "
                           << std::fixed << std::setw(10) << std::right <<
                           (*this->transDen_)(xIA_Alpha,iSt) << endl;
    if(this->iMeth_ == RPA){
      if(absYIA_Alpha > printTol)
          this->fileio_->out << "    "
                             << alphaOccOrb << "A <- " << alphaVirOrb << "A   "
                             << std::fixed << std::setw(10) << std::right <<
                             (*this->transDen_)(yIA_Alpha,iSt) << endl;
    }
  }
  for(auto ia = this->nOAVA_; ia < this->nOAVA_ + this->nOBVB_; ia++){
    auto xIA_Beta = ia; auto yIA_Beta = ia + this->nSingleDim_/2;

    auto betaOccOrb = ((xIA_Beta - this->nOAVA_) % this->nOB_) + 1;
    auto betaVirOrb = ((xIA_Beta - this->nOAVA_) / this->nOB_) + this->nOB_ + 1;

    double absXIA_Beta, absYIA_Beta;

    absXIA_Beta = std::abs((*this->transDen_)(xIA_Beta,iSt));
    if(this->iMeth_ == RPA)
      absYIA_Beta = std::abs((*this->transDen_)(yIA_Beta,iSt));

    if(absXIA_Beta > printTol)
        this->fileio_->out << "    "
                           << betaOccOrb << "B -> " << betaVirOrb << "B   "
                           << std::fixed << std::setw(10) << std::right <<
                           (*this->transDen_)(xIA_Beta,iSt) << endl;
    if(this->iMeth_ == RPA){
      if(absYIA_Beta > printTol)
          this->fileio_->out << "    "
                             << betaOccOrb << "B <- " << betaVirOrb << "B   "
                             << std::fixed << std::setw(10) << std::right <<
                             (*this->transDen_)(yIA_Beta,iSt) << endl;
    }
  }
  this->fileio_->out << bannerMid << endl << endl;
}

void SDResponse::incorePPRPA(){
  enum{a,b,c,d,i,j,k,l,mu,nu,lam,sig};
  Tensor<double> LocMoAO(this->nBasis_,this->nOA_);
  Tensor<double> LocMoAV(this->nBasis_,this->nVA_);
  Tensor<double> LocMoBO(this->nBasis_,this->nOB_);
  Tensor<double> LocMoBV(this->nBasis_,this->nVB_);  


  for(auto ii = 0; ii < this->nBasis_; ii++) {
    for(auto jj = 0; jj < this->nOA_; jj++) {
      LocMoAO(ii,jj) = (*this->singleSlater_->moA())(ii,jj);
    }
    for(auto kk = this->nOA_; kk < this->nBasis_; kk++) {
      LocMoAV(ii,kk-this->nOA_) = (*this->singleSlater_->moA())(ii,kk);
    }
  }

  Tensor<double> IanlsA(this->nVA_,this->nBasis_,this->nBasis_,this->nBasis_); // ( a(A) nu   |  lam   sig  )
  Tensor<double> IablsA(this->nVA_,this->nVA_   ,this->nBasis_,this->nBasis_); // ( a(A) b(A) |  lam   sig  )
  Tensor<double> IabcsA(this->nVA_,this->nVA_   ,this->nVA_   ,this->nBasis_); // ( a(A) b(A) |  c(A)  sig  )
  Tensor<double> SabcdA(this->nVA_,this->nVA_   ,this->nVA_   ,this->nVA_   ); // ( a(A) b(A) |  c(A)  d(A) )
  Tensor<double> dabcdA(this->nVA_,this->nVA_   ,this->nVA_   ,this->nVA_   ); // < a(A) b(A) |  c(A)  d(A) >
  Tensor<double> DabcdA(this->nVA_,this->nVA_   ,this->nVA_   ,this->nVA_   ); // < a(A) b(A) || c(A)  d(A) >

  Tensor<double> IinlsA(this->nOA_,this->nBasis_,this->nBasis_,this->nBasis_); // ( i(A) nu   |  lam   sig  )
  Tensor<double> IijlsA(this->nOA_,this->nOA_   ,this->nBasis_,this->nBasis_); // ( i(A) j(A) |  lam   sig  )
  Tensor<double> IijksA(this->nOA_,this->nOA_   ,this->nOA_   ,this->nBasis_); // ( i(A) j(A) |  k(A)  sig  )
  Tensor<double> SijklA(this->nOA_,this->nOA_   ,this->nOA_   ,this->nOA_   ); // ( i(A) j(A) |  k(A)  l(A) )
  Tensor<double> dijklA(this->nOA_,this->nOA_   ,this->nOA_   ,this->nOA_   ); // < i(A) j(A) |  k(A)  l(A) >
  Tensor<double> DijklA(this->nOA_,this->nOA_   ,this->nOA_   ,this->nOA_   ); // < i(A) j(A) || k(A)  l(A) >

  Tensor<double> IailsA(this->nVA_,this->nOA_   ,this->nBasis_,this->nBasis_); // ( a(A) i(A) |  lam   sig  )
  Tensor<double> IaibsA(this->nVA_,this->nOA_   ,this->nVA_   ,this->nBasis_); // ( a(A) i(A) |  b(A)  sig  )
  Tensor<double> SaibjA(this->nVA_,this->nOA_   ,this->nVA_   ,this->nOA_   ); // ( a(A) i(A) |  b(A)  j(A) )
  Tensor<double> dabijA(this->nVA_,this->nVA_   ,this->nOA_   ,this->nOA_   ); // < a(A) b(A) |  i(A)  j(A) >
  Tensor<double> DabijA(this->nVA_,this->nVA_   ,this->nOA_   ,this->nOA_   ); // < a(A) b(A) || i(A)  j(A) >

  Tensor<double> IabcsAB(this->nVA_,this->nVA_   ,this->nVB_   ,this->nBasis_); // ( a(A) b(A) |  c(B)  sig  )
  Tensor<double> SabcdAB(this->nVA_,this->nVA_   ,this->nVB_   ,this->nVB_   ); // ( a(A) b(A) |  c(B)  d(B) )
  Tensor<double> dabcdAB(this->nVA_,this->nVB_   ,this->nVA_   ,this->nVB_   ); // < a(A) b(B) |  c(A)  d(B) >

  Tensor<double> IijksAB(this->nOA_,this->nOA_   ,this->nOB_   ,this->nBasis_); // ( i(A) j(A) |  k(B)  sig  )
  Tensor<double> SijklAB(this->nOA_,this->nOA_   ,this->nOB_   ,this->nOB_   ); // ( i(A) j(A) |  k(B)  l(B) )
  Tensor<double> dijklAB(this->nOA_,this->nOB_   ,this->nOA_   ,this->nOB_   ); // < i(A) j(B) |  k(A)  l(B) >

  Tensor<double> IaibsAB(this->nVA_,this->nOA_   ,this->nVB_   ,this->nBasis_); // ( a(A) i(A) |  b(B)  sig  )
  Tensor<double> SaibjAB(this->nVA_,this->nOA_   ,this->nVB_   ,this->nOB_   ); // ( a(A) i(A) |  b(B)  j(B) )
  Tensor<double> dabijAB(this->nVA_,this->nVB_   ,this->nOA_   ,this->nOB_   ); // < a(A) b(B) |  i(A)  j(B) >

  // Form <AB||CD>
  
  // (ab | cd) AAAA
  contract(1.0,LocMoAV,{mu ,a},(*this->aoERI_),{mu,nu,lam,sig},0.0,IanlsA,{a,nu,lam,sig});
  contract(1.0,LocMoAV,{nu ,b},IanlsA         ,{a ,nu,lam,sig},0.0,IablsA,{a,b ,lam,sig});
  contract(1.0,LocMoAV,{lam,c},IablsA         ,{a ,b ,lam,sig},0.0,IabcsA,{a,b ,c  ,sig});
  contract(1.0,LocMoAV,{sig,d},IabcsA         ,{a ,b ,c  ,sig},0.0,SabcdA,{a,b ,c  ,d  });

  // (ab | cd) AABB
  contract(1.0,LocMoAV,{lam,c},IablsA ,{a,b,lam,sig},0.0,IabcsAB,{a,b,c,sig});
  contract(1.0,LocMoAV,{sig,d},IabcsAB,{a,b,c  ,sig},0.0,SabcdAB,{a,b,c,d  });

  // <ab | cd> = (ac | bd)
  for(auto a = 0; a < this->nVA_; a++)
  for(auto b = 0; b < this->nVA_; b++)
  for(auto c = 0; c < this->nVA_; c++)
  for(auto d = 0; d < this->nVA_; d++){
    dabcdA(a,b,c,d) = SabcdA(a,c,b,d);
    dabcdAB(a,b,c,d) = SabcdAB(a,c,b,d);
  }

  // <ab || cd> = <ab | cd> - <ab | dc>
  for(auto a = 0; a < this->nVA_; a++)
  for(auto b = 0; b < this->nVA_; b++)
  for(auto c = 0; c < this->nVA_; c++)
  for(auto d = 0; d < this->nVA_; d++){
    DabcdA(a,b,c,d) = dabcdA(a,b,c,d) - dabcdA(a,b,d,c);
  }


  // Form <IJ||KL>
  
  // (ij | kl) AAAA
  contract(1.0,LocMoAO,{mu ,i},(*this->aoERI_),{mu,nu,lam,sig},0.0,IinlsA,{i,nu,lam,sig});
  contract(1.0,LocMoAO,{nu ,j},IinlsA         ,{i ,nu,lam,sig},0.0,IijlsA,{i,j ,lam,sig});
  contract(1.0,LocMoAO,{lam,k},IijlsA         ,{i ,j ,lam,sig},0.0,IijksA,{i,j ,k  ,sig});
  contract(1.0,LocMoAO,{sig,l},IijksA         ,{i ,j ,k  ,sig},0.0,SijklA,{i,j ,k  ,l  });

  // (ij | kl) AABB
  contract(1.0,LocMoAO,{lam,k},IijlsA ,{i,j,lam,sig},0.0,IijksAB,{i,j,k,sig});
  contract(1.0,LocMoAO,{sig,l},IijksAB,{i,j,k  ,sig},0.0,SijklAB,{i,j,k,l  });

  // <ij | kl> = (ik | jl)
  for(auto i = 0; i < this->nOA_; i++)
  for(auto j = 0; j < this->nOA_; j++)
  for(auto k = 0; k < this->nOA_; k++)
  for(auto l = 0; l < this->nOA_; l++){
    dijklA(i,j,k,l) = SijklA(i,k,j,l);
    dijklAB(i,j,k,l) = SijklAB(i,k,j,l);
  }

  // <ij || kl> = <ij | kl> - <ij | lk>
  for(auto i = 0; i < this->nOA_; i++)
  for(auto j = 0; j < this->nOA_; j++)
  for(auto k = 0; k < this->nOA_; k++)
  for(auto l = 0; l < this->nOA_; l++){
    DijklA(i,j,k,l) = dijklA(i,j,k,l) - dijklA(i,j,l,k);
  }

  // Form <AB||IJ>

  // (ai | bj) AAAA
  contract(1.0,LocMoAO,{nu ,i},IanlsA,{a,nu,lam,sig},0.0,IailsA,{a,i,lam,sig});
  contract(1.0,LocMoAV,{lam,b},IailsA,{a,i ,lam,sig},0.0,IaibsA,{a,i,b  ,sig});
  contract(1.0,LocMoAO,{sig,j},IaibsA,{a,i ,b  ,sig},0.0,SaibjA,{a,i,b  ,j  });

  // (ai | bj) AABB
  contract(1.0,LocMoAV,{lam,b},IailsA ,{a,i,lam,sig},0.0,IaibsAB,{a,i,b,sig});
  contract(1.0,LocMoAO,{sig,j},IaibsAB,{a,i,b,  sig},0.0,SaibjAB,{a,i,b,j  });

  // <ab | ij> = (ai | bj)
  for(auto a = 0; a < this->nVA_; a++)
  for(auto b = 0; b < this->nVA_; b++)
  for(auto i = 0; i < this->nOA_; i++)
  for(auto j = 0; j < this->nOA_; j++){
    dabijA(a,b,i,j) = SaibjA(a,i,b,j);
    dabijAB(a,b,i,j) = SaibjAB(a,i,b,j);
  }

  // <ab || ij> = <ab | ij> - <ab | ji>
  for(auto a = 0; a < this->nVA_; a++)
  for(auto b = 0; b < this->nVA_; b++)
  for(auto i = 0; i < this->nOA_; i++)
  for(auto j = 0; j < this->nOA_; j++){
    DabijA(a,b,i,j) = dabijA(a,b,i,j) - dabijA(a,b,j,i);
  }


  double Rmu = (*this->singleSlater_->epsA())(this->nOA_-1) + (*this->singleSlater_->epsA())(this->nOA_);
  Rmu /= 2;

/*
  RealMatrix AAA(this->nVA_*(this->nVA_-1)/2,this->nVA_*(this->nVA_-1)/2);
  RealMatrix BAA(this->nVA_*(this->nVA_-1)/2,this->nOA_*(this->nOA_-1)/2);
  RealMatrix CAA(this->nOA_*(this->nOA_-1)/2,this->nOA_*(this->nOA_-1)/2);
  RealMatrix AAB(this->nVA_*this->nVA_,this->nVA_*this->nVA_);
  RealMatrix BAB(this->nVA_*this->nVA_,this->nOA_*this->nOA_);
  RealMatrix CAB(this->nOA_*this->nOA_,this->nOA_*this->nOA_);

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < a         ; b++, ab++){
    for(auto c = 0, cd = 0; c < this->nVA_; c++      )
    for(auto d = 0        ; d < c         ; d++, cd++){
      AAA(ab,cd) = DabcdA(a,b,c,d);
      if(ab == cd) 
        AAA(ab,cd) += (*this->singleSlater_->epsA())(a+this->nOA_) + 
                    (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu;
    }
  }

  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < i         ; j++, ij++){
    for(auto k = 0, kl = 0; k < this->nOA_; k++      )
    for(auto l = 0        ; l < k         ; l++, kl++){
      CAA(ij,kl) = DijklA(i,j,k,l);
      if(ij == kl) 
        CAA(ij,kl) -= (*this->singleSlater_->epsA())(i) + 
                    (*this->singleSlater_->epsA())(j) - 2*Rmu;
    }
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < a         ; b++, ab++){
    for(auto i = 0, ij = 0; i < this->nOA_; i++      )
    for(auto j = 0        ; j < i         ; j++, ij++){
      BAA(ab,ij) = DabijA(a,b,i,j);
    }
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVA_; b++, ab++){
    for(auto c = 0, cd = 0; c < this->nVA_; c++      )
    for(auto d = 0        ; d < this->nVA_; d++, cd++){
      AAB(ab,cd) = dabcdAB(a,b,c,d);
      if(ab == cd) 
        AAB(ab,cd) += (*this->singleSlater_->epsA())(a+this->nOA_) + 
                    (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu;
    }
  }

  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < this->nOA_; j++, ij++){
    for(auto k = 0, kl = 0; k < this->nOA_; k++      )
    for(auto l = 0        ; l < this->nOA_; l++, kl++){
      CAB(ij,kl) = dijklAB(i,j,k,l);
      if(ij == kl) 
        CAB(ij,kl) -= (*this->singleSlater_->epsA())(i) + 
                    (*this->singleSlater_->epsA())(j) - 2*Rmu;
    }
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVA_; b++, ab++){
    for(auto i = 0, ij = 0; i < this->nOA_; i++      )
    for(auto j = 0        ; j < this->nOA_; j++, ij++){
      BAB(ab,ij) = dabijAB(a,b,i,j);
    }
  }
  

//cout << A - A.adjoint() << endl << endl;

  Eigen::SelfAdjointEigenSolver<RealMatrix> ES;
  ES.compute(AAA);
  cout << ES.eigenvalues() << endl << endl;

  Eigen::VectorXd TDAMOV = ES.eigenvectors().col(0);
  RealMatrix TDAMO(this->nBasis_,this->nBasis_);
  RealMatrix TDAAO(this->nBasis_,this->nBasis_);
  Tensor<double> TDAAOTen(this->nBasis_,this->nBasis_);
  Tensor<double> IXAOTDATen(this->nBasis_,this->nBasis_);
  RealMatrix IXAOTDA(this->nBasis_,this->nBasis_);
  RealMatrix IXMO(this->nBasis_,this->nBasis_);

  for(auto a = 0, ab = 0; a < this->nVA_; a++)
  for(auto b = 0; b < a ; b++, ab++){
    TDAMO(this->nOA_+a,this->nOA_+b) = TDAMOV(ab);
    TDAMO(this->nOA_+b,this->nOA_+a) = -TDAMOV(ab);
  }
  cout << bannerMid << endl;
  cout << TDAMO << endl;
  cout << bannerMid << endl;

  TDAAO = (*this->singleSlater_->moA()) * TDAMO * this->singleSlater_->moA()->adjoint();
  std::memcpy(&TDAAOTen.storage()[0],TDAAO.data(),TDAAO.size()*sizeof(double));
  contract(1.0,TDAAOTen,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,IXAOTDATen,{mu,nu});
  std::memcpy(IXAOTDA.data(),&IXAOTDATen.storage()[0],IXAOTDA.size()*sizeof(double));
  IXMO = this->singleSlater_->moA()->adjoint() * IXAOTDA * (*this->singleSlater_->moA());
  for(auto a = 0, ab = 0; a < this->nVA_; a++)
  for(auto b = 0; b < a ; b++, ab++){
    cout << IXMO(this->nOA_+a,this->nOA_+b) + TDAMOV(ab)*((*this->singleSlater_->epsA())(a+this->nOA_) + (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu) << endl;
  }
  cout << endl;
  cout << AAA*TDAMOV << endl << endl;



  ES.compute(CAA);
  cout << ES.eigenvalues() << endl << endl;
  
  
  RealMatrix FullAA(this->nVA_*(this->nVA_-1)/2+this->nOA_*(this->nOA_-1)/2,this->nVA_*(this->nVA_-1)/2+this->nOA_*(this->nOA_-1)/2);

  FullAA.block(0,0,this->nVA_*(this->nVA_-1)/2,this->nVA_*(this->nVA_-1)/2) = AAA;
  FullAA.block(0,this->nVA_*(this->nVA_-1)/2,this->nVA_*(this->nVA_-1)/2,this->nOA_*(this->nOA_-1)/2) = BAA;
  FullAA.block(this->nVA_*(this->nVA_-1)/2,0,this->nOA_*(this->nOA_-1)/2,this->nVA_*(this->nVA_-1)/2) = -BAA.transpose();
  FullAA.block(this->nVA_*(this->nVA_-1)/2,this->nVA_*(this->nVA_-1)/2,this->nOA_*(this->nOA_-1)/2,this->nOA_*(this->nOA_-1)/2) = -CAA;

  Eigen::EigenSolver<RealMatrix> EA;
  EA.compute(FullAA);
  Eigen::VectorXd ER = -EA.eigenvalues().real();
  std::sort(ER.data(),ER.data()+ER.size());
  ER = -ER;
  cout << ER << endl << endl;
  RealMatrix TAA = EA.eigenvectors().col(0).real();

  cout << "PPRPA" << endl;
  ES.compute(AAB);
  cout << ES.eigenvalues() << endl << endl;

  TDAMOV = ES.eigenvectors().col(0);
  TDAMO.setZero();

  for(auto a = 0, ab = 0; a < this->nVA_; a++)
  for(auto b = 0; b  < this->nVA_ ; b++, ab++){
    TDAMO(this->nOA_+a,this->nOA_+b) = TDAMOV(ab);
  }

  TDAAO =  (*this->singleSlater_->moA()) * TDAMO * this->singleSlater_->moA()->adjoint();
//TDAAO -= (*this->singleSlater_->moA()) * TDAMO.adjoint() * this->singleSlater_->moA()->adjoint();
  std::memcpy(&TDAAOTen.storage()[0],TDAAO.data(),TDAAO.size()*sizeof(double));
  contract(1.0,TDAAOTen,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,IXAOTDATen,{mu,nu});
  std::memcpy(IXAOTDA.data(),&IXAOTDATen.storage()[0],IXAOTDA.size()*sizeof(double));
  IXMO = this->singleSlater_->moA()->adjoint() * IXAOTDA * (*this->singleSlater_->moA());
  for(auto a = 0, ab = 0; a < this->nVA_; a++)
  for(auto b = 0; b  < this->nVA_ ; b++, ab++){
    cout << IXMO(this->nOA_+a,this->nOA_+b) + TDAMOV(ab)*((*this->singleSlater_->epsA())(a+this->nOA_) + (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu) << endl;
  }
  cout << endl;
  cout << AAB*TDAMOV << endl << endl;

  ES.compute(CAB);
  cout << ES.eigenvalues() << endl << endl;
  
  RealMatrix FullAB(this->nVA_*this->nVA_+this->nOA_*this->nOA_,this->nVA_*this->nVA_+this->nOA_*this->nOA_);

  FullAB.block(0,0,this->nVA_*this->nVA_,  this->nVA_*this->nVA_) = AAB;
  FullAB.block(0,  this->nVA_*this->nVA_,  this->nVA_*this->nVA_,this->nOA_*this->nOA_) = BAB;
  FullAB.block(    this->nVA_*this->nVA_,0,this->nOA_*this->nOA_,this->nVA_*this->nVA_) = -BAB.transpose();
  FullAB.block(    this->nVA_*this->nVA_,  this->nVA_*this->nVA_,this->nOA_*this->nOA_,this->nOA_*this->nOA_) = -CAB;

  EA.compute(FullAB);
  Eigen::VectorXd ERAB = -EA.eigenvalues().real();
  std::sort(ERAB.data(),ERAB.data()+ERAB.size());
  ERAB = -ERAB;
  cout << ERAB << endl << endl;
  RealMatrix TAB = EA.eigenvectors().col(0).real();

  Eigen::VectorXd TAAX = TAA.block(0,0,this->nVA_*(this->nVA_-1)/2,1);
  Eigen::VectorXd TAAY = TAA.block(this->nVA_*(this->nVA_-1)/2,0,this->nOA_*(this->nOA_-1)/2,1);
  Eigen::VectorXd TABX = TAB.block(0,0,this->nVA_*this->nVA_,1);
  Eigen::VectorXd TABY = TAB.block(0,0,this->nOA_*this->nOA_,1);
  

  RealMatrix TMOAA(this->nOA_+this->nVA_,this->nOA_+this->nVA_);
  RealMatrix TMOAB(this->nOA_+this->nVA_,this->nOB_+this->nVB_);
  RealMatrix TMO(2*this->nBasis_,2*this->nBasis_);

  cout << "HERE" << endl;
  for(auto a = 0, ab = 0; a < this->nVA_; a++)
  for(auto b = 0; b < a ; b++, ab++){
    TMOAA(this->nOA_+a,this->nOA_+b) = TAAX(ab); 
    TMOAA(this->nOA_+b,this->nOA_+a) = -TAAX(ab); 
  }
  for(auto i = 0, ij = 0; i < this->nOA_; i++)
  for(auto j = 0; j < i ; j++, ij++){
    TMOAA(i,j) = TAAY(ij); 
    TMOAA(j,i) = -TAAY(ij); 
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVA_; b++, ab++){
    TMOAB(this->nOA_+a,this->nOA_+b) = TABX(ab); 
  }
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < this->nOA_; j++, ij++){
    TMOAB(i,j) = TABY(ij); 
  }


  RealMatrix TAO = (*this->singleSlater_->moA()) * TMOAA * this->singleSlater_->moA()->adjoint(); 

  Tensor<double> TAOTen(this->nBasis_,this->nBasis_);
  Tensor<double> IXAOTen(this->nBasis_,this->nBasis_);
  RealMatrix IXAO(this->nBasis_,this->nBasis_);
  RealMatrix IXMOFullAA(this->nBasis_,this->nBasis_);
  Eigen::VectorXd IXMOXAA(this->nVA_*(this->nVA_-1)/2);
  Eigen::VectorXd IXMOYAA(this->nVA_*(this->nVA_-1)/2);

  std::memcpy(&TAOTen.storage()[0],TAO.data(),TAO.size()*sizeof(double));
  contract(1.0,TAOTen,{nu,sig},*this->aoERI_,{mu,nu,lam,sig},0.0,IXAOTen,{mu,lam});
  std::memcpy(IXAO.data(),&IXAOTen.storage()[0],IXAOTen.size()*sizeof(double));

  IXMOFullAA = this->singleSlater_->moA()->adjoint() * IXAO * (*this->singleSlater_->moA());
  cout << endl << FullAA * TAA << endl;

  for(auto a = 0, ab = 0; a < this->nVA_; a++)
  for(auto b = 0; b < a ; b++, ab++){
    IXMOXAA(ab)= IXMOFullAA(this->nOA_+a,this->nOA_+b) + TAAX(ab) * ((*this->singleSlater_->epsA())(a+this->nOA_) + 
                    (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu);
  }
  for(auto i = 0, ij = 0; i < this->nOA_; i++)
  for(auto j = 0; j < i ; j++, ij++){
    IXMOYAA(ij)= IXMOFullAA(this->nOA_+a,this->nOA_+b); 
  }

  cout << endl << IXMOXAA << endl;

  
  RealMatrix ASing(this->nVA_*(this->nVA_+1)/2,this->nVA_*(this->nVA_+1)/2);
  RealMatrix BSing(this->nVA_*(this->nVA_+1)/2,this->nOA_*(this->nOA_+1)/2);
  RealMatrix CSing(this->nOA_*(this->nOA_+1)/2,this->nOA_*(this->nOA_+1)/2);
  RealMatrix FullSing(this->nVA_*(this->nVA_+1)/2+this->nOA_*(this->nOA_+1)/2,this->nVA_*(this->nVA_+1)/2+this->nOA_*(this->nOA_+1)/2);

  for(auto a = 0, ab = 0; a < this->nVA_; a++       )
  for(auto b = 0        ; b <= a        ; b++, ab++){
    for(auto c = 0, cd = 0; c < this->nVA_; c++       )
    for(auto d = 0        ; d <= c        ; d++, cd++){
      double fact = 1.0;
      if(a==b) fact *= std::sqrt(0.5);
      if(c==d) fact *= std::sqrt(0.5);
      ASing(ab,cd) = fact*(dabcdA(a,b,c,d) + dabcdA(a,b,d,c));
      if(ab == cd) 
        ASing(ab,cd) += (*this->singleSlater_->epsA())(a+this->nOA_) + 
                    (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu;
    }
  }
   
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j <= i        ; j++, ij++){
    for(auto k = 0, kl = 0; k < this->nOA_; k++      )
    for(auto l = 0        ; l <= k        ; l++, kl++){
      double fact = 1.0;
      if(i==j) fact *= std::sqrt(0.5);
      if(k==l) fact *= std::sqrt(0.5);

      CSing(ij,kl) = fact*(dijklA(i,j,k,l) + dijklA(i,j,l,k));
      if(ij == kl) 
        CSing(ij,kl) -= (*this->singleSlater_->epsA())(i) + 
                    (*this->singleSlater_->epsA())(j) - 2*Rmu;
    }
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b <= a        ; b++, ab++){
    for(auto i = 0, ij = 0; i < this->nOA_; i++      )
    for(auto j = 0        ; j <= i        ; j++, ij++){
      double fact = 1.0;
      if(i==j) fact *= std::sqrt(0.5);
      if(a==b) fact *= std::sqrt(0.5);

      BSing(ab,ij) = fact*(dabijA(a,b,i,j)+dabijA(a,b,j,i));
    }
  }

  ES.compute(ASing);
  cout << endl << "SA A:" << endl << ES.eigenvalues() << endl;
  ES.compute(CSing);
  cout << endl << "SA C:" << endl << ES.eigenvalues() << endl;

  FullSing.block(0,0,this->nVA_*(this->nVA_+1)/2,this->nVA_*(this->nVA_+1)/2) = ASing;
  FullSing.block(0,this->nVA_*(this->nVA_+1)/2,this->nVA_*(this->nVA_+1)/2,this->nOA_*(this->nOA_+1)/2) = BSing;
  FullSing.block(this->nVA_*(this->nVA_+1)/2,0,this->nOA_*(this->nOA_+1)/2,this->nVA_*(this->nVA_+1)/2) = -BSing.transpose();
  FullSing.block(this->nVA_*(this->nVA_+1)/2,this->nVA_*(this->nVA_+1)/2,this->nOA_*(this->nOA_+1)/2,this->nOA_*(this->nOA_+1)/2) = -CSing;

  EA.compute(FullSing);
  Eigen::VectorXd ERSing = -EA.eigenvalues().real();
  std::sort(ERSing.data(),ERSing.data()+ERSing.size());
  ERSing = -ERSing;
  cout << ERSing << endl << endl;
  RealMatrix TSing = EA.eigenvectors().col(0).real();
*/
  

  int VirSqAASLT   = this->nVA_*(this->nVA_-1)/2;
  int OccSqAASLT   = this->nOA_*(this->nOA_-1)/2;
  int VirSqAALT    = this->nVA_*(this->nVA_+1)/2;
  int OccSqAALT    = this->nOA_*(this->nOA_+1)/2;
  int VirSqBBSLT   = this->nVB_*(this->nVB_-1)/2;
  int OccSqBBSLT   = this->nOB_*(this->nOB_-1)/2;
  int VirSqBBLT    = this->nVB_*(this->nVB_+1)/2;
  int OccSqBBLT    = this->nOB_*(this->nOB_+1)/2;
  int VirSqAA      = this->nVA_*this->nVA_;
  int OccSqAA      = this->nOA_*this->nOA_;
  int VirSqAB      = this->nVA_*this->nVB_;
  int OccSqAB      = this->nOA_*this->nOB_;
  int FullAADim    = VirSqAASLT + OccSqAASLT;
  int FullABDim    = VirSqAB    + OccSqAB   ;
  int FullBBDim    = VirSqBBSLT + OccSqBBSLT;
  int FullSingDim  = VirSqAALT  + OccSqAALT;
  

  // Pure Spin (triplet)
  RealMatrix AAA(VirSqAASLT,VirSqAASLT);
  RealMatrix BAA(VirSqAASLT,OccSqAASLT);
  RealMatrix CAA(OccSqAASLT,OccSqAASLT);
  RealMatrix FullAA(FullAADim,FullAADim);

  RealMatrix ABB(VirSqBBSLT,VirSqBBSLT);
  RealMatrix BBB(VirSqBBSLT,OccSqBBSLT);
  RealMatrix CBB(OccSqBBSLT,OccSqBBSLT);
  RealMatrix FullBB(FullBBDim,FullBBDim);

  // Mixed Spin (singlet + triplet)
  RealMatrix AAB(VirSqAB,VirSqAB);
  RealMatrix BAB(VirSqAB,OccSqAB);
  RealMatrix CAB(OccSqAB,OccSqAB);
  RealMatrix FullAB(FullABDim,FullABDim);

  // Spin-Adapted Singlet
  RealMatrix ASing(VirSqAALT,VirSqAALT);
  RealMatrix BSing(VirSqAALT,OccSqAALT);
  RealMatrix CSing(OccSqAALT,OccSqAALT);
  RealMatrix FullSing(FullSingDim,FullSingDim);

  // Eigensolvers
  Eigen::SelfAdjointEigenSolver<RealMatrix> ES;
  Eigen::EigenSolver<RealMatrix> EA;

  // Eigensolution (Full) storage
  // Eigen values
  Eigen::VectorXd ATDAEAA;
  Eigen::VectorXd ATDAEAB;
  Eigen::VectorXd ATDAEBB;
  Eigen::VectorXd ATDAESing;
  Eigen::VectorXd CTDAEAA;
  Eigen::VectorXd CTDAEAB;
  Eigen::VectorXd CTDAEBB;
  Eigen::VectorXd CTDAESing;
  Eigen::VectorXd RPAEAA;
  Eigen::VectorXd RPAEAB;
  Eigen::VectorXd RPAEBB;
  Eigen::VectorXd RPAESing;
  // Eigen Vectors
  Eigen::VectorXd ATDATAA;
  Eigen::VectorXd ATDATAB;
  Eigen::VectorXd ATDATBB;
  Eigen::VectorXd ATDATSing;
  Eigen::VectorXd CTDATAA;
  Eigen::VectorXd CTDATAB;
  Eigen::VectorXd CTDATBB;
  Eigen::VectorXd CTDATSing;
  Eigen::VectorXd RPATAA;
  Eigen::VectorXd RPATAB;
  Eigen::VectorXd RPATBB;
  Eigen::VectorXd RPATSing;

  // AX Contract Storage
  RealMatrix ATDAMOAA(this->nBasis_,this->nBasis_);
  RealMatrix ATDAMOAB(this->nBasis_,this->nBasis_);
  RealMatrix ATDAMOSing(this->nBasis_,this->nBasis_);
  RealMatrix CTDAMOAA(this->nBasis_,this->nBasis_);
  RealMatrix CTDAMOAB(this->nBasis_,this->nBasis_);
  RealMatrix CTDAMOSing(this->nBasis_,this->nBasis_);
  RealMatrix RPAMOAA(this->nBasis_,this->nBasis_);
  RealMatrix RPAMOAB(this->nBasis_,this->nBasis_);
  RealMatrix RPAMOSing(this->nBasis_,this->nBasis_);
  RealMatrix ATDAAOAA(this->nBasis_,this->nBasis_);
  RealMatrix ATDAAOAB(this->nBasis_,this->nBasis_);
  RealMatrix ATDAAOSing(this->nBasis_,this->nBasis_);
  RealMatrix CTDAAOAA(this->nBasis_,this->nBasis_);
  RealMatrix CTDAAOAB(this->nBasis_,this->nBasis_);
  RealMatrix CTDAAOSing(this->nBasis_,this->nBasis_);
  RealMatrix RPAAOAA(this->nBasis_,this->nBasis_);
  RealMatrix RPAAOAB(this->nBasis_,this->nBasis_);
  RealMatrix RPAAOSing(this->nBasis_,this->nBasis_);

  RealMatrix ATDAIXMOAA(this->nBasis_,this->nBasis_);
  RealMatrix ATDAIXMOAB(this->nBasis_,this->nBasis_);
  RealMatrix ATDAIXMOSing(this->nBasis_,this->nBasis_);
  RealMatrix CTDAIXMOAA(this->nBasis_,this->nBasis_);
  RealMatrix CTDAIXMOAB(this->nBasis_,this->nBasis_);
  RealMatrix CTDAIXMOSing(this->nBasis_,this->nBasis_);
  RealMatrix RPAIXMOAA(this->nBasis_,this->nBasis_);
  RealMatrix RPAIXMOAB(this->nBasis_,this->nBasis_);
  RealMatrix RPAIXMOSing(this->nBasis_,this->nBasis_);
  RealMatrix ATDAIXAOAA(this->nBasis_,this->nBasis_);
  RealMatrix ATDAIXAOAB(this->nBasis_,this->nBasis_);
  RealMatrix ATDAIXAOSing(this->nBasis_,this->nBasis_);
  RealMatrix CTDAIXAOAA(this->nBasis_,this->nBasis_);
  RealMatrix CTDAIXAOAB(this->nBasis_,this->nBasis_);
  RealMatrix CTDAIXAOSing(this->nBasis_,this->nBasis_);
  RealMatrix RPAIXAOAA(this->nBasis_,this->nBasis_);
  RealMatrix RPAIXAOAB(this->nBasis_,this->nBasis_);
  RealMatrix RPAIXAOSing(this->nBasis_,this->nBasis_);

  Eigen::VectorXd ATDAAXMOAA(VirSqAASLT);
  Eigen::VectorXd ATDAAXMOAB(VirSqAB);
  Eigen::VectorXd ATDAAXMOSing(VirSqAALT);
  Eigen::VectorXd CTDAAXMOAA(OccSqAASLT);
  Eigen::VectorXd CTDAAXMOAB(OccSqAB);
  Eigen::VectorXd CTDAAXMOSing(OccSqAALT);
  Eigen::VectorXd RPAAXMOAA(VirSqAASLT + OccSqAASLT);
  Eigen::VectorXd RPAAXMOAB(VirSqAB + OccSqAB);
  Eigen::VectorXd RPAAXMOSing(VirSqAALT + OccSqAALT);

  Tensor<double> ATDAAOTenAA(this->nBasis_,this->nBasis_);
  Tensor<double> ATDAAOTenAB(this->nBasis_,this->nBasis_);
  Tensor<double> ATDAAOTenSing(this->nBasis_,this->nBasis_);
  Tensor<double> CTDAAOTenAA(this->nBasis_,this->nBasis_);
  Tensor<double> CTDAAOTenAB(this->nBasis_,this->nBasis_);
  Tensor<double> CTDAAOTenSing(this->nBasis_,this->nBasis_);
  Tensor<double> RPAAOTenAA(this->nBasis_,this->nBasis_);
  Tensor<double> RPAAOTenAB(this->nBasis_,this->nBasis_);
  Tensor<double> RPAAOTenSing(this->nBasis_,this->nBasis_);
  Tensor<double> ATDAIXAOTenAA(this->nBasis_,this->nBasis_);
  Tensor<double> ATDAIXAOTenAB(this->nBasis_,this->nBasis_);
  Tensor<double> ATDAIXAOTenSing(this->nBasis_,this->nBasis_);
  Tensor<double> CTDAIXAOTenAA(this->nBasis_,this->nBasis_);
  Tensor<double> CTDAIXAOTenAB(this->nBasis_,this->nBasis_);
  Tensor<double> CTDAIXAOTenSing(this->nBasis_,this->nBasis_);
  Tensor<double> RPAIXAOTenAA(this->nBasis_,this->nBasis_);
  Tensor<double> RPAIXAOTenAB(this->nBasis_,this->nBasis_);
  Tensor<double> RPAIXAOTenSing(this->nBasis_,this->nBasis_);
  // Build Pure-Spin Matricies
  // A
  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < a         ; b++, ab++){
    for(auto c = 0, cd = 0; c < this->nVA_; c++      )
    for(auto d = 0        ; d < c         ; d++, cd++){
      AAA(ab,cd) = DabcdA(a,b,c,d);
      if(ab == cd) 
        AAA(ab,cd) += (*this->singleSlater_->epsA())(a+this->nOA_) + 
                    (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu;
    }
  }

  // C
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < i         ; j++, ij++){
    for(auto k = 0, kl = 0; k < this->nOA_; k++      )
    for(auto l = 0        ; l < k         ; l++, kl++){
      CAA(ij,kl) = DijklA(i,j,k,l);
      if(ij == kl) 
        CAA(ij,kl) -= (*this->singleSlater_->epsA())(i) + 
                    (*this->singleSlater_->epsA())(j) - 2*Rmu;
    }
  }

  // B
  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < a         ; b++, ab++){
    for(auto i = 0, ij = 0; i < this->nOA_; i++      )
    for(auto j = 0        ; j < i         ; j++, ij++){
      BAA(ab,ij) = DabijA(a,b,i,j);
    }
  }

  FullAA.block(0,0,VirSqAASLT,VirSqAASLT)                   =  AAA;
  FullAA.block(0,VirSqAASLT,VirSqAASLT,OccSqAASLT)          =  BAA; 
  FullAA.block(VirSqAASLT,0,OccSqAASLT,VirSqAASLT)          = -BAA.adjoint();
  FullAA.block(VirSqAASLT,VirSqAASLT,OccSqAASLT,OccSqAASLT) = -CAA;

  // Build Mixed-Spin Matricies
  // A
  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVA_; b++, ab++){
    for(auto c = 0, cd = 0; c < this->nVA_; c++      )
    for(auto d = 0        ; d < this->nVA_; d++, cd++){
      AAB(ab,cd) = dabcdAB(a,b,c,d);
      if(ab == cd) 
        AAB(ab,cd) += (*this->singleSlater_->epsA())(a+this->nOA_) + 
                    (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu;
    }
  }

  // C
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < this->nOA_; j++, ij++){
    for(auto k = 0, kl = 0; k < this->nOA_; k++      )
    for(auto l = 0        ; l < this->nOA_; l++, kl++){
      CAB(ij,kl) = dijklAB(i,j,k,l);
      if(ij == kl) 
        CAB(ij,kl) -= (*this->singleSlater_->epsA())(i) + 
                    (*this->singleSlater_->epsA())(j) - 2*Rmu;
    }
  }

  // B
  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVA_; b++, ab++){
    for(auto i = 0, ij = 0; i < this->nOA_; i++      )
    for(auto j = 0        ; j < this->nOA_; j++, ij++){
      BAB(ab,ij) = dabijAB(a,b,i,j);
    }
  }

  FullAB.block(0,0,VirSqAB,VirSqAB)             =  AAB;
  FullAB.block(0,VirSqAB,VirSqAB,OccSqAB)       =  BAB; 
  FullAB.block(VirSqAB,0,OccSqAB,VirSqAB)       = -BAB.adjoint();
  FullAB.block(VirSqAB,VirSqAB,OccSqAB,OccSqAB) = -CAB;

  // Spin-Adapted Singlet Matricies
  // A
  for(auto a = 0, ab = 0; a < this->nVA_; a++       )
  for(auto b = 0        ; b <= a        ; b++, ab++){
    for(auto c = 0, cd = 0; c < this->nVA_; c++       )
    for(auto d = 0        ; d <= c        ; d++, cd++){
      double fact = 1.0;
      if(a==b) fact *= std::sqrt(0.5);
      if(c==d) fact *= std::sqrt(0.5);
      ASing(ab,cd) = fact*(dabcdA(a,b,c,d) + dabcdA(a,b,d,c));
      if(ab == cd) 
        ASing(ab,cd) += (*this->singleSlater_->epsA())(a+this->nOA_) + 
                    (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu;
    }
  }
   
  // C
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j <= i        ; j++, ij++){
    for(auto k = 0, kl = 0; k < this->nOA_; k++      )
    for(auto l = 0        ; l <= k        ; l++, kl++){
      double fact = 1.0;
      if(i==j) fact *= std::sqrt(0.5);
      if(k==l) fact *= std::sqrt(0.5);

      CSing(ij,kl) = fact*(dijklA(i,j,k,l) + dijklA(i,j,l,k));
      if(ij == kl) 
        CSing(ij,kl) -= (*this->singleSlater_->epsA())(i) + 
                    (*this->singleSlater_->epsA())(j) - 2*Rmu;
    }
  }

  // B
  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b <= a        ; b++, ab++){
    for(auto i = 0, ij = 0; i < this->nOA_; i++      )
    for(auto j = 0        ; j <= i        ; j++, ij++){
      double fact = 1.0;
      if(i==j) fact *= std::sqrt(0.5);
      if(a==b) fact *= std::sqrt(0.5);

      BSing(ab,ij) = fact*(dabijA(a,b,i,j)+dabijA(a,b,j,i));
    }
  }

  FullSing.block(0,0,VirSqAALT,VirSqAALT)                 =  ASing;
  FullSing.block(0,VirSqAALT,VirSqAALT,OccSqAALT)         =  BSing; 
  FullSing.block(VirSqAALT,0,OccSqAALT,VirSqAALT)         = -BSing.adjoint();
  FullSing.block(VirSqAALT,VirSqAALT,OccSqAALT,OccSqAALT) = -CSing;

  // Full Diagonalization
  // Pure Spin
  ES.compute(AAA);
    ATDAEAA = -ES.eigenvalues();
    std::sort(ATDAEAA.data(),ATDAEAA.data()+ATDAEAA.size());
    ATDAEAA = -ATDAEAA;
    ATDATAA = ES.eigenvectors().col(0);
  ES.compute(-CAA);
    CTDAEAA = -ES.eigenvalues();
    std::sort(CTDAEAA.data(),CTDAEAA.data()+CTDAEAA.size());
    CTDAEAA = -CTDAEAA;
    CTDATAA = ES.eigenvectors().col(0);
  EA.compute(FullAA);
    RPAEAA = -EA.eigenvalues().real();
    std::sort(RPAEAA.data(),RPAEAA.data()+RPAEAA.size());
    RPAEAA = -RPAEAA;
    RPATAA = EA.eigenvectors().col(0).real();
    
  // Mixed Spin
  ES.compute(AAB);
    ATDAEAB = -ES.eigenvalues();
    std::sort(ATDAEAB.data(),ATDAEAB.data()+ATDAEAB.size());
    ATDAEAB = -ATDAEAB;
    ATDATAB = ES.eigenvectors().col(0);
  ES.compute(-CAB);
    CTDAEAB = -ES.eigenvalues();
    std::sort(CTDAEAB.data(),CTDAEAB.data()+CTDAEAB.size());
    CTDAEAB = -CTDAEAB;
    CTDATAB = ES.eigenvectors().col(0);
  EA.compute(FullAB);
    RPAEAB = -EA.eigenvalues().real();
    std::sort(RPAEAB.data(),RPAEAB.data()+RPAEAB.size());
    RPAEAB = -RPAEAB;
    RPATAB = EA.eigenvectors().col(0).real();

  // Spin-Adapted Singlet
  ES.compute(ASing);
    ATDAESing = -ES.eigenvalues();
    std::sort(ATDAESing.data(),ATDAESing.data()+ATDAESing.size());
    ATDAESing = -ATDAESing;
    ATDATSing = ES.eigenvectors().col(0);
  ES.compute(-CSing);
    CTDAESing = -ES.eigenvalues();
    std::sort(CTDAESing.data(),CTDAESing.data()+CTDAESing.size());
    CTDAESing = -CTDAESing;
    CTDATSing = ES.eigenvectors().col(0);
  EA.compute(FullSing);
    RPAESing = -EA.eigenvalues().real();
    std::sort(RPAESing.data(),RPAESing.data()+RPAESing.size());
    RPAESing = -RPAESing;
    RPATSing = EA.eigenvectors().col(0).real();

  Eigen::IOFormat HeavyFmt(8);

  cout << "A TDA (AA) Eigenvalues:" << endl;
  cout << ATDAEAA.format(HeavyFmt) << endl << endl;
  cout << "A TDA (AB) Eigenvalues:" << endl;
  cout << ATDAEAB.format(HeavyFmt) << endl << endl;
  cout << "A TDA (SAS) Eigenvalues:" << endl;
  cout << ATDAESing.format(HeavyFmt) << endl << endl;
  cout << "C TDA (AA) Eigenvalues:" << endl;
  cout << CTDAEAA.format(HeavyFmt) << endl << endl;
  cout << "C TDA (AB) Eigenvalues:" << endl;
  cout << CTDAEAB.format(HeavyFmt) << endl << endl;
  cout << "C TDA (SAS) Eigenvalues:" << endl;
  cout << CTDAESing.format(HeavyFmt) << endl << endl;
  cout << "RPA (AA) Eigenvalues:" << endl;
  cout << RPAEAA.format(HeavyFmt)  << endl << endl;
  cout << "RPA (AB) Eigenvalues:" << endl;
  cout << RPAEAB.format(HeavyFmt)  << endl << endl;
  cout << "RPA (AB) Eigenvalues:" << endl;
  cout << RPAEAB.format(HeavyFmt)  << endl << endl;
  cout << "RPA (SAS) Eigenvalues:" << endl;
  cout << RPAESing.format(HeavyFmt)  << endl << endl;


  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < a         ; b++, ab++){
    ATDAMOAA(this->nOA_+a,this->nOA_+b) = ATDATAA(ab); 
    ATDAMOAA(this->nOA_+b,this->nOA_+a) = -ATDATAA(ab); 
  }
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < i         ; j++, ij++){
    CTDAMOAA(i,j) = CTDATAA(ij); 
    CTDAMOAA(j,i) = -CTDATAA(ij); 
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < a         ; b++, ab++){
    RPAMOAA(this->nOA_+a,this->nOA_+b) =  RPATAA(ab); 
    RPAMOAA(this->nOA_+b,this->nOA_+a) = -RPATAA(ab); 
  }
  for(auto i = 0, ij = VirSqAASLT; i < this->nOA_; i++      )
  for(auto j = 0        ; j < i         ; j++, ij++){
    RPAMOAA(i,j) =  RPATAA(ij); 
    RPAMOAA(j,i) = -RPATAA(ij); 
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVB_; b++, ab++){
    ATDAMOAB(this->nOA_+a,this->nOB_+b) = ATDATAB(ab); 
  }
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < this->nOA_; j++, ij++){
    CTDAMOAB(i,j) = CTDATAB(ij); 
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVA_; b++, ab++){
    RPAMOAB(this->nOA_+a,this->nOA_+b) =  RPATAB(ab); 
  }
  for(auto i = 0, ij = VirSqAB; i < this->nOA_; i++      )
  for(auto j = 0        ; j < this->nOA_; j++, ij++){
    RPAMOAB(i,j) =  RPATAB(ij); 
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++ )
  for(auto b = 0; b <= a; b++           , ab++){
    double fact = 1.0;
    if(a==b) fact = std::sqrt(0.5);
    ATDAMOSing(this->nOA_+a,this->nOA_+b) = fact*ATDATSing(ab);
    if(a > b)
      ATDAMOSing(this->nOA_+b,this->nOA_+a) = ATDATSing(ab);
  }

  for(auto i = 0, ij = 0; i < this->nOA_; i++ )
  for(auto j = 0; j <= i; j++           , ij++){
    double fact = 1.0;
    if(i==j) fact = std::sqrt(0.5);
    CTDAMOSing(i,j) = fact*CTDATSing(ij);
    if(a > b)
      CTDAMOSing(j,i) = -CTDATSing(ij);
  }
   
  for(auto a = 0, ab = 0; a < this->nVA_; a++ )
  for(auto b = 0; b <= a; b++           , ab++){
    double fact = 1.0;
    if(a==b) fact = std::sqrt(0.5);
    RPAMOSing(this->nOA_+a,this->nOA_+b) = fact*RPATSing(ab);
    if(a > b)
      RPAMOSing(this->nOA_+b,this->nOA_+a) = -RPATSing(ab);
  }

  for(auto i = 0, ij = VirSqAALT; i < this->nOA_; i++ )
  for(auto j = 0; j <= i; j++           , ij++){
    double fact = 1.0;
    if(i==j) fact = std::sqrt(0.5);
    RPAMOSing(i,j) = fact*RPATSing(ij);
    if(a > b)
      RPAMOSing(j,i) = -RPATSing(ij);
  }


  ATDAAOAA = (*this->singleSlater_->moA())*ATDAMOAA*this->singleSlater_->moA()->adjoint();  
  CTDAAOAA = (*this->singleSlater_->moA())*CTDAMOAA*this->singleSlater_->moA()->adjoint();  
  RPAAOAA =  (*this->singleSlater_->moA())*RPAMOAA*this->singleSlater_->moA()->adjoint();  
  ATDAAOAB = (*this->singleSlater_->moA())*ATDAMOAB*this->singleSlater_->moA()->adjoint();  
  CTDAAOAB = (*this->singleSlater_->moA())*CTDAMOAB*this->singleSlater_->moA()->adjoint();  
  RPAAOAB =  (*this->singleSlater_->moA())*RPAMOAB*this->singleSlater_->moA()->adjoint();  
  ATDAAOSing=(*this->singleSlater_->moA())*ATDAMOSing*this->singleSlater_->moA()->adjoint();
  CTDAAOSing=(*this->singleSlater_->moA())*CTDAMOSing*this->singleSlater_->moA()->adjoint();
  RPAAOSing=(*this->singleSlater_->moA())*RPAMOSing*this->singleSlater_->moA()->adjoint();  

/*
  std::memcpy(&ATDAAOTenAA.storage()[0],ATDAAOAA.data(),ATDAAOAA.size()*sizeof(double));
  std::memcpy(&CTDAAOTenAA.storage()[0],CTDAAOAA.data(),CTDAAOAA.size()*sizeof(double));
  std::memcpy(&RPAAOTenAA.storage()[0],RPAAOAA.data(),RPAAOAA.size()*sizeof(double));
  std::memcpy(&ATDAAOTenAB.storage()[0],ATDAAOAB.data(),ATDAAOAB.size()*sizeof(double));
  std::memcpy(&CTDAAOTenAB.storage()[0],CTDAAOAB.data(),CTDAAOAB.size()*sizeof(double));
  std::memcpy(&RPAAOTenAB.storage()[0],RPAAOAB.data(),RPAAOAB.size()*sizeof(double));
  std::memcpy(&ATDAAOTenSing.storage()[0],ATDAAOSing.data(),ATDAAOSing.size()*sizeof(double));
  std::memcpy(&CTDAAOTenSing.storage()[0],CTDAAOSing.data(),CTDAAOSing.size()*sizeof(double));
  std::memcpy(&RPAAOTenSing.storage()[0],RPAAOSing.data(),RPAAOSing.size()*sizeof(double));

  contract(1.0,ATDAAOTenAA,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,ATDAIXAOTenAA,{mu,nu});
  contract(1.0,CTDAAOTenAA,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,CTDAIXAOTenAA,{mu,nu});
  contract(1.0,RPAAOTenAA,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,RPAIXAOTenAA,{mu,nu});
  contract(1.0,ATDAAOTenAB,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,ATDAIXAOTenAB,{mu,nu});
  contract(1.0,CTDAAOTenAB,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,CTDAIXAOTenAB,{mu,nu});
  contract(1.0,RPAAOTenAB,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,RPAIXAOTenAB,{mu,nu});
  contract(1.0,ATDAAOTenSing,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,ATDAIXAOTenSing,{mu,nu});
  contract(1.0,CTDAAOTenSing,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,CTDAIXAOTenSing,{mu,nu});
  contract(1.0,RPAAOTenSing,{lam,sig},*this->aoERI_,{mu,lam,nu,sig},0.0,RPAIXAOTenSing,{mu,nu});

//for(auto i = 0; i < this->nBasis_; i++)
//for(auto j = 0; j < this->nBasis_; j++)
//  cout << ATDAIXAOAA(i,j) << " " <<ATDAIXAOTenAA(i,j) << "\t" <<
//  ATDAIXAOAA(i,j) /ATDAIXAOTenAA(i,j) << endl;

//CErr();

  std::memcpy(ATDAIXAOAA.data(),&ATDAIXAOTenAA.storage()[0],ATDAIXAOTenAA.size()*sizeof(double));
  std::memcpy(CTDAIXAOAA.data(),&CTDAIXAOTenAA.storage()[0],CTDAIXAOTenAA.size()*sizeof(double));
  std::memcpy(RPAIXAOAA.data(),&RPAIXAOTenAA.storage()[0],RPAIXAOTenAA.size()*sizeof(double));
  std::memcpy(ATDAIXAOAB.data(),&ATDAIXAOTenAB.storage()[0],ATDAIXAOTenAB.size()*sizeof(double));
  std::memcpy(CTDAIXAOAB.data(),&CTDAIXAOTenAB.storage()[0],CTDAIXAOTenAB.size()*sizeof(double));
  std::memcpy(RPAIXAOAB.data(),&RPAIXAOTenAB.storage()[0],RPAIXAOTenAB.size()*sizeof(double));
  std::memcpy(ATDAIXAOSing.data(),&ATDAIXAOTenSing.storage()[0],ATDAIXAOTenSing.size()*sizeof(double));
  std::memcpy(CTDAIXAOSing.data(),&CTDAIXAOTenSing.storage()[0],CTDAIXAOTenSing.size()*sizeof(double));
  std::memcpy(RPAIXAOSing.data(),&RPAIXAOTenSing.storage()[0],RPAIXAOTenSing.size()*sizeof(double));
*/

  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,ATDAAOAA,ATDAIXAOAA,ATDAAOAA,ATDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,CTDAAOAA,CTDAIXAOAA,CTDAAOAA,CTDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,RPAAOAA,RPAIXAOAA,RPAAOAA,RPAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,ATDAAOAB,ATDAIXAOAB,ATDAAOAB,ATDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,CTDAAOAB,CTDAIXAOAB,CTDAAOAB,CTDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,RPAAOAB,RPAIXAOAB,RPAAOAB,RPAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,ATDAAOSing,ATDAIXAOSing,ATDAAOSing,ATDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,CTDAAOSing,CTDAIXAOSing,CTDAAOSing,CTDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,true,RPAAOSing,RPAIXAOSing,RPAAOSing,RPAIXAOAA);

  ATDAIXMOAA = this->singleSlater_->moA()->adjoint() * ATDAIXAOAA * (*this->singleSlater_->moA()); 
  CTDAIXMOAA = this->singleSlater_->moA()->adjoint() * CTDAIXAOAA * (*this->singleSlater_->moA()); 
  RPAIXMOAA = this->singleSlater_->moA()->adjoint() * RPAIXAOAA * (*this->singleSlater_->moA()); 
  ATDAIXMOAB = this->singleSlater_->moA()->adjoint() * ATDAIXAOAB * (*this->singleSlater_->moA()); 
  CTDAIXMOAB = this->singleSlater_->moA()->adjoint() * CTDAIXAOAB * (*this->singleSlater_->moA()); 
  RPAIXMOAB = this->singleSlater_->moA()->adjoint() * RPAIXAOAB * (*this->singleSlater_->moA()); 
  ATDAIXMOSing = this->singleSlater_->moA()->adjoint() * ATDAIXAOSing * (*this->singleSlater_->moA()); 
  CTDAIXMOSing = this->singleSlater_->moA()->adjoint() * CTDAIXAOSing * (*this->singleSlater_->moA()); 
  RPAIXMOSing = this->singleSlater_->moA()->adjoint() * RPAIXAOSing * (*this->singleSlater_->moA()); 

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < a         ; b++, ab++){
    ATDAAXMOAA(ab) = ATDAIXMOAA(this->nOA_+a,this->nOA_+b) + ATDATAA(ab) *
                     ( (*this->singleSlater_->epsA())(a+this->nOA_) +
                       (*this->singleSlater_->epsA())(b+this->nOA_) - 
                       2*Rmu );
  }
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < i         ; j++, ij++){
    CTDAAXMOAA(ij) = CTDAIXMOAA(i,j) - CTDATAA(ij) *
                     ( (*this->singleSlater_->epsA())(i) +
                       (*this->singleSlater_->epsA())(j) - 
                       2*Rmu );
  }
  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < a         ; b++, ab++){
    RPAAXMOAA(ab) = RPAIXMOAA(this->nOA_+a,this->nOA_+b) + RPATAA(ab) *
                     ( (*this->singleSlater_->epsA())(a+this->nOA_) +
                       (*this->singleSlater_->epsA())(b+this->nOA_) - 
                       2*Rmu );
  }
  for(auto i = 0, ij = VirSqAASLT; i < this->nOA_; i++      )
  for(auto j = 0        ; j < i         ; j++, ij++){
    RPAAXMOAA(ij) = -RPAIXMOAA(i,j) + RPATAA(ij) *
                     ( (*this->singleSlater_->epsA())(i) +
                       (*this->singleSlater_->epsA())(j) - 
                       2*Rmu );
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVA_; b++, ab++){
    ATDAAXMOAB(ab) = ATDAIXMOAB(this->nOA_+a,this->nOA_+b) + ATDATAB(ab) *
                     ( (*this->singleSlater_->epsA())(a+this->nOA_) +
                       (*this->singleSlater_->epsA())(b+this->nOA_) - 
                       2*Rmu );
  }
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < this->nOA_; j++, ij++){
    CTDAAXMOAB(ij) = CTDAIXMOAB(i,j) - CTDATAB(ij) *
                     ( (*this->singleSlater_->epsA())(i) +
                       (*this->singleSlater_->epsA())(j) - 
                       2*Rmu );
  }
  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b < this->nVA_; b++, ab++){
    RPAAXMOAB(ab) = RPAIXMOAB(this->nOA_+a,this->nOA_+b) + RPATAB(ab) *
                     ( (*this->singleSlater_->epsA())(a+this->nOA_) +
                       (*this->singleSlater_->epsA())(b+this->nOA_) - 
                       2*Rmu );
  }
  for(auto i = 0, ij = VirSqAB; i < this->nOA_; i++      )
  for(auto j = 0        ; j < this->nOA_; j++, ij++){
    RPAAXMOAB(ij) = -RPAIXMOAB(i,j) + RPATAB(ij) *
                     ( (*this->singleSlater_->epsA())(i) +
                       (*this->singleSlater_->epsA())(j) - 
                       2*Rmu );
  }

  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b <=a         ; b++, ab++){
    double fact = 1.0;
    if(a==b) fact = std::sqrt(0.5);
    ATDAAXMOSing(ab) = fact*ATDAIXMOSing(this->nOA_+a,this->nOA_+b) + ATDATSing(ab) *
                     ( (*this->singleSlater_->epsA())(a+this->nOA_) +
                       (*this->singleSlater_->epsA())(b+this->nOA_) - 
                       2*Rmu );
  }
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j <=i         ; j++, ij++){
    double fact = 1.0;
    if(i==j) fact = std::sqrt(0.5);
    CTDAAXMOSing(ij) = fact*CTDAIXMOSing(i,j) - CTDATSing(ij) *
                     ( (*this->singleSlater_->epsA())(i) +
                       (*this->singleSlater_->epsA())(j) - 
                       2*Rmu );
  }
  for(auto a = 0, ab = 0; a < this->nVA_; a++      )
  for(auto b = 0        ; b <=a         ; b++, ab++){
    double fact = 1.0;
    if(a==b) fact = std::sqrt(0.5);
    RPAAXMOSing(ab) = fact*RPAIXMOSing(this->nOA_+a,this->nOA_+b) + RPATSing(ab) *
                     ( (*this->singleSlater_->epsA())(a+this->nOA_) +
                       (*this->singleSlater_->epsA())(b+this->nOA_) - 
                       2*Rmu );
  }
  for(auto i = 0, ij = VirSqAALT; i < this->nOA_; i++      )
  for(auto j = 0        ; j <=i         ; j++, ij++){
    double fact = 1.0;
    if(i==j) fact = std::sqrt(0.5);
    RPAAXMOSing(ij) = -fact*RPAIXMOSing(i,j) + RPATSing(ij) *
                     ( (*this->singleSlater_->epsA())(i) +
                       (*this->singleSlater_->epsA())(j) - 
                       2*Rmu );
  }
  cout << "Checking A TDA (AA) AX... |AX| = " << ATDAAXMOAA.norm() << " |R| = " << (ATDAAXMOAA - AAA*ATDATAA).norm() << endl;
  cout << "Checking C TDA (AA) AX... |AX| = " << CTDAAXMOAA.norm() << " |R| = " << (CTDAAXMOAA - CAA*CTDATAA).norm() << endl;
  cout << "Checking RPA   (AA) AX... |AX| = " << RPAAXMOAA.norm() << " |R| = " << (RPAAXMOAA - FullAA*RPATAA).norm() << endl;
  cout << "Checking A TDA (AB) AX... |AX| = " << ATDAAXMOAB.norm() << " |R| = " << (ATDAAXMOAB - AAB*ATDATAB).norm() << endl;
  cout << "Checking C TDA (AB) AX... |AX| = " << CTDAAXMOAB.norm() << " |R| = " << (CTDAAXMOAB - CAB*CTDATAB).norm() << endl;
  cout << "Checking RPA   (AB) AX... |AX| = " << RPAAXMOAB.norm() << " |R| = " << (RPAAXMOAB - FullAB*RPATAB).norm() << endl;
  cout << "Checking A TDA (SA) AX... |AX| = " << ATDAAXMOSing.norm() << " |R| = " << (ATDAAXMOSing - ASing*ATDATSing).norm() << endl;
  cout << "Checking C TDA (SA) AX... |AX| = " << CTDAAXMOSing.norm() << " |R| = " << (CTDAAXMOSing - CSing*CTDATSing).norm() << endl;
  cout << "Checking RPA   (SA) AX... |AX| = " << RPAAXMOSing.norm() << " |R| = " << (RPAAXMOSing - FullSing*RPATSing).norm() << endl;
//for(auto i = 0 ; i < RPAAXMOAB.size() ; i++) cout << RPAAXMOAB(i) << " " << (FullAB*RPATAB)(i) << endl;
}

//dbwye
