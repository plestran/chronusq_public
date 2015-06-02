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
  int nOV = singleSlater->nOV();
  this->CISEnergy_ = std::unique_ptr<RealMatrix>(new RealMatrix(2*nOV,1));
  this->CISTransDen_ = std::unique_ptr<RealMatrix>(new RealMatrix(2*nOV,2*nOV));
  this->TransDipole_ = std::unique_ptr<RealMatrix>(new RealMatrix(1,3));
  cout << "finish initialization" << endl;
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
  int nO = nOA;
  int nVA = this->singleSlater_->nVirA();
  int nV  = nVA;
  cout << "Number of Occupied: " << nO << ", Number of Virtual " << nV << " Number of Basis: "<<this->nBasis_<< ".\n";
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

  // Prepare the 2e integrals
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
  int nOV = this->singleSlater_->nOV();
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
    EigO(i) = (*this->singleSlater_->epsA())(i);
    cout << "The " << (i+1) << " eigenvalue in Occupied is: " << EigO(i) << endl;
  }
  for (auto j=0;j<nV;j++){
    EigV(j) = (*this->singleSlater_->epsA())((j+nO));
    cout << "The " << (j+1) << " eigenvalue in Virtual is: " << EigV(j) << endl;
  }

  ia = 0;
  for (auto a=0;a<nV;a++)
  for (auto i=0;i<nO;i++)
  {
    jb =0;
    for (auto b=0;b<nV;b++)
    for (auto j=0;j<nO;j++)
    {
      Ad(ia,jb) = 0.0;
      if ((a==b)&&(i==j)){
	Ad(ia,jb) = EigV(a)-EigO(i);
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
  
  // Build the ABBA matrix

  ABBA.block(0,0,2*nOV,2*nOV) = A;
  ABBA.block(2*nOV,2*nOV,2*nOV,2*nOV) = -A;
  ABBA.block(0,2*nOV,2*nOV,2*nOV) = B;
  ABBA.block(2*nOV,0,2*nOV,2*nOV) = -B;

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
  Eigen::EigenSolver<RealMatrix> TD;
  TD.compute(ABBA);
  TD.eigenvalues();
  TD.eigenvectors();

  // Print the LR-TDHF Excitation Energies
  //for (auto i=0;i<4*nOV;i++){
  //  cout << "The " << (i+1) << " LR-TDHF Exicitation Energy is: "
  //       << (TD.eigenvalues())(i) << endl;
  //}

}

void SDResponse::DavidsonCIS(){
  int nOV = this->singleSlater_->nOV();
  RealMatrix PDiag(2*nOV,1);
  PDiag = ReturnDiag();
  RealMatrix GVec(2*nOV,2*nOV);
  GVec = Guess(PDiag);
  RealMatrix Gpass = GVec.block(0,0,2*nOV,3);
  cout << Gpass << endl;
  Davidson<double> davA(this,Davidson<double>::CIS,3,&Gpass,3,&PDiag);
  davA.run(this->fileio_->out);
  cout << "The lowest 3 eigenvalue solved by Davidson Algorithm:" <<endl;
  cout << *davA.eigenvalues() << endl;

}

RealMatrix SDResponse::formRM2(RealMatrix &XMO){
  int nO = this->singleSlater_->nOccA();
  int nV = this->singleSlater_->nVirA();
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
  int nOV = this->singleSlater_->nOV();
  RealMatrix EigO(nO,1);
  RealMatrix EigV(nV,1);
  for (auto i=0;i<nO;i++){
    EigO(i) = (*this->singleSlater_->epsA())(i);
    //cout << "The " << (i+1) << " eigenvalue in Occupied is: " << EigO(i) << endl;
  }
  for (auto j=0;j<nV;j++){
    EigV(j) = (*this->singleSlater_->epsA())((j+nO));
    //cout << "The " << (j+1) << " eigenvalue in Virtual is: " << EigV(j) << endl;
  }
  enum{a,j,i,b,mu,nu,lam,sig};

  int nCol = XMO.cols();
  cout << "The dimension of input Matrix XMO:  " << 2*nOV << " * "<< nCol << endl;
  RealMatrix AX(2*nOV,nCol);
  Tensor<double> XAAOTsr(this->nBasis_,this->nBasis_);
  Tensor<double> XBAOTsr(this->nBasis_,this->nBasis_);
  RealMatrix X(nOV,1);
  RealMatrix XAAO(this->nBasis_,this->nBasis_);
  RealMatrix XBAO(this->nBasis_,this->nBasis_);
  Tensor<double> IXAO1(this->nBasis_,this->nBasis_);
  Tensor<double> IIXMO1(nV,this->nBasis_);
  Tensor<double> IXMOTsr1(nV,nO);
  Tensor<double> IXAO2(this->nBasis_,this->nBasis_);
  Tensor<double> IIXMO2(nV,this->nBasis_);
  Tensor<double> IXMOTsr2(nV,nO);
  Tensor<double> IXAO3(this->nBasis_,this->nBasis_);
  Tensor<double> IIXMO3(nV,this->nBasis_);
  Tensor<double> IXMOTsr3(nV,nO);
  Tensor<double> IXAO4(this->nBasis_,this->nBasis_);
  Tensor<double> IIXMO4(nV,this->nBasis_);
  Tensor<double> IXMOTsr4(nV,nO);
  RealMatrix IXMO1(nV,nO);
  RealMatrix IXMO2(nV,nO);
  RealMatrix IXMO3(nV,nO);
  RealMatrix IXMO4(nV,nO);
  RealMatrix IXMOA(nV,nO);
  RealMatrix IXMOB(nV,nO);
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
    RealMap XA(X.data(),nV,nO);
    //cout << "Print the XA" <<endl;
    //cout << XA << endl;
    RealMap XB(X.data()+nOV,nV,nO);
    //cout << "Print the XB" <<endl;
    //cout << XB << endl;
    XAAO = this->singleSlater_->moA()->block(0,nO,this->nBasis_,nV)*XA*this->singleSlater_->moA()->block(0,0,this->nBasis_,nO).transpose();
    XBAO = this->singleSlater_->moA()->block(0,nO,this->nBasis_,nV)*XB*this->singleSlater_->moA()->block(0,0,this->nBasis_,nO).transpose();
    // XAAO,XBAO in tensor form
    for (auto i=0;i<this->nBasis_;i++)
    for (auto j=0;j<this->nBasis_;j++)
    {
      XAAOTsr(i,j) = XAAO(i,j);
      XBAOTsr(i,j) = XBAO(i,j);
    }
    RealMatrix IXAO1t(this->nBasis_,this->nBasis_);
    // Contract A_AAAA( <mn||ls> ) with XA
    contract(1.0,XAAOTsr,{sig,nu},Dmnls,{mu,nu,lam,sig},0.0,IXAO1,{mu,lam});
    contract(1.0,LocMoAV,{mu,a},IXAO1,{mu,lam},0.0,IIXMO1,{a,lam});
    contract(1.0,LocMoAO,{lam,i},IIXMO1,{a,lam},0.0,IXMOTsr1,{a,i});
    for (auto a=0;a<nV;a++)
    for (auto i=0;i<nO;i++)
    {
      IXMO1(a,i)=IXMOTsr1(a,i);
      IXMO1(a,i)= IXMO1(a,i) + XA(a,i)*(EigV(a)-EigO(i));
    }
    // Contract A_AABB( <mn|ls> ) with XB
    contract(1.0,XBAOTsr,{sig,nu},dmnls,{mu,nu,lam,sig},0.0,IXAO2,{mu,lam});
    contract(1.0,LocMoAV,{mu,a},IXAO2,{mu,lam},0.0,IIXMO2,{a,lam});
    contract(1.0,LocMoAO,{lam,i},IIXMO2,{a,lam},0.0,IXMOTsr2,{a,i});
    for (auto a=0;a<nV;a++)
    for (auto i=0;i<nO;i++)
    {
      IXMO2(a,i)=IXMOTsr2(a,i);
    }
    // Contract A_BBAA( <mn|ls> ) with XA
    contract(1.0,XAAOTsr,{sig,nu},dmnls,{mu,nu,lam,sig},0.0,IXAO3,{mu,lam});
    contract(1.0,LocMoAV,{mu,a},IXAO3,{mu,lam},0.0,IIXMO3,{a,lam});
    contract(1.0,LocMoAO,{lam,i},IIXMO3,{a,lam},0.0,IXMOTsr3,{a,i});
    for (auto a=0;a<nV;a++)
    for (auto i=0;i<nO;i++)
    {
      IXMO3(a,i)=IXMOTsr3(a,i);
    }
    // Contract A_BBBB( <mn||ls> ) with XB
    contract(1.0,XBAOTsr,{sig,nu},Dmnls,{mu,nu,lam,sig},0.0,IXAO4,{mu,lam});
    contract(1.0,LocMoAV,{mu,a},IXAO4,{mu,lam},0.0,IIXMO4,{a,lam});
    contract(1.0,LocMoAO,{lam,i},IIXMO4,{a,lam},0.0,IXMOTsr4,{a,i});
    for (auto a=0;a<nV;a++)
    for (auto i=0;i<nO;i++)
    {
      IXMO4(a,i) = IXMOTsr4(a,i);
      IXMO4(a,i) = IXMO4(a,i) + XB(a,i)*(EigV(a)-EigO(i));
    }
    
  
    // Get the Final AX matrix
    IXMOA = IXMO1+IXMO2;
    IXMOB = IXMO3+IXMO4;
 
    // Print AX(i)
    //cout << "Print AX(i) (direct build)" <<endl;
    //cout << IXMOA <<endl;
    //cout << IXMOB << endl;
    RealMap IXMOAV(IXMOA.data(),nOV,1);
    RealMap IXMOBV(IXMOB.data(),nOV,1);
    AX.block(0,idx,nOV,1) = IXMOAV;
    AX.block(nOV,idx,nOV,1) = IXMOBV; 
  }

  return AX;
}

RealMatrix SDResponse::ReturnDiag(){
  int nO = this->singleSlater_->nOccA();
  int nV = this->singleSlater_->nVirA();
  int nOV = this->singleSlater_->nOV();
  RealMatrix EigO(nO,1);
  RealMatrix EigV(nV,1);
  for (auto i=0;i<nO;i++){
    EigO(i) = (*this->singleSlater_->epsA())(i);
  }
  for (auto j=0;j<nV;j++){
    EigV(j) = (*this->singleSlater_->epsA())((j+nO));
  }

  RealMatrix PDiag(2*nOV,1);
  for (auto a=0;a<nV;a++)
  for (auto i=0;i<nO;i++)
  {
    PDiag(a*nO+i,0) = EigV(a)-EigO(i);
    PDiag(a*nO+i+nOV,0) = EigV(a)-EigO(i);
  }

  return PDiag;
}

RealMatrix SDResponse::Guess(RealMatrix &PDiag){
  int nO = this->singleSlater_->nOccA();
  int nV = this->singleSlater_->nVirA();
  int nOV = this->singleSlater_->nOV();
  double temp1;
  double temp2;                                                                                                    
  bool isSorted;
  RealMatrix Permt(2*nOV,1);                                                                                       
  for (auto i=0;i<2*nOV;i++)                                                                                       
  { 
    Permt(i,0) = i;
  }   
  for (auto i=0;i<2*nOV-1;i++)
  {   
    isSorted = true;
    for (auto j=0;j<2*nOV-1;j++)                                                                                   
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
  RealMatrix GVec(2*nOV,2*nOV);                                                                                    
  GVec.Zero(2*nOV,2*nOV);
  for (auto i=0;i<2*nOV;i++)                                                                                       
  {   
    GVec(Permt(i,0),i) = 1.0;                                                                                      
  }   
  return GVec;
}

void SDResponse::TransDipole(){
  int nO   = this->singleSlater_->nOccA();
  int nV   = this->singleSlater_->nVirA();
  int nOV  = this->singleSlater_->nOV();
  int NBSq = this->nBasis_*this->nBasis_;
  double transdipole;
  double Tmax=0.0;
  int order=0;
  (*this->CISTransDen_).transposeInPlace();
  RealMap TransDen(this->CISTransDen_->data()+2*nOV,2*nOV,1);
  cout << "The transition density matrix we use: " << endl;
  cout << TransDen << endl;
  Tmax = TransDen(0,0);
  cout << "Tmax = " << Tmax << endl;
  for (auto i=0;i<nOV;i++)
  {
    if (fabs(TransDen(i,0))>fabs(Tmax))
    {
      Tmax = TransDen(i,0);
      cout << "Tmax = " << Tmax << endl;
      order = i;
    }
  }
  cout << "Maximum contribution is " << Tmax << " from " << order%nO+1 << " -> " << order/nO+nO+1 << endl;
  int nSub = TransDen.rows()/nOV;
  for (auto i=0,IOff=0;i<3;i++,IOff+=NBSq)
  { 
    transdipole = 0.0;
    RealMap Dipole(&this->elecDipole_->storage()[IOff],this->nBasis_,this->nBasis_);
    for (auto j=0, Shift=0;j<nSub;j++,Shift+=nOV)
    {
      RealMap TDenMO(TransDen.data()+Shift,nV,nO);
      RealMatrix TDenAO = this->singleSlater_->moA()->block(0,nO,this->nBasis_,nV)*TDenMO*this->singleSlater_->moA()->block(0,0,this->nBasis_,nO).transpose();
      TDenAO = TDenAO.transpose()*Dipole;
      transdipole += TDenAO.trace();
    }
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

