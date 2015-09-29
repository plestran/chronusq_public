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
using ChronusQ::SDResponse;
using ChronusQ::QuasiNewton;

namespace ChronusQ {
template<>
void SDResponse<double>::getDiag(){
  this->rmDiag_ = std::unique_ptr<RealCMMatrix>(new RealCMMatrix(nSingleDim_,1)); 

  if(this->Ref_ == SingleSlater<double>::TCS) {
    if(this->iMeth_ == RPA || this->iMeth_ == CIS || this->iMeth_ == STAB){
      for(auto a = 0; a < this->nV_; a++)
      for(auto i = 0; i < this->nO_; i++){
        auto ia = a*this->nO_ + i;
        auto ei = (*this->singleSlater_->epsA())(i);
        auto ea = (*this->singleSlater_->epsA())(a+this->nO_);
        (*this->rmDiag_)(ia,0) = ea - ei;
      }
      if(this->iMeth_ == RPA || this->iMeth_ == STAB)
        this->rmDiag_->block(nSingleDim_/2,0,nSingleDim_/2,1)
          = this->rmDiag_->block(0,0,nSingleDim_/2,1);
    } else if(this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA){
      bool doA = ( (this->iMeth_ == PPATDA) || (this->iMeth_ == PPRPA) );
      bool doC = ( (this->iMeth_ == PPCTDA) || (this->iMeth_ == PPRPA) );
      this->initRMu();
 
      if(doA){
          for(auto a = 0, ab = 0; a < this->nV_; a++)
          for(auto b = 0; b <  a ; b++, ab++    ){
            double ea = (*this->singleSlater_->epsA())(a + this->nO_);
            double eb = (*this->singleSlater_->epsA())(b + this->nO_);
            (*this->rmDiag_)(ab) = ea + eb - 2*this->rMu_; 
          }
      } //doA
      if(doC){
        int iOff = this->nVV_SLT_;
        double fact = -1.0;
        if(this->iMeth_ == PPRPA) fact *= -1.0;
        
          for(auto i = 0, ij = iOff; i < this->nO_; i++)
          for(auto j = 0; j <  i ;    j++, ij++    ){
            double ei = (*this->singleSlater_->epsA())(i);
            double ej = (*this->singleSlater_->epsA())(j);
            (*this->rmDiag_)(ij) = fact * (ei + ej - 2*this->rMu_); 
          } // loop I < J (I- J-)
      } // doC
      
    } else {
      CErr("Diagonal elements for given iMeth not defined");
    }
  } else {
    if(this->iMeth_ == RPA || this->iMeth_ == CIS || this->iMeth_ == STAB){
      for(auto aAlpha = 0; aAlpha < this->nVA_; aAlpha++)
      for(auto iAlpha = 0; iAlpha < this->nOA_; iAlpha++){
        auto iaAlpha = aAlpha*this->nOA_ + iAlpha; 
        auto eiAlpha = (*this->singleSlater_->epsA())(iAlpha);
        auto eaAlpha = (*this->singleSlater_->epsA())(aAlpha+this->nOA_);
        (*this->rmDiag_)(iaAlpha,0) = eaAlpha - eiAlpha;
        if(this->Ref_ == SingleSlater<double>::RHF) 
          (*this->rmDiag_)(iaAlpha+this->nOAVA_,0) = eaAlpha - eiAlpha;
      }
      if(this->Ref_ != SingleSlater<double>::RHF){
        for(auto aBeta = 0; aBeta < this->nVB_; aBeta++)
        for(auto iBeta = 0; iBeta < this->nOB_; iBeta++){
          auto iaBeta = aBeta*this->nOB_ + iBeta + this->nOAVA_; 
          auto eiBeta = (*this->singleSlater_->epsB())(iBeta);
          auto eaBeta = (*this->singleSlater_->epsB())(aBeta+this->nOB_);
          (*this->rmDiag_)(iaBeta,0) = eaBeta - eiBeta;
        }
      }
      if(this->iMeth_ == RPA || this->iMeth_ == STAB)
        this->rmDiag_->block(nSingleDim_/2,0,nSingleDim_/2,1)
          = this->rmDiag_->block(0,0,nSingleDim_/2,1);
    } else if(this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA){
      bool doA = ( (this->iMeth_ == PPATDA) || (this->iMeth_ == PPRPA) );
      bool doC = ( (this->iMeth_ == PPCTDA) || (this->iMeth_ == PPRPA) );
      this->initRMu();
 
      if(doA){
        if(this->iPPRPA_ == 0){
          for(auto aAlpha = 0, abAA = 0; aAlpha < this->nVA_; aAlpha++)
          for(auto bAlpha = 0; bAlpha <  aAlpha ; bAlpha++, abAA++    ){
            double eaAlpha = (*this->singleSlater_->epsA())(aAlpha + this->nOA_);
            double ebAlpha = (*this->singleSlater_->epsA())(bAlpha + this->nOA_);
            (*this->rmDiag_)(abAA) = eaAlpha + ebAlpha - 2*this->rMu_; 
          }
        } else if(this->iPPRPA_ == 1) {
          for(auto aAlpha = 0, abAB = 0; aAlpha < this->nVA_; aAlpha++)
          for(auto bBeta  = 0; bBeta < this->nVB_ ; bBeta++, abAB++    ){
            double eaAlpha, ebBeta;
            if(this->Ref_ == SingleSlater<double>::RHF){
              eaAlpha = (*this->singleSlater_->epsA())(aAlpha + this->nOA_);
              ebBeta  = (*this->singleSlater_->epsA())(bBeta  + this->nOA_);
            } else {
              eaAlpha = (*this->singleSlater_->epsA())(aAlpha + this->nOA_);
              ebBeta  = (*this->singleSlater_->epsB())(bBeta  + this->nOB_);
            }
            (*this->rmDiag_)(abAB) = eaAlpha + ebBeta - 2*this->rMu_; 
          } 
        } else if(this->iPPRPA_ == 2) {
          for(auto aBeta = 0, abBB = 0; aBeta < this->nVB_; aBeta++)
          for(auto bBeta = 0; bBeta <  aBeta ; bBeta++, abBB++    ){
            double eaBeta = (*this->singleSlater_->epsA())(aBeta + this->nOB_);
            double ebBeta = (*this->singleSlater_->epsA())(bBeta + this->nOB_);
            (*this->rmDiag_)(abBB) = eaBeta + ebBeta - 2*this->rMu_; 
          }
        }
      }
      if(doC){
        int iOff = 0;
        // Offset in the eigenvector for PPRPA Y amplitudes
        if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 0)) iOff = this->nVAVA_SLT_;
        if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 1)) iOff = this->nVAVB_;
        if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 2)) iOff = this->nVBVB_SLT_;
        double fact = -1.0;
        if(this->iMeth_ == PPRPA) fact *= -1.0;
        
        if(this->iPPRPA_ == 0){
          for(auto iAlpha = 0, ijAA = iOff; iAlpha < this->nOA_; iAlpha++)
          for(auto jAlpha = 0; jAlpha <  iAlpha ;    jAlpha++, ijAA++    ){
            double eiAlpha = (*this->singleSlater_->epsA())(iAlpha);
            double ejAlpha = (*this->singleSlater_->epsA())(jAlpha);
            (*this->rmDiag_)(ijAA) = fact * (eiAlpha + ejAlpha - 2*this->rMu_); 
          } // loop I < J (I-Alpha J-Alpha)
        } else if(this->iPPRPA_ == 1) {
          for(auto iAlpha = 0, ijAB = iOff; iAlpha < this->nOA_; iAlpha++)
          for(auto jBeta  = 0; jBeta < this->nOB_ ;  jBeta++, ijAB++     ){
            double eiAlpha, ejBeta;
            if(this->Ref_ == SingleSlater<double>::RHF){
              eiAlpha = (*this->singleSlater_->epsA())(iAlpha);
              ejBeta  = (*this->singleSlater_->epsA())(jBeta );
            } else {
              eiAlpha = (*this->singleSlater_->epsA())(iAlpha);
              ejBeta  = (*this->singleSlater_->epsB())(jBeta );
            }
            (*this->rmDiag_)(ijAB) = fact * (eiAlpha + ejBeta - 2*this->rMu_); 
          } // loop IJ (I-Alpha J-Beta)
        } else if(this->iPPRPA_ == 2) {
          for(auto iBeta = 0, ijBB = iOff; iBeta < this->nOB_; iBeta++)
          for(auto jBeta = 0; jBeta <  iBeta ;    jBeta++, ijBB++    ){
            double eiBeta = (*this->singleSlater_->epsA())(iBeta);
            double ejBeta = (*this->singleSlater_->epsA())(jBeta);
            (*this->rmDiag_)(ijBB) = fact * (eiBeta + ejBeta - 2*this->rMu_); 
          } // loop I < J (I-Beta J-Beta)
        }
      } // doC
      
    } else {
      CErr("Diagonal elements for given iMeth not defined");
    }
  }
  this->haveDag_ = true;
} // getDiag

template<>
void SDResponse<double>::reoptWF(){
  int maxStabIter = 20;
  double small = 1e-10;
  bool stable = false;
  int NTCSxNBASIS = this->nTCS_*this->nBasis_;

  int lenRealScr    = 0;
  int lenComplexScr = 0;
  int lenMat = NTCSxNBASIS*NTCSxNBASIS;
  int lenEig = NTCSxNBASIS;
  int lWork  = 4*NTCSxNBASIS;
  int INFO;
  char UPLO = 'L', JOBZ = 'V';

  lenRealScr += lenMat; // Stability step in MO basis
  lenRealScr += lenMat; // Matrix exponential
  lenRealScr += lenEig; // Eigenvalues
  lenRealScr += std::max(1,3*NTCSxNBASIS-1); // RWORK LAPACK Workspace


  lenComplexScr += lenMat; // Stability step in MO basis
  lenComplexScr += lenMat; // Matrix exponential
  lenComplexScr += lenMat; // BSCR
  lenComplexScr += lWork;  // LAPACK workspace (WORK)



  double * REAL_SCR = new double[lenRealScr];

  double * realStab    = REAL_SCR;
  double * realExpStab = realStab    + lenMat;
  double * W           = realExpStab + lenMat;
  double * RWORK       = W           + lenEig;

  dcomplex * COMPLEX_SCR = new dcomplex[lenComplexScr];

  dcomplex * complexStab    = COMPLEX_SCR;
  dcomplex * complexExpStab = complexStab    + lenMat;
  dcomplex * BSCR           = complexExpStab + lenMat;
  dcomplex * WORK           = BSCR           + lenMat;

  RealCMMap    AReal(realStab   ,NTCSxNBASIS,NTCSxNBASIS);
  RealCMMap ExpAReal(realExpStab,NTCSxNBASIS,NTCSxNBASIS);

  ComplexCMMap    AComplex(complexStab   ,NTCSxNBASIS,NTCSxNBASIS);
  ComplexCMMap    BComplex(BSCR          ,NTCSxNBASIS,NTCSxNBASIS);
  ComplexCMMap ExpAComplex(complexExpStab,NTCSxNBASIS,NTCSxNBASIS);


  for(auto iter = 0; iter < maxStabIter; iter++){
    AReal.setZero();
    ExpAReal.setZero();
    AComplex.setZero();
    BComplex.setZero();
    ExpAComplex.setZero();

    bool isPositive = (*this->omega_)(0) > 0.0;
    bool isSmall    = std::abs((*this->omega_)(0)) < small;
    if(isPositive || isSmall) {
      stable = true;
      break;
    }
    for(auto a = this->nO_, ia = 0; a < NTCSxNBASIS; a++)
    for(auto i = 0         ; i < this->nO_; i++, ia++){
      AComplex(a,i) = dcomplex(0.0, 0.9*(*this->transDen_)(ia,0));
      AComplex(i,a) = dcomplex(0.0,-0.9*(*this->transDen_)(ia,0));
    }
// cout << AComplex.imag() << endl;

    zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,complexStab,&NTCSxNBASIS,W,WORK,&lWork,RWORK,&INFO);
    std::memcpy(BSCR,complexStab,lenMat*sizeof(dcomplex));
    for(auto i = 0; i < NTCSxNBASIS; i++){
      dcomplex scal = std::exp(dcomplex(0.0,-2.0*W[i]));
      BComplex.col(i) *= scal;
    }
   ExpAComplex = BComplex * AComplex.adjoint();
   ExpAReal = ExpAComplex.real();
// prettyPrint(cout,ExpAReal,"Exp(J)");
// prettyPrint(cout,ExpAReal*ExpAReal.adjoint(),"Exp(J)");

   (*this->singleSlater_->moA()) *= ExpAReal;
   this->singleSlater_->formDensity();
   this->singleSlater_->formFock();
   this->singleSlater_->SCF();
   QuasiNewton<double> dav(this);
   dav.run(this->fileio_->out);
 //CErr();
  } // loop iter
  if(stable){
    this->singleSlater_->computeEnergy();
    this->singleSlater_->computeMultipole();
  } else {
    CErr("Stability failed to Re-Optimize Wavefunction",this->fileio_->out);
  }
} // reoptWF

template<>
void SDResponse<double>::formGuess(){
  this->checkValid();
  if(!this->haveDag_) this->getDiag();
  this->davGuess_ = 
    std::unique_ptr<RealMatrix>(
      new RealMatrix(this->nSingleDim_,this->nGuess_)
    ); 
  int nRPA = 1;
  if(this->iMeth_==RPA || this->iMeth_ == STAB) nRPA *= 2;
  int nCPY;
  if(this->iMeth_ == CIS || this->iMeth_ == RPA || this->iMeth_ == STAB) 
    nCPY = this->nSingleDim_ / nRPA;
  else if(this->iMeth_ == PPRPA){
    if(this->Ref_ == SingleSlater<double>::TCS)
      nCPY = this->nVV_SLT_;
    else {
      if(this->iPPRPA_ == 0)
        nCPY = this->nVAVA_SLT_;
      else if(this->iPPRPA_ = 1)
        nCPY = this->nVAVB_;
      else if(this->iPPRPA_ = 2)
        nCPY = this->nVBVB_SLT_;
    }
  } else {
    nCPY = this->nSingleDim_;
  }
  RealMatrix dagCpy(nCPY,1);
  std::memcpy(dagCpy.data(),this->rmDiag_->data(),dagCpy.size()*sizeof(double));
  std::sort(dagCpy.data(),dagCpy.data()+dagCpy.size());
  std::vector<int> alreadyAdded; 
  for(auto i = 0; i < this->nGuess_; i++){
    int indx;
    for(auto k = 0; k < dagCpy.size(); k++){
      auto it = std::find(alreadyAdded.begin(),alreadyAdded.end(),k);
      if((dagCpy(i % nCPY,0) == (*this->rmDiag_)(k,0)) && 
          it == alreadyAdded.end()){
        indx = k;
        alreadyAdded.push_back(indx);
        break;
      }
    }
    (*this->davGuess_)(indx,i) = 1.0;
  }
} // formGuess


template<>
void SDResponse<double>::IterativeRPA(){
  bool hasProp = ((this->iMeth_==CIS || this->iMeth_==RPA || this->iMeth_==PPATDA));
  this->formGuess();
//if(this->iMeth_ == PPATDA) CErr();
  QuasiNewton<double> davA(this);
  davA.run(this->fileio_->out);
  if(hasProp){
    this->formTransDipole();
    this->formOscStrength();
  }
//prettyPrint(this->fileio_->out,*this->transDen_,"Trans den");
  this->printExcitedStateEnergies();
  
  if(this->iMeth_ == STAB) this->reoptWF();
} // IterativeRPA

template<>
void SDResponse<double>::formRM3(RealCMMap &XMO, RealCMMap &Sigma, RealCMMap &Rho){
/*
 *  Forms sigma (and possibly rho) for the linear transfomation of E^(2)
 *  (and possibly S^(2) ) onto trial vectors (or any vector in general)
 *
 *  Adapted from Helgaker, et al. JCP 113, 8908 (2000)
 *
 *  DBWY (2015)
 */

  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  std::vector<RealMatrix> CommA(XMO.cols(),RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS));
  std::vector<RealMatrix> GCommA(XMO.cols(),RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS));
  RealMatrix XAAO(NTCSxNBASIS,NTCSxNBASIS);
  RealMatrix SigAOA(NTCSxNBASIS,NTCSxNBASIS);

  std::vector<RealMatrix> CommB,GCommB;
  RealMatrix XBAO,SigAOB, RhoAOA, RhoAOB;
  if(this->Ref_ != SingleSlater<double>::TCS){
    CommB  = std::vector<RealMatrix>(XMO.cols(),RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS));
    GCommB = std::vector<RealMatrix>(XMO.cols(),RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS));
    XBAO   = RealMatrix (NTCSxNBASIS,NTCSxNBASIS);
    SigAOB = RealMatrix (NTCSxNBASIS,NTCSxNBASIS);
  }
  if(this->iMeth_ == RPA){
    RhoAOA = RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);
    if(this->Ref_ != SingleSlater<double>::TCS)
      RhoAOB = RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);
  }

  RealMatrix SDA,SDB;
/*
  if(this->Ref_ == SingleSlater<double>::TCS){
    SDA = RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);
    GS  = RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);
    
    for(auto I = 0, i = 0; i < this->nBasis_; I += 2, i++)
    for(auto J = 0, j = 0; j < this->nBasis_; J += 2, j++) {
      GS(I,J)     = (*this->singleSlater_->aointegrals()->overlap_)(i,j);
      GS(I+1,J+1) = (*this->singleSlater_->aointegrals()->overlap_)(i,j);
    }
    SDA = GS * (*this->singleSlater_->densityA());
  } else {
    SDA = 
      (*this->singleSlater_->aointegrals()->overlap_) * (*this->singleSlater_->densityA());
    if(this->Ref_ != SingleSlater<double>::RHF) 
      SDB = 
        (*this->singleSlater_->aointegrals()->overlap_) * (*this->singleSlater_->densityB());
  }
*/
  SDA = (*this->singleSlater_->aointegrals()->overlap_) * (*this->singleSlater_->densityA());
  if(!this->singleSlater_->isClosedShell && this->Ref_ != SingleSlater<double>::TCS)
    SDB = 
      (*this->singleSlater_->aointegrals()->overlap_) * (*this->singleSlater_->densityB());
  

  double fact = 1.0;
  if(this->Ref_ == SingleSlater<double>::RHF) fact = 0.5;
  int iOff = this->nOAVA_ + this->nOBVB_;
  if(this->Ref_ == SingleSlater<double>::TCS) iOff = this->nOV_;

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
    if(this->Ref_ == SingleSlater<double>::RHF){ 
      CommB[idx] =  fact * XBAO * SDA;
      CommB[idx] += fact * SDA.adjoint() * XBAO;
    } else if(this->Ref_ != SingleSlater<double>::TCS) {
      CommB[idx] =  fact * XBAO * SDB;
      CommB[idx] += fact * SDB.adjoint() * XBAO;
    }
  }

//this->singleSlater_->aointegrals()->multTwoEContractDirect(XMO.cols(),false,false,false,
//  (this->nTCS_==2),CommA,GCommA,CommB,GCommB);
  for(auto idx = 0; idx < XMO.cols(); idx++)
    this->singleSlater_->aointegrals()->twoEContractN4(false,false,true,false,(this->nTCS_==2),CommA[idx],
      GCommA[idx],CommB[idx],GCommB[idx]);


  for(auto idx = 0; idx < XMO.cols(); idx++){
/*
    if(this->Ref_ == SingleSlater<double>::TCS){
      SigAOA =  (*this->singleSlater_->fockA()) * CommA[idx] * GS;
      SigAOA -= GS * CommA[idx] * (*this->singleSlater_->fockA());
    } else {
      SigAOA =  
        (*this->singleSlater_->fockA()) * CommA[idx] * 
        (*this->singleSlater_->aointegrals()->overlap_);
      SigAOA -= 
        (*this->singleSlater_->aointegrals()->overlap_) * CommA[idx] * 
        (*this->singleSlater_->fockA());
    }
*/
    SigAOA =  
      (*this->singleSlater_->fockA()) * CommA[idx] * 
      (*this->singleSlater_->aointegrals()->overlap_);
    SigAOA -= 
      (*this->singleSlater_->aointegrals()->overlap_) * CommA[idx] * 
      (*this->singleSlater_->fockA());
    SigAOA += fact * GCommA[idx] * SDA.adjoint();
    SigAOA -= fact * SDA * GCommA[idx];

    if(this->Ref_ == SingleSlater<double>::RHF) {
      SigAOB =  
        (*this->singleSlater_->fockA()) * CommB[idx] * 
        (*this->singleSlater_->aointegrals()->overlap_);
      SigAOB -= 
        (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * 
        (*this->singleSlater_->fockA());
      SigAOB += fact * GCommB[idx] * SDA.adjoint();
      SigAOB -= fact * SDA * GCommB[idx];
    } else if(this->Ref_ != SingleSlater<double>::TCS) {
      SigAOB =  
        (*this->singleSlater_->fockB()) * CommB[idx] * 
        (*this->singleSlater_->aointegrals()->overlap_);
      SigAOB -= 
        (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * 
        (*this->singleSlater_->fockB());
      SigAOB += fact * GCommB[idx] * SDB.adjoint();
      SigAOB -= fact * SDB * GCommB[idx];
    }

    RealVecMap SVec(Sigma.data()+idx*this->nSingleDim_,this->nSingleDim_);
    this->formMOTDen(SVec,SigAOA,SigAOB);

    if(this->iMeth_ == RPA){
/*
      if(this->Ref_ == SingleSlater<double>::TCS)
        RhoAOA =  GS * CommA[idx] * GS;
      else {
        RhoAOA =  
          (*this->singleSlater_->aointegrals()->overlap_) * CommA[idx] * 
          (*this->singleSlater_->aointegrals()->overlap_);
        RhoAOB =  
          (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * 
          (*this->singleSlater_->aointegrals()->overlap_);
      }
*/
      RhoAOA =  
        (*this->singleSlater_->aointegrals()->overlap_) * CommA[idx] * 
        (*this->singleSlater_->aointegrals()->overlap_);
      if(this->Ref_ != SingleSlater<double>::TCS)
        RhoAOB =  
          (*this->singleSlater_->aointegrals()->overlap_) * CommB[idx] * 
          (*this->singleSlater_->aointegrals()->overlap_);
      RealVecMap RVec(Rho.data()+idx*this->nSingleDim_,this->nSingleDim_);
      this->formMOTDen(RVec,RhoAOA,RhoAOB);
    }
  }
} // formRM3

template<>
void SDResponse<double>::formRM4(RealCMMap& XMO, RealCMMap &Sigma, RealCMMap &Rho){

  bool doA = ( (this->iMeth_ == PPATDA) || (this->iMeth_ == PPRPA) );
  bool doC = ( (this->iMeth_ == PPCTDA) || (this->iMeth_ == PPRPA) );

  std::vector<RealMatrix> 
    XAO(XMO.cols(),RealMatrix::Zero(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  std::vector<RealMatrix> 
    IXAO(XMO.cols(),RealMatrix::Zero(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  RealMatrix IXMO(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);

  for(auto idx = 0; idx < XMO.cols(); idx++){
    RealVecMap X(XMO.data()+idx*this->nSingleDim_,this->nSingleDim_);
    this->formAOTDen(X,XAO[idx],XAO[idx]);
  }

  bool RHF   = this->Ref_ != SingleSlater<double>::TCS;
  bool doTCS = this->Ref_ == SingleSlater<double>::TCS;
//this->singleSlater_->aointegrals()->multTwoEContractDirect(XMO.cols(),RHF,false,true,
//  doTCS,XAO,IXAO,XAO,IXAO);
//cout << "HERE" << endl;
  for(auto idx = 0; idx < XMO.cols(); idx++)
    this->singleSlater_->aointegrals()->twoEContractN4(false,false,false,true,doTCS,XAO[idx],
      IXAO[idx],XAO[idx],IXAO[idx]);
//cout << "HERE" << endl;

  for(auto idx = 0; idx < XMO.cols(); idx++){
    RealVecMap X(XMO.data()+idx*this->nSingleDim_,this->nSingleDim_);
    RealVecMap SVec(Sigma.data()+idx*this->nSingleDim_,this->nSingleDim_);
    this->formMOTDen(SVec,IXAO[idx],IXAO[idx]);
    this->scaleDagPPRPA(true,X,SVec);
  }

} // formRM4

} // namespace ChronusQ
