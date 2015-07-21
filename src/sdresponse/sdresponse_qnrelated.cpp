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

void SDResponse::IterativeRPA(){
  this->formGuess();
  QuasiNewton<double> davA(this);
  davA.run(this->fileio_->out);
  this->formTransDipole();
  this->formOscStrength();
  this->printExcitedStateEnergies();
} // IterativeRPA

void SDResponse::formGuess(){
  this->checkValid();
  if(!this->haveDag_) this->getDiag();
  this->davGuess_ = 
    std::unique_ptr<RealMatrix>(
      new RealMatrix(this->nSingleDim_,this->nGuess_)
    ); 
  int nRHF = 1;
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
} // formGuess

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
} // getDiag

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

  this->singleSlater_->aointegrals()->multTwoEContractDirect(XMO.cols(),false, false, false,CommA,GCommA,CommB,GCommB);


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
} // formRM3

void SDResponse::formRM4(RealCMMap& XMO, RealCMMap &Sigma, RealCMMap &Rho){
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

  bool doA = ( (this->iMeth_ == PPATDA) || (this->iMeth_ == PPRPA) );
  bool doC = ( (this->iMeth_ == PPCTDA) || (this->iMeth_ == PPRPA) );

  std::vector<RealMatrix> XAO(XMO.cols(),RealMatrix::Zero(this->nBasis_,this->nBasis_)); 
  std::vector<RealMatrix> IXAO(XMO.cols(),RealMatrix::Zero(this->nBasis_,this->nBasis_)); 
  RealMatrix IXMO(this->nBasis_,this->nBasis_);

  for(auto idx = 0; idx < XMO.cols(); idx++){
    RealVecMap X(XMO.data()+idx*this->nSingleDim_,this->nSingleDim_);
    this->formAOTDen(X,XAO[idx],XAO[idx]);
  }

  this->singleSlater_->aointegrals()->multTwoEContractDirect(XMO.cols(),true,false,true,XAO,IXAO,XAO,IXAO);


  for(auto idx = 0; idx < XMO.cols(); idx++){
    RealVecMap X(XMO.data()+idx*this->nSingleDim_,this->nSingleDim_);
    RealVecMap SVec(Sigma.data()+idx*this->nSingleDim_,this->nSingleDim_);
    this->formMOTDen(SVec,IXAO[idx],IXAO[idx]);
    this->scaleDagPPRPA(true,X,SVec);
  }

} // formRM4

