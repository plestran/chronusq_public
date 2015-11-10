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
void SDResponse<double>::placeVOOV(const RealVecMap &TMOV, RealMatrix &TMOA, RealMatrix &TMOB){
  int iOff = this->nOAVA_ + this->nOBVB_;
  if(this->Ref_ == SingleSlater<double>::TCS) iOff = this->nOV_;

  if(this->Ref_ == SingleSlater<double>::TCS){
    for(auto a = this->nO_, ia = 0; a < this->nTCS_*this->nBasis_; a++)
    for(auto i = 0         ; i < this->nO_; i++, ia++){
      TMOA(a,i) = TMOV(ia);
      if(this->iMeth_ == RPA || this->iMeth_ == STAB) TMOA(i,a) = TMOV(ia+iOff);
    }
  } else {
    for(auto a = this->nOA_, ia = 0; a < this->nBasis_; a++)
    for(auto i = 0         ; i < this->nOA_; i++, ia++){
      TMOA(a,i) = TMOV(ia);
      if(this->iMeth_ == RPA || this->iMeth_ == STAB) TMOA(i,a) = TMOV(ia+iOff);
    }
    for(auto a = this->nOB_, ia = this->nOAVA_; a < this->nBasis_; a++)
    for(auto i = 0         ; i < this->nOB_; i++, ia++){
      TMOB(a,i) = TMOV(ia);
      if(this->iMeth_ == RPA || this->iMeth_ == STAB) TMOB(i,a) = TMOV(ia+iOff);
    }
  }
} //placeVOOV

template<>
void SDResponse<double>::placeVVOO(const RealVecMap &TMOV, RealMatrix &TMO){
  bool doX = ( (this->iMeth_ == PPATDA) || (this->iMeth_ == PPRPA) );
  bool doY = ( (this->iMeth_ == PPCTDA) || (this->iMeth_ == PPRPA) );

  if(this->Ref_ == SingleSlater<double>::TCS){
    if(doX){
      for(auto a = 0, ab = 0; a < this->nV_; a++)
      for(auto b = 0; b  < a;           b++, ab++){
        TMO(this->nO_+a,this->nO_+b) =  TMOV(ab); 
        TMO(this->nO_+b,this->nO_+a) = -TMOV(ab); 
      } // loop AB A < B (A-Alpha B-Alpha)
    } // doX
    if(doY){
      int iOff = this->nVV_SLT_;
      for(auto i = 0, ij = iOff; i < this->nO_; i++)
      for(auto j = 0; j  < i   ;           j++, ij++){
        TMO(i,j) =  TMOV(ij);
        TMO(j,i) = -TMOV(ij);
      } // loop IJ I < J (I-Alpha J-Alpha)
    } // doY
  } else {
    if(doX){
      if(this->iPPRPA_ == 0){ // AA block
        for(auto a = 0, ab = 0; a < this->nVA_; a++)
        for(auto b = 0; b  < a;           b++, ab++){
          TMO(this->nOA_+a,this->nOA_+b) =  TMOV(ab); 
          TMO(this->nOA_+b,this->nOA_+a) = -TMOV(ab); 
        } // loop AB A < B (A-Alpha B-Alpha)
      } else if(this->iPPRPA_ == 1) { // AB block
        for(auto a = 0, ab = 0; a < this->nVA_; a++)
        for(auto b = 0;  b < this->nVB_;  b++, ab++){
          TMO(this->nOA_+a,this->nOB_+b) = TMOV(ab); 
        } // loop AB (A-Alpha B-Beta)
      } else if(this->iPPRPA_ == 2) { // BB block
        for(auto a = 0, ab = 0; a < this->nVB_; a++)
        for(auto b = 0; b  < a;           b++, ab++){
          TMO(this->nOB_+a,this->nOB_+b) =  TMOV(ab); 
          TMO(this->nOB_+b,this->nOB_+a) = -TMOV(ab); 
        } // loop AB A < B (A-Alpha B-Alpha)
      }
    } // doX
    if(doY){
      int iOff = 0;
      // Offset in the eigenvector for PPRPA Y amplitudes
      if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 0)) iOff = this->nVAVA_SLT_;
      if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 1)) iOff = this->nVAVB_;
      if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 2)) iOff = this->nVBVB_SLT_;
 
      if(this->iPPRPA_ == 0) { // AA Block
        for(auto i = 0, ij = iOff; i < this->nOA_; i++)
        for(auto j = 0; j  < i   ;           j++, ij++){
          TMO(i,j) =  TMOV(ij);
          TMO(j,i) = -TMOV(ij);
        } // loop IJ I < J (I-Alpha J-Alpha)
      } else if(this->iPPRPA_ == 1) { // AB Block
        for(auto i = 0, ij = iOff; i < this->nOA_; i++)
        for(auto j = 0;  j < this->nOB_;     j++, ij++){
          TMO(i,j) =  TMOV(ij);
        } // loop IJ (I-Alpha J-Beta)
      } else if(this->iPPRPA_ == 2) { // BB Block
        for(auto i = 0, ij = iOff; i < this->nOB_; i++)
        for(auto j = 0; j  < i   ;           j++, ij++){
          TMO(i,j) =  TMOV(ij);
          TMO(j,i) = -TMOV(ij);
        } // loop IJ I < J (I-Beta J-Beta)
      }
    } // doY
  }
} // placeVVOO

template<>
void SDResponse<double>::retrvVOOV(RealVecMap &TMOV, const RealMatrix &TMOA, const RealMatrix &TMOB){
  int iOff = this->nOAVA_ + this->nOBVB_;
  if(this->Ref_ == SingleSlater<double>::TCS) iOff = this->nOV_;

  if(this->Ref_ == SingleSlater<double>::TCS) {
    for(auto a = this->nO_, ia = 0; a < this->nTCS_*this->nBasis_; a++)
    for(auto i = 0         ; i < this->nO_; i++, ia++){
      TMOV(ia) = TMOA(a,i);
      if(this->iMeth_ == RPA || this->iMeth_ == STAB) TMOV(ia+iOff) = -TMOA(i,a);
    }
  } else {
    for(auto a = this->nOA_, ia = 0; a < this->nBasis_; a++)
    for(auto i = 0         ; i < this->nOA_; i++, ia++){
      TMOV(ia) = TMOA(a,i);
      if(this->iMeth_ == RPA || this->iMeth_ == STAB) TMOV(ia+iOff) = -TMOA(i,a);
    }
    for(auto a = this->nOB_, ia = this->nOAVA_; a < this->nBasis_; a++)
    for(auto i = 0         ; i < this->nOB_; i++, ia++){
      TMOV(ia) = TMOB(a,i);
      if(this->iMeth_ == RPA || this->iMeth_ == STAB) TMOV(ia+iOff) = -TMOB(i,a);
    }
  }
} // retrvVOOV

template<>
void SDResponse<double>::retrvVVOO(RealVecMap &TMOV, const RealMatrix &TMO){
  bool doX = ( (this->iMeth_ == PPATDA) || (this->iMeth_ == PPRPA) );
  bool doY = ( (this->iMeth_ == PPCTDA) || (this->iMeth_ == PPRPA) );

  if(this->Ref_ == SingleSlater<double>::TCS){
    if(doX){
      for(auto a = 0, ab = 0; a < this->nV_; a++)
      for(auto b = 0; b  < a;           b++, ab++){
        TMOV(ab) = TMO(this->nO_+a,this->nO_+b);
      } // loop AB A < B (A-Alpha B-Alpha)
    }
    if(doY){
      auto iOff = this->nVV_SLT_;
      double fact = 1.0;
      if(this->iMeth_ == PPRPA) fact *= -1.0;
      for(auto i = 0, ij = iOff; i < this->nO_; i++)
      for(auto j = 0; j  < i   ;           j++, ij++){
        TMOV(ij) = fact*TMO(i,j);
      } // loop IJ I < J (I-Alpha J-Alpha)
    }
  } else {
    if(doX){
      if(this->iPPRPA_ == 0){ // AA block
        for(auto a = 0, ab = 0; a < this->nVA_; a++)
        for(auto b = 0; b  < a;           b++, ab++){
          TMOV(ab) = TMO(this->nOA_+a,this->nOA_+b);
        } // loop AB A < B (A-Alpha B-Alpha)
      } else if(this->iPPRPA_ == 1) { // AB block
        for(auto a = 0, ab = 0; a < this->nVA_; a++)
        for(auto b = 0;  b < this->nVB_;  b++, ab++){
          TMOV(ab) = TMO(this->nOA_+a,this->nOB_+b);
        } // loop AB (A-Alpha B-Beta)
      } else if(this->iPPRPA_ == 2) { // BB block
        for(auto a = 0, ab = 0; a < this->nVB_; a++)
        for(auto b = 0; b  < a;           b++, ab++){
          TMOV(ab) = TMO(this->nOA_+a,this->nOB_+b);
        } // loop AB A < B (A-Alpha B-Alpha)
      }
    } // doX
    if(doY){
      int iOff = 0;
      // Offset in the eigenvector for PPRPA Y amplitudes
      if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 0)) iOff = this->nVAVA_SLT_;
      if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 1)) iOff = this->nVAVB_;
      if((this->iMeth_ == PPRPA) && (this->iPPRPA_ == 2)) iOff = this->nVBVB_SLT_;
      double fact = 1.0;
      if(this->iMeth_ == PPRPA) fact *= -1.0;
 
      if(this->iPPRPA_ == 0) { // AA Block
        for(auto i = 0, ij = iOff; i < this->nOA_; i++)
        for(auto j = 0; j  < i   ;           j++, ij++){
          TMOV(ij) = fact*TMO(i,j);
        } // loop IJ I < J (I-Alpha J-Alpha)
      } else if(this->iPPRPA_ == 1) { // AB Block
        for(auto i = 0, ij = iOff; i < this->nOA_; i++)
        for(auto j = 0;  j < this->nOB_;     j++, ij++){
          TMOV(ij) =  fact*TMO(i,j);
        } // loop IJ (I-Alpha J-Beta)
      } else if(this->iPPRPA_ == 2) { // BB Block
        for(auto i = 0, ij = iOff; i < this->nOB_; i++)
        for(auto j = 0; j  < i   ;           j++, ij++){
          TMOV(ij) = fact*TMO(i,j);
        } // loop IJ I < J (I-Beta J-Beta)
      }
    } // doY
  }
} //retrvVVOO

template<>
void SDResponse<double>::formAOTDen(const RealVecMap &TMOV, RealMatrix &TAOA, RealMatrix &TAOB){
  bool doVOOV = (this->iMeth_ == CIS || this->iMeth_ == RPA || this->iMeth_ == STAB); 
  bool doVVOO = (this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA);

  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  RealMatrix TMOA,TMOB;
  TMOA = RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);
  if(!doVVOO && this->Ref_ != SingleSlater<double>::TCS) 
    TMOB = RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);

  if(doVOOV) this->placeVOOV(TMOV,TMOA,TMOB);
  else if(doVVOO) this->placeVVOO(TMOV,TMOA);

  TAOA = (*this->singleSlater_->moA()) * TMOA * this->singleSlater_->moA()->adjoint();
  if(!doVVOO && this->Ref_ != SingleSlater<double>::TCS){
    if(this->Ref_ == SingleSlater<double>::RHF)
      TAOB = (*this->singleSlater_->moA()) * TMOB * this->singleSlater_->moA()->adjoint();
    else
      TAOB = (*this->singleSlater_->moB()) * TMOB * this->singleSlater_->moB()->adjoint();
  }
} //formAOTDen

template<>
void SDResponse<double>::formMOTDen(RealVecMap &TMOV, const RealMatrix &TAOA, const RealMatrix &TAOB){
  bool doVOOV = (this->iMeth_ == CIS || this->iMeth_ == RPA || this->iMeth_ == STAB); 
  bool doVVOO = (this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA);
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;

  RealMatrix TMOA,TMOB;
  TMOA = RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);
  if(!doVVOO && this->Ref_ != SingleSlater<double>::TCS) 
    TMOB = RealMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);

  TMOA = this->singleSlater_->moA()->adjoint() * TAOA * (*this->singleSlater_->moA());
  if(!doVVOO && this->Ref_ != SingleSlater<double>::TCS) {
    if(this->Ref_ == SingleSlater<double>::RHF)
      TMOB = this->singleSlater_->moA()->adjoint() * TAOB * (*this->singleSlater_->moA());
    else 
      TMOB = this->singleSlater_->moB()->adjoint() * TAOB * (*this->singleSlater_->moB());
  }
  if(doVOOV) this->retrvVOOV(TMOV,TMOA,TMOB);
  else if(doVVOO) this->retrvVVOO(TMOV,TMOA);
  
}// formMOTDen

template<>
void SDResponse<double>::initRMu(){
  // RMu = [ e(HOMO) + e(LUMO) ] / 2
  if(this->Ref_ == SingleSlater<double>::TCS){
    this->rMu_ = ( (*this->singleSlater_->epsA())(this->nO_-1) + 
                   (*this->singleSlater_->epsA())(this->nO_)    ) / 2.0;
  } else {
    if(this->Ref_ == SingleSlater<double>::RHF || this->nOB_ == 0)
      this->rMu_ = ( (*this->singleSlater_->epsA())(this->nOA_-1) + 
                     (*this->singleSlater_->epsA())(this->nOA_)    ) / 2.0;
    else if(this->nOB_ > 0)
      this->rMu_ = ( std::max( (*this->singleSlater_->epsA())(this->nOA_-1), 
                               (*this->singleSlater_->epsB())(this->nOB_-1) ) +
                     std::max( (*this->singleSlater_->epsA())(this->nOA_  ), 
                               (*this->singleSlater_->epsB())(this->nOB_  ) ) ) / 2.0;
  }
} // initRMu

template<>
void SDResponse<double>::scaleDagPPRPA(bool inplace, RealVecMap &T, RealVecMap &IX, RealVecMap *AX){
  // inplace triggers whether or not AX is populated (or touched for that matter)
  if(!this->haveDag_) this->getDiag();
  if(inplace) IX    = IX + T.cwiseProduct(*this->rmDiag_);
  else        (*AX) = IX + T.cwiseProduct(*this->rmDiag_);
} // scaleDagPPRPA

template<>
void SDResponse<double>::initMeth(){
  if(this->nSek_ == 0) 
    CErr("Must set NSek before initializing a PSCF method",this->fileio_->out);
  if((this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA)
      && !(this->iPPRPA_ >= 0 && this->iPPRPA_ <= 2))
    CErr("Invalid iPPRPA_ in SDResponse::initMeth()",this->fileio_->out);
/*
  if((this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA)
      && this->Ref_ == SingleSlater<double>::TCS)
    CErr("PPRPA/PPTDA NYI for Spinor Reference");
*/

  if(this->iMeth_ == CIS)                              this->nSingleDimCIS();
  else if(this->iMeth_ == RPA || this->iMeth_ == STAB) this->nSingleDimFOPP();
  else if(this->iMeth_ == PPRPA)                       this->nSingleDimPPRPA();
  else if(this->iMeth_ == PPATDA)                      this->nSingleDimPPATDA();
  else if(this->iMeth_ == PPCTDA)                      this->nSingleDimPPCTDA();
  else 
    CErr("PSCF Method " + std::to_string(this->iMeth_) + " NYI",this->fileio_->out);
} // initMeth

template<>
void SDResponse<double>::checkValid(){
  if(this->nSek_ == 0)
    CErr("Specification of zero desired roots is not acceptable",
         this->fileio_->out);
  if(this->nGuess_ == 0)
    CErr("Specification of guess vectors is not acceptable",
         this->fileio_->out);
  if(this->iMeth_ == 0)
    CErr("Invalid Method: SDResponse::iMeth_ = " + std::to_string(this->iMeth_),
         this->fileio_->out);
  if(this->nSingleDim_ == 0)
    CErr("Leading Dimenstion not defined for SDResponse::iMeth_ = " +
         std::to_string(this->iMeth_),this->fileio_->out);
  bool gtNSD   = this->nSingleDim_     < this->nGuess_;
  bool gtNSDd2 = (this->nSingleDim_/2) < this->nGuess_;
  bool gtADim;
  if(this->Ref_ == SingleSlater<double>::TCS)
    gtADim = this->nVV_SLT_ < this->nGuess_;
  else {
    if(this->iPPRPA_ == 0) gtADim = nVAVA_SLT_ < this->nGuess_;
    if(this->iPPRPA_ == 1) gtADim = nVAVB_     < this->nGuess_;
    if(this->iPPRPA_ == 2) gtADim = nVBVB_SLT_ < this->nGuess_;
  }

  if(this->nGuess_ < this->nSek_)
    CErr("Must specify more guess vectors than desired roots",this->fileio_->out);

  if((this->iMeth_ == PPRPA && gtADim ) ||
     (this->iMeth_ == RPA   && gtNSDd2) ||
     gtNSD)
    CErr("Number of guess roots exceeds number of number of possible roots");
} //checkValid


template<>
void SDResponse<double>::mpiSend(int toID,int tag) {
  //OOMPI_COMM_WORLD[toID].Send(this->nAtoms_,tag);
  //OOMPI_COMM_WORLD[toID].Send(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiSend(toID,tag);
}; //mpiSend
template<>
void SDResponse<double>::mpiRecv(int fromID,int tag) {
  //OOMPI_COMM_WORLD[fromID].Recv(this->nAtoms_,tag);
  //this->index_=new int[this->nAtoms_];
  //this->cart_ =new Matrix(3, this->nAtoms_, "Molecule");
  //OOMPI_COMM_WORLD[fromID].Recv(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiRecv(fromID,tag);
}; //mpiRecv

template<>
void SDResponse<double>::incorePPRPA(){
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


/*
  double Rmu = (*this->singleSlater_->epsA())(this->nOA_-1) + (*this->singleSlater_->epsA())(this->nOA_);
  Rmu /= 2;
  this->rMu_ = Rmu;
*/
  this->initRMu();
  double Rmu = this->rMu_;

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
  cout << "CORRECT" << endl;
  for(auto i = 0, ij = 0; i < this->nOA_; i++      )
  for(auto j = 0        ; j < i         ; j++, ij++){
    for(auto k = 0, kl = 0; k < this->nOA_; k++      )
    for(auto l = 0        ; l < k         ; l++, kl++){
      cout << "Single " << dijklA(i,j,k,l) << endl; 
      cout << "Double " << DijklA(i,j,k,l) << endl; 
      if(ij==kl) {
        cout << "EPS " << (*this->singleSlater_->epsA())(i) + (*this->singleSlater_->epsA())(j) - 2*Rmu << endl;;
      }
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
    cout << CAA << endl << endl;
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

  RealVecMap ATDATAAMap(ATDATAA.data(),VirSqAASLT);
  RealVecMap CTDATAAMap(CTDATAA.data(),OccSqAASLT);
  RealVecMap RPATAAMap(RPATAA.data(),VirSqAASLT + OccSqAASLT);
  RealVecMap ATDATABMap(ATDATAB.data(),VirSqAB);
  RealVecMap CTDATABMap(CTDATAB.data(),OccSqAB);
  RealVecMap RPATABMap(RPATAB.data(),VirSqAB + OccSqAB);
  RealMatrix ATDAAOAATmp(this->nBasis_,this->nBasis_);
  RealMatrix CTDAAOAATmp(this->nBasis_,this->nBasis_);
  RealMatrix RPAAOAATmp(this->nBasis_,this->nBasis_);
  RealMatrix ATDAAOABTmp(this->nBasis_,this->nBasis_);
  RealMatrix CTDAAOABTmp(this->nBasis_,this->nBasis_);
  RealMatrix RPAAOABTmp(this->nBasis_,this->nBasis_);

  this->iPPRPA_ = 0;
  this->iMeth_  = PPATDA;
  this->formAOTDen(ATDATAAMap,ATDAAOAATmp,ATDAAOAATmp);
  this->iMeth_  = PPCTDA;
  this->formAOTDen(CTDATAAMap,CTDAAOAATmp,CTDAAOAATmp);
  this->iMeth_  = PPRPA;
  this->formAOTDen(RPATAAMap,RPAAOAATmp,RPAAOAATmp);

  this->iPPRPA_ = 1;
  this->iMeth_  = PPATDA;
  this->formAOTDen(ATDATABMap,ATDAAOABTmp,ATDAAOABTmp);
  this->iMeth_  = PPCTDA;
  this->formAOTDen(CTDATABMap,CTDAAOABTmp,CTDAAOABTmp);
  this->iMeth_  = PPRPA;
  this->formAOTDen(RPATABMap,RPAAOABTmp,RPAAOABTmp);

  cout << "Checking AO Trans function A TDA (AA)... |T| = " << ATDAAOAATmp.norm() << " |R| = " << (ATDAAOAATmp - ATDAAOAA).norm() << endl;
  cout << "Checking AO Trans function C TDA (AA)... |T| = " << CTDAAOAATmp.norm() << " |R| = " << (CTDAAOAATmp - CTDAAOAA).norm() << endl;
  cout << "Checking AO Trans function RPA   (AA)... |T| = " << RPAAOAATmp.norm() << " |R| = " << (RPAAOAATmp - RPAAOAA).norm() << endl;
  cout << "Checking AO Trans function A TDA (AB)... |T| = " << ATDAAOABTmp.norm() << " |R| = " << (ATDAAOABTmp - ATDAAOAB).norm() << endl;
  cout << "Checking AO Trans function C TDA (AB)... |T| = " << CTDAAOABTmp.norm() << " |R| = " << (CTDAAOABTmp - CTDAAOAB).norm() << endl;
  cout << "Checking AO Trans function RPA   (AB)... |T| = " << RPAAOABTmp.norm() << " |R| = " << (RPAAOABTmp - RPAAOAB).norm() << endl;

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

  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,ATDAAOAA,ATDAIXAOAA,ATDAAOAA,ATDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,CTDAAOAA,CTDAIXAOAA,CTDAAOAA,CTDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,RPAAOAA,RPAIXAOAA,RPAAOAA,RPAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,ATDAAOAB,ATDAIXAOAB,ATDAAOAB,ATDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,CTDAAOAB,CTDAIXAOAB,CTDAAOAB,CTDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,RPAAOAB,RPAIXAOAB,RPAAOAB,RPAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,ATDAAOSing,ATDAIXAOSing,ATDAAOSing,ATDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,CTDAAOSing,CTDAIXAOSing,CTDAAOSing,CTDAIXAOAA);
  this->singleSlater_->aointegrals()->twoEContractDirect(true,false,false,true,false,RPAAOSing,RPAIXAOSing,RPAAOSing,RPAIXAOAA);


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
    ATDAAXMOAA(ab) = ATDAIXMOAA(this->nOA_+a,this->nOA_+b)
                     + ATDATAA(ab) *
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
//cout << "Checking A TDA (SA) AX... |AX| = " << ATDAAXMOSing.norm() << " |R| = " << (ATDAAXMOSing - ASing*ATDATSing).norm() << endl;
//cout << "Checking C TDA (SA) AX... |AX| = " << CTDAAXMOSing.norm() << " |R| = " << (CTDAAXMOSing - CSing*CTDATSing).norm() << endl;
//cout << "Checking RPA   (SA) AX... |AX| = " << RPAAXMOSing.norm() << " |R| = " << (RPAAXMOSing - FullSing*RPATSing).norm() << endl;
//for(auto i = 0 ; i < RPAAXMOAB.size() ; i++) cout << RPAAXMOAB(i) << " " << (FullAB*RPATAB)(i) << endl;
//


  Eigen::VectorXd ATDAAXAAMOTmp(VirSqAASLT); 
  Eigen::VectorXd CTDAAXAAMOTmp(OccSqAASLT); 
  Eigen::VectorXd RPAAXAAMOTmp(VirSqAASLT + OccSqAASLT); 
  Eigen::VectorXd ATDAAXABMOTmp(VirSqAB); 
  Eigen::VectorXd CTDAAXABMOTmp(OccSqAB); 
  Eigen::VectorXd RPAAXABMOTmp(VirSqAB + OccSqAB); 

  RealVecMap ATDAAXAAMOTmpMap(ATDAAXAAMOTmp.data(),VirSqAASLT); 
  RealVecMap CTDAAXAAMOTmpMap(CTDAAXAAMOTmp.data(),OccSqAASLT); 
  RealVecMap RPAAXAAMOTmpMap(RPAAXAAMOTmp.data(),VirSqAASLT + OccSqAASLT); 
  RealVecMap ATDAAXABMOTmpMap(ATDAAXABMOTmp.data(),VirSqAB); 
  RealVecMap CTDAAXABMOTmpMap(CTDAAXABMOTmp.data(),OccSqAB); 
  RealVecMap RPAAXABMOTmpMap(RPAAXABMOTmp.data(),VirSqAB + OccSqAB); 
/*
  this->iPPRPA_ = 0;
  this->iMeth_  = PPATDA;
  this->formMOTDen(ATDAAXAAMOTmpMap,ATDAIXAOAA,ATDAIXAOAA);
  this->scaleDagPPRPA(true,ATDATAAMap,ATDAAXAAMOTmpMap);
  this->iMeth_  = PPCTDA;
  this->formMOTDen(CTDAAXAAMOTmpMap,CTDAIXAOAA,CTDAIXAOAA);
  this->scaleDagPPRPA(true,CTDATAAMap,CTDAAXAAMOTmpMap);
  this->iMeth_  = PPRPA;
  this->formMOTDen(RPAAXAAMOTmpMap,RPAIXAOAA,RPAIXAOAA);
  this->scaleDagPPRPA(true,RPATAAMap,RPAAXAAMOTmpMap);

  this->iPPRPA_ = 1;
  this->iMeth_  = PPATDA;
  this->formMOTDen(ATDAAXABMOTmpMap,ATDAIXAOAB,ATDAIXAOAB);
  this->scaleDagPPRPA(true,ATDATABMap,ATDAAXABMOTmpMap);
  this->iMeth_  = PPCTDA;
  this->formMOTDen(CTDAAXABMOTmpMap,CTDAIXAOAB,CTDAIXAOAB);
  this->scaleDagPPRPA(true,CTDATABMap,CTDAAXABMOTmpMap);
  this->iMeth_  = PPRPA;
  this->formMOTDen(RPAAXABMOTmpMap,RPAIXAOAB,RPAIXAOAB);
  this->scaleDagPPRPA(true,RPATABMap,RPAAXABMOTmpMap);

  cout << "Checking MO Trans function A TDA (AA)... |AX| = " << ATDAAXAAMOTmp.norm() << " |R| = " << (ATDAAXMOAA - ATDAAXAAMOTmp).norm() << endl;
  cout << "Checking MO Trans function C TDA (AA)... |AX| = " << CTDAAXAAMOTmp.norm() << " |R| = " << (CTDAAXMOAA - CTDAAXAAMOTmp).norm() << endl;
  cout << "Checking MO Trans function RPA   (AA)... |AX| = " << RPAAXAAMOTmp.norm() << " |R| = " << (RPAAXMOAA - RPAAXAAMOTmp).norm() << endl;
  cout << "Checking MO Trans function A TDA (AB)... |AX| = " << ATDAAXABMOTmp.norm() << " |R| = " << (ATDAAXMOAB - ATDAAXABMOTmp).norm() << endl;
  cout << "Checking MO Trans function C TDA (AB)... |AX| = " << CTDAAXABMOTmp.norm() << " |R| = " << (CTDAAXMOAB - CTDAAXABMOTmp).norm() << endl;
  cout << "Checking MO Trans function RPA   (AB)... |AX| = " << RPAAXABMOTmp.norm() << " |R| = " << (RPAAXMOAB - RPAAXABMOTmp).norm() << endl;
*/

  Eigen::VectorXd ATDAAXAAMOTmp2(VirSqAASLT); 
  Eigen::VectorXd CTDAAXAAMOTmp2(OccSqAASLT); 
  Eigen::VectorXd RPAAXAAMOTmp2(VirSqAASLT + OccSqAASLT); 
  Eigen::VectorXd ATDAAXABMOTmp2(VirSqAB); 
  Eigen::VectorXd CTDAAXABMOTmp2(OccSqAB); 
  Eigen::VectorXd RPAAXABMOTmp2(VirSqAB + OccSqAB); 

  RealCMMap ATDAAXAAMOTmp2Map(ATDAAXAAMOTmp2.data(),VirSqAASLT,1); 
  RealCMMap CTDAAXAAMOTmp2Map(CTDAAXAAMOTmp2.data(),OccSqAASLT,1); 
  RealCMMap RPAAXAAMOTmp2Map(RPAAXAAMOTmp2.data(),VirSqAASLT + OccSqAASLT,1); 
  RealCMMap ATDAAXABMOTmp2Map(ATDAAXABMOTmp2.data(),VirSqAB,1); 
  RealCMMap CTDAAXABMOTmp2Map(CTDAAXABMOTmp2.data(),OccSqAB,1); 
  RealCMMap RPAAXABMOTmp2Map(RPAAXABMOTmp2.data(),VirSqAB + OccSqAB,1); 
  RealCMMap ATDATAAMap2(ATDATAA.data(),VirSqAASLT,1);
  RealCMMap CTDATAAMap2(CTDATAA.data(),OccSqAASLT,1);
  RealCMMap RPATAAMap2(RPATAA.data(),VirSqAASLT + OccSqAASLT,1);
  RealCMMap ATDATABMap2(ATDATAB.data(),VirSqAB,1);
  RealCMMap CTDATABMap2(CTDATAB.data(),OccSqAB,1);
  RealCMMap RPATABMap2(RPATAB.data(),VirSqAB + OccSqAB,1);

  this->iPPRPA_ = 0;
  this->iMeth_  = PPATDA;
  this->nSingleDim_ = VirSqAASLT;
  this->haveDag_ = false;
  this->formGuess();
  this->formRM4(ATDATAAMap2,ATDAAXAAMOTmp2Map,ATDAAXAAMOTmp2Map);
  this->iMeth_  = PPCTDA;
  this->nSingleDim_ = OccSqAASLT;
  this->haveDag_ = false;
  this->formRM4(CTDATAAMap2,CTDAAXAAMOTmp2Map,CTDAAXAAMOTmp2Map);
  this->iMeth_  = PPRPA;
  this->nSingleDim_ = VirSqAASLT + OccSqAASLT;
  this->haveDag_ = false;
  this->formRM4(RPATAAMap2,RPAAXAAMOTmp2Map,RPAAXAAMOTmp2Map);

  this->iPPRPA_ = 1;
  this->iMeth_  = PPATDA;
  this->nSingleDim_ = VirSqAB;
  this->haveDag_ = false;
  this->formRM4(ATDATABMap2,ATDAAXABMOTmp2Map,ATDAAXABMOTmp2Map);
  this->iMeth_  = PPCTDA;
  this->nSingleDim_ = OccSqAB;
  this->haveDag_ = false;
  this->formRM4(CTDATABMap2,CTDAAXABMOTmp2Map,CTDAAXABMOTmp2Map);
  this->iMeth_  = PPRPA;
  this->nSingleDim_ = VirSqAB + OccSqAB;
  this->haveDag_ = false;
  this->formRM4(RPATABMap2,RPAAXABMOTmp2Map,RPAAXABMOTmp2Map);

  cout << "Checking Full AX function A TDA (AA)... |AX| = " << ATDAAXAAMOTmp2.norm() << " |R| = " << (ATDAAXMOAA - ATDAAXAAMOTmp2).norm() << endl;
  cout << "Checking Full AX function A TDA (AA)... |AX| = " << CTDAAXAAMOTmp2.norm() << " |R| = " << (CTDAAXMOAA - CTDAAXAAMOTmp2).norm() << endl;
  cout << "Checking Full AX function A TDA (AA)... |AX| = " << RPAAXAAMOTmp2.norm() << " |R| = " << (RPAAXMOAA - RPAAXAAMOTmp2).norm() << endl;
  cout << "Checking Full AX function A TDA (AB)... |AX| = " << ATDAAXABMOTmp2.norm() << " |R| = " << (ATDAAXMOAB - ATDAAXABMOTmp2).norm() << endl;
  cout << "Checking Full AX function A TDA (AB)... |AX| = " << CTDAAXABMOTmp2.norm() << " |R| = " << (CTDAAXMOAB - CTDAAXABMOTmp2).norm() << endl;
  cout << "Checking Full AX function A TDA (AB)... |AX| = " << RPAAXABMOTmp2.norm() << " |R| = " << (RPAAXMOAB - RPAAXABMOTmp2).norm() << endl;
} //incorePPRPA

template<>
void SDResponse<double>::formRM(){
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
  if (this->Ref_ == SingleSlater<double>::RHF)
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

  if (this->Ref_ != SingleSlater<double>::RHF)
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
  cout << Aud << endl << endl << Auod << endl;
  // Build the ABBA matrix
  ABBA.block(0,0,nOVA+nOVB,nOVA+nOVB) = A;
  ABBA.block(nOVA+nOVB,nOVA+nOVB,nOVA+nOVB,nOVA+nOVB) = -A;
  ABBA.block(0,nOVA+nOVB,nOVA+nOVB,nOVA+nOVB) = B;
  ABBA.block(nOVA+nOVB,0,nOVA+nOVB,nOVA+nOVB) = -B;

  // CIS routine
  // Diagonalize the A matrix
  Eigen::SelfAdjointEigenSolver<RealMatrix> CIS;
  CIS.compute(A);
  cout << CIS.eigenvalues() << endl << endl;;
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
//cout <<"ABBA" << ABBA << endl;
/*
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
*/
} //formRM

template<>
RealMatrix SDResponse<double>::formRM2(RealMatrix &XMO){
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
    if (this->Ref_ == SingleSlater<double>::RHF)
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
    if(this->singleSlater_->aointegrals()->integralAlgorithm == AOIntegrals::DIRECT)
      this->singleSlater_->aointegrals()->twoEContractDirect(false,false,false,false,false,XAAO,IXA,XBAO,IXB);
    else if(this->singleSlater_->aointegrals()->integralAlgorithm == AOIntegrals::INCORE)
      this->singleSlater_->aointegrals()->twoEContractN4(false,false,true,false,false,XAAO,IXA,XBAO,IXB);
    

    /*
     *  AX(s)   += IXMO(s)
     *  IXMO(s) =  Cv(s)**H * IXAO(s) * Co(s)
     */  
    AXA = 
      this->singleSlater_->moA()->block(0,this->nOA_,this->nBasis_,this->nVA_).adjoint()*
      IXA*
      this->singleSlater_->moA()->block(0,0,this->nBasis_,this->nOA_);
    if (this->Ref_ == SingleSlater<double>::RHF)
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
} //formRM2

template<>
void SDResponse<double>::incoreCIS(){
  RealMatrix A;
  if(this->Ref_ == SingleSlater<double>::TCS) A = RealMatrix(this->nOV_,this->nOV_);
  else A = RealMatrix(this->nOAVA_ + this->nOBVB_, this->nOAVA_ + this->nOBVB_);
  this->mointegrals_->formIAJB(false); 
  this->mointegrals_->formIJAB(false); 


  if(this->Ref_ == SingleSlater<double>::TCS){
    for(auto a = 0, ia = 0        ; a < this->nV_; a++)
    for(auto i = 0; i < this->nO_; i++, ia++      )
    for(auto b = 0, jb = 0        ; b < this->nV_; b++)
    for(auto j = 0; j < this->nO_; j++, jb++      ){
      A(ia,jb) = this->mointegrals_->IAJB(i,a,j,b) - this->mointegrals_->IJAB(i,j,a,b);
      if(ia==jb) A(ia,jb) += (*this->singleSlater_->epsA())(a+this->nO_) -
                             (*this->singleSlater_->epsA())(i);
    }
  }
 
  Eigen::SelfAdjointEigenSolver<RealMatrix> ES(A);
  cout << ES.eigenvalues()<< endl;
//RealMatrix Vec = ES.eigenvectors().col(0);
//RealMatrix AX(this->nOV_,1);
//RealCMMap VMap(Vec.data(),this->nOV_,1);
//RealCMMap AXMap(AX.data(),this->nOV_,1);
//cout << "AX Before" << endl << AX << endl;
//this->formRM3(VMap,AXMap,AXMap);
//cout << "NOV " << this->nOV_ << endl;
//cout << "VEC" << endl << Vec << endl;
//cout << "VEC An" << endl << A*Vec << endl;
//cout << "VEC An" << endl << A*VMap << endl;
//cout << "AX Gen" << endl << AX << endl;
/*
  RealVecMap VcMap(Vec.data(),this->nOV_);
  RealMatrix TAO(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  RealMatrix SDA = (*this->singleSlater_->aointegrals()->overlap_) * (*this->singleSlater_->densityA());
  this->formAOTDen(VcMap,TAO,TAO);
  RealMatrix Comm = TAO*SDA + SDA.adjoint() * TAO;
  RealMatrix GComm(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  this->singleSlater_->aointegrals()->multTwoEContractDirect(1,false,false,false,true,Comm,
    GComm,Comm,GComm);
  RealMatrix AXAO = (*this->singleSlater_->fockA()) * Comm * (*this->singleSlater_->aointegrals()->overlap_) - (*this->singleSlater_->aointegrals()->overlap_) * Comm * (*this->singleSlater_->fockA()) + GComm * SDA - SDA * GComm;
  RealVecMap AXMap(AX.data(),this->nOV_);
  this->formMOTDen(AXMap,AXAO,AXAO);
  cout << "AX Gen 2" << endl << AXMap << endl;
*/

}
template<>
void SDResponse<double>::incoreRPA(){
  RealMatrix A,B,ABBA;
  if(this->Ref_ == SingleSlater<double>::TCS) {
    A = RealMatrix(this->nOV_,this->nOV_);
    B = RealMatrix(this->nOV_,this->nOV_);
    ABBA = RealMatrix(2*this->nOV_,2*this->nOV_);
  } else {
    A = RealMatrix(this->nOAVA_ + this->nOBVB_, this->nOAVA_ + this->nOBVB_);
    B = RealMatrix(this->nOAVA_ + this->nOBVB_, this->nOAVA_ + this->nOBVB_);
    ABBA = RealMatrix(2*(this->nOAVA_ + this->nOBVB_), 2*(this->nOAVA_ + this->nOBVB_));
  }
  this->mointegrals_->formIAJB(false); 
  this->mointegrals_->formIJAB(false); 


  if(this->Ref_ == SingleSlater<double>::TCS){
    for(auto a = 0, ia = 0        ; a < this->nV_; a++)
    for(auto i = 0; i < this->nO_; i++, ia++      )
    for(auto b = 0, jb = 0        ; b < this->nV_; b++)
    for(auto j = 0; j < this->nO_; j++, jb++      ){
      A(ia,jb) = this->mointegrals_->IAJB(i,a,j,b) - this->mointegrals_->IJAB(i,j,a,b);
      if(ia==jb) A(ia,jb) += (*this->singleSlater_->epsA())(a+this->nO_) -
                             (*this->singleSlater_->epsA())(i);
      B(ia,jb) = this->mointegrals_->IAJB(i,a,j,b) - this->mointegrals_->IAJB(i,b,j,a);
    }
  }
  ABBA.block(0,0,this->nOV_,this->nOV_) = A;
  ABBA.block(this->nOV_,0,this->nOV_,this->nOV_) = B;
  ABBA.block(0,this->nOV_,this->nOV_,this->nOV_) = B;
  ABBA.block(this->nOV_,this->nOV_,this->nOV_,this->nOV_) = A;
 
  Eigen::SelfAdjointEigenSolver<RealMatrix> ES(ABBA);
  cout << std::fixed << std::setprecision(12);
  cout << ES.eigenvalues()<< endl;

  RealCMMatrix TStab = ES.eigenvectors().real();
  RealCMMap TStabMap(TStab.data(),2*this->nOV_,2*this->nOV_);
  RealCMMatrix ATStab(2*this->nOV_,2*this->nOV_);
  RealCMMap ATStabMap(ATStab.data(),2*this->nOV_,2*this->nOV_);

//this->formRM3(TStabMap,ATStabMap,ATStabMap);
//prettyPrint(cout,ABBA*TStab - ATStabMap,"DIFF");
  

  ABBA.block(0,0,this->nOV_,this->nOV_) = A;
  ABBA.block(this->nOV_,0,this->nOV_,this->nOV_) = -B;
  ABBA.block(0,this->nOV_,this->nOV_,this->nOV_) = B;
  ABBA.block(this->nOV_,this->nOV_,this->nOV_,this->nOV_) = -A;
  Eigen::EigenSolver<RealMatrix> EA(ABBA);
  cout << "LR EIG" << endl;
  cout << endl << EA.eigenvalues() << endl;
  RealMatrix AmB = A-B;
  RealMatrix ApB = A+B;

//EA.compute(ApB);
//cout << endl << EA.eigenvalues() << endl;
//EA.compute(AmB);
//cout << endl << EA.eigenvalues() << endl;


  
  ABBA.block(0,0,this->nOV_,this->nOV_) = A;
  ABBA.block(this->nOV_,0,this->nOV_,this->nOV_) = B;
  ABBA.block(0,this->nOV_,this->nOV_,this->nOV_) = B;
  ABBA.block(this->nOV_,this->nOV_,this->nOV_,this->nOV_) = A;
    if(!this->haveDag_) this->getDiag();
    RealCMMatrix Vec(2*this->nOV_,3);
    VectorXd Eig(3);
    RealCMMatrix ACM = ABBA;
    QuasiNewton<double> davA(false,true,3,&ACM,this->rmDiag_.get(),&Vec,&Eig);
    davA.run(cout);


//RealMatrix Vec = ES.eigenvectors().col(0);
//RealMatrix AX(this->nOV_,1);
//RealCMMap VMap(Vec.data(),this->nOV_,1);
//RealCMMap AXMap(AX.data(),this->nOV_,1);
//cout << "AX Before" << endl << AX << endl;
//this->formRM3(VMap,AXMap,AXMap);
//cout << "NOV " << this->nOV_ << endl;
//cout << "VEC" << endl << Vec << endl;
//cout << "VEC An" << endl << A*Vec << endl;
//cout << "VEC An" << endl << A*VMap << endl;
//cout << "AX Gen" << endl << AX << endl;
/*
  RealVecMap VcMap(Vec.data(),this->nOV_);
  RealMatrix TAO(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  RealMatrix SDA = (*this->singleSlater_->aointegrals()->overlap_) * (*this->singleSlater_->densityA());
  this->formAOTDen(VcMap,TAO,TAO);
  RealMatrix Comm = TAO*SDA + SDA.adjoint() * TAO;
  RealMatrix GComm(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  this->singleSlater_->aointegrals()->multTwoEContractDirect(1,false,false,false,true,Comm,
    GComm,Comm,GComm);
  RealMatrix AXAO = (*this->singleSlater_->fockA()) * Comm * (*this->singleSlater_->aointegrals()->overlap_) - (*this->singleSlater_->aointegrals()->overlap_) * Comm * (*this->singleSlater_->fockA()) + GComm * SDA - SDA * GComm;
  RealVecMap AXMap(AX.data(),this->nOV_);
  this->formMOTDen(AXMap,AXAO,AXAO);
  cout << "AX Gen 2" << endl << AXMap << endl;
*/

}

//void SDResponse::incorePPRPAnew(){
//  cout << "BEFORE ABCD" << endl;
//  this->mointegrals_->formABCD(false); 
////cout << "BEFORE IJKL" << endl;
////this->mointegrals_->formIJKL(true); 
//  RealMatrix A  ,B  ,C;
//  RealMatrix AAA,AAB,ABB;
//  RealMatrix CAA,CAB,CBB;
//  RealMatrix BAA,BAB,BBB;
//  // Eigensolvers
//  Eigen::SelfAdjointEigenSolver<RealMatrix> ES;
//  Eigen::EigenSolver<RealMatrix> EA;
//  // Eigen values
//  Eigen::VectorXd ATDAEAA;
//  Eigen::VectorXd ATDAEAB;
//  Eigen::VectorXd ATDAEBB;
//  Eigen::VectorXd ATDAESing;
//  Eigen::VectorXd CTDAEAA;
//  Eigen::VectorXd CTDAEAB;
//  Eigen::VectorXd CTDAEBB;
//  Eigen::VectorXd CTDAESing;
//  Eigen::VectorXd RPAEAA;
//  Eigen::VectorXd RPAEAB;
//  Eigen::VectorXd RPAEBB;
//  Eigen::VectorXd RPAESing;
//  // Eigen Vectors
//  Eigen::VectorXd ATDATAA;
//  Eigen::VectorXd ATDATAB;
//  Eigen::VectorXd ATDATBB;
//  Eigen::VectorXd ATDATSing;
//  Eigen::VectorXd CTDATAA;
//  Eigen::VectorXd CTDATAB;
//  Eigen::VectorXd CTDATBB;
//  Eigen::VectorXd CTDATSing;
//  Eigen::VectorXd RPATAA;
//  Eigen::VectorXd RPATAB;
//  Eigen::VectorXd RPATBB;
//  Eigen::VectorXd RPATSing;
//
//  AAA = RealMatrix(this->nVAVA_SLT_,this->nVAVA_SLT_);
//  AAB = RealMatrix(this->nVAVB_,this->nVAVB_);
////CAA = RealMatrix(this->nOAOA_SLT_,this->nOAOA_SLT_);
////CAB = RealMatrix(this->nOAOB_,this->nOAOB_);
// 
//
//  this->initRMu();
//  double Rmu = 0.0;
//
//  if(this->Ref_ != SingleSlater<double>::TCS) {
////  for(auto a = 0, ab = 0; a < this->nV_; a++      )
////  for(auto b = 0        ; b < a        ; b++, ab++){
////    for(auto c = 0, cd = 0; c < this->nV_; c++      )
////    for(auto d = 0        ; d < c        ; d++, cd++){
////      A(ab,cd) = this->mointegrals_->ABCD(a,c,b,d) -
////                   this->mointegrals_->ABCD(a,d,b,c);
////      if(ab == cd) 
////        A(ab,cd) += (*this->singleSlater_->epsA())(a+this->nO_) + 
////                      (*this->singleSlater_->epsA())(b+this->nO_) - 
////                      2*Rmu;
////    }
////  }
//
//  } else {
//    for(auto a = 0, ab = 0; a < this->nVA_; a++      )
//    for(auto b = 0        ; b < a         ; b++, ab++){
//      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
//      for(auto d = 0        ; d < c         ; d++, cd++){
//        AAA(ab,cd) = this->mointegrals_->ABCD(a,c,b,d,"AAAA") -
//                     this->mointegrals_->ABCD(a,d,b,c,"AAAA");
//        if(ab == cd) 
//          AAA(ab,cd) += (*this->singleSlater_->epsA())(a+this->nOA_) + 
//                        (*this->singleSlater_->epsA())(b+this->nOA_) - 
//                        2*Rmu;
//      }
//    }
//    for(auto a = 0, ab = 0; a < this->nVA_; a++      )
//    for(auto b = 0        ; b < this->nVB_; b++, ab++){
//      for(auto c = 0, cd = 0; c < this->nVA_; c++      )
//      for(auto d = 0        ; d < this->nVB_; d++, cd++){
//        AAB(ab,cd) = this->mointegrals_->ABCD(a,c,b,d,"AABB") -
//                     this->mointegrals_->ABCD(a,d,b,c,"AABB");
//        if(ab == cd) 
//          AAB(ab,cd) += (*this->singleSlater_->epsA())(a+this->nOA_) + 
//                        (*this->singleSlater_->epsA())(b+this->nOA_) - 
//                        2*Rmu;
//      }
//    }
//    // Full Diagonalization
//    // Pure Spin
//    ES.compute(AAA);
//      ATDAEAA = -ES.eigenvalues();
//      std::sort(ATDAEAA.data(),ATDAEAA.data()+ATDAEAA.size());
//      ATDAEAA = -ATDAEAA;
//      ATDATAA = ES.eigenvectors().col(0);
//  /*
//    ES.compute(-CAA);
//      CTDAEAA = -ES.eigenvalues();
//      std::sort(CTDAEAA.data(),CTDAEAA.data()+CTDAEAA.size());
//      CTDAEAA = -CTDAEAA;
//      CTDATAA = ES.eigenvectors().col(0);
//    EA.compute(FullAA);
//      RPAEAA = -EA.eigenvalues().real();
//      std::sort(RPAEAA.data(),RPAEAA.data()+RPAEAA.size());
//      RPAEAA = -RPAEAA;
//      RPATAA = EA.eigenvectors().col(0).real();
//  */
//      
//    // Mixed Spin
//    ES.compute(AAB);
//      ATDAEAB = -ES.eigenvalues();
//      std::sort(ATDAEAB.data(),ATDAEAB.data()+ATDAEAB.size());
//      ATDAEAB = -ATDAEAB;
//      ATDATAB = ES.eigenvectors().col(0);
//      cout << ATDATAB << endl;
//  /*
//    ES.compute(-CAB);
//      CTDAEAB = -ES.eigenvalues();
//      std::sort(CTDAEAB.data(),CTDAEAB.data()+CTDAEAB.size());
//      CTDAEAB = -CTDAEAB;
//      CTDATAB = ES.eigenvectors().col(0);
//    EA.compute(FullAB);
//      RPAEAB = -EA.eigenvalues().real();
//      std::sort(RPAEAB.data(),RPAEAB.data()+RPAEAB.size());
//      RPAEAB = -RPAEAB;
//      RPATAB = EA.eigenvectors().col(0).real();
//  */
//    Eigen::IOFormat HeavyFmt(8);
//    cout << "A TDA (AA) Eigenvalues:" << endl;
//    cout << (ATDAEAA*phys.eVPerHartree).format(HeavyFmt) << endl << endl;
//    cout << "A TDA (AB) Eigenvalues:" << endl;
//    cout << (ATDAEAB).format(HeavyFmt) << endl << endl;
//  /*
//    cout << "A TDA (SAS) Eigenvalues:" << endl;
//    cout << ATDAESing.format(HeavyFmt) << endl << endl;
//    cout << "C TDA (AA) Eigenvalues:" << endl;
//    cout << CTDAEAA.format(HeavyFmt) << endl << endl;
//    cout << "C TDA (AB) Eigenvalues:" << endl;
//    cout << CTDAEAB.format(HeavyFmt) << endl << endl;
//    cout << "C TDA (SAS) Eigenvalues:" << endl;
//    cout << CTDAESing.format(HeavyFmt) << endl << endl;
//    cout << "RPA (AA) Eigenvalues:" << endl;
//    cout << RPAEAA.format(HeavyFmt)  << endl << endl;
//    cout << "RPA (AB) Eigenvalues:" << endl;
//    cout << RPAEAB.format(HeavyFmt)  << endl << endl;
//    cout << "RPA (AB) Eigenvalues:" << endl;
//    cout << RPAEAB.format(HeavyFmt)  << endl << endl;
//    cout << "RPA (SAS) Eigenvalues:" << endl;
//    cout << RPAESing.format(HeavyFmt)  << endl << endl;
//  */
//  }
//  
//}
template<>
void SDResponse<double>::incorePPRPAnew(){
  RealMatrix A,B,C;
  RealMatrix AAA,BAA,CAA;
  RealMatrix AAB,BAB,CAB;
  RealMatrix ABB,BBB,CBB;
  RealMatrix Full;
  RealMatrix FullAA, FullAB, FullBB;

  if(this->Ref_ == SingleSlater<double>::TCS){
    A = RealMatrix(this->nVV_SLT_,this->nVV_SLT_);
    B = RealMatrix(this->nVV_SLT_,this->nOO_SLT_);
    C = RealMatrix(this->nOO_SLT_,this->nOO_SLT_);
    Full = RealMatrix(this->nVV_SLT_+this->nOO_SLT_,this->nVV_SLT_+this->nOO_SLT_);
  } else {
    AAA = RealMatrix(this->nVAVA_SLT_,this->nVAVA_SLT_);
    BAA = RealMatrix(this->nVAVA_SLT_,this->nOAOA_SLT_);
    CAA = RealMatrix(this->nOAOA_SLT_,this->nOAOA_SLT_);
    AAB = RealMatrix(this->nVAVB_,this->nVAVB_);
    BAB = RealMatrix(this->nVAVB_,this->nOAOB_);
    CAB = RealMatrix(this->nOAOB_,this->nOAOB_);
    FullAA =RealMatrix(this->nVAVA_SLT_+this->nOAOA_SLT_,this->nVAVA_SLT_+this->nOAOA_SLT_);
    FullAB = RealMatrix(this->nVAVB_+this->nOAOB_,this->nVAVB_+this->nOAOB_);
    if(!this->singleSlater_->isClosedShell){
      ABB = RealMatrix(this->nVBVB_SLT_,this->nVBVB_SLT_);
      BBB = RealMatrix(this->nVBVB_SLT_,this->nOBOB_SLT_);
      CBB = RealMatrix(this->nOBOB_SLT_,this->nOBOB_SLT_);
      FullBB =
        RealMatrix(this->nVBVB_SLT_+this->nOBOB_SLT_,this->nVBVB_SLT_+this->nOBOB_SLT_);
    }
  }
  this->mointegrals_->formABCD(false);
  this->mointegrals_->formIAJB(false);
  this->mointegrals_->formIJKL(false);

  this->initRMu();
  double Rmu = this->rMu_;

  if(this->Ref_ == SingleSlater<double>::TCS) {
    for(auto a = 0, ab = 0; a < this->nV_; a++      )
    for(auto b = 0        ; b < a        ; b++, ab++)
    for(auto c = 0, cd = 0; c < this->nV_; c++      )
    for(auto d = 0        ; d < c        ; d++, cd++){
      A(ab,cd) = this->mointegrals_->ABCD(a,c,b,d) - this->mointegrals_->ABCD(a,d,b,c);
      if(ab == cd) A(ab,cd) +=
        (*this->singleSlater_->epsA())(a+this->nO_) + 
        (*this->singleSlater_->epsA())(b+this->nO_) - 2*Rmu; 
    }

    for(auto i = 0, ij = 0; i < this->nO_; i++      )
    for(auto j = 0        ; j < i        ; j++, ij++)
    for(auto k = 0, kl = 0; k < this->nO_; k++      )
    for(auto l = 0        ; l < k        ; l++, kl++){
      C(ij,kl) = this->mointegrals_->IJKL(i,k,j,l) - this->mointegrals_->IJKL(i,l,j,k);
      if(ij == kl) C(ij,kl) -=
        (*this->singleSlater_->epsA())(i) + 
        (*this->singleSlater_->epsA())(j) - 2*Rmu; 
    }

    for(auto a = 0, ab = 0; a < this->nV_; a++      )
    for(auto b = 0        ; b < a        ; b++, ab++)
    for(auto i = 0, ij = 0; i < this->nO_; i++      )
    for(auto j = 0        ; j < i        ; j++, ij++){
      B(ab,ij) = this->mointegrals_->IAJB(i,a,j,b) - this->mointegrals_->IAJB(i,b,j,a);
    }
 
    Eigen::SelfAdjointEigenSolver<RealMatrix> ES;
    ES.compute(A);
    Eigen::VectorXd EATDA = ES.eigenvalues().real();

    double valATDA = EATDA(0);
    EATDA = -EATDA;
    std::sort(EATDA.data(),EATDA.data()+EATDA.size());
    EATDA = -EATDA;
    VectorXd LowATDA = EATDA;
    for(auto i = 0; i < LowATDA.size(); i++) LowATDA(i) = valATDA;
    cout << std::fixed << std::setprecision(12);
    cout << EATDA-LowATDA << endl;

    if(!this->haveDag_) this->getDiag();
    RealCMMatrix TATDA = ES.eigenvectors().real();
    RealCMMap TATDAMap(TATDA.data(),this->nVV_SLT_,this->nVV_SLT_);
    RealCMMatrix ATATDA(this->nVV_SLT_,this->nVV_SLT_);
    RealCMMap ATATDAMap(ATATDA.data(),this->nVV_SLT_,this->nVV_SLT_);
/*
    for(auto iSt = 0; iSt < this->nVV_SLT_;iSt++)
    for(auto a = 0, ab = 0; a < this->nV_; a++      )
    for(auto b = 0        ; b < a        ; b++, ab++)
    for(auto c = 0, cd = 0; c < this->nV_; c++      )
    for(auto d = 0        ; d < c        ; d++, cd++){
      ATATDA(ab,iSt) += 
        (this->mointegrals_->ABCD(a,c,b,d) - this->mointegrals_->ABCD(a,d,b,c))*
        TATDA(cd,iSt);
      if(ab==cd)
        ATATDA(ab,iSt) += (*this->rmDiag_)(ab)*TATDA(ab,iSt);
    }
*/
   
    this->formRM4(TATDAMap,ATATDAMap,ATATDAMap);
//  prettyPrint(cout,A*TATDA - ATATDAMap,"DIFF");
/*
    if(!this->haveDag_) this->getDiag();
    RealCMMatrix Vec(this->nVV_SLT_,8);
    VectorXd Eig(8,1);
    RealCMMatrix ACM = A;
    QuasiNewton<double> davA(8,&ACM,this->rmDiag_.get(),&Vec,&Eig);
    davA.run(cout);
*/
    
    
    cout << endl << endl;
    ES.compute(-C);
    Eigen::VectorXd ECTDA = ES.eigenvalues().real();
    ECTDA = -ECTDA;
    std::sort(ECTDA.data(),ECTDA.data()+ECTDA.size());
    ECTDA = -ECTDA;
    cout << ECTDA << endl;

    Full.block(0,0,this->nVV_SLT_,this->nVV_SLT_) = A;
    Full.block(this->nVV_SLT_,this->nVV_SLT_,this->nOO_SLT_,this->nOO_SLT_) = -C;
    Full.block(0,this->nVV_SLT_,this->nVV_SLT_,this->nOO_SLT_) = B;
    Full.block(this->nVV_SLT_,0,this->nOO_SLT_,this->nVV_SLT_) = -B.adjoint();
    Eigen::EigenSolver<RealMatrix> EA;
    EA.compute(Full);
    Eigen::VectorXd ERPA = EA.eigenvalues().real();
    ERPA = -ERPA;
    std::sort(ERPA.data(),ERPA.data()+ERPA.size());
    ERPA = -ERPA;
    double valRPA = ERPA(this->nVV_SLT_-1);
    VectorXd LowRPA = ERPA;
    for(auto i = 0; i < LowRPA.size(); i++) LowRPA(i) = valRPA;
    cout << endl << endl << ERPA-LowRPA << endl;
    
  } else {
    for(auto a = 0, ab = 0; a < this->nVA_; a++      )
    for(auto b = 0        ; b < a        ; b++, ab++)
    for(auto c = 0, cd = 0; c < this->nVA_; c++      )
    for(auto d = 0        ; d < c        ; d++, cd++){
      AAA(ab,cd) = this->mointegrals_->ABCD(a,c,b,d,"AAAA") - this->mointegrals_->ABCD(a,d,b,c,"AAAA");
      if(ab == cd) AAA(ab,cd) +=
        (*this->singleSlater_->epsA())(a+this->nOA_) + 
        (*this->singleSlater_->epsA())(b+this->nOA_) - 2*Rmu; 
    }
    for(auto i = 0, ij = 0; i < this->nOA_; i++      )
    for(auto j = 0        ; j < i        ; j++, ij++)
    for(auto k = 0, kl = 0; k < this->nOA_; k++      )
    for(auto l = 0        ; l < k        ; l++, kl++){
      CAA(ij,kl) = this->mointegrals_->IJKL(i,k,j,l,"AAAA") - this->mointegrals_->IJKL(i,l,j,k,"AAAA") ;;
      if(ij == kl) CAA(ij,kl) -=
        (*this->singleSlater_->epsA())(i) + 
        (*this->singleSlater_->epsA())(j) - 2*Rmu; 
    }
    for(auto a = 0, ab = 0; a < this->nVA_; a++      )
    for(auto b = 0        ; b < a        ; b++, ab++)
    for(auto i = 0, ij = 0; i < this->nOA_; i++      )
    for(auto j = 0        ; j < i        ; j++, ij++){
      BAA(ab,ij) = this->mointegrals_->IAJB(i,a,j,b,"AAAA") - this->mointegrals_->IAJB(i,b,j,a,"AAAA");
    }
 
    Eigen::SelfAdjointEigenSolver<RealMatrix> ES;
    ES.compute(AAA);
    Eigen::VectorXd EATDA = ES.eigenvalues().real();
    EATDA = -EATDA;
    std::sort(EATDA.data(),EATDA.data()+EATDA.size());
    EATDA = -EATDA; 
    cout << EATDA << endl;
 
    cout << endl << endl;
    ES.compute(-CAA);
    Eigen::VectorXd ECTDA = ES.eigenvalues().real();
    ECTDA = -ECTDA;
    std::sort(ECTDA.data(),ECTDA.data()+ECTDA.size());
    ECTDA = -ECTDA;
    cout << ECTDA << endl;
  }
  
}
} // namespace ChronusQ
