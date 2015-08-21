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
void SDResponse<dcomplex>::placeVOOV(const ComplexVecMap &TMOV, ComplexMatrix &TMOA, ComplexMatrix &TMOB){
  int iOff = this->nOAVA_ + this->nOBVB_;
  if(this->Ref_ == SingleSlater<dcomplex>::TCS) iOff = this->nOV_;

  if(this->Ref_ == SingleSlater<dcomplex>::TCS){
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
void SDResponse<dcomplex>::placeVVOO(const ComplexVecMap &TMOV, ComplexMatrix &TMO){
  bool doX = ( (this->iMeth_ == PPATDA) || (this->iMeth_ == PPRPA) );
  bool doY = ( (this->iMeth_ == PPCTDA) || (this->iMeth_ == PPRPA) );

  if(this->Ref_ == SingleSlater<dcomplex>::TCS){
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
void SDResponse<dcomplex>::retrvVOOV(ComplexVecMap &TMOV, const ComplexMatrix &TMOA, const ComplexMatrix &TMOB){
  int iOff = this->nOAVA_ + this->nOBVB_;
  if(this->Ref_ == SingleSlater<dcomplex>::TCS) iOff = this->nOV_;

  if(this->Ref_ == SingleSlater<dcomplex>::TCS) {
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
void SDResponse<dcomplex>::retrvVVOO(ComplexVecMap &TMOV, const ComplexMatrix &TMO){
  bool doX = ( (this->iMeth_ == PPATDA) || (this->iMeth_ == PPRPA) );
  bool doY = ( (this->iMeth_ == PPCTDA) || (this->iMeth_ == PPRPA) );

  if(this->Ref_ == SingleSlater<dcomplex>::TCS){
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
void SDResponse<dcomplex>::formAOTDen(const ComplexVecMap &TMOV, ComplexMatrix &TAOA, ComplexMatrix &TAOB){
  bool doVOOV = (this->iMeth_ == CIS || this->iMeth_ == RPA || this->iMeth_ == STAB); 
  bool doVVOO = (this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA);

  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  ComplexMatrix TMOA,TMOB;
  TMOA = ComplexMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);
  if(!doVVOO && this->Ref_ != SingleSlater<dcomplex>::TCS) 
    TMOB = ComplexMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);

  if(doVOOV) this->placeVOOV(TMOV,TMOA,TMOB);
  else if(doVVOO) this->placeVVOO(TMOV,TMOA);

  TAOA = (*this->singleSlater_->moA()) * TMOA * this->singleSlater_->moA()->adjoint();
  if(!doVVOO && this->Ref_ != SingleSlater<dcomplex>::TCS){
    if(this->Ref_ == SingleSlater<dcomplex>::RHF)
      TAOB = (*this->singleSlater_->moA()) * TMOB * this->singleSlater_->moA()->adjoint();
    else
      TAOB = (*this->singleSlater_->moB()) * TMOB * this->singleSlater_->moB()->adjoint();
  }
} //formAOTDen

template<>
void SDResponse<dcomplex>::formMOTDen(ComplexVecMap &TMOV, const ComplexMatrix &TAOA, const ComplexMatrix &TAOB){
  bool doVOOV = (this->iMeth_ == CIS || this->iMeth_ == RPA || this->iMeth_ == STAB); 
  bool doVVOO = (this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA);
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;

  ComplexMatrix TMOA,TMOB;
  TMOA = ComplexMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);
  if(!doVVOO && this->Ref_ != SingleSlater<dcomplex>::TCS) 
    TMOB = ComplexMatrix::Zero(NTCSxNBASIS,NTCSxNBASIS);

  TMOA = this->singleSlater_->moA()->adjoint() * TAOA * (*this->singleSlater_->moA());
  if(!doVVOO && this->Ref_ != SingleSlater<dcomplex>::TCS) {
    if(this->Ref_ == SingleSlater<dcomplex>::RHF)
      TMOB = this->singleSlater_->moA()->adjoint() * TAOB * (*this->singleSlater_->moA());
    else 
      TMOB = this->singleSlater_->moB()->adjoint() * TAOB * (*this->singleSlater_->moB());
  }
  if(doVOOV) this->retrvVOOV(TMOV,TMOA,TMOB);
  else if(doVVOO) this->retrvVVOO(TMOV,TMOA);
  
}// formMOTDen

template<>
void SDResponse<dcomplex>::initRMu(){
  // RMu = [ e(HOMO) + e(LUMO) ] / 2
  if(this->Ref_ == SingleSlater<dcomplex>::TCS){
    this->rMu_ = ( (*this->singleSlater_->epsA())(this->nO_-1) + 
                   (*this->singleSlater_->epsA())(this->nO_)    ) / 2.0;
  } else {
    if(this->Ref_ == SingleSlater<dcomplex>::RHF || this->nOB_ == 0)
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
void SDResponse<dcomplex>::scaleDagPPRPA(bool inplace, ComplexVecMap &T, ComplexVecMap &IX, ComplexVecMap *AX){
  if(!this->haveDag_) this->getDiag();
  if(inplace) IX    = IX + T.cwiseProduct(*this->rmDiag_);
  else        (*AX) = IX + T.cwiseProduct(*this->rmDiag_);
} // scaleDagPPRPA

template<>
void SDResponse<dcomplex>::initMeth(){
  if(this->nSek_ == 0) 
    CErr("Must set NSek before initializing a PSCF method",this->fileio_->out);
  if((this->iMeth_ == PPRPA || this->iMeth_ == PPATDA || this->iMeth_ == PPCTDA)
      && !(this->iPPRPA_ >= 0 && this->iPPRPA_ <= 2))
    CErr("Invalid iPPRPA_ in SDResponse::initMeth()",this->fileio_->out);

  if(this->iMeth_ == CIS){
    /******************/
    /* CIS Single Dim */
    if(this->Ref_ == SingleSlater<dcomplex>::TCS)
      this->nSingleDim_ = this->nOV_;
    else 
      this->nSingleDim_ = this->nOAVA_ + this->nOBVB_;
    /******************/
  } else if(this->iMeth_ == RPA || this->iMeth_ == STAB){
    /******************/
    /* RPA Single Dim */
    if(this->Ref_ == SingleSlater<dcomplex>::TCS)
      this->nSingleDim_ = 2*this->nOV_;
    else
      this->nSingleDim_ = 2*(this->nOAVA_ + this->nOBVB_);
    /******************/
  } else if(this->iMeth_ == PPRPA){
    /*********************/
    /* pp-RPA Single Dim */
    if(this->Ref_ == SingleSlater<dcomplex>::TCS)
      this->nSingleDim_ = this->nVV_SLT_ + this->nOO_SLT_;
    else {
      if(this->iPPRPA_ == 0)      this->nSingleDim_ = this->nVAVA_SLT_ + this->nOAOA_SLT_;
      else if(this->iPPRPA_ == 1) this->nSingleDim_ = this->nVAVB_     + this->nOAOB_;
      else if(this->iPPRPA_ == 2) this->nSingleDim_ = this->nVBVB_SLT_ + this->nOBOB_SLT_;
    }
    /*********************/
  } else if(this->iMeth_ == PPATDA){
    /*************************/
    /* pp-TDA (A) Single Dim */
    if(this->Ref_ == SingleSlater<dcomplex>::TCS)
      this->nSingleDim_ = this->nVV_SLT_;
    else {
      if(this->iPPRPA_ == 0)      this->nSingleDim_ = this->nVAVA_SLT_;
      else if(this->iPPRPA_ == 1) this->nSingleDim_ = this->nVAVB_    ; 
      else if(this->iPPRPA_ == 2) this->nSingleDim_ = this->nVBVB_SLT_;
    }
    /*************************/
  } else if(this->iMeth_ == PPCTDA){
    /*************************/
    /* pp-TDA (C) Single Dim */
    if(this->Ref_ == SingleSlater<dcomplex>::TCS)
      this->nSingleDim_ = this->nOO_SLT_;
    else {
      if(this->iPPRPA_ == 0)      this->nSingleDim_ = this->nOAOA_SLT_;
      else if(this->iPPRPA_ == 1) this->nSingleDim_ = this->nOAOB_;
      else if(this->iPPRPA_ == 2) this->nSingleDim_ = this->nOBOB_SLT_;
    }
    /*************************/
  } else {
    CErr("PSCF Method " + std::to_string(this->iMeth_) + " NYI",this->fileio_->out);
  }
} // initMeth

template<>
void SDResponse<dcomplex>::checkValid(){
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
  if(this->Ref_ == SingleSlater<dcomplex>::TCS)
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
void SDResponse<dcomplex>::mpiSend(int toID,int tag) {
  //OOMPI_COMM_WORLD[toID].Send(this->nAtoms_,tag);
  //OOMPI_COMM_WORLD[toID].Send(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiSend(toID,tag);
}; //mpiSend
template<>
void SDResponse<dcomplex>::mpiRecv(int fromID,int tag) {
  //OOMPI_COMM_WORLD[fromID].Recv(this->nAtoms_,tag);
  //this->index_=new int[this->nAtoms_];
  //this->cart_ =new Matrix(3, this->nAtoms_, "Molecule");
  //OOMPI_COMM_WORLD[fromID].Recv(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiRecv(fromID,tag);
}; //mpiRecv

} // namespace ChronusQ
