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
#include <realtime.h>
#include <aointegrals.h>
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;
using ChronusQ::RealTime;

namespace ChronusQ {

template<>
void RealTime<double>::iniDensity() {
  bool inOrthoBas;
  bool idempotent;

  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
  // Set up Eigen Maps
  ComplexMap oTrans1(this->oTrans1Mem_,NTCSxNBASIS,NTCSxNBASIS); 
  ComplexMap oTrans2(this->oTrans2Mem_,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap POA    (this->POAMem_    ,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap POAsav (this->POAsavMem_ ,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap FOA    (this->FOAMem_    ,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap initMOA(this->initMOAMem_,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap scratch(this->scratchMem_,NTCSxNBASIS,NTCSxNBASIS);

  ComplexMap POB    (this->POBMem_    ,0,0);
  ComplexMap POBsav (this->POBsavMem_ ,0,0);
  ComplexMap FOB    (this->FOBMem_    ,0,0);
  ComplexMap initMOB(this->initMOBMem_,0,0);

  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS){
    new (&POB    ) ComplexMap(this->POBMem_    ,NTCSxNBASIS,NTCSxNBASIS);
    new (&POBsav ) ComplexMap(this->POBsavMem_ ,NTCSxNBASIS,NTCSxNBASIS);
    new (&FOB    ) ComplexMap(this->FOBMem_    ,NTCSxNBASIS,NTCSxNBASIS);
    new (&initMOB) ComplexMap(this->initMOBMem_,NTCSxNBASIS,NTCSxNBASIS);
  }

// Form the orthonormal transformation matrices
  if (this->typeOrtho_ == 1) {  
   // Lowdin transformation 
   // V1 = S^(-1/2)
   // V2 = S^(1/2)

/*
    oTrans1.real() = (*this->aointegrals_->overlap_).pow(-0.5);
    oTrans2.real() = (*this->aointegrals_->overlap_).pow(0.5);
*/

    char JOBZ = 'V';
    char UPLO = 'L';
    int INFO;

    double *A =    this->REAL_LAPACK_SCR;
    double *W =    A + NTCSxNBASIS * NTCSxNBASIS;
    double *WORK = W + NTCSxNBASIS;

    RealVecMap E(W,NTCSxNBASIS);
    RealMap    V(A,NTCSxNBASIS,NTCSxNBASIS);
    RealMap    S(WORK,NTCSxNBASIS,NTCSxNBASIS); // Requires WORK to be NBSq

    E.setZero();
    V.setZero();
    S.setZero();

    std::memcpy(A,this->aointegrals_->overlap_->data(),
      NTCSxNBASIS*NTCSxNBASIS*sizeof(double));

    dsyev_(&JOBZ,&UPLO,&NTCSxNBASIS,A,&NTCSxNBASIS,W,WORK,&this->lWORK,&INFO);

    V.transposeInPlace(); // BC Col major
    std::memcpy(WORK,A,NTCSxNBASIS*NTCSxNBASIS*sizeof(double));

    for(auto i = 0; i < NTCSxNBASIS; i++){ S.col(i) *= std::sqrt(W[i]); }
    oTrans2.real() = S * V.adjoint();

    for(auto i = 0; i < NTCSxNBASIS; i++){ S.col(i) /= W[i]; }
    oTrans1.real() = S * V.adjoint();

    if(this->controls_->printLevel>3) {
      prettyPrintComplex(this->fileio_->out,oTrans1,"S^(-1/2)");
      prettyPrintComplex(this->fileio_->out,oTrans2,"S^(1/2)");
    }
  }
  else if (this->typeOrtho_ == 2) {  
  // Cholesky transformation
    CErr("Cholesky orthogonalization NYI",this->fileio_->out);
  }
  else if (this->typeOrtho_ == 3) {  	
    CErr("Canonical orthogonalization NYI",this->fileio_->out);
  // Canonical orthogonalization
  // V1 = U*s^(-1/2)
  // V2 = S*V1
  }

// Form the initial density
  if (this->initDensity_ == 0) { 
// Use converged ground-state density
    inOrthoBas = false;
    idempotent = true;
  }
  else if (this->initDensity_ == 1) { 
// Form the initial density by swaping MOs
    inOrthoBas = false;
    idempotent = true;
    if (this->swapMOA_ != 0) {
      // MOs to swap
      int iA = ((this->swapMOA_)/1000); 
      int jA = ((this->swapMOA_)%1000); 

      this->fileio_->out << endl << "Alpha MOs swapped: "
                         << iA << " <-> " << jA << endl;

      if(this->controls_->printLevel > 3) {
        prettyPrint(this->fileio_->out,
                    (*this->ssPropagator_->moA()),"Initial Alpha MO");
      }
      this->ssPropagator_->moA()->col(jA-1).swap(
        this->ssPropagator_->moA()->col(iA-1)
        );
    }
    if (this->swapMOB_ != 0 && 
        !this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      // MOs to swap
      int iB = (this->swapMOB_/1000); 
      int jB = (this->swapMOB_%1000); 
      this->fileio_->out << endl << "Beta MOs swapped: "
                         << iB << " <-> " << jB << endl;

      if(this->controls_->printLevel > 3) {
        prettyPrint(this->fileio_->out,
                    (*this->ssPropagator_->moB()),"Initial Beta MO");
      }
      this->ssPropagator_->moB()->col(jB-1).swap(
        this->ssPropagator_->moB()->col(iB-1)
        );
    }
    this->ssPropagator_->formDensity();

    // Transform the ground state MO to orthonormal basis
    initMOA.setZero();
    initMOA.real() = *this->groundState_->moA();
    initMOA = oTrans2 * initMOA;
    if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      initMOB.setZero();
      initMOB.real() = *this->groundState_->moB();
      initMOB = oTrans2 * initMOB;
    }
  }
  else if (this->initDensity_ == 2) { 
// Read in the AO density from checkpoint file
    CErr("Read in the AO density from checkpint file NYI",this->fileio_->out);
  }
  else if (this->initDensity_ == 3) { 
// Read in the orthonormal density from checkpoint file
    CErr("Read in the orthonormal density from checkpoint file NYI",
         this->fileio_->out);
  }

  if (!inOrthoBas) { 
// Transform density from AO to orthonormal basis
    POA    = oTrans2 * (*this->ssPropagator_->densityA()) * oTrans2;

    POAsav = POA;
    if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      POB    = oTrans2 * (*this->ssPropagator_->densityB()) * oTrans2;
      POBsav = POB;
    }
  }
  else { 
// Transform density from orthonormal to AO basis
    (*this->ssPropagator_->densityA()) = oTrans1 * POAsav * oTrans1;

    if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) 
      (*this->ssPropagator_->densityB()) = oTrans1 * POB * oTrans1;
  }
};


template<>
void RealTime<double>::formUTrans() {  
//
// Form the unitary transformation matrix: 
// U = exp(-i*dT*F)
//

  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
  // Set up Eigen Maps
  ComplexMap uTransA(this->uTransAMem_,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap scratch(this->scratchMem_,NTCSxNBASIS,NTCSxNBASIS);

  ComplexMap uTransB(this->uTransBMem_,0,0);

  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS){
    new (&uTransB) ComplexMap(this->uTransBMem_,NTCSxNBASIS,NTCSxNBASIS);
  }
  // FIXME: Eigen's Eigensolver is terrible, replace with LAPACK routines
  if (this->methFormU_ == 1) { 
   //  Eigen-decomposition

/*     
    ComplexMatrix EVec(NTCSxNBASIS,NTCSxNBASIS);
    RealMatrix 	  EVal(NTCSxNBASIS,1);

    Eigen::SelfAdjointEigenSolver<ComplexMatrix> 
      sys(*this->ssPropagator_->fockA());

    EVec = sys.eigenvectors();
    EVal = sys.eigenvalues();
    
    uTransA.setZero();
    for (int i = 0; i < NTCSxNBASIS; i++) {
      uTransA(i,i) = 
        dcomplex( cos(deltaT_ * EVal(i,0)), -sin(deltaT_ * EVal(i,0)) );
    }

    uTransA = EVec * uTransA * EVec.adjoint();
    if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      Eigen::SelfAdjointEigenSolver<ComplexMatrix> 
        sys(*this->ssPropagator_->fockB());

      EVec = sys.eigenvectors();
      EVal = sys.eigenvalues();
      uTransB.setZero();
      for (int i = 0; i < NTCSxNBASIS; i++) {
        uTransB(i,i) = 
          dcomplex( cos(deltaT_ * EVal(i,0)), -sin(deltaT_ * EVal(i,0)) );
      }
      uTransB = EVec * uTransB * EVec.adjoint();
    }
*/
    char JOBZ = 'V';
    char UPLO = 'L';
    int INFO;

    dcomplex *A     = this->scratchMem_;
    double   *W     = this->REAL_LAPACK_SCR;
    double   *RWORK = W + std::max(1,3*NTCSxNBASIS-2);
    dcomplex *WORK  = this->CMPLX_LAPACK_SCR;

    RealVecMap E(W,NTCSxNBASIS);
    ComplexMap V(A,NTCSxNBASIS,NTCSxNBASIS);
    ComplexMap S(WORK,NTCSxNBASIS,NTCSxNBASIS);

    E.setZero();
    V.setZero();
    S.setZero();

    std::memcpy(A,this->ssPropagator_->fockA()->data(),
      NTCSxNBASIS*NTCSxNBASIS*sizeof(dcomplex));

    V.transposeInPlace(); // BC Col major
    zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,A,&NTCSxNBASIS,W,WORK,&this->lWORK,RWORK,
      &INFO);

    V.transposeInPlace(); // BC Col major
    std::memcpy(WORK,A,NTCSxNBASIS*NTCSxNBASIS*sizeof(dcomplex));


    for(auto i = 0; i < NTCSxNBASIS; i++){ 
      S.col(i) *= dcomplex(std::cos(this->deltaT_ * W[i]),
                          -std::sin(this->deltaT_ * W[i])); 
    }

    uTransA = S * V.adjoint();

    if(!this->isClosedShell_ && this->Ref_ != SingleSlater<dcomplex>::TCS) {
      E.setZero();
      V.setZero();
      S.setZero();
     
      std::memcpy(A,this->ssPropagator_->fockB()->data(),
        NTCSxNBASIS*NTCSxNBASIS*sizeof(dcomplex));
     
      V.transposeInPlace(); // BC Col major
      zheev_(&JOBZ,&UPLO,&NTCSxNBASIS,A,&NTCSxNBASIS,W,WORK,&this->lWORK,RWORK,
        &INFO);
     
      V.transposeInPlace(); // BC Col major
      std::memcpy(WORK,A,NTCSxNBASIS*NTCSxNBASIS*sizeof(dcomplex));
     
     
      for(auto i = 0; i < NTCSxNBASIS; i++){ 
        S.col(i) *= dcomplex(std::cos(this->deltaT_ * W[i]),
                            -std::sin(this->deltaT_ * W[i])); 
      }
      uTransB = S * V.adjoint();
    }
  }
  else if (this->methFormU_ == 2) { 
  // Taylor expansion
 
/*  This is not taylor and breaks with the new memory scheme
 
    scratch = -math.ii * deltaT_ * (*this->ssPropagator_->fockA());
    uTransA = scratch.exp(); // FIXME
    if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      scratch = -math.ii * deltaT_ * (*this->ssPropagator_->fockB());
      uTransB = scratch.exp(); // FIXME
    }
*/
  }
//    prettyPrint(this->fileio_->out,(*this->uTransA_),"uTransA");
//    if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) prettyPrint(this->fileio_->out,(*this->uTransB_),"uTransB");
};
  
template<>
void RealTime<double>::doPropagation() {
  long int iStep;
  bool checkFP = false;

  currentTime_ = 0.0;

  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
  // Set up Eigen Maps
  ComplexMap oTrans1(this->oTrans1Mem_,NTCSxNBASIS,NTCSxNBASIS); 
  ComplexMap oTrans2(this->oTrans2Mem_,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap POA    (this->POAMem_    ,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap POAsav (this->POAsavMem_ ,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap FOA    (this->FOAMem_    ,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap initMOA(this->initMOAMem_,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap uTransA(this->uTransAMem_,NTCSxNBASIS,NTCSxNBASIS);
  ComplexMap scratch(this->scratchMem_,NTCSxNBASIS,NTCSxNBASIS);

  ComplexMap POB    (this->POBMem_    ,0,0);
  ComplexMap POBsav (this->POBsavMem_ ,0,0);
  ComplexMap FOB    (this->FOBMem_    ,0,0);
  ComplexMap initMOB(this->initMOBMem_,0,0);
  ComplexMap uTransB(this->uTransBMem_,0,0);

  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS){
    new (&POB    ) ComplexMap(this->POBMem_    ,NTCSxNBASIS,NTCSxNBASIS);
    new (&POBsav ) ComplexMap(this->POBsavMem_ ,NTCSxNBASIS,NTCSxNBASIS);
    new (&FOB    ) ComplexMap(this->FOBMem_    ,NTCSxNBASIS,NTCSxNBASIS);
    new (&initMOB) ComplexMap(this->initMOBMem_,NTCSxNBASIS,NTCSxNBASIS);
    new (&uTransB) ComplexMap(this->uTransBMem_,NTCSxNBASIS,NTCSxNBASIS);
  }

  for (iStep = 0; iStep <= this->maxSteps_; iStep++) {
    //this->fileio_->out<<"\nStep "<<iStep<<":\n"<<endl;
    //this->fileio_->out<<std::right<<std::setw(20)<<"Time = "<<std::setw(15)<<currentTime_<<std::setw(5)<<" a.u. "<<endl;

    if (iStep == 0) deltaT_ = this->stepSize_;
    else            deltaT_ = 2.0 * (this->stepSize_);

    scratch = POA;
    POA     = POAsav;
    POAsav  = scratch;
    if (!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      scratch = POB;
      POB     = POBsav;
      POBsav  = scratch;
    }

//  Print 
    if(this->controls_->printLevel >= 1) {
      //prettyPrintComplex(this->fileio_->out,(*this->ssPropagator_->densityA()),"Alpha AO Density");
      //if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) prettyPrintComplex(this->fileio_->out,(*this->ssPropagator_->densityB()),"Beta AO Density");

//  Form AO Fock matrix
    this->formEDField();
    this->ssPropagator_->setField(this->EDField_);
    this->ssPropagator_->formFock();
    this->ssPropagator_->computeEnergy();
    this->ssPropagator_->computeMultipole();
    this->printRT();

//  Transform Fock from AO to orthonormal basis
    scratch = (*this->ssPropagator_->fockA());
    (*this->ssPropagator_->fockA()) = oTrans1 * scratch * oTrans1;

    if (!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      scratch = (*this->ssPropagator_->fockB());
      (*this->ssPropagator_->fockB()) = oTrans1 * scratch * oTrans1;
    }

//  Form the unitary propagation matrix
    this->formUTrans();

    if ((this->initDensity_ == 0) && checkFP) {
//    Check [F,P] for converged density, should equal to zero
      scratch =  (*this->ssPropagator_->fockA()) * POAsav; 
      scratch -= POAsav * (*this->ssPropagator_->fockA());

      prettyPrint(this->fileio_->out,scratch,"[FOA,POA]");

      if (!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
        scratch =  (*this->ssPropagator_->fockB()) * POBsav; 
        scratch -= POBsav * (*this->ssPropagator_->fockB());

        prettyPrint(this->fileio_->out,scratch,"[FOB,POB]");
      }
    }

//  Propagate the density matrix
//
    scratch = POA;
    POA     = uTransA * scratch * uTransA.adjoint();

    if (!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      scratch = POB;
      POB     = uTransB * scratch * uTransB.adjoint();
    }

//  Transform density matrix from orthonormal to AO basis
    (*this->ssPropagator_->densityA()) = oTrans1 * POA * oTrans1;

    if (!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS) {
      (*this->ssPropagator_->densityB()) = oTrans1 * POB * oTrans1;
    }

//  Advance step
    currentTime_ += this->stepSize_;
    };
  }
  delete [] this->SCR;
  delete [] this->REAL_LAPACK_SCR;
};

} // namespace ChronusQ
