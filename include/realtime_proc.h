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
template<typename T>
void RealTime<T>::doPropagation(){
  long int iStep;
  bool checkFP = false;

  currentTime_ = 0.0;

  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
  // Set up Eigen Maps
  ComplexMap oTrans1(this->oTrans1Mem_,0,0); 
  ComplexMap oTrans2(this->oTrans2Mem_,0,0);
  ComplexMap POA    (this->POAMem_    ,0,0);
  ComplexMap POAsav (this->POAsavMem_ ,0,0);
  ComplexMap FOA    (this->FOAMem_    ,0,0);
  ComplexMap initMOA(this->initMOAMem_,0,0);
  ComplexMap uTransA(this->uTransAMem_,0,0);
  ComplexMap scratch(this->scratchMem_,0,0);

  ComplexMap POB    (this->POBMem_    ,0,0);
  ComplexMap POBsav (this->POBsavMem_ ,0,0);
  ComplexMap FOB    (this->FOBMem_    ,0,0);
  ComplexMap initMOB(this->initMOBMem_,0,0);
  ComplexMap uTransB(this->uTransBMem_,0,0);

  if(getRank() == 0) {
    new (&oTrans1) ComplexMap(this->oTrans1Mem_,NTCSxNBASIS,NTCSxNBASIS); 
    new (&oTrans2) ComplexMap(this->oTrans2Mem_,NTCSxNBASIS,NTCSxNBASIS);
    new (&POA    ) ComplexMap(this->POAMem_    ,NTCSxNBASIS,NTCSxNBASIS);
    new (&POAsav ) ComplexMap(this->POAsavMem_ ,NTCSxNBASIS,NTCSxNBASIS);
    new (&FOA    ) ComplexMap(this->FOAMem_    ,NTCSxNBASIS,NTCSxNBASIS);
    new (&initMOA) ComplexMap(this->initMOAMem_,NTCSxNBASIS,NTCSxNBASIS);
    new (&uTransA) ComplexMap(this->uTransAMem_,NTCSxNBASIS,NTCSxNBASIS);
    new (&scratch) ComplexMap(this->scratchMem_,NTCSxNBASIS,NTCSxNBASIS);
 
    if(!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS){
      new (&POB    ) ComplexMap(this->POBMem_    ,NTCSxNBASIS,NTCSxNBASIS);
      new (&POBsav ) ComplexMap(this->POBsavMem_ ,NTCSxNBASIS,NTCSxNBASIS);
      new (&FOB    ) ComplexMap(this->FOBMem_    ,NTCSxNBASIS,NTCSxNBASIS);
      new (&initMOB) ComplexMap(this->initMOBMem_,NTCSxNBASIS,NTCSxNBASIS);
      new (&uTransB) ComplexMap(this->uTransBMem_,NTCSxNBASIS,NTCSxNBASIS);
    }
    this->initCSV();
  }
  for (iStep = 0; iStep <= this->maxSteps_; iStep++) {

    if (iStep == 0) deltaT_ = this->stepSize_;
    else            deltaT_ = 2.0 * (this->stepSize_);

    if(getRank() == 0) {
      scratch = POA;
      POA     = POAsav;
      POAsav  = scratch;
      if (!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS) {
        scratch = POB;
        POB     = POBsav;
        POBsav  = scratch;
      }
    }

//  Form AO Fock matrix
    if(getRank() == 0) {
      this->formEDField();
      this->ssPropagator_->setField(this->EDField_);
    }

#ifdef CQ_ENABLE_MPI
    MPI_Bcast(&this->EDField_[0],3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

    this->ssPropagator_->mpiBCastDensity();
    this->ssPropagator_->formFock();
    this->ssPropagator_->computeEnergy();
    this->ssPropagator_->computeProperties();
    if(this->printLevel_ > 3)
      this->ssPropagator_->printProperties();

    if(getRank() == 0) {
      this->ssPropagator_->mullikenPop();
      this->printRT();
    }

    if(getRank() == 0) {
//    Transform Fock from AO to orthonormal basis
      scratch = (*this->ssPropagator_->fockA());
      (*this->ssPropagator_->fockA()) = oTrans1 * scratch * oTrans1.adjoint();
      if (!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS) {
        scratch = (*this->ssPropagator_->fockB());
        (*this->ssPropagator_->fockB()) = oTrans1 * scratch * oTrans1.adjoint();
      }

//    Form the unitary propagation matrix
      this->formUTrans();
     
      if ((this->initDensity_ == 0) && checkFP) {
//      Check [F,P] for converged density, should equal to zero
        scratch =  (*this->ssPropagator_->fockA()) * POAsav; 
        scratch -= POAsav * (*this->ssPropagator_->fockA());
     
        prettyPrint(this->fileio_->out,scratch,"[FOA,POA]");
     
        if (!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS) {
          scratch =  (*this->ssPropagator_->fockB()) * POBsav; 
          scratch -= POBsav * (*this->ssPropagator_->fockB());
     
          prettyPrint(this->fileio_->out,scratch,"[FOB,POB]");
        }
      }

    } // end serial code
  //  Populate information for printing.
    PropInfo rec;
      rec.timeStep  = this->currentTime_;
      rec.energy    = this->ssPropagator_->totalEnergy;
      rec.dipole[0] = (*this->ssPropagator_->dipole())(0)/phys.debye;
      rec.dipole[1] = (*this->ssPropagator_->dipole())(1)/phys.debye;
      rec.dipole[2] = (*this->ssPropagator_->dipole())(2)/phys.debye;
      rec.dipole[3] = std::sqrt( std::pow(rec.dipole[0],2.0) +
                                 std::pow(rec.dipole[1],2.0) +
                                 std::pow(rec.dipole[2],2.0) );
      rec.appliedfield[0] = this->EDField_[0];
      rec.appliedfield[1] = this->EDField_[1];
      rec.appliedfield[2] = this->EDField_[2];
      rec.appliedfield[3] = std::sqrt( std::pow(rec.appliedfield[0],2.0) +
                                 std::pow(rec.appliedfield[1],2.0) +
                                 std::pow(rec.appliedfield[2],2.0) );
    if(getRank() == 0) {
      rec.mullPop    = (this->ssPropagator_->mullPop());
      scratch = (initMOA.adjoint() * POA * initMOA);
      for(auto idx = 0; idx != NTCSxNBASIS; idx++) {
        rec.orbitalOccA.push_back(scratch(idx,idx).real());
      }
      if(!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS){
        scratch = (initMOB.adjoint() * POB * initMOB);
        for(auto idx = 0; idx != NTCSxNBASIS; idx++) {
          rec.orbitalOccB.push_back(scratch(idx,idx).real());
        }
      }
    }
    if(getRank() == 0) {
      this->writeDipoleCSV(rec,iStep);
      this->writeAppliedFieldCSV(rec,iStep);
      this->writeMullikenCSV(rec,iStep);
      this->writeOrbitalCSV(rec,iStep);
    }
     
    this->propInfo.push_back(rec);
     
    if(getRank() == 0){
//    Propagate the density matrix
      scratch = POA;
      POA     = uTransA * scratch * uTransA.adjoint();
     
      if (!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS) {
        scratch = POB;
        POB     = uTransB * scratch * uTransB.adjoint();
      }
     
//    Transform density matrix from orthonormal to AO basis
      (*this->ssPropagator_->densityA()) = oTrans1.adjoint() * POA * oTrans1;
     
      if (!this->isClosedShell_ && this->Ref_ != SingleSlater<T>::TCS) 
        (*this->ssPropagator_->densityB()) = oTrans1.adjoint() * POB * oTrans1;
    } // end serial code

//  Advance step
    currentTime_ += this->stepSize_;
#ifdef CQ_ENABLE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  if(getRank() == 0) {
    delete [] this->SCR;
    delete [] this->REAL_LAPACK_SCR;

    if(this->tarCSVs) this->tarCSVFiles();
  }
#ifdef CQ_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}; // doPropagation
