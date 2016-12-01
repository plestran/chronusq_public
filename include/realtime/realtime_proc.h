/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
template <typename T>
void RealTime<T>::doPropagation() {

  auto NBT = ssPropagator_->nBasis() * ssPropagator_->nTCS();

  bool Start(false); // Start the MMUT iterations
  bool FinMM(false); // Wrap up the MMUT iterations

  long double currentTime(0.);
  PropagationStep currentStep;

  initCSV();

  // Logic for Delta pulse
//if(this->tOff_ != 0 and (this->tOff_ - this->tOn_) < this->stepSize_){
//  this->tOff_ = this->stepSize_; 
//}
  printRTHeader();

  for(auto iStep = 0ul; iStep <= maxSteps_; iStep++) {

    if(iScheme_ == MMUT) {
      Start = (iStep == 0);
      Start = Start or FinMM;
      if(iRstrt_ > 0) Start = Start or (iStep % iRstrt_) == 0;
      
      FinMM = (iStep == maxSteps_);
      FinMM = FinMM or (tOff_ != 0.0 and currentTime > tOff_ and 
                        currentTime <= tOff_ + stepSize_);

      if(iRstrt_ > 0) FinMM = FinMM or ( (iStep + 1) % iRstrt_ == 0 );

      if(Start or FinMM) currentStep = iRstScheme_;
      else               currentStep = ModifiedMidpoint;

      if(Start or FinMM) 
        fileio_->out << "  *** Performing MMUT Restart ***" << endl;
    } else if(iScheme_ == ExpMagnus2) {
      currentStep = ExplicitMagnus2;
    } else if(iScheme_ == ExpMagnus3) {
      currentStep = ExplicitMagnus3;
    }

    if(currentStep == ModifiedMidpoint) {

      // Swap POSav densities with those from ssPropagator
      // for leapfrog step
      //
      // POSav(k) <-> PO(k)
      //
      // POSav(k) = PO(k)
      // PO(k)    = PO(k-1)
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        POSav_[iODen]->swap(*ssPropagator_->onePDMOrtho()[iODen]);
      
      // Double the step size for MMUT Step
      deltaT_ = 2*stepSize_;

    } else {

      // Copy the orthonormal density from ssPropagator to POSav
      // POSav(k) = PO(k)
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        (*POSav_[iODen]) = (*ssPropagator_->onePDMOrtho()[iODen]);

      deltaT_ = stepSize_;

    }

    // Obtain field value for current time point
    formField(currentTime);

    // Form AO Fock Matrix and compute properties
    // F(k) = H + G[P(k)]
    ssPropagator_->formFock();
    ssPropagator_->computeEnergy();
    ssPropagator_->computeProperties();

    // Print line in output file
    printRTStep(currentTime);

    // Orthonormalize the AO Fock
    // F(k) -> FO(k)
    ssPropagator_->orthoFock3();

    // Copy Fock's over for EM2 / EM3
    if(currentStep == ExplicitMagnus2 or currentStep == ExplicitMagnus3) {
      // Copy the orthonormal fock from ssPropagator to FOSav
      // FOSav(k) = FO(k)
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++){
        (*FOSav_[iOFock]) = (*ssPropagator_->fockOrtho()[iOFock]);
      }

    }

    // Scale Fock's for EM3
    if(currentStep == ExplicitMagnus3) {
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++){
        //(*FOSav_[iOFock]) *= 0.5;
        (*ssPropagator_->fockOrtho()[iOFock]) *= 0.5;
      }
    }

    // Form the unitary propagation matrix
    // U**H(k) = exp(-i * dt * F(k))
    formUTrans();

    // Propagate the (orthonormal) Density Matrix
    // PO (in ssPropagator) now stores PO(k+1)
    // PO(k+1) = U**H(k) * PO(k-1) * U(k)
    propDen();

    // Add a time point record onto the list
    addRecord(currentTime);

    // Write CSV Files
    writeCSVs();

    // Unorthonormalize the density for next Fock build
    // PO(k+1) -> P(k+1)
    ssPropagator_->unOrthoDen3();

    // Further propagate using Explicit Magnus 2
    if(currentStep == ExplicitMagnus2) {
      // Form and orthonormalize Fock at t + dt using FE density
      formField(currentTime + stepSize_);
      ssPropagator_->formFock();
      ssPropagator_->orthoFock3();

      // FO = 0.5 * (FO + FOSav)
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++) {
        ssPropagator_->fockOrtho()[iOFock]->noalias() += (*FOSav_[iOFock]);
        (*ssPropagator_->fockOrtho()[iOFock]) *= dcomplex(0.5);
      }

      // Copy POSav to PO
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        (*ssPropagator_->onePDMOrtho()[iODen]) = (*POSav_[iODen]); 

      // Propagate using new "averaged" fock
      formUTrans();
      propDen();
      ssPropagator_->unOrthoDen3();
    } else if(currentStep == ExplicitMagnus3) {
      // Form and orthonormalize Fock at t + 0.5*dt using new density
      formField(currentTime + 0.5*stepSize_);
      ssPropagator_->formFock();
      ssPropagator_->orthoFock3();

      // Form FOSav2 and FOSav3
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++){
        (*FOSav2_[iOFock]) = (*ssPropagator_->fockOrtho()[iOFock]);

        FOSav3_[iOFock]->noalias() = 
          0.25 * ((*FOSav_[iOFock]) + (*FOSav2_[iOFock]));
      }
     
      // update ssPropagator
      //  FO = FOSav3
      //  PO = POSav 
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++) {
        ssPropagator_->fockOrtho()[iOFock]->noalias() = (*FOSav3_[iOFock]);
      }

      // Copy POSav to PO
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        (*ssPropagator_->onePDMOrtho()[iODen]) = (*POSav_[iODen]); 

      // Propagate and update Fock
      formUTrans();
      propDen();
      ssPropagator_->unOrthoDen3();
      formField(currentTime + 0.5*stepSize_);
      ssPropagator_->formFock();
      ssPropagator_->orthoFock3();

      // Form FOSav4
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++){
        (*FOSav4_[iOFock]) = (*ssPropagator_->fockOrtho()[iOFock]);
      }
  
      //update ssPropagator 
      //  FO = FOSav2
      //  PO = POSav 
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++) {
        ssPropagator_->fockOrtho()[iOFock]->noalias() = (*FOSav2_[iOFock]);
      }

      // Copy POSav to PO
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        (*ssPropagator_->onePDMOrtho()[iODen]) = (*POSav_[iODen]); 

      // Propagate and update Fock
      formUTrans();
      propDen();
      ssPropagator_->unOrthoDen3();
      formField(currentTime + stepSize_);
      ssPropagator_->formFock();
      ssPropagator_->orthoFock3();

      // Form FOSav5
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++){
        (*FOSav5_[iOFock]) = (*ssPropagator_->fockOrtho()[iOFock]);
      }

      // Final fock update
      // FO =   (1/6)  * (FOSav + 4*FOSav4 + FOSav5)
      //      - (i/3)  * [FOSav3,FOSav4]
      //      - (i/12) * [FOSav2,FOSav5]
      //
      // PO = POSav
      dcomplex h = -dcomplex(0,1)*stepSize_;

      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++){
        ssPropagator_->fockOrtho()[iOFock]->noalias() = (1./6.) * 
          ((*FOSav_[iOFock]) + 4*(*FOSav4_[iOFock]) + (*FOSav5_[iOFock]));
      }

/*
      *ssPropagator_->fockOrtho()[0] -= 
        (h/3.)*( (*FOSav3_[0])*(*FOSav4_[0]) - (*FOSav4_[0])*(*FOSav3_[0]));

      *ssPropagator_->fockOrtho()[0] -= 
        (h/12.)*( (*FOSav2_[0])*(*FOSav5_[0]) - (*FOSav5_[0])*(*FOSav2_[0]));
*/

      // Store [F3,F4] in F (not being used)
      incoreHComm((-h/3.), ssPropagator_->fock(),FOSav3_,FOSav4_);

      // Increment FO by F
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++)
        ssPropagator_->fockOrtho()[iOFock]->noalias() +=
          (*ssPropagator_->fock()[iOFock]);
        
      // Store [F2,F5] in F (not being used)
      incoreHComm((-h/12.),ssPropagator_->fock(),FOSav2_,FOSav5_);

      // Increment FO by F
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++)
        ssPropagator_->fockOrtho()[iOFock]->noalias() +=
          (*ssPropagator_->fock()[iOFock]);


      // Copy POSav to PO
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        (*ssPropagator_->onePDMOrtho()[iODen]) = (*POSav_[iODen]); 

      // Propagate and update Fock
      formUTrans();
      propDen();
      ssPropagator_->unOrthoDen3();
    }

    // Increment the current time
    currentTime += stepSize_;
  }
  if(tarCSVs)  tarCSVFiles();
};

template <typename T>
void RealTime<T>::addRecord(const long double currentTime) {
  auto NB  = ssPropagator_->basisset()->nBasis();
  auto NBT = NB * ssPropagator_->nTCS();

  PropInfo rec;
  rec.timeStep  = currentTime;
  rec.energy    = ssPropagator_->totalEnergy();
  rec.dipole[0] = ssPropagator_->elecDipole()[0]/phys.debye;
  rec.dipole[1] = ssPropagator_->elecDipole()[1]/phys.debye;
  rec.dipole[2] = ssPropagator_->elecDipole()[2]/phys.debye;
  rec.dipole[3] = std::sqrt( std::pow(rec.dipole[0],2.0) +
                             std::pow(rec.dipole[1],2.0) +
                             std::pow(rec.dipole[2],2.0) );
  rec.appliedfield[0] = EDField_[0];
  rec.appliedfield[1] = EDField_[1];
  rec.appliedfield[2] = EDField_[2];
  rec.appliedfield[3] = std::sqrt( std::pow(rec.appliedfield[0],2.0) +
                             std::pow(rec.appliedfield[1],2.0) +
                             std::pow(rec.appliedfield[2],2.0) );

  rec.mullPop = ssPropagator_->mullPop();
  rec.lowPop  = ssPropagator_->lowPop();

  orbPop(rec.orbitalOccA,rec.orbitalOccB);

  propInfo.push_back(rec);
}

template <typename T>
void RealTime<T>::orbPop(std::vector<double> &a, std::vector<double> &b){
  auto NB  = ssPropagator_->basisset()->nBasis();
  auto NBT = NB * ssPropagator_->nTCS();

  a.resize(NBT);
  if(ssPropagator_->nTCS() == 1 and !ssPropagator_->isClosedShell)
    b.resize(NBT);

  
  ComplexMap S(NBTSqScratch_,NBT,NBT);
  ComplexMap S2(NBTSqScratch2_,NBT,NBT);

  // Scatter the orthonormal densities into the MO storage
  if(ssPropagator_->nTCS() == 1 and ssPropagator_->isClosedShell)
    (*ssPropagator_->moA()) = 0.5 * (*ssPropagator_->onePDMOrtho()[0]);
  else if(ssPropagator_->nTCS() == 1 and !ssPropagator_->isClosedShell) {
    (*ssPropagator_->moA()) = 
      0.5 * ((*ssPropagator_->onePDMOrtho()[0]) + 
             (*ssPropagator_->onePDMOrtho()[1]));
    (*ssPropagator_->moB()) = 
      0.5 * ((*ssPropagator_->onePDMOrtho()[0]) - 
             (*ssPropagator_->onePDMOrtho()[1]));
  } else {
    std::vector<std::reference_wrapper<ComplexMap>> scattered;
    int I = 0;
    for(auto iF = ssPropagator_->onePDMOrtho().begin(); 
        iF != ssPropagator_->onePDMOrtho().end(); iF++){
      scattered.emplace_back(*(*iF));
    }
    Quantum<dcomplex>::spinGather(*ssPropagator_->moA(),scattered);
  }

  // S2 = C(0)**H * PO(t) * C(0)
  S.noalias() = 
    groundState_->moA()->template cast<dcomplex>().adjoint() *
    (*ssPropagator_->moA());
  S2.noalias() = S * groundState_->moA()->template cast<dcomplex>();

  for(auto i = 0; i < NBT; i++) a[i] = std::real(S2(i,i));

  if(ssPropagator_->nTCS() == 1 and !ssPropagator_->isClosedShell){
    S.noalias() = 
      groundState_->moB()->template cast<dcomplex>().adjoint() *
      (*ssPropagator_->moB());
    S2.noalias() = S * groundState_->moB()->template cast<dcomplex>();

    for(auto i = 0; i < NBT; i++) b[i] = std::real(S2(i,i));
  }
  
}

template <typename T>
void RealTime<T>::incoreHComm(dcomplex fact,
  std::vector<ComplexMap*> &COMM, std::vector<ComplexMap*> &A, 
  std::vector<ComplexMap*> &B) {

  for( auto x : COMM) x->setZero();

  // Scalar Part
  // AB(S) = (A(S) * B(S) + A(k) * B(k))
  for(auto i = 0; i < COMM.size(); i++)
    COMM[0]->noalias() += (*A[i]) * (*B[i]);

  
  if(COMM.size() > 1) {
    // Non cross product part
    // AB(k) = (A(S) * B(k) + A(k) * B(S))
    for(auto i = 1; i < COMM.size(); i++){
      COMM[i]->noalias() += (*A[0]) * (*B[i]);
      COMM[i]->noalias() += (*A[i]) * (*B[0]);
    }
  }

  if(COMM.size() > 2) {
    // Z Cross product part
    // AB(z) += i * (A(x) * B(y) - A(y) * B(x))
    COMM[3]->noalias() += dcomplex(0,1) * (*A[1]) * (*B[2]);
    COMM[3]->noalias() -= dcomplex(0,1) * (*A[2]) * (*B[1]);

    // X Cross product part
    // AB(x) += i * (A(y) * B(z) - A(z) * B(y))
    COMM[1]->noalias() += dcomplex(0,1) * (*A[2]) * (*B[3]);
    COMM[1]->noalias() -= dcomplex(0,1) * (*A[3]) * (*B[2]);

    // Y Cross product part
    // AB(y) += i * (A(z) * B(x) - A(x) * B(z))
    COMM[2]->noalias() += dcomplex(0,1) * (*A[3]) * (*B[1]);
    COMM[2]->noalias() -= dcomplex(0,1) * (*A[1]) * (*B[3]);
  }


  auto NB  = ssPropagator_->basisset()->nBasis();
  ComplexMap S(NBTSqScratch_,NB,NB);

  // Because A and B are hermetian, the commutator
  // fact * [A,B] = fact * ( AB - (AB)**H )
  for(auto i = 0; i < COMM.size(); i++){
    S.noalias()        = (*COMM[i]);
    COMM[i]->noalias() = 0.5 * fact * ( S - S.adjoint() );   
  }
};
