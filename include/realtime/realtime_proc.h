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

/*
    Start = (iStep == 0);
    Start = Start or FinMM;
    if(iRstrt_ > 0) Start = Start or (iStep % iRstrt_) == 0;
    
    FinMM = (iStep == maxSteps_);
    FinMM = FinMM or (tOff_ != 0.0 and currentTime_ > tOff_ and currentTime_ <= tOff_ + stepSize_);
    if(iRstrt_ > 0) FinMM = FinMM or ( (iStep + 1) % iRstrt_ == 0 );


    if(Start) this->fileio_->out << "  *** STARTING MMUT ***" << endl;
    if(FinMM) this->fileio_->out << "  *** FINISHING MMUT ***" << endl;

    // Initial entry into MMUT
    if(Start or FinMM) {
      // Copy the orthonormal density from ssPropagator to POSav
      // POSav(k) = PO(k)
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        (*POSav_[iODen]) = (*ssPropagator_->onePDMOrtho()[iODen]);

      deltaT_ = stepSize_;
    // Subsequent MMUT iterations
    } else {
      // Swap POSav densities with those from ssPropagator
      // POSav(k) <-> PO(k)
      //
      // POSav(k) = PO(k)
      // PO(k)    = PO(k-1)
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        POSav_[iODen]->swap(*ssPropagator_->onePDMOrtho()[iODen]);
    }
*/

    //JJG forcing Magnus3
    iScheme_ = ExpMagnus3;
    // JJG alloc scratch for Magnus3
    ComplexMatrix POSav1(NBT, NBT); 
    ComplexMatrix FOSav1(NBT, NBT); 
    ComplexMatrix FOSav2(NBT, NBT); 
    ComplexMatrix FOSav3(NBT, NBT); 
    ComplexMatrix FOSav4(NBT, NBT); 
    ComplexMatrix FOSav5(NBT, NBT); 

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

    if(currentStep == ExplicitMagnus2) {
      // Copy the orthonormal fock from ssPropagator to FOSav
      // FOSav(k) = FO(k)
      for(auto iOFock = 0; iOFock < FOSav_.size(); iOFock++){
        (*FOSav_[iOFock]) = (*ssPropagator_->fockOrtho()[iOFock]);
      }
    }

    if (currentStep == ExplicitMagnus3) {
      // Copy the orthonormal density from ssPropagator to POSav1
      // POSav1(k) = PO(k)
      POSav1 = *ssPropagator_->onePDMOrthoScalar();
      // Copy the orthonormal fock from ssPropagator to FOSav1
      // FOSav1(k) = FO(k)
      FOSav1 = *ssPropagator_->fockOrtho()[0];
      // FO = 0.5 * FO 
      *ssPropagator_->fockOrtho()[0] *= dcomplex(0.5);
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
      FOSav2 = *ssPropagator_->fockOrtho()[0];
      FOSav3 = 0.25*(FOSav1 + FOSav2); 
     
      // update ssPropagator 
      *ssPropagator_->fockOrtho()[0] = FOSav3;
      *ssPropagator_->onePDMOrthoScalar() = POSav1;
      // Propagate and update Fock
      formUTrans();
      propDen();
      ssPropagator_->unOrthoDen3();
      formField(currentTime + 0.5*stepSize_);
      ssPropagator_->formFock();
      ssPropagator_->orthoFock3();

      // Form FOSav4
      FOSav4 = *ssPropagator_->fockOrtho()[0];
  
      //update ssPropagator 
      *ssPropagator_->fockOrtho()[0] = FOSav2;
      *ssPropagator_->onePDMOrthoScalar() = POSav1;
      // Propagate and update Fock
      formUTrans();
      propDen();
      ssPropagator_->unOrthoDen3();
      formField(currentTime + stepSize_);
      ssPropagator_->formFock();
      ssPropagator_->orthoFock3();

      // Form FOSav5
      FOSav5 = *ssPropagator_->fockOrtho()[0];

      // Final fock update
      dcomplex h = -dcomplex(0,1)*stepSize_;
      *ssPropagator_->fockOrtho()[0]  = (1/6.)*(FOSav1 + 4*FOSav4 + FOSav5);
      *ssPropagator_->fockOrtho()[0] -= (h/3.)*(FOSav3*FOSav4 - FOSav4*FOSav3);
      *ssPropagator_->fockOrtho()[0] -= (h/12.)*(FOSav2*FOSav5 - FOSav5*FOSav2);
      *ssPropagator_->onePDMOrthoScalar() = POSav1;
      // Propagate and update Fock
      formUTrans();
      propDen();
      ssPropagator_->unOrthoDen3();
      //formField(currentTime + stepSize_);
      //ssPropagator_->formFock();
      //ssPropagator_->orthoFock3();
    }

    // Increment the current time
    currentTime += stepSize_;
  }
  if(tarCSVs)  tarCSVFiles();
};

template <typename T>
void RealTime<T>::addRecord(const long double currentTime) {
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

  propInfo.push_back(rec);
}
