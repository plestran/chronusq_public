template <typename T>
void RealTime<T>::doPropagation() {

  bool Start; // Start the MMUT iterations
  bool FinMM; // Wrap up the MMUT iterations

  currentTime_ = 0.0;

  initCSV();

  for(auto iStep = 0ul; iStep < maxSteps_; iStep++) {

    // Logic for MMUT restart
    if(iRstrt_ > 0) {
      Start = (iStep == 0) or (iStep % iRstrt_) == 0;
      FinMM = (iStep + 1) % iRstrt_ == 0;
    } else {
      Start = (iStep == 0);
      FinMM = (iStep == maxSteps_);
    }

    // Initial entry for MMUT
    if(Start) {
      deltaT_ = stepSize_;
    // Subsequent MMUT iterations
    } else {
      deltaT_ = FinMM ? stepSize_ : 2*stepSize_;
    }

    // Obtain field value for current time point
    formField();

    // Form AO Fock Matrix and compute properties
    ssPropagator_->formFock();
    ssPropagator_->computeEnergy();
    ssPropagator_->computeProperties();
//  prettyPrintSmart(cout,*ssPropagator_->onePDMScalar(),"PS");
//  prettyPrintSmart(cout,*ssPropagator_->fockScalar(),"FS");
    this->printRTStep();

    // Orthonormalize the AO Fock
    ssPropagator_->orthoFock3();

    // Form the unitary propagation matrix
    formUTrans();

    // Propagate the (orthonormal) Density Matrix
    propDen();

    // Add a time point record onto the list
    addRecord();

    // Write CSV Files
    writeCSVs();

    // Unorthonormalize the density ofr next Fock build
//  ssPropagator_->cpyAOtoOrthoDen();
    ssPropagator_->unOrthoDen3();

    // Increment the current time
    currentTime_ += stepSize_;
  }
  if(tarCSVs)  tarCSVFiles();
};

template <typename T>
void RealTime<T>::addRecord() {
  PropInfo rec;
  rec.timeStep  = currentTime_;
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
