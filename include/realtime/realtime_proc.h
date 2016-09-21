template <typename T>
void RealTime<T>::doPropagation() {

  bool Start; // Start the MMUT iterations
  bool FinMM; // Wrap up the MMUT iterations

  currentTime_ = 0.0;

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
    prettyPrintSmart(cout,*ssPropagator_->onePDMScalar(),"PS");
//  prettyPrintSmart(cout,*ssPropagator_->fockScalar(),"FS");
    cout << std::setprecision(13) << ssPropagator_->totalEnergy() << endl;

    // Orthonormalize the AO Fock
    ssPropagator_->orthoFock3();

    // Form the unitary propagation matrix
    formUTrans();

    // Propagate the (orthonormal) Density Matrix
    propDen();

    // Unorthonormalize the density ofr next Fock build
//  ssPropagator_->cpyAOtoOrthoDen();
    ssPropagator_->unOrthoDen3();

    // Increment the current time
    currentTime_ += stepSize_;
  }
};
