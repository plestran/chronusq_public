template <typename T>
void RealTime<T>::doPropagation() {

  bool Start(false); // Start the MMUT iterations
  bool FinMM(false); // Wrap up the MMUT iterations

  currentTime_ = 0.0;

  initCSV();

  // Logic for Delta pulse
//if(this->tOff_ != 0 and (this->tOff_ - this->tOn_) < this->stepSize_){
//  this->tOff_ = this->stepSize_; 
//}
  this->printRTHeader();

  for(auto iStep = 0ul; iStep <= maxSteps_; iStep++) {

    // Logic for MMUT restart
    /*
    if(iRstrt_ > 0) {
      Start = (iStep == 0) or (iStep % iRstrt_) == 0;
      FinMM = (iStep + 1) % iRstrt_ == 0;
    } else {
      Start = (iStep == 0);
      FinMM = (iStep == maxSteps_);
    }
    */
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

      deltaT_ = FinMM ? stepSize_ : 2*stepSize_;
    }

    // Obtain field value for current time point
    formField();

    // Form AO Fock Matrix and compute properties
    // F(k) = H + G[P(k)]
    ssPropagator_->formFock();
    ssPropagator_->computeEnergy();
    ssPropagator_->computeProperties();
    this->printRTStep();
//  cout << iStep << " " << Start << " " << FinMM << " " << deltaT_ << 
//       " " << ssPropagator_->elecDipole()[0] << endl;

    // Orthonormalize the AO Fock
    // F(k) -> FO(k)
    ssPropagator_->orthoFock3();

    // Form the unitary propagation matrix
    // U**H(k) = exp(-i * dt * F(k))
    formUTrans();

    // Propagate the (orthonormal) Density Matrix
    // PO (in ssPropagator) now stores PO(k+1)
    // PO(k+1) = U**H(k) * PO(k-1) * U(k)
    propDen();

    // Add a time point record onto the list
    addRecord();

    // Write CSV Files
    writeCSVs();

    // Logic for MMUT restart setup
    if(false) {
      // PO = 0.5 * (PO + POSav)
      // PO(k+1) = 0.5 * (PO(k+1) + PO(k-1))
      for(auto iODen = 0; iODen < POSav_.size(); iODen++){
/*
        (*ssPropagator_->onePDMOrtho()[iODen]) *= 0.5;
        ssPropagator_->onePDMOrtho()[iODen]->noalias() = 
          (*ssPropagator_->onePDMOrtho()[iODen]) +
          0.5 * (*POSav_[iODen]) ;
*/
        (*ssPropagator_->onePDMOrtho()[iODen]) = 0.5 * (
          (*ssPropagator_->onePDMOrtho()[iODen]) + (*POSav_[iODen])
        );
      }
        
      prettyPrintSmart(cout,(*ssPropagator_->onePDMOrtho()[0])
      -0.5*(*ssPropagator_->onePDMOrtho()[0])*
       (*ssPropagator_->onePDMOrtho()[0]),"DIFF");

      // FIXME: Need to add McWeeny Purification
      ssPropagator_->McWeeny(ssPropagator_->onePDMOrtho(),10);

      // POSav = PO
      // POSav now stores the next step's PO
      for(auto iODen = 0; iODen < POSav_.size(); iODen++)
        (*POSav_[iODen]) = (*ssPropagator_->onePDMOrtho()[iODen]);
    }

    // Unorthonormalize the density for next Fock build
    // PO(k+1) -> P(k+1)
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
