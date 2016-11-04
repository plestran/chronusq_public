template <typename T>
void RealTime<T>::printRTHeader() {
  std::ostream &output = this->fileio_->out;

  output << BannerTop << endl;
  output << "Real-Time Propagation Settings:" << endl << endl;

  output << std::setw(38) << std::left << "  * Simulation Parameters:" 
         << endl; 
  output << std::setw(38) << std::left << "  Number of Steps:" 
                          << this->maxSteps_ << endl;
  output << std::setw(38) << std::left << "  Step Size:" 
                          << std::setprecision(5) 
                          << this->stepSize_ << " \u0127 / Eh" << endl;
  output << std::setw(38) << std::left << " " 
                          << std::setprecision(5) 
                          << this->stepSize_*phys.AuToFs << " fs" << endl;
  if(this->iRstrt_ > 0) {
    output << std::setw(38) << std::left << "  Restarting MMUT every:" 
                            << this->iRstrt_ << " steps" << endl;
  }

  output << endl;
  output << std::setw(38) << std::left << "  * Integration Parameters:" 
         << endl; 

  output << std::setw(38) << std::left << "  Electronic Integration:" 
                          << "Modified Midpoint Unitary Transformation (MMUT)" 
                          << endl;

  output << endl;
  output << std::setw(38) << std::left << "  * External Field Parameters:" 
         << endl; 

  output << std::setw(38) << std::left << "  Field Envelope:";
  if (iEnvlp_ == Constant)      output << "Constant";
  else if (iEnvlp_ == LinRamp)  output << "Linear Ramp";
  else if (iEnvlp_ == Gaussian) output << "Gaussian";
  else if (iEnvlp_ == Step)     output << "Step Function";
  output << endl;

  output << std::setw(38) << std::left << "  Field Components Included:";
  output << "Electric Dipole Only" << endl;

  output << std::setw(38) << std::left << "  Gauge:";
  output << "Length" << endl;

  output << std::setw(38) << std::left << "  Electric Field Amplitude (Dipole):"
         << "{" << this->staticEDAmp_[0] << ", " << this->staticEDAmp_[1] << ", "
         << this->staticEDAmp_[2] << "} Eh / (e bohr)" << endl;

  output << std::setw(38) << std::left << "  Frequency:" 
         << this->freq_  << " Eh / \u0127 " << endl; 

  output << std::setw(38) << std::left << "  Time On:" 
         << this->tOn_ << endl;
  output << std::setw(38) << std::left << "  Time Off:" 
         << this->tOff_ << endl;

  if(this->iEnvlp_ == Gaussian) {
    output << std::setw(38) << std::left << "  FWHM:" 
           << this->sigma_ << endl;
  }

  output << BannerEnd << endl;

  output << endl;
  output << std::right;
  output << std::setw(11) << "Time (a.u.)" << " ";
  output << std::setw(16) << "Energy (Eh)" << " ";
  output << std::setw(12) << "Ex (a.u.)" << " ";
  output << std::setw(12) << "Ey (a.u.)" << " ";
  output << std::setw(12) << "Ez (a.u.)" << " ";
  output << std::setw(12) << "|F| (a.u.)" << " ";
  output << endl << bannerTop << endl;

}

template <typename T>
void RealTime<T>::printRTStep() {
  fileio_->out << std::setw(11) << std::right << std::setprecision(4) 
               << currentTime_  << std::setw(1)  << " ";

  fileio_->out << std::setw(16) << std::right << std::setprecision(10) 
               << this->ssPropagator_->totalEnergy() << std::setw(1)  << " ";

/*
  fileio_->out << std::setw(12) << std::right << std::setprecision(6) 
               << EDField_[0] << std::setw(1)  << " ";

  fileio_->out << std::setw(12) << std::right << std::setprecision(6) 
               << EDField_[1] << std::setw(1)  << " ";

  fileio_->out << std::setw(12) << std::right << std::setprecision(6) 
               << EDField_[2] << std::setw(1)  << " ";

  fileio_->out << std::setw(12) << std::right << std::setprecision(6) <<
                  std::sqrt( std::pow(EDField_[0],2.0) +
                             std::pow(EDField_[1],2.0) +
                             std::pow(EDField_[2],2.0));
*/
  fileio_->out << std::setw(12) << std::right << std::setprecision(6) 
               << this->ssPropagator_->elecDipole()[0] << std::setw(1)  << " ";

  fileio_->out << std::setw(12) << std::right << std::setprecision(6) 
               << this->ssPropagator_->elecDipole()[1] << std::setw(1)  << " ";

  fileio_->out << std::setw(12) << std::right << std::setprecision(6) 
               << this->ssPropagator_->elecDipole()[2] << std::setw(1)  << " ";

/*
  fileio_->out << std::setw(12) << std::right << std::setprecision(6) <<
                  std::sqrt( std::pow(this->ssPropagator_->elecDipole()[0],2.0) +
                             std::pow(this->ssPropagator_->elecDipole()[1],2.0) +
                             std::pow(this->ssPropagator_->elecDipole()[2],2.0));
*/
  fileio_->out << std::setw(12) << std::right << std::setprecision(6) <<
                  std::sqrt( std::pow(EDField_[0],2.0) +
                             std::pow(EDField_[1],2.0) +
                             std::pow(EDField_[2],2.0));

  fileio_->out << std::endl;
};
