template <typename T>
void RealTime<T>::printRTStep() {
  fileio_->out << std::setw(11) << std::right << std::setprecision(4) 
               << currentTime_  << std::setw(1)  << " ";

  fileio_->out << std::setw(16) << std::right << std::setprecision(10) 
               << this->ssPropagator_->totalEnergy() << std::setw(1)  << " ";

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

  fileio_->out << std::endl;
};
