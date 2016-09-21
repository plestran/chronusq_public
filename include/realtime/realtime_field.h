template <typename T>
void RealTime<T>::formField() {
  double TOn   = tOn_   / phys.AuToFs;
  double TOff  = tOff_  / phys.AuToFs;
  double Omega = freq_  / phys.eVPerHartree;
  double Sigma = sigma_ / phys.eVPerHartree;

  double omegaT = Omega * (currentTime_ - TOn) + phase_;
  omegaT = std::cos(omegaT);

  // Constant Envelope (Plane Wave)
  if( iEnvlp_ == Constant ) {
    if( currentTime_ >= TOn and currentTime_ <= TOff ) {
      EDField_ = {{ staticEDAmp_[0] * omegaT ,
                    staticEDAmp_[1] * omegaT ,
                    staticEDAmp_[2] * omegaT }};
    } else {
      EDField_ = {{0.0,0.0,0.0}};
    }
  } else if( iEnvlp_ == LinRamp ) {
    double TMax = 2.0 * math.pi / Omega;
    if( currentTime_ <= (TOn + TMax) ) {
      EDField_ = {{ staticEDAmp_[0] * ((currentTime_-TOn)/TMax) * omegaT ,
                   staticEDAmp_[1] * ((currentTime_-TOn)/TMax) * omegaT ,
                   staticEDAmp_[2] * ((currentTime_-TOn)/TMax) * omegaT }};
    } else if( currentTime_ > (TOn + TMax) and currentTime_ < (TOff - TMax) ) {
      EDField_ = {{ staticEDAmp_[0] * omegaT ,
                   staticEDAmp_[1] * omegaT ,
                   staticEDAmp_[2] * omegaT }};
    } else if( currentTime_ > (TOff - TMax)) {
      EDField_ = {{ staticEDAmp_[0] * ((TOff-currentTime_)/TMax) * omegaT ,
                   staticEDAmp_[1] * ((TOff-currentTime_)/TMax) * omegaT ,
                   staticEDAmp_[2] * ((TOff-currentTime_)/TMax) * omegaT }};
    } else {
      EDField_ = {{0.0,0.0,0.0}};
    } 
  } else if( iEnvlp_ == Gaussian ) {
    if( currentTime_ >= TOn and currentTime_ <= TOff ) {
      double tCntr  = std::sqrt(std::log(1.0e3)) / Sigma;
      double gauFac = std::exp(-std::pow(Sigma * (currentTime_ - tCntr),2.0));
      EDField_ = {{ staticEDAmp_[0] * omegaT * gauFac,
                   staticEDAmp_[1] * omegaT * gauFac,
                   staticEDAmp_[2] * omegaT * gauFac}};
    } else {
      EDField_ = {{0.0,0.0,0.0}};
    }
  } else if( iEnvlp_ == Step ) {
    if( currentTime_ >= TOn and currentTime_ <= TOff ) {
      EDField_ = {{ staticEDAmp_[0] ,
                   staticEDAmp_[1] ,
                   staticEDAmp_[2] }};
    } else {
      EDField_ = {{0.0,0.0,0.0}};
    }
  } else if( iEnvlp_ == Elliptic ) {
    EDField_ = {{0.0,0.0,0.0}};
    omegaT = std::acos(omegaT);
    double cosOT = std::cos(omegaT);
    double sinOT = std::sin(omegaT);

    if( currentTime_ >= TOn and currentTime_ <= TOff ) {
      if( iEllPol_ == RXY ) {
        EDField_[0] = staticEDAmp_[0] * cosOT;
        EDField_[1] = staticEDAmp_[1] * sinOT;
      } else if( iEllPol_ == RXZ ) {
        EDField_[0] = staticEDAmp_[0] * cosOT;
        EDField_[2] = staticEDAmp_[2] * sinOT;
      } else if( iEllPol_ == RYZ ) {
        EDField_[1] = staticEDAmp_[1] * cosOT;
        EDField_[2] = staticEDAmp_[2] * sinOT;
      } else if( iEllPol_ == LXY ) {
        EDField_[0] = staticEDAmp_[0] * sinOT;
        EDField_[1] = staticEDAmp_[1] * cosOT;
      } else if( iEllPol_ == LXZ ) {
        EDField_[0] = staticEDAmp_[0] * sinOT;
        EDField_[2] = staticEDAmp_[2] * cosOT;
      } else if( iEllPol_ == LYZ ) {
        EDField_[1] = staticEDAmp_[1] * sinOT;
        EDField_[2] = staticEDAmp_[2] * cosOT;
      }
    } 
  }

  // Copy field over to SS Propagator
  ssPropagator_->setField(EDField_);
};
