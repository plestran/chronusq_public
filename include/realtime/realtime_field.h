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
void RealTime<T>::formField(const long double currentTime) {
/*
  double TOn   = tOn_   / phys.AuToFs;
  double TOff  = tOff_  / phys.AuToFs;
  double Omega = freq_  / phys.eVPerHartree;
  double Sigma = sigma_ / phys.eVPerHartree;
*/
  double TOn   = tOn_  ; 
  double TOff  = tOff_ ; 
  double Omega = freq_ ; 
  double Sigma = sigma_; 

  double omegaT = Omega * (currentTime - TOn) + phase_;
  omegaT = std::cos(omegaT);

  // Constant Envelope (Plane Wave)
  if( iEnvlp_ == Constant ) {
    if( currentTime >= TOn and currentTime <= TOff ) {
      EDField_ = {{ staticEDAmp_[0] * omegaT ,
                    staticEDAmp_[1] * omegaT ,
                    staticEDAmp_[2] * omegaT }};
    } else {
      EDField_ = {{0.0,0.0,0.0}};
    }
  } else if( iEnvlp_ == LinRamp ) {
    double TMax = 2.0 * math.pi / Omega;
    if( currentTime <= (TOn + TMax) ) {
      EDField_ = {{ staticEDAmp_[0] * ((currentTime-TOn)/TMax) * omegaT ,
                    staticEDAmp_[1] * ((currentTime-TOn)/TMax) * omegaT ,
                    staticEDAmp_[2] * ((currentTime-TOn)/TMax) * omegaT }};
    } else if( currentTime > (TOn + TMax) and currentTime < (TOff - TMax) ) {
      EDField_ = {{ staticEDAmp_[0] * omegaT ,
                    staticEDAmp_[1] * omegaT ,
                    staticEDAmp_[2] * omegaT }};
    } else if( currentTime > (TOff - TMax)) {
      EDField_ = {{ staticEDAmp_[0] * ((TOff-currentTime)/TMax) * omegaT ,
                    staticEDAmp_[1] * ((TOff-currentTime)/TMax) * omegaT ,
                    staticEDAmp_[2] * ((TOff-currentTime)/TMax) * omegaT }};
    } else {
      EDField_ = {{0.0,0.0,0.0}};
    } 
  } else if( iEnvlp_ == Gaussian ) {
    if( currentTime >= TOn and currentTime <= TOff ) {
      double tCntr  = std::sqrt(std::log(1.0e3)) / Sigma;
      double gauFac = std::exp(-std::pow(Sigma * (currentTime - tCntr),2.0));
      EDField_ = {{ staticEDAmp_[0] * omegaT * gauFac,
                    staticEDAmp_[1] * omegaT * gauFac,
                    staticEDAmp_[2] * omegaT * gauFac}};
    } else {
      EDField_ = {{0.0,0.0,0.0}};
    }
  } else if( iEnvlp_ == Step ) {
    if( currentTime >= TOn and currentTime <= TOff ) {
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

    if( currentTime >= TOn and currentTime <= TOff ) {
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
