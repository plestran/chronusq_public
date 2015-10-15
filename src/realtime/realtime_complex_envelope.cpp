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
#include <realtime.h>
#include <aointegrals.h>
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;
using ChronusQ::RealTime;

namespace ChronusQ {

template<>
void RealTime<dcomplex>::formEDField() {  
  int IEnvlp   = this->IEnvlp_;
  double Ex    = this->Ex_;
  double Ey    = this->Ey_;
  double Ez    = this->Ez_;
  double TOn   = (this->TOn_)/phys.AuToFs;
  double TOff  = (this->TOff_)/phys.AuToFs;
  double Omega = (this->Freq_)/phys.eVPerHartree;
  double Phase = this->Phase_;
  double Sigma = (this->Sigma_)/phys.eVPerHartree;
  double Time  = currentTime_;
  double OmegT;
  double TMax;
  double TCntr;

  if (IEnvlp == 1) { 
  //   Constant envelope (plane wave)
    if(Time >= TOn && Time <= TOff) {
      OmegT = Omega*(Time-TOn) + Phase;
      this->EDField_[0] = Ex*std::cos(OmegT);
      this->EDField_[1] = Ey*std::cos(OmegT);
      this->EDField_[2] = Ez*std::cos(OmegT);
    } else {
      this->EDField_[0] = 0.0;
      this->EDField_[1] = 0.0;
      this->EDField_[2] = 0.0;
    }
  } 
  else if (IEnvlp == 2) { 
  /*   
       Linearly ramping up to the maximum in the first cycle
       Constant envelope afterwards
       Linearly ramping off to zero
  */
    if(Time >= TOn && Time <= TOff) {
      OmegT = Omega * (Time-TOn) + Phase;
      TMax = 2.0*math.pi/Omega;
      if(Time <= (TOn+TMax)) {
        this->EDField_[0] = Ex*((Time-TOn)/TMax)*std::cos(OmegT);
        this->EDField_[1] = Ey*((Time-TOn)/TMax)*std::cos(OmegT);
        this->EDField_[2] = Ez*((Time-TOn)/TMax)*std::cos(OmegT);
      }
      else if(Time > (TOn+TMax) && Time < (TOff-TMax)) {
        this->EDField_[0] = Ex*std::cos(OmegT);
        this->EDField_[1] = Ey*std::cos(OmegT);
        this->EDField_[2] = Ez*std::cos(OmegT);
      }
      else if(Time > (TOff-TMax)) {
        this->EDField_[0] = Ex*((TOff-Time)/TMax)*std::cos(OmegT);
        this->EDField_[1] = Ey*((TOff-Time)/TMax)*std::cos(OmegT);
        this->EDField_[2] = Ez*((TOff-Time)/TMax)*std::cos(OmegT);
      }
    } else {
      this->EDField_[0] = 0.0;
      this->EDField_[1] = 0.0;
      this->EDField_[2] = 0.0;
    }

  }
  else if (IEnvlp == 3) { 
  /*
       Gaussian envelope: E(t) = E * exp(-(a*(t-t0))^2) * Sin(wt)
           a  = the range of frequency (FWHM) 
           t0 = the time when the amplitude reaches maximum 
                  the default for t0 = sqrt(LN(1000))/a 
                  (at t=0, the amplitude is 1/1000 times maximum. 
                  This ensures a smooth turning-on of the field)             
        Dipole length approximation only
  */
    if(Time >= TOn && Time <= TOff) {
      OmegT = Omega * (Time-TOn) + Phase;
      TCntr = std::sqrt(std::log(1000.0))/Sigma;
      this->EDField_[0] = Ex*std::cos(OmegT)*std::exp(-std::pow(Sigma*(Time-TCntr),2.0));
      this->EDField_[1] = Ey*std::cos(OmegT)*std::exp(-std::pow(Sigma*(Time-TCntr),2.0));
      this->EDField_[2] = Ez*std::cos(OmegT)*std::exp(-std::pow(Sigma*(Time-TCntr),2.0));
    } else {
      this->EDField_[0] = 0.0;
      this->EDField_[1] = 0.0;
      this->EDField_[2] = 0.0;
    }

  }
  else if (IEnvlp == 4) { 
  /*
        Step function
  */
    if(Time >= TOn && Time <= TOff) {
      this->EDField_[0] = Ex;
      this->EDField_[1] = Ey;
      this->EDField_[2] = Ez;
    } else {
      this->EDField_[0] = 0.0;
      this->EDField_[1] = 0.0;
      this->EDField_[2] = 0.0;
    }

  }
  else if (IEnvlp == 5) { 
  /*
        Sine-square envelope: 
          A(r,t) = (E0/w)*(sin(Pi*t/T))^2*sin(k*r-w*t) 
          E(r,t) = E0*(sin(Pi*t/T))^2*cos(k*r-w*t)
                    -(E0*Pi/w*T)*sin(2*Pi*t/T)*sin(k*r-w*t)
        only dipole length/velocity and quardrupole velocity
  */
    if(Time >= TOn && Time <= TOff) {
      OmegT = Omega * (Time-TOn) + Phase;
    } else {
      this->EDField_[0] = 0.0;
      this->EDField_[1] = 0.0;
      this->EDField_[2] = 0.0;
    }
  }
};

} // namespace ChronusQ
