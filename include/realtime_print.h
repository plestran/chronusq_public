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
template<typename T>
void RealTime<T>::printRT() {
  double EDx = (*this->ssPropagator_->dipole())(0)/phys.debye;
  double EDy = (*this->ssPropagator_->dipole())(1)/phys.debye;
  double EDz = (*this->ssPropagator_->dipole())(2)/phys.debye;

  // FIXME: fix print so that all of the fields are consistent
  // i.e. Electric Dipole Onl and Modified mid-point
  if (currentTime_ <= 0.0) {
    this->fileio_->out << std::left << "\n *** Real-Time Time-Dependent Hartree-Fock ***" << endl; 
    this->fileio_->out << std::left << "\n * General Parameters:" << endl; 
    this->fileio_->out << std::setw(33) << std::right << " Number of steps =" 
                       << std::setw(14) << this->maxSteps_ << endl;
    this->fileio_->out << std::setw(33) << std::right << " Step size =" 
                       << std::setw(14) << std::setprecision(5) 
                       << this->stepSize_ << std::setw(5) << " a.u."
                       << std::setw(14) << std::setprecision(7) 
                       << this->stepSize_*phys.AuToFs 
                       << std::setw(3) << " fs" << endl;
    this->fileio_->out << std::setw(33) << std::right << " Orthonormalization method =" 
                       << std::setw(14) << std::right ;

    if (this->typeOrtho_ == Lowdin) this->fileio_->out << " L\u00F6wdin" << endl; 
    if (this->typeOrtho_ == Cholesky) this->fileio_->out << " Cholesky" << endl; 
    if (this->typeOrtho_ == Canonical) this->fileio_->out << " Canonical (NYI)" << endl; 

    this->fileio_->out << std::left << "\n * Integration Parameters:" << endl; 
    this->fileio_->out << std::right << std::setw(33) 
                       << " Electronic integration =" << std::setw(14) << std::right; 
    this->fileio_->out << " Modified mid-point" << endl;

    this->fileio_->out << std::left << "\n * External Field Parameters:" << endl; 
    this->fileio_->out << std::right << std::setw(33) 
                       << " Field type =" << std::setw(14) << std::right; 
    if (IEnvlp_ == Constant) this->fileio_->out << " Constant" << endl;
    if (IEnvlp_ == LinRamp) this->fileio_->out << " Linear ramp" << endl;
    if (IEnvlp_ == Gaussian) this->fileio_->out << " Gaussian" << endl;
    if (IEnvlp_ == Step) this->fileio_->out << " Step function" << endl;
    if (IEnvlp_ == SinSq) this->fileio_->out << " Sine-square (NYI)" << endl;
    this->fileio_->out << std::right << std::setw(33) 
                       << " Terms Included =" << std::setw(14) << std::right; 
    this->fileio_->out << " Electric Dipole Only" << endl;
    this->fileio_->out << std::right << std::setw(33) 
                       << " Gauge =" << std::setw(14) << std::right; 
    this->fileio_->out << " Length" << endl;
    this->fileio_->out << std::right << std::setw(33) 
                       << " Update Method =" << std::setw(14) << std::right; 
    this->fileio_->out << " MMUT" << endl;
 
    this->fileio_->out << std::setw(33) << std::right << " Ex =" 
                       << std::setw(14) << std::setprecision(5) 
                       << this->Ex_ << std::setw(5) << " a.u."
                       << endl; 
    this->fileio_->out << std::setw(33) << std::right << " Ey =" 
                       << std::setw(14) << std::setprecision(5) 
                       << this->Ey_ << std::setw(5) << " a.u."
                       << endl; 
    this->fileio_->out << std::setw(33) << std::right << " Ez =" 
                       << std::setw(14) << std::setprecision(5) 
                       << this->Ez_ << std::setw(5) << " a.u."
                       << endl; 
    if (this->IEnvlp_ != Step) {
      this->fileio_->out << std::setw(33) << std::right << " Frequency =" 
                         << std::setw(14) << std::setprecision(5) 
                         << this->Freq_ << std::setw(3) << " eV"
                         << endl;
      this->fileio_->out << std::setw(33) << std::right << " Phase =" 
                         << std::setw(14) << std::setprecision(5) 
                         << this->Phase_ << std::setw(4) << " rad"
                         << endl; 
    }
    this->fileio_->out << std::setw(33) << std::right << " Time On =" 
                       << std::setw(14) << std::setprecision(5) 
                       << this->TOn_ << std::setw(3) << " fs"
                       << endl; 
    this->fileio_->out << std::setw(33) << std::right << " Time Off =" 
                       << std::setw(14) << std::setprecision(5) 
                       << this->TOff_ << std::setw(3) << " fs"
                       << endl; 
    if (this->IEnvlp_ == Gaussian) {
      this->fileio_->out << std::setw(33) << std::right << " FWHM =" 
                         << std::setw(14) << std::setprecision(5) 
                         << this->Sigma_ << std::setw(3) << " eV"
                         << endl; 
    }
/*
  int           rtInitDensity;  // Initial Density for RT-TDSCF
  int           rtSwapMOA;      // which alpha MOs to swap for RT-TDSCF
  int           rtSwapMOB;      // which beta MOs to swap for RT-TDSCF
  int           rtMethFormU;    // Which method to form U propagator for RT-TDSCF
*/
                        
    this->fileio_->out << endl;
    this->fileio_->out << bannerTop << endl;
    this->fileio_->out << std::setw(11) << std::right << "Time" 
                       << std::setw(1)  << " "
                       << std::setw(16) << std::right << "Energy"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << "Ex Dipole"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " Ey Dipole" 
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " Ez Dipole"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " Etot Dipole"
                       << std::endl;
    this->fileio_->out << std::setw(11) << std::right << "    (au)" 
                       << std::setw(1)  << " "
                       << std::setw(16) << std::right << "    (Eh)"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " (debye)"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " (debye)" 
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " (debye)"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " (debye)"
                       << std::endl;
    this->fileio_->out << bannerMid << endl;
  }
  this->fileio_->out << std::setw(11) << std::right << std::setprecision(4) 
                     << currentTime_  << std::setw(1)  << " "
                     << std::setw(16) << std::right << std::setprecision(10) 
                     << this->ssPropagator_->totalEnergy << std::setw(1)  << " "
                     << std::setw(12) << std::right << std::setprecision(6) 
                     << EDx << std::setw(1)  << " "
                     << std::setw(12) << std::right << std::setprecision(6) 
                     << EDy << std::setw(1)  << " "
                     << std::setw(12) << std::right << std::setprecision(6) 
                     << EDz << std::setw(1)  << " "
                     << std::setw(12) << std::right << std::setprecision(6) <<
                        std::sqrt( std::pow(EDx,2.0) +
                                   std::pow(EDy,2.0) +
                                   std::pow(EDz,2.0))
                     << std::endl;
};


