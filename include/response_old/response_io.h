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
template<typename T>
void Response<T>::printInfo(){
  this->fileio_->out << bannerTop << endl;
  
  this->fileio_->out << "Response Module Settings:" << endl;
  this->fileio_->out << endl;

  this->fileio_->out << std::setw(35) << std::left <<"  Response Method (Common):"; 
  if(this->iMeth_ == RESPONSE_TYPE::CIS)
    this->fileio_->out << "Configuration Interaction Singles (CIS)";
  else if(this->iMeth_ == RESPONSE_TYPE::RPA)
    this->fileio_->out << "Random Phase Approximation (RPA)";
  else if(this->iMeth_ == RESPONSE_TYPE::STAB)
    this->fileio_->out << "Stability";
  else if(this->iMeth_ == RESPONSE_TYPE::PPRPA)
    this->fileio_->out << "Particle-Particle Random Phase Approximation";
  else if(this->iMeth_ == RESPONSE_TYPE::PPTDA)
    this->fileio_->out << "Particle-Particle Tamm-Dancoff Approximation";
  this->fileio_->out << endl;


  this->fileio_->out << std::setw(35) << std::left <<"  Response Function:"; 
  if(this->iClass_ == RESPONSE_CLASS::FOPPA)
    this->fileio_->out << "First-Order Polarization Propagator";
  else if(this->iClass_ == RESPONSE_CLASS::SOPPA)
    this->fileio_->out << "Second-Order Polarization Propagator";
  else if(this->iClass_ == RESPONSE_CLASS::TOPPA)
    this->fileio_->out << "Third-Order Polarization Propagator";
  else if(this->iClass_ == RESPONSE_CLASS::PPPA)
    this->fileio_->out << "Particle-Particle Polarization Propagator";
  this->fileio_->out << endl;

  this->fileio_->out << std::setw(35) << std::left <<"  Tamm-Dancoff Approximation:";
  if(this->doTDA_) this->fileio_->out << "Yes";
  else             this->fileio_->out << "No";
  this->fileio_->out << endl;

  this->fileio_->out << std::setw(35) << std::left <<"  Partitioning Scheme:";
  if(this->iPart_ == RESPONSE_MATRIX_PARTITIONING::SPIN_SEPARATED)
    this->fileio_->out << "Spin Separated";
  else if(this->iPart_ == RESPONSE_MATRIX_PARTITIONING::SPIN_ADAPTED)
    this->fileio_->out << "Spin Adapted";
  this->fileio_->out << endl;

  this->fileio_->out << std::setw(35) << std::left <<"  Number of Response Matricies:" << iMatIter_.size() << endl;

  this->fileio_->out << std::setw(35) << std::left <<"  Job Type:";
  if(this->iJob_ == RESPONSE_JOB_TYPE::EIGEN)
    this->fileio_->out << "Eigenspectrum";
  else if(this->iJob_ == RESPONSE_JOB_TYPE::DYNAMIC)
    this->fileio_->out << "Dynamic Response";
  this->fileio_->out << endl;

  this->fileio_->out << endl;
  this->fileio_->out << bannerEnd << endl;

} // printInfo

