/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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

template<typename Op>
void computeElecDipole(const std::vector<Op> &DipoleOperators){
  std::vector<double> dple = 
    this-> template computeProperty<double,TOTAL>(DipoleOperators);
  for(std::size_t iXYZ = 0; iXYZ < 3; iXYZ++){
    this->elecDipole_[iXYZ] = -dple[iXYZ];
  }
}

template<typename Op>
void computeElecQuadpole(const std::vector<Op> &QuadpoleOperators){
  std::vector<double> qdple = 
    this-> template computeProperty<double,TOTAL>(QuadpoleOperators);
  for(std::size_t iXYZ = 0, iX = 0; iXYZ < 3; iXYZ++)
  for(std::size_t jXYZ = iXYZ     ; jXYZ < 3; jXYZ++, iX++){
    this->elecQuadpole_[iXYZ][jXYZ] = -qdple[iX];
  }

  RealMap QP( &this->elecQuadpole_[0][0],3,3);
  QP = QP.template selfadjointView<Upper>();

  RealMap TQP(&this->elecTracelessQuadpole_[0][0],3,3);
  
  TQP = QP - (RealMatrix::Identity(3,3) * QP.trace() / 3); 
}

template<typename Op>
void computeElecOctpole(const std::vector<Op> &OctpoleOperators){
  std::vector<double> ople = 
    this-> template computeProperty<double,TOTAL>(OctpoleOperators);
  for(std::size_t iXYZ = 0, iX = 0; iXYZ < 3; iXYZ++)
  for(std::size_t jXYZ = iXYZ     ; jXYZ < 3; jXYZ++)
  for(std::size_t kXYZ = jXYZ     ; kXYZ < 3; kXYZ++, iX++){
    this->elecOctpole_[iXYZ][jXYZ][kXYZ] = -ople[iX];
    this->elecOctpole_[iXYZ][kXYZ][jXYZ] = -ople[iX];
    this->elecOctpole_[jXYZ][iXYZ][kXYZ] = -ople[iX];
    this->elecOctpole_[jXYZ][kXYZ][iXYZ] = -ople[iX];
    this->elecOctpole_[kXYZ][iXYZ][jXYZ] = -ople[iX];
    this->elecOctpole_[kXYZ][jXYZ][iXYZ] = -ople[iX];
  }

}

template<typename Op>
void computeElecMultipoles(const std::vector<Op> &DipoleOps,
    const std::vector<Op> &QuadpoleOps,
    const std::vector<Op> &OctpoleOps){

  bool needToScatter = !this->isScattered_;
  this->scatterDensity();

  if(this->maxMultipole_ >= 1)
    this->template computeElecDipole(DipoleOps);
  if(this->maxMultipole_ >= 2)
    this->template computeElecQuadpole(QuadpoleOps);
  if(this->maxMultipole_ >= 3)
    this->template computeElecOctpole(OctpoleOps);

  if(needToScatter) this->gatherDensity();

}
