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
void Quantum<T>::rotateDensities(const std::array<double,3> &n, double theta){
  if(this->nTCS_ != 2) return;

  double sint = std::sin(theta);
  double cost = std::cos(theta);
  double omcost = 1 - cost;
  for(auto i = 0; i < this->onePDMScalar_->rows(); i++)
  for(auto j = 0; j < this->onePDMScalar_->rows(); j++) {
    T x = (*this->onePDMMx_)(i,j);
    T y = (*this->onePDMMy_)(i,j);
    T z = (*this->onePDMMz_)(i,j);

    (*this->onePDMMx_)(i,j)  = (cost + n[0]*n[0] * omcost) * x;
    (*this->onePDMMx_)(i,j) += (n[0]*n[1] * omcost - n[2] * sint) * y;
    (*this->onePDMMx_)(i,j) += (n[0]*n[2] * omcost + n[1] * sint) * z;
    
    (*this->onePDMMy_)(i,j)  = (n[0]*n[1] * omcost + n[2] * sint) * x;
    (*this->onePDMMy_)(i,j) += (cost + n[1]*n[1] * omcost) * y;
    (*this->onePDMMy_)(i,j) += (n[1]*n[2] * omcost - n[0] * sint) * z;

    (*this->onePDMMz_)(i,j)  = (n[0]*n[2] * omcost - n[1] * sint) * x;
    (*this->onePDMMz_)(i,j) += (n[1]*n[2] * omcost + n[0] * sint) * y;
    (*this->onePDMMz_)(i,j) += (cost + n[2]*n[2] * omcost) * z;
  }

};
