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
void Response<T>::formDensity(){
  if(this->iDen_ == TRANSITION)      this->formTransitionDensity();
  else if(this->iDen_ == DIFFERENCE) this->formDifferenceDensity();
};

template<typename T>
void Response<T>::formTransitionDensity(){
  if(this->iClass_ == RESPONSE_CLASS::FOPPA) 
    this->formTransitionDensityFOPPA();
  else if(this->iClass_ == RESPONSE_CLASS::PPPA)
    this->formTransitionDensityPPRPA();
};

template<typename T>
void Response<T>::formDifferenceDensity(){
  if(this->iClass_ == RESPONSE_CLASS::FOPPA) 
    this->formDifferenceDensityFOPPA();
  else if(this->iClass_ == RESPONSE_CLASS::PPPA)
    this->formDifferenceDensityPPRPA();
};
