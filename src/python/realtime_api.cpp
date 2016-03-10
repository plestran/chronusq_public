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
#include <workers.h>
#include <global.h>

using ChronusQ::RealTime;

namespace ChronusQ {

/*
  template<>
  boost::python::list RealTime<double>::Wrapper_recs(){
    boost::python::list result;
    for(auto i = this->propInfo.begin(); i != this->propInfo.end(); i++){
      Wrapper_PropInfo tmpRec(*i);
      result.append(tmpRec);
    }
    return result;
  }

  template<>
  boost::python::list RealTime<dcomplex>::Wrapper_recs(){
    boost::python::list result;
    for(auto i = this->propInfo.begin(); i != this->propInfo.end(); i++){
      Wrapper_PropInfo tmpRec(*i);
      result.append(tmpRec);
    }
    return result;
  }
*/

  template<>
  boost::python::list RealTime<double>::lastDipole(){
    boost::python::list result;
    for(auto i = 0; i < 4; i++)
      result.append(this->propInfo[this->propInfo.size()-1].dipole[i]);
    return result;
  }

  template<>
  double RealTime<double>::lastEnergy(){
    return this->propInfo[this->propInfo.size()-1].energy;
  }

  template<>
  double RealTime<double>::getTimeStep(){
    return this->propInfo[this->propInfo.size()-1].timeStep;
  }
  template<>
  boost::python::list RealTime<dcomplex>::lastDipole(){
    boost::python::list result;
    for(auto i = 0; i < 4; i++)
      result.append(this->propInfo[this->propInfo.size()-1].dipole[i]);
    return result;
  }

  template<>
  double RealTime<dcomplex>::lastEnergy(){
    return this->propInfo[this->propInfo.size()-1].energy;
  }

  template<>
  double RealTime<dcomplex>::getTimeStep(){
    return this->propInfo[this->propInfo.size()-1].timeStep;
  }
  



}; // namespace ChronusQ

