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

using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::AOIntegrals;

namespace ChronusQ {
  template<>
  void SingleSlater<double>::Wrapper_iniSingleSlater(Molecule &mol, BasisSet &basis,
         AOIntegrals &ints, FileIO &fileio, Controls &controls) {
     this->iniSingleSlater( &mol, &basis, &ints, &fileio, &controls); 
  }
  template<>
  void SingleSlater<dcomplex>::Wrapper_iniSingleSlater(Molecule &mol, BasisSet &basis,
         AOIntegrals &ints, FileIO &fileio, Controls &controls) {
     this->iniSingleSlater( &mol, &basis, &ints, &fileio, &controls); 
  }

  template<>
  boost::python::list SingleSlater<double>::Wrapper_dipole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++)
      result.append(this->elecDipole_[i]);
    return result;
  }

  template<>
  boost::python::list SingleSlater<double>::Wrapper_quadrupole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++) {
      boost::python::list oneDim;
      for(auto j = 0; j < 3; j++) 
        oneDim.append(this->elecQuadpole_[i][j]);
      result.append(oneDim);
    }
    return result;
  }

  template<>
  boost::python::list SingleSlater<double>::Wrapper_octupole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++) {
      boost::python::list twoDim;
      for(auto j = 0; j < 3; j++) {
        boost::python::list oneDim;
        for(auto k = 0; k < 3; k++)
          oneDim.append(this->elecOctpole_[i][j][k]);
        twoDim.append(oneDim);
      }
      result.append(twoDim);
    }
    return result;
  }

  template<>
  boost::python::list SingleSlater<dcomplex>::Wrapper_dipole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++)
      result.append(this->elecDipole_[i]);
    return result;
  }

  template<>
  boost::python::list SingleSlater<dcomplex>::Wrapper_quadrupole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++) {
      boost::python::list oneDim;
      for(auto j = 0; j < 3; j++) 
        oneDim.append(this->elecQuadpole_[i][j]);
      result.append(oneDim);
    }
    return result;
  }

  template<>
  boost::python::list SingleSlater<dcomplex>::Wrapper_octupole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++) {
      boost::python::list twoDim;
      for(auto j = 0; j < 3; j++) {
        boost::python::list oneDim;
        for(auto k = 0; k < 3; k++)
          oneDim.append(this->elecOctpole_[i][j][k]);
        twoDim.append(oneDim);
      }
      result.append(twoDim);
    }
    return result;
  }
};


