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
#include <quantum.h>
namespace ChronusQ {
  #include <quantum/quantum_scattergather.h>
  #include <quantum/quantum_misc.h>
};

// Explicitly instantiate Quantum Classes
template class ChronusQ::Quantum<double>;
template class ChronusQ::Quantum<dcomplex>;

// Instantiate Scatter / Gather
template void ChronusQ::Quantum<double>::spinScatter(RealMatrix&,RealMatrix&,
  std::vector<std::reference_wrapper<RealMatrix>>&);
template void ChronusQ::Quantum<double>::spinScatter(RealMap&,RealMap&,
  std::vector<std::reference_wrapper<RealMap>>&);
template void ChronusQ::Quantum<double>::spinScatter(RealMatrix&,
  std::vector<std::reference_wrapper<RealMatrix>>&);
template void ChronusQ::Quantum<double>::spinScatter(RealMap&,
  std::vector<std::reference_wrapper<RealMap>>&);

template void ChronusQ::Quantum<dcomplex>::spinScatter(ComplexMatrix&,
  ComplexMatrix&,std::vector<std::reference_wrapper<ComplexMatrix>>&);
template void ChronusQ::Quantum<dcomplex>::spinScatter(ComplexMap&,ComplexMap&,
  std::vector<std::reference_wrapper<ComplexMap>>&);
template void ChronusQ::Quantum<dcomplex>::spinScatter(ComplexMatrix&,
  std::vector<std::reference_wrapper<ComplexMatrix>>&);
template void ChronusQ::Quantum<dcomplex>::spinScatter(ComplexMap&,
  std::vector<std::reference_wrapper<ComplexMap>>&);

template void ChronusQ::Quantum<double>::spinGather(RealMatrix&,RealMatrix&,
  std::vector<std::reference_wrapper<RealMatrix>>&);
template void ChronusQ::Quantum<double>::spinGather(RealMap&,RealMap&,
  std::vector<std::reference_wrapper<RealMap>>&);
template void ChronusQ::Quantum<double>::spinGather(RealMatrix&,
  std::vector<std::reference_wrapper<RealMatrix>>&);
template void ChronusQ::Quantum<double>::spinGather(RealMap&,
  std::vector<std::reference_wrapper<RealMap>>&);

template void ChronusQ::Quantum<dcomplex>::spinGather(ComplexMatrix&,
  ComplexMatrix&,std::vector<std::reference_wrapper<ComplexMatrix>>&);
template void ChronusQ::Quantum<dcomplex>::spinGather(ComplexMap&,ComplexMap&,
  std::vector<std::reference_wrapper<ComplexMap>>&);
template void ChronusQ::Quantum<dcomplex>::spinGather(ComplexMatrix&,
  std::vector<std::reference_wrapper<ComplexMatrix>>&);
template void ChronusQ::Quantum<dcomplex>::spinGather(ComplexMap&,
  std::vector<std::reference_wrapper<ComplexMap>>&);


// Instantiate Misc functions
//template void ChronusQ::Quantum<double>::rotateDensities(
//  const std::array<double,3>&, double);
//template void ChronusQ::Quantum<dcomplex>::rotateDensities(
//  const std::array<double,3>&, double);
