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
#include <aointegrals/aointegrals_impl.h>

template void ChronusQ::AOIntegrals::newTwoEContractDirect(
      const std::vector<std::reference_wrapper<RealMatrix>>&,
      std::vector<std::reference_wrapper<RealMatrix>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<double>&) ;
template void ChronusQ::AOIntegrals::newTwoEContractDirect(
      const std::vector<std::reference_wrapper<ComplexMatrix>>&,
      std::vector<std::reference_wrapper<ComplexMatrix>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<double>&) ;
template void ChronusQ::AOIntegrals::newTwoEContractDirect(
      const std::vector<std::reference_wrapper<RealMap>>&,
      std::vector<std::reference_wrapper<RealMap>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<double>&) ;
template void ChronusQ::AOIntegrals::newTwoEContractDirect(
      const std::vector<std::reference_wrapper<ComplexMap>>&,
      std::vector<std::reference_wrapper<ComplexMap>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<double>&) ;
template void ChronusQ::AOIntegrals::newTwoEContractIncore(
      const std::vector<std::reference_wrapper<RealMatrix>>&,
      std::vector<std::reference_wrapper<RealMatrix>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<double>&) ;
template void ChronusQ::AOIntegrals::newTwoEContractIncore(
      const std::vector<std::reference_wrapper<ComplexMatrix>>&,
      std::vector<std::reference_wrapper<ComplexMatrix>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<double>&) ;
template void ChronusQ::AOIntegrals::newTwoEContractIncore(
      const std::vector<std::reference_wrapper<RealMap>>&,
      std::vector<std::reference_wrapper<RealMap>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<double>&) ;
template void ChronusQ::AOIntegrals::newTwoEContractIncore(
      const std::vector<std::reference_wrapper<ComplexMap>>&,
      std::vector<std::reference_wrapper<ComplexMap>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<double>&) ;

template void ChronusQ::AOIntegrals::Ortho1Trans(RealMatrix&,RealMatrix&);
template void ChronusQ::AOIntegrals::Ortho1Trans(ComplexMatrix&,ComplexMatrix&);
template void ChronusQ::AOIntegrals::Ortho1Trans(RealMap&,RealMap&);
template void ChronusQ::AOIntegrals::Ortho1Trans(ComplexMap&,ComplexMap&);
template void ChronusQ::AOIntegrals::Ortho1TransT(RealMatrix&,RealMatrix&);
template void ChronusQ::AOIntegrals::Ortho1TransT(ComplexMatrix&,ComplexMatrix&);
template void ChronusQ::AOIntegrals::Ortho1TransT(RealMap&,RealMap&);
template void ChronusQ::AOIntegrals::Ortho1TransT(ComplexMap&,ComplexMap&);
template void ChronusQ::AOIntegrals::Ortho2Trans(RealMatrix&,RealMatrix&);
template void ChronusQ::AOIntegrals::Ortho2Trans(ComplexMatrix&,ComplexMatrix&);
template void ChronusQ::AOIntegrals::Ortho2Trans(RealMap&,RealMap&);
template void ChronusQ::AOIntegrals::Ortho2Trans(ComplexMap&,ComplexMap&);
//template void ChronusQ::AOIntegrals::Ortho2TransT(RealMatrix&,RealMatrix&);
//template void ChronusQ::AOIntegrals::Ortho2TransT(ComplexMatrix&,ComplexMatrix&);
//template void ChronusQ::AOIntegrals::Ortho2TransT(RealMap&,RealMap&);
//template void ChronusQ::AOIntegrals::Ortho2TransT(ComplexMap&,ComplexMap&);
