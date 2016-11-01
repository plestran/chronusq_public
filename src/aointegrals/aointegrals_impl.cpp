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
