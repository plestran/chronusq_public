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
