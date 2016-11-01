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
