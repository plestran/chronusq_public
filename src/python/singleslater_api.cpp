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
};


