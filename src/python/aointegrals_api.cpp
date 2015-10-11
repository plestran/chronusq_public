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
  void AOIntegrals::Wrapper_iniAOIntegrals(Molecule &mol, BasisSet &basis,
         FileIO &fileio, Controls &controls) {
     this->iniAOIntegrals( &mol, &basis, &fileio, &controls); 
  }
};



