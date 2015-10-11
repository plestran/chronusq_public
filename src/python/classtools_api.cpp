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
  void Wrapper_readInput(FileIO &fileio, Molecule &mol, BasisSet &basis, 
    Controls &controls, BasisSet &dfBasis) {
    readInput(&fileio,&mol,&basis,&controls,&dfBasis);
  }
};

