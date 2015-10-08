
#include <workers.h>
#include <global.h>

using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::AOIntegrals;

namespace ChronusQ{
  void BasisSet::Wrapper_constructLocal(Molecule &mol){
    this->constructLocal(&mol);
  }
  void BasisSet::Wrapper_makeMaps(int nTCS, Molecule &mol){
    this->makeMaps(nTCS,&mol);
  }
}
