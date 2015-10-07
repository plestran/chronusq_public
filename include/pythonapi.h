#include <workers.h>
#include <global.h>
#include <atoms.h>
#include <boost/python.hpp>

using namespace boost::python;
using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::AOIntegrals;

namespace ChronusQ{
  void Wrapper_readInput(FileIO&,Molecule&,BasisSet&,Controls&,BasisSet&);
  int getAtomicNumber(int);
}
