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
using ChronusQ::RealTime;
using ChronusQ::CErr;

namespace ChronusQ{
  void Wrapper_readInput(FileIO&,Molecule&,BasisSet&,Controls&,BasisSet&);

  void Wrapper_CErr_Default(FileIO &);
  void Wrapper_CErr_Message(FileIO &, std::string &);
  int getAtomicNumber(int);
}
