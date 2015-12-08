#include <workers.h>
#include <global.h>
#include <atoms.h>

using namespace boost::python;
using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::AOIntegrals;
using ChronusQ::MOIntegrals;
using ChronusQ::RealTime;
using ChronusQ::CErr;
using ChronusQ::SDResponse;

namespace ChronusQ{
  void Wrapper_readInput(FileIO&,Molecule&,BasisSet&,Controls&,BasisSet&);

  void Wrapper_CErr_Default(FileIO &);
  void Wrapper_CErr_Message(FileIO &, std::string );
  inline void CQSetNumThreads(int n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#endif
     int h = 0;
  };
  int getAtomicNumber(int);
  void Wrapper_initCQ(int,boost::python::list);
}
