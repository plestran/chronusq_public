#include <boost/python.hpp>
#include <singleslater.h>
using namespace boost::python;
using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;

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

BOOST_PYTHON_MODULE(libpythonapi){
  class_<SingleSlater<double>,boost::noncopyable>("SingleSlater_double",init<>())
    .def("iniSingleSlater", &SingleSlater<double>::Wrapper_iniSingleSlater)
    .def("printInfo"      , &SingleSlater<double>::printInfo              )
    .def("SCF"            , &SingleSlater<double>::SCF                    )
  ;

  class_<Molecule,boost::noncopyable>("Molecule")
    .def(init<int,FileIO>())
    .def(init<Atoms,FileIO>())
    .def("printInfo", &Molecule::printInfo)
  ;
};

