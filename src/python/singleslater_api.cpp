#include <boost/python.hpp>
#include <workers.h>
#include <global.h>

using namespace boost::python;
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
  void Wrapper_readInput(FileIO &fileio, Molecule &mol, BasisSet &basis, Controls &controls,
    BasisSet &dfBasis) {
    readInput(&fileio,&mol,&basis,&controls,&dfBasis);
  }

  void FileIO::write(std::string str){
    this->out << str << std::endl; 
  }
  void AOIntegrals::Wrapper_iniAOIntegrals(Molecule &mol, BasisSet &basis,
         FileIO &fileio, Controls &controls) {
     this->iniAOIntegrals( &mol, &basis, &fileio, &controls); 
  }
  
  void Molecule::Wrapper_printInfo(FileIO& fileio, Controls &cont){
    this->printInfo(&fileio, &cont);
  }
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(printSettings_Overload,
  Controls::printSettings,0,1);

BOOST_PYTHON_MODULE(libpythonapi){
  class_<SingleSlater<double>,boost::noncopyable>("SingleSlater_double",init<>())
    .def("iniSingleSlater" , &SingleSlater<double>::Wrapper_iniSingleSlater)
    .def("printInfo"       , &SingleSlater<double>::printInfo              )
    .def("formGuess"       , &SingleSlater<double>::formGuess              )
    .def("formFock"        , &SingleSlater<double>::formFock               )
    .def("computeEnergy"   , &SingleSlater<double>::computeEnergy          )
    .def("computeMultipole", &SingleSlater<double>::computeMultipole       )
    .def("printMultipole"  , &SingleSlater<double>::printMultipole         )
    .def("SCF"             , &SingleSlater<double>::SCF                    )
  ;

  class_<Molecule,boost::noncopyable>("Molecule",init<>())
    .def("printInfo", &Molecule::Wrapper_printInfo)
  ;

  class_<BasisSet,boost::noncopyable>("BasisSet",init<>())
    .def("printInfo", &BasisSet::printInfo)
  ;

  class_<Controls,boost::noncopyable>("Controls",init<>())
    .def("iniControls"  , &Controls::iniControls  )
    .def("printSettings", &Controls::printSettings,
      printSettings_Overload(args("x"),"printSettings docstring")
      )
  ;

  class_<FileIO,boost::noncopyable>("FileIO",init<std::string>())
    .def("write", &FileIO::write)
  ;

  
  class_<AOIntegrals,boost::noncopyable>("AOIntegrals",init<>())
    .def("iniAOIntegrals", &AOIntegrals::Wrapper_iniAOIntegrals)
    .def("printTimings"  , &AOIntegrals::printTimings          )
  ;

  def("readInput", ChronusQ::Wrapper_readInput);
};

