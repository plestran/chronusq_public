#include <pythonapi.h>
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(printSettings_Overload,
  Controls::printSettings,0,1);

BOOST_PYTHON_MODULE(libpythonapi){
  class_<SingleSlater<double>,boost::noncopyable>("SingleSlater_double",
    init<>())
    .def("iniSingleSlater" , &SingleSlater<double>::Wrapper_iniSingleSlater)
    .def("printInfo"       , &SingleSlater<double>::printInfo              )
    .def("formGuess"       , &SingleSlater<double>::formGuess              )
    .def("formFock"        , &SingleSlater<double>::formFock               )
    .def("computeEnergy"   , &SingleSlater<double>::computeEnergy          )
    .def("computeMultipole", &SingleSlater<double>::computeMultipole       )
    .def("printMultipole"  , &SingleSlater<double>::printMultipole         )
    .def("SCF"             , &SingleSlater<double>::SCF                    )
    .def("communicate"     , &SingleSlater<double>::communicate            )
    .def("initMeta"        , &SingleSlater<double>::initMeta               )
    .def("alloc"           , &SingleSlater<double>::alloc                  )
  ;

  class_<Molecule,boost::noncopyable>("Molecule",init<>())
    .def("printInfo", &Molecule::Wrapper_printInfo)
    .def("alloc"    , &Molecule::Wrapper_alloc    )
    .def("setCharge", &Molecule::setCharge        )
    .def("setMultip", &Molecule::setMultip        )
    .def("setNAtoms", &Molecule::setNAtoms        )
    .def("setIndex",  &Molecule::setIndex         )
    .def("setCart",   &Molecule::setCart          )
    .def("computeRij",&Molecule::computeRij       )
    .def("toCOM",     &Molecule::toCOM            )
    .def("computeI",  &Molecule::computeI         )
    .def("setNTotalE",&Molecule::setNTotalE       )
    .def("convBohr",  &Molecule::convBohr         )
    .def("computeNucRep", &Molecule::computeNucRep)
    .def("energyNuclei",  &Molecule::energyNuclei )
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

  def("readInput",       ChronusQ::Wrapper_readInput);
  def("HashAtom",        ChronusQ::HashAtom         );
  def("getAtomicNumber", ChronusQ::getAtomicNumber  );
};

