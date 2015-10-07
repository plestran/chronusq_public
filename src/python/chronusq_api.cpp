#include <pythonapi.h>
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(printSettings_Overload,
  Controls::printSettings,0,1);

BOOST_PYTHON_MODULE(libpythonapi){
  class_<SingleSlater<double>,boost::noncopyable>("SingleSlater_double",
    init<>())
    .def("iniSingleSlater" , &SingleSlater<double>::Wrapper_iniSingleSlater)
    .def("printInfo"       , &SingleSlater<double>::printInfo              )
    .def("printDensity"    , &SingleSlater<double>::printDensity           )
    .def("formGuess"       , &SingleSlater<double>::formGuess              )
    .def("formFock"        , &SingleSlater<double>::formFock               )
    .def("computeEnergy"   , &SingleSlater<double>::computeEnergy          )
    .def("computeMultipole", &SingleSlater<double>::computeMultipole       )
    .def("printMultipole"  , &SingleSlater<double>::printMultipole         )
    .def("SCF"             , &SingleSlater<double>::SCF                    )
    .def("communicate"     , &SingleSlater<double>::communicate            )
    .def("initMeta"        , &SingleSlater<double>::initMeta               )
    .def("alloc"           , &SingleSlater<double>::alloc                  )
    .def("genMethString"   , &SingleSlater<double>::genMethString          )
    .def("setRef"          , &SingleSlater<double>::setRef                 )
    .def("setNTCS"         , &SingleSlater<double>::setNTCS                )

    .def("Ref"             , &SingleSlater<double>::Ref                    )
    .def("nTCS"            , &SingleSlater<double>::nTCS                   ) 

    .def_readwrite("isClosedShell", &SingleSlater<double>::isClosedShell   )
  ;

  enum_<SingleSlater<double>::REFERENCE>("Reference")
    .value("_INVALID"      , SingleSlater<double>::_INVALID               )
    .value("RHF"           , SingleSlater<double>::RHF                    )
    .value("UHF"           , SingleSlater<double>::UHF                    )
    .value("CUHF"          , SingleSlater<double>::CUHF                   )
    .value("TCS"           , SingleSlater<double>::TCS                    )
    .value("RKS"           , SingleSlater<double>::RKS                    )
    .value("UKS"           , SingleSlater<double>::UKS                    )
    .value("CUKS"          , SingleSlater<double>::CUKS                   )
    .value("GKS"           , SingleSlater<double>::GKS                    )
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
    .def("multip"      ,  &Molecule::multip       )
  ;

  class_<BasisSet,boost::noncopyable>("BasisSet",init<>())
    .def("communicate",    &BasisSet::communicate            )
    .def("printInfo",      &BasisSet::printInfo              )
    .def("findBasisFile",  &BasisSet::findBasisFile          )
    .def("parseGlobal",    &BasisSet::parseGlobal            )
    .def("constructLocal", &BasisSet::Wrapper_constructLocal ) 
    .def("makeMaps"     ,  &BasisSet::Wrapper_makeMaps       )
    .def("renormShells" ,  &BasisSet::renormShells           )
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

