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
    .def("setMaxMultipole" , &SingleSlater<double>::setMaxMultipole        )
    .def("printLevel"      , &SingleSlater<double>::printLevel             )
    .def("setPrintLevel"   , &SingleSlater<double>::setPrintLevel          )

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
    .def("printInfo",     &Molecule::Wrapper_printInfo)
    .def("alloc"    ,     &Molecule::Wrapper_alloc    )
    .def("setCharge",     &Molecule::setCharge        )
    .def("setMultip",     &Molecule::setMultip        )
    .def("setNAtoms",     &Molecule::setNAtoms        )
    .def("setIndex",      &Molecule::setIndex         )
    .def("setCart",       &Molecule::setCart          )
    .def("setPrintLevel", &Molecule::setPrintLevel    )
    .def("setNTotalE",    &Molecule::setNTotalE       )
    .def("printLevel",    &Molecule::printLevel       )
    .def("energyNuclei",  &Molecule::energyNuclei     )
    .def("multip"      ,  &Molecule::multip           )
    .def("toCOM",         &Molecule::toCOM            )
    .def("convBohr",      &Molecule::convBohr         )
    .def("computeNucRep", &Molecule::computeNucRep    )
    .def("computeI",      &Molecule::computeI         )
    .def("computeRij",    &Molecule::computeRij       )
  ;

  class_<BasisSet,boost::noncopyable>("BasisSet",init<>())
    .def("communicate",    &BasisSet::communicate            )
    .def("printInfo",      &BasisSet::printInfo              )
    .def("findBasisFile",  &BasisSet::findBasisFile          )
    .def("parseGlobal",    &BasisSet::parseGlobal            )
    .def("constructLocal", &BasisSet::Wrapper_constructLocal ) 
    .def("makeMaps"     ,  &BasisSet::Wrapper_makeMaps       )
    .def("renormShells" ,  &BasisSet::renormShells           )
    .def("printLevel"   ,  &BasisSet::printLevel             )
    .def("setPrintLevel",  &BasisSet::setPrintLevel          )
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
    .def("iniAOIntegrals" , &AOIntegrals::Wrapper_iniAOIntegrals)
    .def("printTimings"   , &AOIntegrals::printTimings          )
    .def("communicate"    , &AOIntegrals::communicate           )
    .def("initMeta"       , &AOIntegrals::initMeta              )
    .def("alloc"          , &AOIntegrals::alloc                 )
    .def("nTCS"           , &AOIntegrals::nTCS                  )
    .def("setNTCS"        , &AOIntegrals::setNTCS               )
    .def("setMaxMultipole", &AOIntegrals::setMaxMultipole       )
    
    .def_readwrite("allocERI", &AOIntegrals::allocERI           )
    .def_readwrite("doDF"    , &AOIntegrals::doDF               )
  ;

  class_<RealTime<double>,boost::noncopyable>("RealTime_double",init<>())
    .def("communicate"  , &RealTime<double>::communicate  )
    .def("initMeta"     , &RealTime<double>::initMeta     )
    .def("alloc"        , &RealTime<double>::alloc        )
    .def("iniDensity"   , &RealTime<double>::iniDensity   )
    .def("doPropagation", &RealTime<double>::doPropagation)
    .def("setMaxSteps"  , &RealTime<double>::setMaxSteps  ) 
    .def("setStepSize"  , &RealTime<double>::setStepSize  )
    .def("setOrthoTyp"  , &RealTime<double>::setOrthoTyp  )
    .def("setInitDen"   , &RealTime<double>::setInitDen   )
    .def("setSwapMOA"   , &RealTime<double>::setSwapMOA   )
    .def("setSwapMOB"   , &RealTime<double>::setSwapMOB   )
    .def("setFormU"     , &RealTime<double>::setFormU     )
    .def("setEnvelope"  , &RealTime<double>::setEnvelope  )
    .def("setFieldAmp"  , &RealTime<double>::setFieldAmp  )
    .def("setTOn"       , &RealTime<double>::setTOn       )
    .def("setTOff"      , &RealTime<double>::setTOff      )
    .def("setFreq"      , &RealTime<double>::setFreq      )
    .def("setPhase"     , &RealTime<double>::setPhase     )
    .def("setSigma"     , &RealTime<double>::setSigma     )
  ;

  enum_<RealTime<double>::ORTHO>("RealTime_ORTHO"   )
    .value("Lowdin"   , RealTime<double>::Lowdin    )
    .value("Cholesky" , RealTime<double>::Cholesky  )
    .value("Canonical", RealTime<double>::Canonical )
  ;

  enum_<RealTime<double>::FORM_U>("RealTime_FORM_U"    )
    .value("EigenDecomp", RealTime<double>::EigenDecomp)
    .value("Taylor"     , RealTime<double>::Taylor     )
  ;

  enum_<RealTime<double>::ENVELOPE>("RealTime_ENVELOPE")
    .value("Constant", RealTime<double>::Constant      )
    .value("LinRamp" , RealTime<double>::LinRamp       )
    .value("Gaussian", RealTime<double>::Gaussian      )
    .value("Step"    , RealTime<double>::Step          )
    .value("SinSq"   , RealTime<double>::SinSq         )
  ;

  def("readInput",       ChronusQ::Wrapper_readInput);
  def("HashAtom",        ChronusQ::HashAtom         );
  def("getAtomicNumber", ChronusQ::getAtomicNumber  );
};

