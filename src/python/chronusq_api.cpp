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
    .def("setSCFDenTol"    , &SingleSlater<double>::setSCFDenTol           )
    .def("setSCFEneTol"    , &SingleSlater<double>::setSCFEneTol           )
    .def("setSCFMaxIter"   , &SingleSlater<double>::setSCFMaxIter          )
    .def("setField"        , &SingleSlater<double>::Wrapper_setField       )
    .def("setGuess"        , &SingleSlater<double>::setGuess               )

    .def("dipole"          , &SingleSlater<double>::Wrapper_dipole         )
    .def("quadrupole"      , &SingleSlater<double>::Wrapper_quadrupole     )
    .def("octupole"        , &SingleSlater<double>::Wrapper_octupole       )

    .def("Ref"             , &SingleSlater<double>::Ref                    )
    .def("nTCS"            , &SingleSlater<double>::nTCS                   ) 

    .def_readwrite("isClosedShell", &SingleSlater<double>::isClosedShell   )
    .def_readonly("totalEnergy"   , &SingleSlater<double>::totalEnergy     )
    .def_readonly("nSCFIter"      , &SingleSlater<double>::nSCFIter        )
  ;

  class_<SingleSlater<dcomplex>,boost::noncopyable>("SingleSlater_complex",
    init<>())
    .def("iniSingleSlater" , &SingleSlater<dcomplex>::Wrapper_iniSingleSlater)
    .def("printInfo"       , &SingleSlater<dcomplex>::printInfo              )
    .def("printDensity"    , &SingleSlater<dcomplex>::printDensity           )
    .def("formGuess"       , &SingleSlater<dcomplex>::formGuess              )
    .def("formFock"        , &SingleSlater<dcomplex>::formFock               )
    .def("computeEnergy"   , &SingleSlater<dcomplex>::computeEnergy          )
    .def("computeMultipole", &SingleSlater<dcomplex>::computeMultipole       )
    .def("printMultipole"  , &SingleSlater<dcomplex>::printMultipole         )
    .def("SCF"             , &SingleSlater<dcomplex>::SCF                    )
    .def("communicate"     , &SingleSlater<dcomplex>::communicate            )
    .def("initMeta"        , &SingleSlater<dcomplex>::initMeta               )
    .def("alloc"           , &SingleSlater<dcomplex>::alloc                  )
    .def("genMethString"   , &SingleSlater<dcomplex>::genMethString          )
    .def("setRef"          , &SingleSlater<dcomplex>::setRef                 )
    .def("setNTCS"         , &SingleSlater<dcomplex>::setNTCS                )
    .def("setMaxMultipole" , &SingleSlater<dcomplex>::setMaxMultipole        )
    .def("printLevel"      , &SingleSlater<dcomplex>::printLevel             )
    .def("setPrintLevel"   , &SingleSlater<dcomplex>::setPrintLevel          )
    .def("setSCFDenTol"    , &SingleSlater<dcomplex>::setSCFDenTol           )
    .def("setSCFEneTol"    , &SingleSlater<dcomplex>::setSCFEneTol           )
    .def("setSCFMaxIter"   , &SingleSlater<dcomplex>::setSCFMaxIter          )
    .def("setField"        , &SingleSlater<dcomplex>::Wrapper_setField       )
    .def("setGuess"        , &SingleSlater<dcomplex>::setGuess               )

    .def("dipole"          , &SingleSlater<dcomplex>::Wrapper_dipole         )
    .def("quadrupole"      , &SingleSlater<dcomplex>::Wrapper_quadrupole     )
    .def("octupole"        , &SingleSlater<dcomplex>::Wrapper_octupole       )

    .def("Ref"             , &SingleSlater<dcomplex>::Ref                    )
    .def("nTCS"            , &SingleSlater<dcomplex>::nTCS                   ) 

    .def_readwrite("isClosedShell", &SingleSlater<dcomplex>::isClosedShell   )
    .def_readonly("totalEnergy"   , &SingleSlater<dcomplex>::totalEnergy     )
    .def_readonly("nSCFIter"      , &SingleSlater<dcomplex>::nSCFIter        )
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

  enum_<SingleSlater<double>::GUESS>("Guess")
    .value("SAD"  , SingleSlater<double>::SAD  )
    .value("CORE" , SingleSlater<double>::CORE )
    .value("READ" , SingleSlater<double>::READ )
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
    .def("write"      , &FileIO::write         )
    .def("iniH5Files" , &FileIO::iniH5Files    )
    .def("iniStdGroups", &FileIO::iniStdGroups )
    .def("iniStdSCFFilesDouble", &FileIO::iniStdSCFFilesDouble)
    .def("iniStdSCFFilesComplex", &FileIO::iniStdSCFFilesComplex)
    .def_readwrite("doRestart", &FileIO::doRestart   )
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
    .def("setAlgorithm"   , &AOIntegrals::setAlgorithm          )
    
//  .def_readwrite("allocERI", &AOIntegrals::allocERI           )
//  .def_readwrite("doDF"    , &AOIntegrals::doDF               )
  ;

  enum_<AOIntegrals::INTEGRAL_ALGORITHM>("AOIntegrals_INTEGRAL_ALGORITHM")
    .value("DIRECT"   , AOIntegrals::DIRECT)
    .value("INCORE"   , AOIntegrals::INCORE)
    .value("DENFIT"   , AOIntegrals::DENFIT)
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
    .def("setFieldAmp"  , &RealTime<double>::Wrapper_setFieldAmp  )
    .def("setTOn"       , &RealTime<double>::setTOn       )
    .def("setTOff"      , &RealTime<double>::setTOff      )
    .def("setFreq"      , &RealTime<double>::setFreq      )
    .def("setPhase"     , &RealTime<double>::setPhase     )
    .def("setSigma"     , &RealTime<double>::setSigma     )
    .def("printRT"      , &RealTime<double>::printRT      )
  ;

  class_<RealTime<dcomplex>,boost::noncopyable>("RealTime_complex",init<>())
    .def("communicate"  , &RealTime<dcomplex>::communicate  )
    .def("initMeta"     , &RealTime<dcomplex>::initMeta     )
    .def("alloc"        , &RealTime<dcomplex>::alloc        )
    .def("iniDensity"   , &RealTime<dcomplex>::iniDensity   )
    .def("doPropagation", &RealTime<dcomplex>::doPropagation)
    .def("setMaxSteps"  , &RealTime<dcomplex>::setMaxSteps  ) 
    .def("setStepSize"  , &RealTime<dcomplex>::setStepSize  )
    .def("setOrthoTyp"  , &RealTime<dcomplex>::setOrthoTyp  )
    .def("setInitDen"   , &RealTime<dcomplex>::setInitDen   )
    .def("setSwapMOA"   , &RealTime<dcomplex>::setSwapMOA   )
    .def("setSwapMOB"   , &RealTime<dcomplex>::setSwapMOB   )
    .def("setFormU"     , &RealTime<dcomplex>::setFormU     )
    .def("setEnvelope"  , &RealTime<dcomplex>::setEnvelope  )
    .def("setFieldAmp"  , &RealTime<dcomplex>::Wrapper_setFieldAmp  )
    .def("setTOn"       , &RealTime<dcomplex>::setTOn       )
    .def("setTOff"      , &RealTime<dcomplex>::setTOff      )
    .def("setFreq"      , &RealTime<dcomplex>::setFreq      )
    .def("setPhase"     , &RealTime<dcomplex>::setPhase     )
    .def("setSigma"     , &RealTime<dcomplex>::setSigma     )
    .def("printRT"      , &RealTime<dcomplex>::printRT      )
  ;

  enum_<RealTime<double>::ORTHO>("RealTime_ORTHO"   )
    .value("Lowdin"   , RealTime<double>::Lowdin    )
    .value("Cholesky" , RealTime<double>::Cholesky  )
    .value("Canonical", RealTime<double>::Canonical )
  ;

  enum_<RealTime<dcomplex>::ORTHO>("RealTime_ORTHO"   )
    .value("Lowdin"   , RealTime<dcomplex>::Lowdin    )
    .value("Cholesky" , RealTime<dcomplex>::Cholesky  )
    .value("Canonical", RealTime<dcomplex>::Canonical )
  ;

  enum_<RealTime<double>::FORM_U>("RealTime_FORM_U"    )
    .value("EigenDecomp", RealTime<double>::EigenDecomp)
    .value("Taylor"     , RealTime<double>::Taylor     )
  ;

  enum_<RealTime<dcomplex>::FORM_U>("RealTime_FORM_U"    )
    .value("EigenDecomp", RealTime<dcomplex>::EigenDecomp)
    .value("Taylor"     , RealTime<dcomplex>::Taylor     )
  ;

  enum_<RealTime<double>::ENVELOPE>("RealTime_ENVELOPE")
    .value("Constant", RealTime<double>::Constant      )
    .value("LinRamp" , RealTime<double>::LinRamp       )
    .value("Gaussian", RealTime<double>::Gaussian      )
    .value("Step"    , RealTime<double>::Step          )
    .value("SinSq"   , RealTime<double>::SinSq         )
  ;

  enum_<RealTime<dcomplex>::ENVELOPE>("RealTime_ENVELOPE")
    .value("Constant", RealTime<dcomplex>::Constant      )
    .value("LinRamp" , RealTime<dcomplex>::LinRamp       )
    .value("Gaussian", RealTime<dcomplex>::Gaussian      )
    .value("Step"    , RealTime<dcomplex>::Step          )
    .value("SinSq"   , RealTime<dcomplex>::SinSq         )
  ;

  class_<MOIntegrals<double>,boost::noncopyable>("MOIntegrals_double",init<>())
    .def("communicate" , &MOIntegrals<double>::communicate)
    .def("initMeta"    , &MOIntegrals<double>::initMeta   )
  ;
  class_<MOIntegrals<dcomplex>,boost::noncopyable>("MOIntegrals_complex",init<>())
    .def("communicate" , &MOIntegrals<dcomplex>::communicate)
    .def("initMeta"    , &MOIntegrals<dcomplex>::initMeta   )
  ;

  class_<SDResponse<double>,boost::noncopyable>("SDResponse_double",init<>())
    .def("communicate"       , &SDResponse<double>::communicate               )
    .def("initMeta"          , &SDResponse<double>::initMeta                  )
    .def("alloc"             , &SDResponse<double>::alloc                     )
    .def("setNSek"           , &SDResponse<double>::setNSek                   )
    .def("setMeth"           , &SDResponse<double>::setMeth                   )
    .def("initMeth"          , &SDResponse<double>::initMeth                  )
    .def("IterativeRPA"      , &SDResponse<double>::IterativeRPA              )
    .def("excitationEnergies", &SDResponse<double>::Wrapper_excitationEnergies)
    .def("oscStrengths"      , &SDResponse<double>::Wrapper_oscStrengths      )
    
    .def_readonly("nIter"    ,&SDResponse<double>::nQNIter                   )
    
    
  ;

  class_<SDResponse<dcomplex>,boost::noncopyable>("SDResponse_complex",init<>())
    .def("communicate"       , &SDResponse<dcomplex>::communicate              )
    .def("initMeta"          , &SDResponse<dcomplex>::initMeta                 )
    .def("alloc"             , &SDResponse<dcomplex>::alloc                    )
    .def("setNSek"           , &SDResponse<dcomplex>::setNSek                  )
    .def("setMeth"           , &SDResponse<dcomplex>::setMeth                  )
    .def("initMeth"          , &SDResponse<dcomplex>::initMeth                 )
    .def("IterativeRPA"      , &SDResponse<dcomplex>::IterativeRPA             )
    .def("excitationEnergies",&SDResponse<dcomplex>::Wrapper_excitationEnergies)
    .def("oscStrengths"      ,&SDResponse<dcomplex>::Wrapper_oscStrengths      )
    
    .def_readonly("nIter"    ,&SDResponse<dcomplex>::nQNIter                   )
  ;

  enum_<SDResponse<double>::METHOD>("SDResponse_METHOD")
    .value("INVALID", SDResponse<double>::__invalid    )
    .value("CIS"    , SDResponse<double>::CIS          )
    .value("RPA"    , SDResponse<double>::RPA          )
    .value("PPRPA"  , SDResponse<double>::PPRPA        )
    .value("PPATDA" , SDResponse<double>::PPATDA       )
    .value("PPCTDA" , SDResponse<double>::PPCTDA       )
    .value("STAB"   , SDResponse<double>::STAB         )
  ;


  def("readInput",       ChronusQ::Wrapper_readInput    );
  def("HashAtom",        ChronusQ::HashAtom             );
  def("getAtomicNumber", ChronusQ::getAtomicNumber      );
  def("CErr",            ChronusQ::Wrapper_CErr_Default );
  def("CErrMsg",         ChronusQ::Wrapper_CErr_Message );
  def("initCQ",          ChronusQ::initCQ               );
  def("finalizeCQ",      ChronusQ::finalizeCQ           );
  def("CQSetNumThreads", ChronusQ::CQSetNumThreads      );
};

