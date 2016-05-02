/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
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
    .def("computeProperties",&SingleSlater<double>::computeProperties      )
    .def("printMultipole"  , &SingleSlater<double>::printMultipole         )
    .def("printProperties" , &SingleSlater<double>::printProperties        )
    .def("SCF"             , &SingleSlater<double>::SCF2                   )
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
    .def("setCorrKernel"   , &SingleSlater<double>::setCorrKernel          )
    .def("setExchKernel"   , &SingleSlater<double>::setExchKernel          )
    .def("setDFTKernel"    , &SingleSlater<double>::setDFTKernel           )
    .def("setDFTWeightScheme", &SingleSlater<double>::setDFTWeightScheme)
    .def("setDFTGrid"      , &SingleSlater<double>::setDFTGrid             ) 
    .def("setDFTNGridPts"  , &SingleSlater<double>::setDFTNGridPts         )
    .def("setDFTNRad"  , &SingleSlater<double>::setDFTNRad         )
    .def("setDFTNAng"  , &SingleSlater<double>::setDFTNAng         )
    .def("setDFTScreenTol" , &SingleSlater<double>::setDFTScreenTol        )

    .def("checkDFTType", &SingleSlater<double>::checkDFTType)
  
    .def("turnOffDFTScreening", &SingleSlater<double>::turnOffDFTScreening ) 

    .def("dipole"          , &SingleSlater<double>::Wrapper_dipole         )
    .def("quadrupole"      , &SingleSlater<double>::Wrapper_quadrupole     )
    .def("octupole"        , &SingleSlater<double>::Wrapper_octupole       )

    .def("Ref"             , &SingleSlater<double>::Ref                    )
    .def("nTCS"            , &SingleSlater<double>::nTCS                   ) 

    .def_readwrite("isClosedShell", &SingleSlater<double>::isClosedShell   )
    .def_readonly("totalEnergy"   , &SingleSlater<double>::totalEnergy     )
    .def_readonly("nSCFIter"      , &SingleSlater<double>::nSCFIter        )
    .def_readwrite("isDFT"        , &SingleSlater<double>::isDFT           )
    .def_readwrite("isGGA"        , &SingleSlater<double>::isGGA           )
    .def_readwrite("isHF"         , &SingleSlater<double>::isHF            )
    .def_readwrite("doDIIS"       , &SingleSlater<double>::doDIIS          )
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
    .def("computeProperties",&SingleSlater<dcomplex>::computeProperties      )
    .def("printMultipole"  , &SingleSlater<dcomplex>::printMultipole         )
    .def("printProperties" , &SingleSlater<dcomplex>::printProperties        )
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
    .def("setCorrKernel"   , &SingleSlater<dcomplex>::setCorrKernel          )
    .def("setExchKernel"   , &SingleSlater<dcomplex>::setExchKernel          )
    .def("setDFTKernel"    , &SingleSlater<dcomplex>::setDFTKernel           )
    .def("setDFTWeightScheme", &SingleSlater<dcomplex>::setDFTWeightScheme)
    .def("setDFTGrid"      , &SingleSlater<dcomplex>::setDFTGrid             ) 
    .def("setDFTNGridPts"  , &SingleSlater<dcomplex>::setDFTNGridPts         )
    .def("setDFTNRad"  , &SingleSlater<dcomplex>::setDFTNRad         )
    .def("setDFTNAng"  , &SingleSlater<dcomplex>::setDFTNAng         )
    .def("setDFTScreenTol" , &SingleSlater<dcomplex>::setDFTScreenTol        )

    .def("checkDFTType", &SingleSlater<dcomplex>::checkDFTType)
  
    .def("turnOffDFTScreening", &SingleSlater<dcomplex>::turnOffDFTScreening ) 

    .def("dipole"          , &SingleSlater<dcomplex>::Wrapper_dipole         )
    .def("quadrupole"      , &SingleSlater<dcomplex>::Wrapper_quadrupole     )
    .def("octupole"        , &SingleSlater<dcomplex>::Wrapper_octupole       )

    .def("Ref"             , &SingleSlater<dcomplex>::Ref                    )
    .def("nTCS"            , &SingleSlater<dcomplex>::nTCS                   ) 

    .def_readwrite("isClosedShell", &SingleSlater<dcomplex>::isClosedShell   )
    .def_readonly("totalEnergy"   , &SingleSlater<dcomplex>::totalEnergy     )
    .def_readonly("nSCFIter"      , &SingleSlater<dcomplex>::nSCFIter        )
    .def_readwrite("isDFT"        , &SingleSlater<dcomplex>::isDFT           )
    .def_readwrite("isGGA"        , &SingleSlater<dcomplex>::isGGA           )
    .def_readwrite("isHF"         , &SingleSlater<dcomplex>::isHF            )
    .def_readwrite("doDIIS"       , &SingleSlater<dcomplex>::doDIIS          )
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

  enum_<SingleSlater<double>::EXCH>("EXCH")
    .value("NOEXCH"  , SingleSlater<double>::NOEXCH) 
    .value("EXACT"   , SingleSlater<double>::EXACT )
    .value("SLATER"  , SingleSlater<double>::SLATER)
    .value("B88"     , SingleSlater<double>::B88   )
  ;
  enum_<SingleSlater<double>::CORR>("CORR")
    .value("NOCORR", SingleSlater<double>::NOCORR) 
    .value("VWN3"  , SingleSlater<double>::VWN3  )
    .value("VWN5"  , SingleSlater<double>::VWN5  )
    .value("LYP"   , SingleSlater<double>::LYP   )
  ;
  enum_<SingleSlater<double>::DFT>("DFT")
    .value("NODFT"      , SingleSlater<double>::NODFT      )
    .value("USERDEFINED", SingleSlater<double>::USERDEFINED)
    .value("LSDA"       , SingleSlater<double>::LSDA       )
  ;

  enum_<SingleSlater<double>::DFT_GRID>("DFT_GRID")
    .value("EULERMACL", SingleSlater<double>::EULERMACL)
    .value("GAUSSCHEB", SingleSlater<double>::GAUSSCHEB)
  ;
  enum_<SingleSlater<double>::DFT_WEIGHT_SCHEME>("DFT_WEIGHT_SCHEME")
    .value("BECKE",  SingleSlater<double>::BECKE )
    .value("FRISCH", SingleSlater<double>::FRISCH)
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
    .def("cart",          &Molecule::Wrapper_cart     )
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
    .def("forceCart"    ,  &BasisSet::forceCart              )
  ;

  class_<Controls,boost::noncopyable>("Controls",init<>())
    .def("iniControls"  , &Controls::iniControls  )
    .def("printSettings", &Controls::printSettings,
      printSettings_Overload(args("x"),"printSettings docstring")
      )
  ;

  class_<FileIO,boost::noncopyable>("FileIO",init<std::string>())
    .def(init<std::string,std::string>())
    .def(init<std::string,std::string,std::string>())
    .def("write"      , &FileIO::write         )
    .def("iniH5Files" , &FileIO::iniH5Files    )
    .def("iniStdGroups", &FileIO::iniStdGroups )
    .def("iniStdSCFFilesDouble", &FileIO::iniStdSCFFilesDouble)
    .def("iniStdSCFFilesComplex", &FileIO::iniStdSCFFilesComplex)
    .def_readwrite("doRestart", &FileIO::doRestart   )
  ;

  
  class_<AOIntegrals,boost::noncopyable>("AOIntegrals",init<>())
    .def("printTimings"   , &AOIntegrals::printTimings          )
    .def("communicate"    , &AOIntegrals::communicate           )
    .def("initMeta"       , &AOIntegrals::initMeta              )
    .def("alloc"          , &AOIntegrals::alloc                 )
    .def("setMaxMultipole", &AOIntegrals::setMaxMultipole       )
    .def("setAlgorithm"   , &AOIntegrals::setAlgorithm          )
    .def("setPrintLevel"   , &AOIntegrals::setPrintLevel        )
    
    .def_readonly("integralAlgorithm", &AOIntegrals::integralAlgorithm)
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
    .def("setEllPol"    , &RealTime<double>::setEllPol    )
    .def("setFieldAmp"  , &RealTime<double>::Wrapper_setFieldAmp  )
    .def("setTOn"       , &RealTime<double>::setTOn       )
    .def("setTOff"      , &RealTime<double>::setTOff      )
    .def("setFreq"      , &RealTime<double>::setFreq      )
    .def("setPhase"     , &RealTime<double>::setPhase     )
    .def("setSigma"     , &RealTime<double>::setSigma     )
    .def("setPrintLevel", &RealTime<double>::setPrintLevel)
    .def("printRT"      , &RealTime<double>::printRT      )
    //.def("recs"         , &RealTime<double>::Wrapper_recs )
    .def("lastDipole"   , &RealTime<double>::lastDipole   )
    .def("lastEnergy"   , &RealTime<double>::lastEnergy   )
    .def("getTimeStep"  , &RealTime<double>::getTimeStep  )
    .def("doNotTarCSV"  , &RealTime<double>::doNotTarCSV  )
    .def("setIRstrt"       , &RealTime<double>::setIRstrt       )
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
    .def("setEllPol"    , &RealTime<dcomplex>::setEllPol    )
    .def("setFieldAmp"  , &RealTime<dcomplex>::Wrapper_setFieldAmp  )
    .def("setTOn"       , &RealTime<dcomplex>::setTOn       )
    .def("setTOff"      , &RealTime<dcomplex>::setTOff      )
    .def("setFreq"      , &RealTime<dcomplex>::setFreq      )
    .def("setPhase"     , &RealTime<dcomplex>::setPhase     )
    .def("setSigma"     , &RealTime<dcomplex>::setSigma     )
    .def("setPrintLevel", &RealTime<dcomplex>::setPrintLevel)
    .def("printRT"      , &RealTime<dcomplex>::printRT      )
//    .def("recs"         , &RealTime<dcomplex>::Wrapper_recs )
    .def("lastDipole"   , &RealTime<dcomplex>::lastDipole   )
    .def("lastEnergy"   , &RealTime<dcomplex>::lastEnergy   )
    .def("getTimeStep"  , &RealTime<dcomplex>::getTimeStep  )
    .def("doNotTarCSV"  , &RealTime<dcomplex>::doNotTarCSV  )
    .def("setIRstrt"       , &RealTime<dcomplex>::setIRstrt       )
  ;

/*
  class_<RealTime<double>::Wrapper_PropInfo,boost::noncopyable>("RTD_PropInfo",init<>())
    .def_readwrite("timeStep"    , &RealTime<double>::Wrapper_PropInfo::timeStep)
    .def_readwrite("energy"      , &RealTime<double>::Wrapper_PropInfo::energy)
    .def("dipole"      , &RealTime<double>::Wrapper_PropInfo::dipole)
    .def_readwrite("mullPop"     , &RealTime<double>::Wrapper_PropInfo::mullPop)
    .def_readwrite("orbitalOccA" , &RealTime<double>::Wrapper_PropInfo::orbitalOccA)
    .def_readwrite("orbitalOccB" , &RealTime<double>::Wrapper_PropInfo::orbitalOccB)
  ;

  class_<RealTime<dcomplex>::Wrapper_PropInfo,boost::noncopyable>("RTC_PropInfo",init<>())
    .def_readwrite("timeStep"    , &RealTime<dcomplex>::Wrapper_PropInfo::timeStep)
    .def_readwrite("energy"      , &RealTime<dcomplex>::Wrapper_PropInfo::energy)
    .def_readwrite("dipole"      , &RealTime<dcomplex>::Wrapper_PropInfo::dipole)
    .def_readwrite("mullPop"     , &RealTime<dcomplex>::Wrapper_PropInfo::mullPop)
    .def_readwrite("orbitalOccA" , &RealTime<dcomplex>::Wrapper_PropInfo::orbitalOccA)
    .def_readwrite("orbitalOccB" , &RealTime<dcomplex>::Wrapper_PropInfo::orbitalOccB)
  ;
*/

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
    .value("Elliptic", RealTime<double>::Elliptic      )
  ;

  enum_<RealTime<dcomplex>::ENVELOPE>("RealTime_ENVELOPE")
    .value("Constant", RealTime<dcomplex>::Constant      )
    .value("LinRamp" , RealTime<dcomplex>::LinRamp       )
    .value("Gaussian", RealTime<dcomplex>::Gaussian      )
    .value("Step"    , RealTime<dcomplex>::Step          )
    .value("SinSq"   , RealTime<dcomplex>::SinSq         )
    .value("Elliptic", RealTime<dcomplex>::Elliptic      )
  ;

  enum_<RealTime<double>::ELL_POL>("RealTime_ELL_POL")
    .value("LXY", RealTime<double>::LXY)
    .value("LXZ", RealTime<double>::LXZ)
    .value("LYZ", RealTime<double>::LYZ)
    .value("RXY", RealTime<double>::RXY)
    .value("RXZ", RealTime<double>::RXZ)
    .value("RYZ", RealTime<double>::RYZ)
  ;

  enum_<RealTime<dcomplex>::ELL_POL>("RealTime_ELL_POL")
    .value("LXY", RealTime<dcomplex>::LXY)
    .value("LXZ", RealTime<dcomplex>::LXZ)
    .value("LYZ", RealTime<dcomplex>::LYZ)
    .value("RXY", RealTime<dcomplex>::RXY)
    .value("RXZ", RealTime<dcomplex>::RXZ)
    .value("RYZ", RealTime<dcomplex>::RYZ)
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
    .def("IterativeRESP"      , &SDResponse<double>::IterativeRESP              )
    .def("excitationEnergies", &SDResponse<double>::Wrapper_excitationEnergies)
    .def("oscStrengths"      , &SDResponse<double>::Wrapper_oscStrengths      )
    .def("setPPRPA"          , &SDResponse<double>::setPPRPA                  )
    
    .def_readonly("nIter"    ,&SDResponse<double>::nQNIter                   )
    
    
  ;

  class_<SDResponse<dcomplex>,boost::noncopyable>("SDResponse_complex",init<>())
    .def("communicate"       , &SDResponse<dcomplex>::communicate              )
    .def("initMeta"          , &SDResponse<dcomplex>::initMeta                 )
    .def("alloc"             , &SDResponse<dcomplex>::alloc                    )
    .def("setNSek"           , &SDResponse<dcomplex>::setNSek                  )
    .def("setMeth"           , &SDResponse<dcomplex>::setMeth                  )
    .def("initMeth"          , &SDResponse<dcomplex>::initMeth                 )
    .def("IterativeRESP"      , &SDResponse<dcomplex>::IterativeRESP             )
    .def("excitationEnergies",&SDResponse<dcomplex>::Wrapper_excitationEnergies)
    .def("oscStrengths"      ,&SDResponse<dcomplex>::Wrapper_oscStrengths      )
    .def("setPPRPA"          , &SDResponse<dcomplex>::setPPRPA                  )
    
    .def_readonly("nIter"    ,&SDResponse<dcomplex>::nQNIter                   )
  ;

  class_<NumericalDifferentiation<double>,boost::noncopyable>(
    "NumericalDifferentiationDouble",init<>())
    .def("differentiate"  ,&NumericalDifferentiation<double>::differentiate)
    .def("setSingleSlater",&NumericalDifferentiation<double>::setSingleSlater)
    .def("generateESObjs"  ,&NumericalDifferentiation<double>::generateESObjs)
    .def("setRespNRoots"  ,&NumericalDifferentiation<double>::setRespNRoots)
    .def("setRespRoot"  ,&NumericalDifferentiation<double>::setRespRoot)
    .def("setRespType"  ,&NumericalDifferentiation<double>::setRespType)
    .def("GSEnergy" , &NumericalDifferentiation<double>::Wrapper_GSEnergy )
    .def("GSGrad"   , &NumericalDifferentiation<double>::Wrapper_GSGrad   ) 
    .def("ESEnergy" , &NumericalDifferentiation<double>::Wrapper_ESEnergy )
    .def("ESGrad"   , &NumericalDifferentiation<double>::Wrapper_ESGrad   )
    .def("ESGSNACME", &NumericalDifferentiation<double>::Wrapper_ESGSNACME)
    .def("ESESNACME", &NumericalDifferentiation<double>::Wrapper_ESESNACME)

    .def_readwrite("computeGSGradient", 
                   &NumericalDifferentiation<double>::computeGSGradient)
    .def_readwrite("computeESGradient", 
                   &NumericalDifferentiation<double>::computeESGradient)
    .def_readwrite("computeES2GSNACME", 
                   &NumericalDifferentiation<double>::computeES2GSNACME)
    .def_readwrite("computeES2ESNACME", 
                   &NumericalDifferentiation<double>::computeES2ESNACME)
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

  enum_<RESPONSE_TYPE>("RESPONSE_TYPE")
    .value("NOMETHOD", RESPONSE_TYPE::NOMETHOD    )
    .value("CIS"    , RESPONSE_TYPE::CIS          )
    .value("RPA"    , RESPONSE_TYPE::RPA          )
    .value("PPRPA"  , RESPONSE_TYPE::PPRPA        )
    .value("PPTDA" , RESPONSE_TYPE::PPTDA       )
    .value("STAB"   , RESPONSE_TYPE::STAB         )
  ;

  

  def("readInput",       ChronusQ::Wrapper_readInput    );
  def("HashAtom",        ChronusQ::HashAtom             );
  def("HashZ",           ChronusQ::HashZ                );
  def("getAtomicNumber", ChronusQ::getAtomicNumber      );
  def("CErr",            ChronusQ::Wrapper_CErr_Default );
  def("CErrMsg",         ChronusQ::Wrapper_CErr_Message );
  def("initCQ",          ChronusQ::Wrapper_initCQ       );
  def("finalizeCQ",      ChronusQ::finalizeCQ           );
  def("CQSetNumThreads", ChronusQ::CQSetNumThreads      );
  def("getRank",         ChronusQ::getRank              );
  def("getSize",         ChronusQ::getSize              );
  def("mpiBarrier",      ChronusQ::mpiBarrier           );

  
};

