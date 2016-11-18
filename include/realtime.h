/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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
#ifndef INCLUDED_REALTIME
#define INCLUDED_REALTIME
#include <global.h>
#include <cerr.h>
#include <singleslater.h>

namespace ChronusQ {

enum PropagatorFormation {
  EigenDecomp,
  Taylor,
  Chebyshev
};

enum Envelope {
  Constant, 
  LinRamp,
  Gaussian,
  Step,
  SinSq,
  Elliptic
};

enum EllipticalPolarization {
  LXY,
  LXZ,
  LYZ,
  RXY,
  RXZ,
  RYZ
};

enum PropagationStep {
  ForwardEuler,
  ModifiedMidpoint,
  ExplicitMagnus2,
  ExplicitMagnus3,
  ExplicitMagnus4
};


enum IntegrationScheme {
  MMUT,
  ExpMagnus2,
  ExpMagnus3
};

struct PropInfo {
  double timeStep;
  double energy;
  std::array<double,4> dipole;
  std::array<double,4> appliedfield;
  std::vector<double> mullPop;
  std::vector<double> lowPop;
  std::vector<double> orbitalOccA;
  std::vector<double> orbitalOccB;
};

template <typename T>
class RealTime {
  typedef Eigen::Matrix<T,Dynamic,Dynamic> TMat; 
  typedef Eigen::Map<TMat> TMap;


  // Other classes in the ChronusQ namespace
  FileIO *        fileio_;
  SingleSlater<T> *  groundState_;
  CQMemManager *  memManager_;

  // SS object to propagate
  std::unique_ptr<SingleSlater<dcomplex>> ssPropagator_;


  dcomplex * NBTSqScratch_;
  dcomplex * NBTSqScratch2_;
  dcomplex * NBTSqScratch3_;
  dcomplex * UTransScalar_;
  dcomplex * UTransMz_;
  dcomplex * UTransMy_;
  dcomplex * UTransMx_;

  std::unique_ptr<ComplexMap> POScalarSav_;
  std::unique_ptr<ComplexMap> POMzSav_;
  std::unique_ptr<ComplexMap> POMySav_;
  std::unique_ptr<ComplexMap> POMxSav_;
  std::vector<ComplexMap *> POSav_;

  std::unique_ptr<ComplexMap> FOScalarSav_;
  std::unique_ptr<ComplexMap> FOMzSav_;
  std::unique_ptr<ComplexMap> FOMySav_;
  std::unique_ptr<ComplexMap> FOMxSav_;
  std::vector<ComplexMap *> FOSav_;

  // RT Simulation Control
  int maxSteps_;
  int nSkip_;
  int iRstrt_;

  // Polynomial Expansion of the propagator parameters
  size_t nPolyExpMax_;
  double polyEps_;

  double tOn_;
  double tOff_;
  double stepSize_;
  double deltaT_;
  double freq_;
  double phase_;
  double sigma_;

//  long double currentTime_;

  IntegrationScheme      iScheme_;
  PropagationStep        iRstScheme_;
  PropagatorFormation    iMethFormU_;
  Envelope               iEnvlp_;
  EllipticalPolarization iEllPol_;

  std::array<double,3> staticEDAmp_;
  std::array<double,3> EDField_;


  // Output
  int printLevel_;
  
  std::vector<std::ofstream*> csvs_;
  std::map<std::ofstream*,std::string> csvFiles_;
  bool tarCSVs;
  
public:
  
  std::vector<PropInfo> propInfo;

  RealTime() : fileio_(NULL), groundState_(NULL), memManager_(NULL),
    ssPropagator_(nullptr),
    deltaT_(0.0), 
//  currentTime_(0.0), 
    freq_(0.0), phase_(0.0), sigma_(0.0),
    nSkip_(0),
    tOn_(0.0), tOff_(1.0e4), maxSteps_(10), stepSize_(0.05), iRstrt_(50),
    iScheme_(MMUT), iRstScheme_(ExplicitMagnus2),
    iMethFormU_(EigenDecomp), iEnvlp_(Constant), iEllPol_(LXZ),
    nPolyExpMax_(100), polyEps_(1e-15),
    EDField_({0.0,0.0,0.0}),
    tarCSVs(true){ };

  ~RealTime() {
    for(auto CSV : this->csvs_) {
      CSV->close(); delete CSV;
    }
  }

  void alloc();


  inline void communicate(SingleSlater<T> &groundState) {
    this->groundState_ = &groundState;
    this->memManager_  = groundState.memManager();
    this->fileio_      = groundState.fileio();
  };


  // RT Propagation
  void doPropagation();
  void formUTrans();
  void formField(const long double);
  void propDen();
  void addRecord(const long double);

  // Extra Properties
  void orbPop(std::vector<double>&,std::vector<double>&);

  // Print Functions
  void printRTStep(const long double);
  void printRTHeader();

  // CSVs
  void initCSV();
  void writeDipoleCSV();
  void writeAppliedFieldCSV();
  void writeMullikenCSV();
  void writeOrbitalCSV();
  void tarCSVFiles();

  inline void writeCSVs() {
    this->writeDipoleCSV();
    this->writeAppliedFieldCSV();
    this->writeMullikenCSV();
    this->writeOrbitalCSV();
  }

  inline void resetPropagation() {
    for(auto iD = 0; iD < POSav_.size(); iD++){
     (*ssPropagator_->onePDM()[iD]) = 
       groundState_->onePDM()[iD]->template cast<dcomplex>();
     (*ssPropagator_->onePDMOrtho()[iD]) = 
       groundState_->onePDMOrtho()[iD]->template cast<dcomplex>();
    }
    propInfo.clear();
  }

  // Setters
    
  void setIntScheme(const std::string&); 
  void setMMUTRstScheme(const std::string&); 
  void setEnvlp(const std::string&);
  void setPropMeth(const std::string&);
    
  inline void setMaxSteps(int x)  { maxSteps_ = x; };
  inline void setNSkip(int x)     { nSkip_    = x; };
  inline void setIRstrt(int x)    { iRstrt_   = x; };
  inline void setPrintLevel(int x){ printLevel_ = x;};
  inline void setNPolyExpMax(size_t x){ nPolyExpMax_ = x; };
             
  inline void setTOn(double x)      { tOn_      = x; };
  inline void setTOff(double x)     { tOff_     = x; };
  inline void setStepSize(double x) { stepSize_ = x; };
  inline void setDeltaT(double x)   { deltaT_   = x; };
  inline void setFreq(double x)     { freq_     = x; };
  inline void setPhase(double x)    { phase_    = x; };
  inline void setSigma(double x)    { sigma_    = x; };
  inline void setPolyEps(double x)  { polyEps_  = x; };

  inline void setEDFieldAmp(std::array<double,3> x) { staticEDAmp_ = x; };

  inline void doNotTarCSV() { tarCSVs = false;};

/*
  void setIMethFormU(PropagatorFormation x ) { this->iMethFormU_ = x;};
  void setIEnvlp(Envelope x                ) { this->iEnvlp_     = x;};
  void setIEllPol(EllipticalPolarization x ) { this->iEllPol_    = x;};
*/

/*
  // Place holder for when we implement more matrix exponentials
  void setFormU(const std::string &x) {

  }
*/

  // Python API
  inline void Wrapper_setFieldAmp(double x, double y, double z){
    this->setEDFieldAmp({{x,y,z}});
  }

  boost::python::list Python_propEnergy();
  boost::python::list Python_propDipole();
};

}; // namespace ChronusQ

#endif
