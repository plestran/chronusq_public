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


  dcomplex * NBSqScratch_;
  dcomplex * NBTSqScratch_;
  dcomplex * NBTSqScratch2_;
  dcomplex * UTransScalar_;
  dcomplex * UTransMz_;
  dcomplex * UTransMy_;
  dcomplex * UTransMx_;

  // RT Simulation Control
  int maxSteps_;
  int nSkip_;
  int iRstrt_;

  double tOn_;
  double tOff_;
  double stepSize_;
  double deltaT_;
  double freq_;
  double phase_;
  double sigma_;

  long double currentTime_;

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
  
  struct PropInfo {
    double timeStep;
    double energy;
    std::array<double,4> dipole;
    std::array<double,4> appliedfield;
    std::vector<double> mullPop;
    std::vector<double> orbitalOccA;
    std::vector<double> orbitalOccB;
  };
  std::vector<PropInfo> propInfo;

  RealTime() : fileio_(NULL), groundState_(NULL), memManager_(NULL),
    ssPropagator_(nullptr),
    deltaT_(0.0), currentTime_(0.0), freq_(0.0), phase_(0.0), sigma_(0.0),
    nSkip_(0),
    tOn_(0.0), tOff_(1.0e4), maxSteps_(10), stepSize_(0.05), iRstrt_(50),
    iMethFormU_(EigenDecomp), iEnvlp_(Constant), iEllPol_(LXZ),
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


  void doPropagation();
  void formUTrans();
  void formField();
  void propDen();
};

#include <realtime/realtime_alloc.h>
#include <realtime/realtime_field.h>
#include <realtime/realtime_propagator.h>
#include <realtime/realtime_proc.h>
}; // namespace ChronusQ

#endif
