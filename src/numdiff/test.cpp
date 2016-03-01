#define EIGEN_DONT_PARALLELIZE
#include <response.h>
#include <workers.h>
#include <pythonapi.h>

template<typename T>
std::vector<T> genScan(T min,T max,T step){
  std::vector<T> scan;
  for(T i = min; i <= max; i += step){
    scan.push_back(i);
  }
  return scan;
};

using namespace ChronusQ;
void twoDScan(std::vector<double>&,std::vector<double>&,std::string&,RESPONSE_TYPE,int,int,int,int,
std::string&);
void oneDScan(std::vector<double>&,std::string&,RESPONSE_TYPE,int,int,int,std::string&);

void cartDiff(std::string&,std::string&,RESPONSE_TYPE,int,int);

int main(int argc, char*argv[]){
  
  std::string scan1MinStrx,scan1MaxStrx,scan1StepStrx,scan2MinStrx,
    scan2MaxStrx,scan2StepStrx,basisName,respTypeStr,chargeStr,
    nFreqStr,pxStr,pyStr,PREFIX,XYZFName; 

  bool do2D;
  bool doScan;
  // Input parameters
  if(argc == 14) {
    do2D = true;
    doScan = true;
    scan1MinStrx  = argv[1];
    scan1MaxStrx  = argv[2];
    scan1StepStrx = argv[3];
    scan2MinStrx  = argv[4];
    scan2MaxStrx  = argv[5];
    scan2StepStrx = argv[6];

    basisName   = argv[7];
    respTypeStr = argv[8];
    chargeStr   = argv[9];
    nFreqStr    = argv[10];
    pxStr       = argv[11];
    pyStr       = argv[12];

    PREFIX      = argv[13];
  } else if(argc == 10) {
    do2D = false;
    doScan = true;
    scan1MinStrx  = argv[1];
    scan1MaxStrx  = argv[2];
    scan1StepStrx = argv[3];

    basisName   = argv[4];
    respTypeStr = argv[5];
    chargeStr   = argv[6];
    nFreqStr    = argv[7];
    pxStr       = argv[8];

    PREFIX      = argv[9];
  } else if(argc == 6) {
    doScan = false;
    XYZFName  = argv[1];
    basisName = argv[2];  
    respTypeStr = argv[3];
    chargeStr = argv[4];
    nFreqStr = argv[5];
  } else {
    cout << "Incorrect arguements" << endl;
    exit(EXIT_FAILURE);
  }

  double scan1Min , scan1Max , scan1Step, scan2Min , scan2Max , scan2Step; 

  if(doScan) {
    // Parse inout parameters
    scan1Min  = std::atof(scan1MinStrx.c_str()); 
    scan1Max  = std::atof(scan1MaxStrx.c_str()); 
    scan1Step = std::atof(scan1StepStrx.c_str());
    if(do2D){
      scan2Min  = std::atof(scan2MinStrx.c_str()); 
      scan2Max  = std::atof(scan2MaxStrx.c_str()); 
      scan2Step = std::atof(scan2StepStrx.c_str());
    }
  }

  RESPONSE_TYPE respType;
  if(!respTypeStr.compare("CIS")) respType = RESPONSE_TYPE::CIS;
  if(!respTypeStr.compare("PPTDA")) respType = RESPONSE_TYPE::PPTDA;

  int charge = std::atoi(chargeStr.c_str());
  int nFreq  = std::atoi(nFreqStr.c_str());

  int px,py;

  if(doScan) {
    px     = std::atoi(pxStr.c_str());
    if(do2D)
       py     = std::atoi(pyStr.c_str()); 
  }

  std::vector<double> scan1,scan2;
  if(doScan) {
    scan1 = genScan(scan1Min,scan1Max,scan1Step);
    if(do2D)
      scan2 = genScan(scan2Min,scan2Max,scan2Step);
  }
  

  if(doScan) {
    if(do2D)
      twoDScan(scan1,scan2,basisName,respType,charge,
        nFreq,px,py,PREFIX);
    else
      oneDScan(scan1,basisName,respType,charge,nFreq,px,PREFIX);
  } else {
    cartDiff(XYZFName,basisName,respType,charge,nFreq);
  }

  return 0;
}
