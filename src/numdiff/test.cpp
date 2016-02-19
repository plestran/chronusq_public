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


int main(int argc, char*argv[]){
  
  // Input parameters
  std::string scan1MinStrx  = argv[1];
  std::string scan1MaxStrx  = argv[2];
  std::string scan1StepStrx = argv[3];
  std::string scan2MinStrx  = argv[4];
  std::string scan2MaxStrx  = argv[5];
  std::string scan2StepStrx = argv[6];

  std::string basisName   = argv[7];
  std::string respTypeStr = argv[8];
  std::string chargeStr   = argv[9];
  std::string nFreqStr    = argv[10];
  std::string pxStr       = argv[11];
  std::string pyStr       = argv[12];

  std::string PREFIX      = argv[13];

  // Parse inout parameters
  double scan1Min  = std::atof(scan1MinStrx.c_str()); 
  double scan1Max  = std::atof(scan1MaxStrx.c_str()); 
  double scan1Step = std::atof(scan1StepStrx.c_str());
  double scan2Min  = std::atof(scan2MinStrx.c_str()); 
  double scan2Max  = std::atof(scan2MaxStrx.c_str()); 
  double scan2Step = std::atof(scan2StepStrx.c_str());

  RESPONSE_TYPE respType;
  if(!respTypeStr.compare("CIS")) respType = RESPONSE_TYPE::CIS;

  int charge = std::atoi(chargeStr.c_str());
  int nFreq  = std::atoi(nFreqStr.c_str());
  int px     = std::atoi(pxStr.c_str());
  int py     = std::atoi(pyStr.c_str()); 

  std::vector<double> scan1 = genScan(scan1Min,scan1Max,scan1Step);
  std::vector<double> scan2 = genScan(scan2Min,scan2Max,scan2Step);
  
  twoDScan(scan1,scan2,basisName,respType,charge,
    nFreq,px,py,PREFIX);

  return 0;
}
