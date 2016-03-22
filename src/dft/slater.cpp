#include<dft.h>

SlaterExchange::SlaterExchange(){
  this->CxVx  = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));  
  this->small = 1.0e-12; 
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;
};

DFTFunctional::DFTInfo SlaterExchange::eval(double rhoA, double rhoB){
  DFTFunctional::DFTInfo info;
  double rhoT     = rhoA + rhoB;
  double spindensity     = (rhoA - rhoB) / rhoT;
  info.eps        = std::pow(rhoT,d1over3);      
  if (spindensity > this->small){
  info.ddrhoA     = this->d4over3*info.eps*std::pow((1.0+spindensity),this->d1over3);       
  info.ddrhoB     = this->d4over3*info.eps*std::pow((1.0-spindensity),this->d1over3);   
  } else{   
  info.ddrhoA     = this->d4over3*info.eps;       
  }
  return info;
};

DFTFunctional::DFTInfo SlaterExchange::eval(double rhoA, double rhoB, double gammaAA, double gammaAB){
};

DFTFunctional::DFTInfo SlaterExchange::eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB){
};

