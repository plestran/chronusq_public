#include<dft.h>

SlaterExchange::SlaterExchange(){
// Memo Factor to be added at the end for numerical stability
  this->CxVx  = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));  
  this->small = 1.0e-12; 
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;
  cout << "Created Slater Object" << endl;
};

DFTFunctional::DFTInfo SlaterExchange::eval(double rhoA, double rhoB){
  DFTFunctional::DFTInfo info;
  this->rhoT            = rhoA + rhoB;
  this->spindensity     = (rhoA - rhoB) / this->rhoT;
  info.eps        = std::pow(this->rhoT,this->d1over3)*this->CxVx;      
  if (this->spindensity > this->small){
// Open Shell
    info.ddrhoA = 
      this->d4over3*info.eps*std::pow((1.0+this->spindensity),this->d1over3); 
    info.ddrhoB = 
      this->d4over3*info.eps*std::pow((1.0-this->spindensity),this->d1over3);   
  } else {   
// Closed Shell
    info.ddrhoA     = this->d4over3*info.eps;       
  }
  return info;
};

DFTFunctional::DFTInfo SlaterExchange::eval(double rhoA, double rhoB, double gammaAA, double gammaAB){
};

DFTFunctional::DFTInfo SlaterExchange::eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB){
};

