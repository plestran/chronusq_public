#include <dft.h>

SlaterExchange::SlaterExchange(double X, double eps):
DFTFunctional(X,eps){
// Memo Factor to be added at the end for numerical stability
  this->CxVx  = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));  
  this->small = 1.0e-16; 
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;

  this->name = "Slater";
#ifdef CQ_ENABLE_LIBXC
  xc_func_init(&this->func,XC_LDA_X,XC_POLARIZED);
#endif
};

DFTFunctional::DFTInfo SlaterExchange::eval(const double &rhoA, const double &rhoB){

};

DFTFunctional::DFTInfo SlaterExchange::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB){
};

DFTFunctional::DFTInfo SlaterExchange::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){

  DFTFunctional::DFTInfo info;
  double rhoT            = rhoA + rhoB;
#ifndef CQ_ENABLE_LIBXC
  double spindensity     = (rhoA - rhoB) / rhoT;
  double eps_spin;
  double tmp1Sppow;
  double tmp2Sppow;
  info.eps        = std::pow(rhoT,this->d1over3)*this->CxVx;      
  if (std::abs(spindensity) > this->small){
// Open Shell
    tmp1Sppow  = std::pow((1.0+spindensity),this->d1over3);
    tmp2Sppow  = std::pow((1.0-spindensity),this->d1over3);
    eps_spin   = tmp1Sppow*(1.0+spindensity) ;      
    eps_spin  += tmp2Sppow*(1.0-spindensity) ;      
    eps_spin   /= 2.0;      
    info.ddrhoA = 
      this->d4over3*info.eps*tmp1Sppow; 
    info.ddrhoB = 
      this->d4over3*info.eps*tmp2Sppow;   
    info.eps *= eps_spin;
  } else {   
// Closed Shell
    info.ddrhoA     = this->d4over3*info.eps;       
    info.ddrhoB     = info.ddrhoA;       
  }
#else
  std::array<double,2> rho = {rhoA,rhoB};
  std::array<double,2> vxc;
  
  xc_lda_exc_vxc(&this->func,1,&rho[0],&info.eps,&vxc[0]);
  info.ddrhoA = vxc[0];
  info.ddrhoB = vxc[1];
#endif
  info.eps  *= rhoT; 
  return info;
};

