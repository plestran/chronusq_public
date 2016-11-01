#include<dft.h>

PBE::PBE(double X, double eps):
DFTFunctional(X,eps){
// Memo Factor to be added at the end for numerical stability
  this->CxVx  =  0.930525736349100;  // (3/2)*((3/(4*pi))^(1/3)) ;  
//  this->small = 1.0e-18; 
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;
  this-> beta  =  std::pow((36*math.pi),(5.0/3.0));
  this-> beta *=  5;
  this-> b     =  0.0042;
//  cout << "B88 object created " <<endl;

  this->name = "PBE";
#ifdef CQ_ENABLE_LIBXC
  xc_func_init(&this->func,XC_GGA_X_B88,XC_POLARIZED);
#endif
};


DFTFunctional::DFTInfo PBE::eval(const double &rhoA, const double &rhoB){
};

DFTFunctional::DFTInfo PBE::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
  DFTFunctional::DFTInfo info;
#ifndef CQ_ENABLE_LIBXC
  double spindensity   = (rhoA - rhoB); 
  spindensity         /= (rhoA + rhoB);
  double rhoA1ov3 = std::pow(rhoA,this->d1over3);
  double rhoB1ov3 ;
//rhoA4ov3 = std::pow(rhoA,this->d4over3);
  double rhoA4ov3 = rhoA1ov3 * rhoA; 
  double rhoB4ov3 ; 
  if(std::abs(spindensity) > this->small) {
    rhoB1ov3 = std::pow(rhoB,this->d1over3);
    rhoB4ov3 = rhoB1ov3 * rhoB; 
  }
// Note that in Eq A5 xA   = gammaAA / rhoA4ov3; 
// but actually they meants xA = sqrt(gammaAA) /rhoA4ov3
// and also eq A6 rho is rho^(1/3) instead of rho^(4/3)
#else
  std::array<double,2> rho = {rhoA,rhoB};
  std::array<double,3> sigma = {gammaAA,gammaAB,gammaBB};
  std::array<double,2> vrho;
  std::array<double,3> vsigma;
  xc_gga_exc_vxc(&this->func,1,&rho[0],&sigma[0],&info.eps,&vrho[0],
    &vsigma[0]);

  info.ddrhoA = vrho[0];
  info.ddrhoB = vrho[1];
  info.ddgammaAA = vsigma[0];
  info.ddgammaAB = vsigma[1];
  info.ddgammaBB = vsigma[2];

  info.eps *= (rhoA + rhoB); //?
#endif
  info *= this->scalingFactor; 
  return info;
};

DFTFunctional::DFTInfo PBE::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB){
};

