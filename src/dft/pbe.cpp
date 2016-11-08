#include<dft.h>

PBE::PBE(double X, double eps):
DFTFunctional(X,eps){
// Memo Factor to be added at the end for numerical stability
  this->CxVx      = -(3.0/2.0)*(std::pow((3.0/(4.0*math.pi)),(1.0/3.0)));  
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;
  this-> mu       =  0.2195102405736;
  this-> kappa    =  0.8040004238;
  this-> Ckf      =  std::pow((math.pi*math.pi*6.0),(this->d1over3));
//  this-> Ckf       =  std::pow((math.pi*math.pi*3.0),(this->d1over3));
//  cout << "B88 object created " <<endl;

  this->name = "PBE";
#ifdef CQ_ENABLE_LIBXC
  xc_func_init(&this->func,XC_GGA_X_PBE,XC_POLARIZED);
#endif
};

double PBE::g0 (double x,double &den){
  double gx;
  den  = this->mu * x * x;
  den /= this->kappa;
  den += 1.0; 
  gx = -this->kappa;
  gx /= den;
  gx += 1.0 + this->kappa;
  return this->CxVx*gx;
};  //End Form g(x) or F9in Barone)

double PBE::g1 (double x, double &den){
  double gx;
  gx = 2.0 * this->mu * x ;
  gx /= (den*den);
  return this->CxVx*gx;
};  //End Form g function for B88 Exchange

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
  double xA   = std::sqrt(gammaAA) / (2.0 * this->Ckf*rhoA4ov3); 
  double xB ;
  double den;
  double g0 = this->g0(xA,den);
  double g1 = this->g1(xA,den);
  info.eps        = rhoA4ov3*g0;
  info.ddrhoA     = g0 - xA*g1;
  info.ddrhoA    *= this->d4over3*rhoA1ov3;
  info.ddgammaAA  = g1/(4.0*this->Ckf*std::sqrt(gammaAA));
  if(std::abs(spindensity) > this->small) {
    //Open Shell
    xB   = std::sqrt(gammaBB) / (2.0 * this->Ckf*rhoB4ov3); 
    g0 = this->g0(xB,den);
    g1 = this->g1(xB,den);
    info.eps       += rhoB4ov3*g0;
    info.ddrhoB     = g0 - xB*g1;
    info.ddrhoB    *= this->d4over3*rhoB1ov3;
    info.ddgammaBB  = g1/(4.0*this->Ckf*std::sqrt(gammaBB));

  } else {
    //Closed Shell
    info.eps *= 2.0;
    info.ddrhoB     = info.ddrhoA;
    info.ddgammaBB  = info.ddgammaAA;
  }
  info.ddgammaAB  = 0.0;
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

