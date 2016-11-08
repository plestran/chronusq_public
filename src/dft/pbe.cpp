#include<dft.h>

PBE::PBE(double X, double eps):
DFTFunctional(X,eps){
// Memo Factor to be added at the end for numerical stability
  this->CxVx  =  0.930525736349100;  // (3/2)*((3/(4*pi))^(1/3)) ;  
//  this->small = 1.0e-18; 
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;
  this-> mu      =  0.2195102405736;
  this-> kappa    =  0.2195102405736;
  this-> Ckf       =  std::pow((math.pi*math.pi*6.0),(this->d1over3));
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
  return  this->CxVx*gx;
};  //End Form g(x) or F9in Barone)

double PBE::g1 (double x, double &den){
  double gx;
  gx = 2.0 * this->mu * x ;
  gx /= den;
  gx /= den;
//  return this->CxVx*gx;
  return gx;
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
/*
  cout << "Ax " <<this->CxVx << endl;
  cout << "b "  <<this->b << endl;
  cout << "c "  <<this->c << endl;
  cout << "d "  <<this->d << endl;
  cout << "beta "  <<this->beta << endl;
  cout << "xA " <<xA << " g0 " << g0 << " g1 " << g1 <<endl;
*/
  info.eps        = rhoA4ov3*g0;
  info.ddrhoA     = g0 - xA*g1;
  info.ddrhoA    *= this->d4over3*rhoA1ov3;
  info.ddgammaAA  = g1/(4.0*this->Ckf*std::sqrt(gammaAA));
  if(std::abs(spindensity) > this->small) {
    //Open Shell
  //Paper  xB   = gammaBB / rhoA4ov3; 
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
//  cout <<  "eps " <<info.eps << endl ;
//  cout <<  "ddrhoA " <<info.ddrhoA << endl ;
//  cout <<  "ddrhoB " <<info.ddrhoB << endl ;
//  cout <<  "ddgammaAA" <<info.ddgammaAA << endl ;
//  cout <<  "ddgammaBB" <<info.ddgammaBB << endl ;
  info.ddgammaAB  = 0.0;
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

