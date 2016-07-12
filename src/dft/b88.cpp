#include<dft.h>

BEightEight::BEightEight(){
// Memo Factor to be added at the end for numerical stability
  this->CxVx  =  0.930525736349100;  // (3/2)*((3/(4*pi))^(1/3)) ;  
  this->small = 1.0e-12; 
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;
  this-> beta =  0.0042;
//  cout << "B88 object created " <<endl;
};

double BEightEight::g0B88 (double x){
  double gx;
  double bx = this->beta * x;
//  Eq A4 Pople, J. Chem. Phys. 5612, (1992) nder =0
 gx  = -bx*x;
 gx /= (1.0 + 6.0 *bx * boost::math::asinh(x)) ;
 gx -= this->CxVx;
 return gx; 
};  //End Form g function for B88 Exchange

double BEightEight::g1B88 (double x){
  double gx;
  double asx = boost::math::asinh(x);
  double bx = this->beta * x;
  double denom = (1.0 + 6.0 * bx * asx);

//  Eq A8 Pople, J. Chem. Phys. 5612, (1992)
/*
 gx  = x / (std::sqrt(x*x+1.0));
 gx -= asx ;
 gx *= 6.0 * bx * bx;
 gx -= 2.0 * bx;
 gx /= denom * denom; 
*/
 gx  = x / (std::sqrt(x*x+1.0));
 gx -= asx ;
 gx *= 3.0*bx;
 gx -= 1.0;
 gx *= 2.0*bx;
 gx /= denom * denom; 
 return gx; 
};  //End Form g function for B88 Exchange

DFTFunctional::DFTInfo BEightEight::eval(const double &rhoA, const double &rhoB){
};

DFTFunctional::DFTInfo BEightEight::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
  DFTFunctional::DFTInfo info;
  this->rhoT          = rhoA + rhoB;
  this->spindensity   = (rhoA - rhoB) / this->rhoT;
  rhoA1ov3 = std::pow(rhoA,this->d1over3);
//rhoA4ov3 = std::pow(rhoA,this->d4over3);
  rhoA4ov3 = rhoA1ov3 * rhoA; 
  if(std::abs(spindensity) > this->small) {
    rhoB1ov3 = std::pow(rhoB,this->d1over3);
    rhoB4ov3 = rhoB1ov3 * rhoB; 
  }
// Note that in Eq A5 xA   = gammaAA / rhoA4ov3; 
// but actually they meants xA = sqrt(gammaAA) /rhoA4ov3
// and also eq A6 rho is rho^(1/3) instead of rho^(4/3)
  this->xA   = std::sqrt(gammaAA) / rhoA4ov3; 
  info.eps        = rhoA4ov3*this->g0B88(xA);
  info.ddrhoA     = this->g0B88(xA) - xA*this->g1B88(xA);
  info.ddrhoA    *= this->d4over3*rhoA1ov3;
  info.ddgammaAA  = 0.5*(this->g1B88(xA))/(std::sqrt(gammaAA));
  if(std::abs(this->spindensity) > this->small) {
    //Open Shell
  //Paper  xB   = gammaBB / rhoA4ov3; 
    this->xB   = std::sqrt(gammaBB) / rhoB4ov3; 
    info.eps       += rhoB4ov3*this->g0B88(xB);
    info.ddrhoB     = this->g0B88(xB) - xB*this->g1B88(xB);
    info.ddrhoB    *= this->d4over3*rhoB1ov3;
    info.ddgammaBB  = 0.5*(this->g1B88(xB))/(std::sqrt(gammaBB));

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
  return info;
};

DFTFunctional::DFTInfo BEightEight::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB){
};

