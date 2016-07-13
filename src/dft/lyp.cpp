#include<dft.h>

lyp::lyp(){
// Memo Factor to be added at the end for numerical stability
  this-> Cfact = 36.4623989787648;  //-(2.0^(11.0/3.0))*(3.0/10.0)*((3.0*pi*pi)^(2.0/3.0))
  this-> small = 1.0e-10; 
  this-> d1over3  = 1.0/3.0;
  this-> a = 0.04918;
  this-> b = 0.132;
  this-> c = 0.2533;
  this-> d = 0.349;
};
  std::vector<double> denspow(const double &rho)  {
    std::vector<double> rhopow;
    return rhopow;
    }; 

void lyp::popLYPdens(double rhoA, double rhoB){
  this->rhoT          = rhoA + rhoB;
//  this->spindensity   = (rhoA - rhoB) / this->rhoT;
//  if( this->rhoT <  this-> small) {return ;}
//  cout << "rhoT " << this->rhoT <<endl;
//  cout << "rhoMz " << this->spindensity <<endl;
  this->rho1over3 = std::pow(this->rhoT,this->d1over3);
  this->rhoA8over3 = std::pow(rhoA,(8.0/3.0));
  this->rhoB8over3 = std::pow(rhoB,(8.0/3.0));
//function used in LYP Correlation (Eq. A26, A32 Ref.LYP2)
  this->omega0  =  std::exp(-this->c / this->rho1over3);
  this->omega0 /=  1.0 + (this->d / this->rho1over3);
  this->omega0 *= std::pow(this->rhoT,(-11.0/3.0));
  this->omega1  = this->omega0*(11.0*this->rho1over3 -this->c - 
    (this->d / (1 +this->d/this->rho1over3) )  );
  this->omega1 /= -3.0*this->rhoT*this->rho1over3 ;
//function used in LYP Correlation (Eq. A27, A33 Ref.LYP2)
  this->delta0  = this->d / this->rho1over3;
  this->delta0 /= (1.0 + this->d / this->rho1over3);
  this->delta0 += this->c / this->rho1over3;
  this->delta1  =  this->delta0/(-3.0*this->rhoT);
  this->delta1 += (this->d *this->d / std::pow(this->rhoT,(5.0/3.0))) 
    /( 3.0 * (1.0 + this->d / this->rho1over3) 
    * (1.0 + this->d / this->rho1over3) );
  // Eq A23(delLYP/delgammaAA) and eq A25(delLYP/delgammaBB) (remind for A25 call the function inverting rhoA with rhoB)
     this->dLYPdgAA   = 11.0 - this->delta0;
     this->dLYPdgAA  *= rhoA/this->rhoT;
     this->dLYPdgAA  += 1.0 - 3.0*this->delta0;
     this->dLYPdgAA  *= rhoA*rhoB/9.0;
     this->dLYPdgAA  -= rhoB*rhoB;
     this->dLYPdgAA  *= -this->a * this->b * this->omega0;
     this->dLYPdgBB   = 11.0 - this->delta0;
     this->dLYPdgBB  *= rhoB/this->rhoT;
     this->dLYPdgBB  += 1.0 - 3.0*this->delta0;
     this->dLYPdgBB  *= rhoB*rhoA/9.0;
     this->dLYPdgBB  -= rhoA*rhoA;
     this->dLYPdgBB  *= -this->a * this->b * this->omega0;
//  Eq A24 (delLYP/delgammaAB)
      this->dLYPdgAB  = 47.0 - 7.0 *  this->delta0;
      this->dLYPdgAB *= rhoA * rhoB / 9.0;
      this->dLYPdgAB -= 4.0 *  this->rhoT *  this->rhoT / 3.0;
      this->dLYPdgAB *= - this->a *  this->b * this->omega0;
//  Eq A29 (del^2 LYP / (del rho_X del gammaXX))
      this->d2LYPdrhgAA =   - (rhoA*rhoB/9.0) 
             *(  ( 3.0 + rhoA/this->rhoT ) * 
              this->delta1 + rhoB*( this->delta0 - 11.0)/
              (this->rhoT*this->rhoT)  ) 
           + (rhoB/9.0)
             *( 1.0 - 3.0*this->delta0 -rhoA*(this->delta0 - 11)/
             this->rhoT );
    this->d2LYPdrhgAA *= -this->a *  this->b * this->omega0;
    this->d2LYPdrhgAA += this->dLYPdgAA*this->omega1/this->omega0;
//  Eq A30 (del^2 LYP / (del rho_X del gammaXY))
    this->d2LYPdrhgAB  =  - 8.0*this->rhoT/3.0;
    this->d2LYPdrhgAB +=  -  this->delta1*(7.0*rhoA*rhoB/9.0);
    this->d2LYPdrhgAB +=  rhoB*(47.0-7.0*this->delta0)/9.0;
    this->d2LYPdrhgAB *= -this->a *  this->b * this->omega0;
    this->d2LYPdrhgAB += this->dLYPdgAB*this->omega1/this->omega0;
//  Eq A31 (del^2 LYP / (del rho_X del gammaYY)) (debugged alread)
    this->d2LYPdrhgBB =   - (rhoA*rhoB/9.0) 
             *(  ( 3.0 + rhoB/ this->rhoT ) *  this->delta1 - rhoB*
           (this->delta0 - 11.0)/( this->rhoT* this->rhoT)  ) 
           + (rhoB/9.0)
           * ( 1.0 - 3.0* this->delta0 
           - rhoB*(this->delta0 - 11)/ this->rhoT )
           - 2.0 * rhoA;
    this->d2LYPdrhgBB *= -this->a *  this->b * this->omega0;
    this->d2LYPdrhgBB += this->dLYPdgBB*this->omega1/this->omega0;
};  //End Form denisuty related  LYP Corr


DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB){
};

DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB){
};

DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
  DFTFunctional::DFTInfo info;
//  if( rhoA <  this-> small) {return info;}
//  if( rhoB <  this-> small) {return info;}
//  if( (rhoA+rhoB) <  this-> small) {return info;}
  this->popLYPdens(rhoA, rhoB);
//  if( gammaAA <  this-> small) {return info ;}
//  if( gammaAB <  this-> small) {return info ;}
//  if( gammaBB <  this-> small) {return info ;}
//  Eq. A23  dLYP/dgammaAA (debugged)
  info.ddgammaAA  = this->dLYPdgAA ;
//  Eq. A23* dLYP/dgammaBB (debugged)
  info.ddgammaBB  = this->dLYPdgBB ;
//  Eq. A24  dLYP/dgammaAB (debugged)
  info.ddgammaAB  = this->dLYPdgAB ;
//  Eq. A28  dLYP/dRhoA (debugged)
  info.ddrhoA = this->dLYPdgAA ;
    info.ddrhoA   = - 4.0 *  this->a * rhoA * rhoB 
      / ( ( this->rhoT)*(1.0 +  this->d/this->rho1over3 ) );
    info.ddrhoA  *= ( (1.0/rhoA)
                  -(1.0/this->rhoT) 
                  +((this->d/3.0) 
                  * (std::pow((this->rhoT),(-4.0/3.0)))/ (1.0 + this->d 
                  * std::pow(this->rhoT,(-1.0/3.0))))
                 );
    info.ddrhoA  += - this->Cfact * this->a * this->b 
             *( 
              (this->omega1 * rhoA * rhoB 
               * (this->rhoA8over3 + this->rhoB8over3))
              +(this->omega0 * rhoB * ( (11.0*this->rhoA8over3/3.0) 
               + this->rhoB8over3))
             );
    info.ddrhoA  += gammaAA* this->d2LYPdrhgAA;
    info.ddrhoA  += gammaAB* this->d2LYPdrhgAB;
    info.ddrhoA  += gammaBB* this->d2LYPdrhgBB;
//  Eq. A28* dLYP/dRhoB (debugged)
    info.ddrhoB   = - 4.0 *  this->a * rhoA * rhoB 
      / ( ( this->rhoT)*(1.0 +  this->d/this->rho1over3 ) );
    info.ddrhoB  *= ( (1.0/rhoB)
                  -(1.0/this->rhoT) 
                  +((this->d/3.0) 
                  * (std::pow((this->rhoT),(-4.0/3.0)))/ (1.0 + this->d 
                  * std::pow(this->rhoT,(-1.0/3.0))))
                 );
    info.ddrhoB  += - this->Cfact * this->a * this->b 
             *( 
              (this->omega1 * rhoA * rhoB 
               * (this->rhoA8over3 + this->rhoB8over3))
              +(this->omega0 * rhoA * ( (11.0*this->rhoB8over3/3.0) 
               + this->rhoA8over3))
             );
    info.ddrhoB  += gammaBB* this->d2LYPdrhgAA;
    info.ddrhoB  += gammaAB* this->d2LYPdrhgAB;
    info.ddrhoB  += gammaAA* this->d2LYPdrhgBB;
    info.eps  = - 4.0 * this->a * rhoA * rhoB 
      / ( this->rhoT*(1.0 + this->d / std::pow(this->rhoT,this-> d1over3) ) );
    info.eps += - this->Cfact * this->a * this->b 
      * this->omega0 * rhoA * rhoB * (this->rhoA8over3 + this->rhoB8over3); 
    info.eps += info.ddgammaAA * gammaAA;
    info.eps += info.ddgammaBB * gammaBB;
    info.eps += info.ddgammaAB * gammaAB;
/*
    cout << "eps " << info.eps <<endl;
    cout << "ddrhoA " <<  info.ddrhoA <<endl;
    cout << "ddrhoB " <<  info.ddrhoB <<endl;
    cout << "ddgammaAA " <<  info.ddgammaAA <<endl;
    cout << "ddgammaAB " <<  info.ddgammaAB <<endl;
    cout << "ddgammaBB " <<  info.ddgammaBB <<endl;
*/
    return info;
};

