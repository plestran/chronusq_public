#include<dft.h>

lyp::lyp(double X, double eps):
DFTFunctional(X,eps){
// Memo Factor to be added at the end for numerical stability
  this-> Cfact = 36.4623989787648;  //-(2.0^(11.0/3.0))*(3.0/10.0)*((3.0*pi*pi)^(2.0/3.0))
  this-> small = 1.0e-18; 
  this-> d1over3  = 1.0/3.0;
  this-> a = 0.04918;
  this-> b = 0.132;
  this-> c = 0.2533;
  this-> d = 0.349;

  this->name = "LYP";
};


    void lyp::popDensPow(const double &rhoA, const double &rhoB, denspow &denquant)  {
    double tmp = 8.0/3.0;
    denquant.rhoT        = rhoA + rhoB;
    denquant.spindensity = (rhoA - rhoB) / denquant.rhoT;
    denquant.rho1over3   = std::pow(denquant.rhoT,this->d1over3);
//    if(denquant.rhoT < this->small) {
//      denquant.zeta = 1.0;
//    }else{
      denquant.zeta        = 1.0/(1.0 + this->d / denquant.rho1over3);
//    }
    denquant.rhoA8over3   = std::pow(rhoA,tmp);
//    if(std::abs(spindensity) > this->small) 
    denquant.rhoB8over3   = std::pow(rhoB,tmp);
    denquant.rho4over3  = (denquant.rho1over3*denquant.rhoT);
    denquant.rho2        = denquant.rhoT*denquant.rhoT;
//    denquant.rho5over3   = denquant.rho1over3*denquant.rho4over3;
//    cout << denquant.rho1over3 << endl ;
//    cout << denquant.rho2 << endl ;
//    cout << denquant.rho5over3 << endl ;
//    cout << denquant.rhom11over3 << endl ;
    return ;
    }; 

void lyp::popLYPdens(const double &rhoA, const double &rhoB, denspow &RhoTQuant){
//  this->spindensity   = (rhoA - rhoB) / this->rhoT;
//function used in LYP Correlation (Eq. A26, A32 Ref.LYP2)
//A  if(RhoTQuant.rhoT < this->small) {
    RhoTQuant.omega0 = 0.0;
    RhoTQuant.omega1 = 0.0;
//A  }else{
    RhoTQuant.omega0  =  std::exp(-this->c / RhoTQuant.rho1over3);
    RhoTQuant.omega0 *=  RhoTQuant.zeta;
    RhoTQuant.omega0 /=  RhoTQuant.rho2*RhoTQuant.rho4over3*RhoTQuant.rho1over3; //^(11/3)
    RhoTQuant.omega1  = this->c ;
    RhoTQuant.omega1 += this->d*RhoTQuant.zeta;
    RhoTQuant.omega1 -= 11.0*RhoTQuant.rho1over3;
    RhoTQuant.omega1 /= 3.0 *RhoTQuant.rho4over3;
    RhoTQuant.omega1 *= RhoTQuant.omega0;
//A  }
//A  if(RhoTQuant.rhoT < this->small) {
    RhoTQuant.delta0 = 0.0;
    RhoTQuant.delta1 = 0.0;
  //function used in LYP Correlation (Eq. A27, A33 Ref.LYP2)
//A  }else{ 
    RhoTQuant.delta0   = this->c + this->d*RhoTQuant.zeta;
    RhoTQuant.delta0  /= RhoTQuant.rho1over3;
    RhoTQuant.delta1 = (RhoTQuant.zeta* RhoTQuant.zeta* this->d* this->d/ (RhoTQuant.rho1over3*RhoTQuant.rho1over3) ); 
    RhoTQuant.delta1 -= RhoTQuant.delta0;
    RhoTQuant.delta1 /= (3.0*RhoTQuant.rhoT);
//A  }
  // Eq A23(delLYP/delgammaAA) and eq A25(delLYP/delgammaBB) (remind for A25 call the function inverting rhoA with rhoB)
//     if(std::abs(spindensity) > this->small) {
//     UKS
       RhoTQuant.dLYPdgAA   = 11.0 - RhoTQuant.delta0;
       RhoTQuant.dLYPdgAA  *= rhoA/RhoTQuant.rhoT;
       RhoTQuant.dLYPdgAA  += 1.0 - 3.0*RhoTQuant.delta0;
       RhoTQuant.dLYPdgAA  *= rhoA*rhoB/9.0;
       RhoTQuant.dLYPdgAA  -= rhoB*rhoB;
       RhoTQuant.dLYPdgAA  *= -this->a * this->b * RhoTQuant.omega0;
       RhoTQuant.dLYPdgBB   = 11.0 - RhoTQuant.delta0;
       RhoTQuant.dLYPdgBB  *= rhoB/RhoTQuant.rhoT;
       RhoTQuant.dLYPdgBB  += 1.0 - 3.0*RhoTQuant.delta0;
       RhoTQuant.dLYPdgBB  *= rhoB*rhoA/9.0;
       RhoTQuant.dLYPdgBB  -= rhoA*rhoA;
       RhoTQuant.dLYPdgBB  *= -this->a * this->b * RhoTQuant.omega0;
//     } else {




//      }
//  Eq A24 (delLYP/delgammaAB)
      RhoTQuant.dLYPdgAB  = 47.0 - 7.0 *  RhoTQuant.delta0;
      RhoTQuant.dLYPdgAB *= rhoA * rhoB / 9.0;
      RhoTQuant.dLYPdgAB -= 4.0 *  RhoTQuant.rho2/ 3.0;
      RhoTQuant.dLYPdgAB *= - this->a *  this->b * RhoTQuant.omega0;
//  Eq A29 (del^2 LYP / (del rho_X del gammaXX))
      RhoTQuant.d2LYPdrhgAA =   - (rhoA*rhoB/9.0) 
             *(  ( 3.0 + rhoA/RhoTQuant.rhoT ) * 
              RhoTQuant.delta1 + rhoB*( RhoTQuant.delta0 - 11.0)/
              (RhoTQuant.rho2)  ) 
           + (rhoB/9.0)
             *( 1.0 - 3.0*RhoTQuant.delta0 -rhoA*(RhoTQuant.delta0 - 11)/
             RhoTQuant.rhoT );
    RhoTQuant.d2LYPdrhgAA *= -this->a *  this->b * RhoTQuant.omega0;
    RhoTQuant.d2LYPdrhgAA += RhoTQuant.dLYPdgAA*RhoTQuant.omega1/RhoTQuant.omega0;
//  Eq A30 (del^2 LYP / (del rho_X del gammaXY))
    RhoTQuant.d2LYPdrhgAB  =  - 8.0*RhoTQuant.rhoT/3.0;
    RhoTQuant.d2LYPdrhgAB +=  -  RhoTQuant.delta1*(7.0*rhoA*rhoB/9.0);
    RhoTQuant.d2LYPdrhgAB +=  rhoB*(47.0-7.0*RhoTQuant.delta0)/9.0;
    RhoTQuant.d2LYPdrhgAB *= -this->a *  this->b * RhoTQuant.omega0;
    RhoTQuant.d2LYPdrhgAB += RhoTQuant.dLYPdgAB*RhoTQuant.omega1/RhoTQuant.omega0;
//  Eq A31 (del^2 LYP / (del rho_X del gammaYY)) (debugged alread)
    RhoTQuant.d2LYPdrhgBB =   - (rhoA*rhoB/9.0) 
             *(  ( 3.0 + rhoB/ RhoTQuant.rhoT ) *  RhoTQuant.delta1 - rhoB*
           (RhoTQuant.delta0 - 11.0)/( RhoTQuant.rho2)  ) 
           + (rhoB/9.0)
           * ( 1.0 - 3.0* RhoTQuant.delta0 
           - rhoB*(RhoTQuant.delta0 - 11)/ RhoTQuant.rhoT )
           - 2.0 * rhoA;
    RhoTQuant.d2LYPdrhgBB *= -this->a *  this->b * RhoTQuant.omega0;
    RhoTQuant.d2LYPdrhgBB += RhoTQuant.dLYPdgBB*RhoTQuant.omega1/RhoTQuant.omega0;
};  //End Form denisuty related  LYP Corr


DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB){
};

DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB){
};

DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB, bool gammasmall){
  DFTFunctional::DFTInfo info;
  denspow RhoTQuant;
  this->popDensPow(rhoA, rhoB, RhoTQuant);
  this->popLYPdens(rhoA, rhoB, RhoTQuant);
  info.ddgammaAA  = RhoTQuant.dLYPdgAA ;
//  Eq. A23* dLYP/dgammaBB (debugged)
  info.ddgammaBB  = RhoTQuant.dLYPdgBB ;
//  Eq. A24  dLYP/dgammaAB (debugged)
  info.ddgammaAB  = RhoTQuant.dLYPdgAB ;
//  Eq. A28  dLYP/dRhoA (debugged)
    info.ddrhoA   = - 4.0 *  this->a * rhoA * rhoB * RhoTQuant.zeta
                    / RhoTQuant.rhoT;
    info.ddrhoA  *= ( (1.0/rhoA)
                    - (1.0/RhoTQuant.rhoT) 
                    + (this->d * RhoTQuant.zeta /(3.0* RhoTQuant.rho4over3) )
                    );
    info.ddrhoA  += - this->Cfact * this->a * this->b 
             *( 
              (RhoTQuant.omega1 * rhoA * rhoB 
               * (RhoTQuant.rhoA8over3 + RhoTQuant.rhoB8over3))
              +(RhoTQuant.omega0 * rhoB * ( (11.0*RhoTQuant.rhoA8over3/3.0) 
               + RhoTQuant.rhoB8over3))
             );
    if(gammaAA > this->small) info.ddrhoA  += gammaAA* RhoTQuant.d2LYPdrhgAA;
    if(gammaAB > this->small) info.ddrhoA  += gammaAB* RhoTQuant.d2LYPdrhgAB;
    if(gammaBB > this->small) info.ddrhoA  += gammaBB* RhoTQuant.d2LYPdrhgBB;
//  Eq. A28* dLYP/dRhoB (debugged)
    info.ddrhoB   = - 4.0 *  this->a * rhoA * rhoB 
      / ( ( RhoTQuant.rhoT)*(1.0 +  this->d/RhoTQuant.rho1over3 ) );
    info.ddrhoB  *= ( (1.0/rhoB)
                  -(1.0/RhoTQuant.rhoT) 
                  +((this->d/3.0) 
                  / RhoTQuant.rho4over3 / (1.0 + this->d 
                  / RhoTQuant.rho1over3  ) )
                 );
    info.ddrhoB  += - this->Cfact * this->a * this->b 
             *( 
              (RhoTQuant.omega1 * rhoA * rhoB 
               * (RhoTQuant.rhoA8over3 + RhoTQuant.rhoB8over3))
              +(RhoTQuant.omega0 * rhoA * ( (11.0*RhoTQuant.rhoB8over3/3.0) 
               + RhoTQuant.rhoA8over3))
             );
    if (gammaBB > this->small) info.ddrhoB  += gammaBB* RhoTQuant.d2LYPdrhgAA;
    if (gammaAB > this->small) info.ddrhoB  += gammaAB* RhoTQuant.d2LYPdrhgAB;
    if (gammaAA > this->small) info.ddrhoB  += gammaAA* RhoTQuant.d2LYPdrhgBB;
    info.eps  = - 4.0 * this->a * rhoA * rhoB 
      / ( RhoTQuant.rhoT*(1.0 + this->d / RhoTQuant.rho1over3 ) );
    info.eps += - this->Cfact * this->a * this->b 
      * RhoTQuant.omega0 * rhoA * rhoB * (RhoTQuant.rhoA8over3 + RhoTQuant.rhoB8over3); 
    if(gammaAA > this->small) info.eps += info.ddgammaAA * gammaAA;
    if(gammaAB > this->small) info.eps += info.ddgammaBB * gammaBB;
    if(gammaBB > this->small) info.eps += info.ddgammaAB * gammaAB;
};


DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
  DFTFunctional::DFTInfo info;
  denspow RhoTQuant;
  this->popDensPow(rhoA, rhoB, RhoTQuant);
  this->popLYPdens(rhoA, rhoB, RhoTQuant);
  if( gammaAA <  this-> small ||  gammaBB <  this-> small ||  gammaAB <  this-> small) {
    this->eval(rhoA,rhoB,gammaAA,gammaAB,gammaBB,true);
  } else {
//  Eq. A23  dLYP/dgammaAA (debugged)
  info.ddgammaAA  = RhoTQuant.dLYPdgAA ;
//  Eq. A23* dLYP/dgammaBB (debugged)
  info.ddgammaBB  = RhoTQuant.dLYPdgBB ;
//  Eq. A24  dLYP/dgammaAB (debugged)
  info.ddgammaAB  = RhoTQuant.dLYPdgAB ;
//  Eq. A28  dLYP/dRhoA (debugged)
    info.ddrhoA   = - 4.0 *  this->a * rhoA * rhoB * RhoTQuant.zeta
                    / RhoTQuant.rhoT;
    info.ddrhoA  *= ( (1.0/rhoA)
                    - (1.0/RhoTQuant.rhoT) 
                    + (this->d * RhoTQuant.zeta /(3.0* RhoTQuant.rho4over3) )
                    );
    info.ddrhoA  += - this->Cfact * this->a * this->b 
             *( 
              (RhoTQuant.omega1 * rhoA * rhoB 
               * (RhoTQuant.rhoA8over3 + RhoTQuant.rhoB8over3))
              +(RhoTQuant.omega0 * rhoB * ( (11.0*RhoTQuant.rhoA8over3/3.0) 
               + RhoTQuant.rhoB8over3))
             );
    info.ddrhoA  += gammaAA* RhoTQuant.d2LYPdrhgAA;
    info.ddrhoA  += gammaAB* RhoTQuant.d2LYPdrhgAB;
    info.ddrhoA  += gammaBB* RhoTQuant.d2LYPdrhgBB;
//  Eq. A28* dLYP/dRhoB (debugged)
    info.ddrhoB   = - 4.0 *  this->a * rhoA * rhoB * RhoTQuant.zeta  
      / RhoTQuant.rhoT;
    info.ddrhoB  *= ( (1.0/rhoB)
                    -(1.0/RhoTQuant.rhoT) 
                    + (this->d * RhoTQuant.zeta /(3.0* RhoTQuant.rho4over3) )
                    );
    info.ddrhoB  += - this->Cfact * this->a * this->b 
             *( 
              (RhoTQuant.omega1 * rhoA * rhoB 
               * (RhoTQuant.rhoA8over3 + RhoTQuant.rhoB8over3))
              +(RhoTQuant.omega0 * rhoA * ( (11.0*RhoTQuant.rhoB8over3/3.0) 
               + RhoTQuant.rhoA8over3))
             );
    info.ddrhoB  += gammaBB* RhoTQuant.d2LYPdrhgAA;
    info.ddrhoB  += gammaAB* RhoTQuant.d2LYPdrhgAB;
    info.ddrhoB  += gammaAA* RhoTQuant.d2LYPdrhgBB;
    info.eps  = - 4.0 * this->a * rhoA * rhoB 
      / ( RhoTQuant.rhoT*(1.0 + this->d / RhoTQuant.rho1over3 ) );
    info.eps += - this->Cfact * this->a * this->b 
      * RhoTQuant.omega0 * rhoA * rhoB * (RhoTQuant.rhoA8over3 + RhoTQuant.rhoB8over3); 
    info.eps += info.ddgammaAA * gammaAA;
    info.eps += info.ddgammaBB * gammaBB;
    info.eps += info.ddgammaAB * gammaAB;
    }
    return info;
};

