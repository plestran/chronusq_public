#include<dft.h>

lyp::lyp(double X, double eps):
DFTFunctional(X,eps){
// Memo Factor to be added at the end for numerical stability
  this-> Cfact = 36.4623989787648;  //-(2.0^(11.0/3.0))*(3.0/10.0)*((3.0*pi*pi)^(2.0/3.0))
//  this-> small = 1.0e-18; 
  this-> d1over3  = 1.0/3.0;
  this-> a = 0.04918;
  this-> b = 0.132;
  this-> c = 0.2533;
  this-> d = 0.349;

  this->name = "LYP";
#ifdef CQ_ENABLE_LIBXC
  xc_func_init(&this->func,XC_GGA_C_LYP,XC_POLARIZED);
#endif
};


    void lyp::popDensPow(const double &rhoA, const double &rhoB, denspow &denquant)  {
    double tmp = 8.0/3.0;
    denquant.rhoT        = rhoA + rhoB;
//    denquant.spindensity = (rhoA - rhoB) / denquant.rhoT;
    denquant.spindensity = (rhoA - rhoB) ;
    denquant.rho1over3   = std::pow(denquant.rhoT,this->d1over3);
//    if(denquant.rho1over3 < this->small) {
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

double lyp::d2LYPdrhoXdgXX(const double &rhoX, const double &rhoY,const double &dLYPdgXX, denspow &denquant){
//     Eq A29 (del^2 LYP / (del rho_X del gammaXX))
       double tmp;
       tmp =   - (rhoX*rhoY/9.0) 
             *(  ( 3.0 + rhoX/denquant.rhoT ) * 
              denquant.delta1 + rhoY*( denquant.delta0 - 11.0)/
              (denquant.rho2)  ) 
           + (rhoY/9.0)
             *( 1.0 - 3.0*denquant.delta0 -rhoX*(denquant.delta0 - 11)/
             denquant.rhoT );
       tmp *= -this->a *  this->b * denquant.omega0;
       tmp += dLYPdgXX*denquant.omegainf;
       return tmp;
};

double lyp::d2LYPdrhoXdgXY(const double &rhoX, const double &rhoY,const double &dLYPdgXY, denspow &denquant){
//     Eq A30 (del^2 LYP / (del rho_X del gammaXY))
       double tmp;
       tmp   =  - 8.0*denquant.rhoT/3.0;
       tmp  +=  -  denquant.delta1*(7.0*rhoX*rhoY/9.0);
       tmp  +=  rhoY*(47.0-7.0*denquant.delta0)/9.0;
       tmp  *= -this->a *  this->b * denquant.omega0;
       tmp  += dLYPdgXY*denquant.omegainf;
       return tmp;
};

double lyp::d2LYPdrhoXdgYY(const double &rhoX, const double &rhoY,const double &dLYPdgYY, denspow &denquant){
//     Eq A31 (del^2 LYP / (del rho_X del gammaYY)) (debugged alread)
       double tmp;
       tmp =   - (rhoX*rhoY/9.0) 
             *(  ( 3.0 + rhoY/ denquant.rhoT ) *  denquant.delta1 - rhoY*
           (denquant.delta0 - 11.0)/( denquant.rho2)  ) 
           + (rhoY/9.0)
           * ( 1.0 - 3.0* denquant.delta0 
           - rhoY*(denquant.delta0 - 11)/ denquant.rhoT )
           - 2.0 * rhoX;
       tmp *= -this->a *  this->b * denquant.omega0;
       tmp += dLYPdgYY*denquant.omegainf;
       return tmp;
};

void lyp::popLYPdens(const double &rhoA, const double &rhoB, denspow &RhoTQuant){
//function used in LYP Correlation (Eq. A26, A32 Ref.LYP2)
  if(RhoTQuant.rhoT < this->small) {
    RhoTQuant.omega0 = 0.0;
    RhoTQuant.omega1 = 0.0;
    RhoTQuant.omegainf  = this->c ;
    RhoTQuant.omegainf += this->d*RhoTQuant.zeta;
    RhoTQuant.omegainf -= 11.0*RhoTQuant.rho1over3;
    RhoTQuant.omegainf /= 3.0 *RhoTQuant.rho4over3;
  }else{
    RhoTQuant.omega0  =  std::exp(-this->c / RhoTQuant.rho1over3);
    RhoTQuant.omega0 *=  RhoTQuant.zeta;
    RhoTQuant.omega0 /=  RhoTQuant.rho2*RhoTQuant.rho4over3*RhoTQuant.rho1over3; //^(11/3)
    RhoTQuant.omega1  = this->c ;
    RhoTQuant.omega1 += this->d*RhoTQuant.zeta;
    RhoTQuant.omega1 -= 11.0*RhoTQuant.rho1over3;
    RhoTQuant.omega1 /= 3.0 *RhoTQuant.rho4over3;
    RhoTQuant.omegainf = RhoTQuant.omega1;
    RhoTQuant.omega1  *= RhoTQuant.omega0;
    
  }
  if(RhoTQuant.rhoT < this->small) {
    RhoTQuant.delta0 = 0.0;
    RhoTQuant.delta1 = 0.0;
  //function used in LYP Correlation (Eq. A27, A33 Ref.LYP2)
  }else{ 
    RhoTQuant.delta0   = this->c + this->d*RhoTQuant.zeta;
    RhoTQuant.delta0  /= RhoTQuant.rho1over3;
    RhoTQuant.delta1 = (RhoTQuant.zeta* RhoTQuant.zeta* this->d* this->d/ (RhoTQuant.rho1over3*RhoTQuant.rho1over3) ); 
    RhoTQuant.delta1 -= RhoTQuant.delta0;
    RhoTQuant.delta1 /= (3.0*RhoTQuant.rhoT);
  }
  // Eq A23(delLYP/delgammaAA) 
  RhoTQuant.dLYPdgAA   = 11.0 - RhoTQuant.delta0;
  RhoTQuant.dLYPdgAA  *= rhoA/RhoTQuant.rhoT;
  RhoTQuant.dLYPdgAA  += 1.0 - 3.0*RhoTQuant.delta0;
  RhoTQuant.dLYPdgAA  *= rhoA*rhoB/9.0;
  RhoTQuant.dLYPdgAA  -= rhoB*rhoB;
  RhoTQuant.dLYPdgAA  *= -this->a * this->b * RhoTQuant.omega0;
//Eq A29 (del^2 LYP / (del rho_X del gammaXX))
  RhoTQuant.d2LYPdrhoAdgAA = d2LYPdrhoXdgXX(rhoA,rhoB,RhoTQuant.dLYPdgAA,RhoTQuant);
/*
  RhoTQuant.d2LYPdrhoAdgAA =   - (rhoA*rhoB/9.0) 
        *(  ( 3.0 + rhoA/RhoTQuant.rhoT ) * 
         RhoTQuant.delta1 + rhoB*( RhoTQuant.delta0 - 11.0)/
         (RhoTQuant.rho2)  ) 
      + (rhoB/9.0)
        *( 1.0 - 3.0*RhoTQuant.delta0 -rhoA*(RhoTQuant.delta0 - 11)/
        RhoTQuant.rhoT );
  RhoTQuant.d2LYPdrhoAdgAA *= -this->a *  this->b * RhoTQuant.omega0;
  RhoTQuant.d2LYPdrhoAdgAA += RhoTQuant.dLYPdgAA*RhoTQuant.omegainf;
*/
  // Eq A25(delLYP/delgammaBB) 
  RhoTQuant.dLYPdgBB   = 11.0 - RhoTQuant.delta0;
  RhoTQuant.dLYPdgBB  *= rhoB/RhoTQuant.rhoT;
  RhoTQuant.dLYPdgBB  += 1.0 - 3.0*RhoTQuant.delta0;
  RhoTQuant.dLYPdgBB  *= rhoB*rhoA/9.0;
  RhoTQuant.dLYPdgBB  -= rhoA*rhoA;
  RhoTQuant.dLYPdgBB  *= -this->a * this->b * RhoTQuant.omega0;
 //UKS 
//Eq A29 (del^2 LYP / (del rho_X del gammaXX))
  RhoTQuant.d2LYPdrhoBdgBB = d2LYPdrhoXdgXX(rhoB,rhoA,RhoTQuant.dLYPdgBB,RhoTQuant);
//Eq A24 (delLYP/delgammaAB)
  RhoTQuant.dLYPdgAB  = 47.0 - 7.0 *  RhoTQuant.delta0;
  RhoTQuant.dLYPdgAB *= rhoA * rhoB / 9.0;
  RhoTQuant.dLYPdgAB -= 4.0 *  RhoTQuant.rho2/ 3.0;
  RhoTQuant.dLYPdgAB *= - this->a *  this->b * RhoTQuant.omega0;
//Eq A30 (del^2 LYP / (del rho_X del gammaXY))
  RhoTQuant.d2LYPdrhoAdgAB = d2LYPdrhoXdgXY(rhoA,rhoB,RhoTQuant.dLYPdgAB,RhoTQuant);
//UKS
  RhoTQuant.d2LYPdrhoBdgAB = d2LYPdrhoXdgXY(rhoB,rhoA,RhoTQuant.dLYPdgAB,RhoTQuant);
/*
  RhoTQuant.d2LYPdrhoAdgAB  =  - 8.0*RhoTQuant.rhoT/3.0;
  RhoTQuant.d2LYPdrhoAdgAB +=  -  RhoTQuant.delta1*(7.0*rhoA*rhoB/9.0);
  RhoTQuant.d2LYPdrhoAdgAB +=  rhoB*(47.0-7.0*RhoTQuant.delta0)/9.0;
  RhoTQuant.d2LYPdrhoAdgAB *= -this->a *  this->b * RhoTQuant.omega0;
  RhoTQuant.d2LYPdrhoAdgAB += RhoTQuant.dLYPdgAB*RhoTQuant.omegainf;
*/
//Eq A31 (del^2 LYP / (del rho_X del gammaYY)) (debugged alread)
  RhoTQuant.d2LYPdrhoAdgBB = d2LYPdrhoXdgYY(rhoA,rhoB,RhoTQuant.dLYPdgBB,RhoTQuant); 
//UKS
  RhoTQuant.d2LYPdrhoBdgAA = d2LYPdrhoXdgYY(rhoB,rhoA,RhoTQuant.dLYPdgAA,RhoTQuant); 
/*
  RhoTQuant.d2LYPdrhoAdgBB =   - (rhoA*rhoB/9.0) 
        *(  ( 3.0 + rhoB/ RhoTQuant.rhoT ) *  RhoTQuant.delta1 - rhoB*
      (RhoTQuant.delta0 - 11.0)/( RhoTQuant.rho2)  ) 
      + (rhoB/9.0)
      * ( 1.0 - 3.0* RhoTQuant.delta0 
      - rhoB*(RhoTQuant.delta0 - 11)/ RhoTQuant.rhoT )
      - 2.0 * rhoA;
*/
/////  RhoTQuant.d2LYPdrhoAdgBB *= -this->a *  this->b * RhoTQuant.omega0;
//  RhoTQuant.d2LYPdrhoAdgBB += RhoTQuant.dLYPdgBB*RhoTQuant.omega1/RhoTQuant.omega0;
/////  RhoTQuant.d2LYPdrhoAdgBB += RhoTQuant.dLYPdgBB*RhoTQuant.omegainf;
///NEW      }
};  //End Form denisuty related  LYP Corr


DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB){
};

DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB){
};

DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB, bool gammasmall, DFTInfo &info){
  cout << "SMALL" <<endl;
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
    if(gammaAA > this->small) {
      info.ddrhoA  += gammaAA* RhoTQuant.d2LYPdrhoAdgAA;
      }
    if(gammaAB > this->small) {
      info.ddrhoA  += gammaAB* RhoTQuant.d2LYPdrhoAdgAB;
      }
    if(gammaBB > this->small) {
      info.ddrhoA  += gammaBB* RhoTQuant.d2LYPdrhoAdgBB;
      }
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
    if (gammaBB > this->small) info.ddrhoB  += gammaBB* RhoTQuant.d2LYPdrhoBdgBB;
    if (gammaAB > this->small) info.ddrhoB  += gammaAB* RhoTQuant.d2LYPdrhoBdgAB;
    if (gammaAA > this->small) info.ddrhoB  += gammaAA* RhoTQuant.d2LYPdrhoBdgAA;
    info.eps  = - 4.0 * this->a * rhoA * rhoB 
      / ( RhoTQuant.rhoT*(1.0 + this->d / RhoTQuant.rho1over3 ) );
    info.eps += - this->Cfact * this->a * this->b 
      * RhoTQuant.omega0 * rhoA * rhoB * (RhoTQuant.rhoA8over3 + RhoTQuant.rhoB8over3); 
    if(gammaAA > this->small) info.eps += info.ddgammaAA * gammaAA;
    if(gammaBB > this->small) info.eps += info.ddgammaBB * gammaBB;
    if(gammaAB > this->small) info.eps += info.ddgammaAB * gammaAB;
};


DFTFunctional::DFTInfo lyp::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
  DFTFunctional::DFTInfo info;
//#ifndef CQ_ENABLE_LIBXC
  denspow RhoTQuant;
  this->popDensPow(rhoA, rhoB, RhoTQuant);
  this->popLYPdens(rhoA, rhoB, RhoTQuant);

  if( gammaAA <  this-> small ||  gammaBB <  this-> small ||  gammaAB <  this-> small) {
    this->eval(rhoA,rhoB,gammaAA,gammaAB,gammaBB,true, info);
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
    info.ddrhoA  += gammaAA* RhoTQuant.d2LYPdrhoAdgAA;
    info.ddrhoA  += gammaAB* RhoTQuant.d2LYPdrhoAdgAB;
    info.ddrhoA  += gammaBB* RhoTQuant.d2LYPdrhoAdgBB;
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
    info.ddrhoB  += gammaBB* RhoTQuant.d2LYPdrhoBdgBB;
    info.ddrhoB  += gammaAB* RhoTQuant.d2LYPdrhoBdgAB;
    info.ddrhoB  += gammaAA* RhoTQuant.d2LYPdrhoBdgAA;
    info.eps  = - 4.0 * this->a * rhoA * rhoB 
      / ( RhoTQuant.rhoT*(1.0 + this->d / RhoTQuant.rho1over3 ) );
    info.eps += - this->Cfact * this->a * this->b 
      * RhoTQuant.omega0 * rhoA * rhoB * (RhoTQuant.rhoA8over3 + RhoTQuant.rhoB8over3); 
    info.eps += info.ddgammaAA * gammaAA;
    info.eps += info.ddgammaBB * gammaBB;
    info.eps += info.ddgammaAB * gammaAB;
    }


// debug

  DFTFunctional::DFTInfo infoZ;
  std::array<double,2> rhoZ = {rhoA,rhoB};
  std::array<double,3> sigmaZ = {gammaAA,gammaAB,gammaBB};
  std::array<double,2> vrhoZ;
  std::array<double,3> vsigmaZ;
  xc_gga_exc_vxc(&this->func,1,&rhoZ[0],&sigmaZ[0],&infoZ.eps,&vrhoZ[0],
    &vsigmaZ[0]);

  infoZ.ddrhoA = vrhoZ[0];
  infoZ.ddrhoB = vrhoZ[1];
  infoZ.ddgammaAA = vsigmaZ[0];
  infoZ.ddgammaAB = vsigmaZ[1];
  infoZ.ddgammaBB = vsigmaZ[2];

  infoZ.eps *= (rhoA + rhoB);

  cout << " eps "<<info.eps       -   infoZ.eps << endl;
  cout << " ddrhoA "<<info.ddrhoA -   infoZ.ddrhoA << endl;
  cout << " ddrhoB "<<info.ddrhoB -   infoZ.ddrhoB << endl;
  cout << " ddgammaAA "<<info.ddgammaAA -  infoZ.ddgammaAA << endl;
  cout << " ddgammaBB "<<info.ddgammaBB -  infoZ.ddgammaBB << endl;
  cout << " ddgammaAB "<<info.ddgammaAB -  infoZ.ddgammaAB << endl;

//#else
/*
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

  info.eps *= (rhoA + rhoB);
*/
//#endif
    return infoZ;
};

