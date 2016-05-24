#include<global.h>
#ifndef INCLUDED_DFT
#define INCLUDED_DFT
class DFTFunctional{
public:
  double scalingFactor;    //< Hybrid Scaling
  double epsScreen;        //< screening 
  DFTFunctional(double X = 1.0){
    this->scalingFactor = X;
    this->epsScreen = 1.0e-10;
  };

  struct DFTInfo {
    double eps;
    double ddrhoA;
    double ddrhoB;
    double ddgammaAA;
    double ddgammaAB;
    double ddgammaBB;

    DFTInfo(){
      this->eps = 0;
      this->ddrhoA = 0;
      this->ddrhoB = 0;
      this->ddgammaAA = 0;
      this->ddgammaAB = 0;
      this->ddgammaBB = 0;
    };

  };

  virtual DFTInfo eval(double rA,double rB) = 0;
  virtual DFTInfo eval(double rA,double rB,double gammaAA, double gammaAB) = 0;
  virtual DFTInfo eval(double rA,double rB,double gammaAA, double gammaAB, double gammaBB) = 0;
};

class SlaterExchange : public DFTFunctional {
  double CxVx;  //TF LDA Prefactor (for Vx)  
  double small;    
  double d1over3 ;
  double d4over3 ;
  double rhoT     ;    
  double spindensity   ;
public:
  SlaterExchange();
  DFTInfo eval(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB);
};

class VWNIII : public DFTFunctional {
// General Constant
public:
  double small;    
  double over2;
  double over3;
  double over4;
  double over6;
  double fourover3;
// Functional Constant
// Parameter for the fit (according Eq 4.4 Vosko Can. J. Phys. 1980
  double A_p  ; 
  double A1  ; 
  double A_f  ; 
  double A_a  ;
  double b_p  ;
  double b_f  ;
  double b_a  ;
  double c_p  ;
  double c_f  ;
  double c_a  ;
  double x0_p ;
  double x0_f ;
  double x0_a ;
  double X_x0p ; 
  double X_x0f ; 
// Density related constant
  double Q        ; 
  double r_s      ;
  double r_s_sqrt ;
  double r_s_32   ;
  double Xp       ; 
  double Xf       ; 
  double rhoT     ;    
  double spindensity   ;
  double spindensity_4 ;
  double spindensity_3 ;
  double f0_spindensity ;
  double f1_spindensity ;
  double df_spindensity ;
//Interm Quantities
  double b1p      ; 
  double b2p      ; 
  double b3p      ; 
  double Qp       ; 
  double b1f      ; 
  double b2f      ; 
  double b3f      ; 
  double Qf       ; 
  double tmp1 ;
  double tmp2 ;
  double tmp3 ;
  double eps_p       ;
  double eps_f       ;
  double delta_eps_1 ;
  double S1          ;
  double S2          ;
  double M3_A        ;
  double M3_B        ;
  double alpha       ;
  double mu_p        ;
  double mu_f        ;
  double beta        ;
  double S3          ;
  double S4          ;
  double S5          ;
  double M1          ;
  double db_dr       ; 
  double delta_eps_etha ;
  VWNIII();
  DFTInfo eval(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB);
  double Eveps0VWN(double &A_x, double &b_x, double &Q, double &X, double &x0_x, double &X_x0);
  double Eveps1VWN(double &A_x, double &b1, double &b2, double &b3);
  double Eveps2VWN(double A_x, double &b_x, double &c_x, double &X, double &x0_x);
  void popVWNconst();
  void popVWNdens(double rhoA, double rhoB);
};

class VWNV : public VWNIII {
public:
  VWNV();
  double alpha;
  double Xa;
  double Qa;
  double X_x0a;
  double beta;
  double df2_spindensity;
  double delta_eps_etha;
  double db_dr;
  DFTInfo eval(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB);
  double Eveps0VWN(double &A_x, double &b_x, double &Q, double &X, double &x0_x, double &X_x0);
  double Eveps1VWN(double &A_x, double &b1, double &b2, double &b3);
  double Eveps2VWN(double A_x, double &b_x, double &c_x, double &X, double &x0_x);
  void popVWNconst();
  void popVWNdens(double rhoA, double rhoB);
};

class BEightEight : public DFTFunctional {
  double  d1over3;
  double  d4over3;
  double CxVx;  //TF LDA Prefactor (for Vx)  
  double beta;    
  double rhoT     ;    
  double spindensity   ;
  double small;    
  double xA;
  double xB;
  double rhoA1ov3 ;
  double rhoA4ov3 ;
  double rhoB1ov3 ;
  double rhoB4ov3 ;
public:
  BEightEight();
  double g0B88 (double x);
  double g1B88 (double x);
  DFTInfo eval(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaBB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB);
};

class lyp : public DFTFunctional {
  double d1over3;
  double Cfact;  //TF LDA Prefactor (for Vx)  
  double small;    
  double a     ;    
  double b     ;    
  double c     ;    
  double d     ;    
  double omega1     ;    
  double omega0     ;    
  double delta1     ;    
  double delta0     ;    
  double rhoT     ;    
  double spindensity   ;
  double rho1over3 ;
  double rhoA8over3 ;
  double rhoB8over3 ;
  double dLYPdgAA   ;
  double dLYPdgBB   ;
  double dLYPdgAB   ;
  double d2LYPdrhgAA   ;
  double d2LYPdrhgAB   ;
  double d2LYPdrhgBB   ;
public:
  lyp();
  void popLYPdens(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaBB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB);
};
#endif
