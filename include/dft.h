#include<global.h>
#ifndef INCLUDED_DFT
#define INCLUDED_DFT
class DFTFunctional{
public:
  double scalingFactor;
  DFTFunctional(double X = 1.0){
    this->scalingFactor = X;
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

class VWN3 : public DFTFunctional {
// General Constant
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
public:
  VWN3();
  DFTInfo eval(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB);
  double Eveps0VWN(double &A_x, double &b_x, double &Q, double &X, double &x0_x, double &X_x0);
  double Eveps1VWN(double &A_x, double &b1, double &b2, double &b3);
  double Eveps2VWN(double A_x, double &b_x, double &c_x, double &X, double &x0_x);
  void popVWNconst();
  void popVWNdens(double rhoA, double rhoB);
};
#endif
