#include <global.h>
#ifdef CQ_ENABLE_LIBXC
#include <xc.h>
#endif
#ifndef INCLUDED_DFT
#define INCLUDED_DFT
class DFTFunctional{
public:
  double scalingFactor;    //< Hybrid Scaling
  double epsScreen;        //< screening 

#ifdef CQ_ENABLE_LIBXC
  xc_func_type func;
#endif

  std::string name;

  DFTFunctional(double X = 1.0, double eps = 1e-10){
#ifdef CQ_ENABLE_LIBXC
  xc_func_init(&this->func,XC_LDA_X,XC_POLARIZED);
#endif
    this->scalingFactor = X;
    this->epsScreen = eps;
  };

  ~DFTFunctional(){
#ifdef CQ_ENABLE_LIBXC
     xc_func_end(&func);
#endif
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
    inline void operator+=(DFTInfo other){
      this->eps       += other.eps      ;
      this->ddrhoA    += other.ddrhoA   ;
      this->ddrhoB    += other.ddrhoB   ;
      this->ddgammaAA += other.ddgammaAA;
      this->ddgammaAB += other.ddgammaAB;
      this->ddgammaBB += other.ddgammaBB;
    }

  };

  virtual DFTInfo eval(const double &rA,const double &rB) = 0;
  virtual DFTInfo eval(const double &rA,const double &rB,const double &gammaAA, const double &gammaAB) = 0;
  virtual DFTInfo eval(const double &rA,const double &rB,const double &gammaAA, const double &gammaAB, const double &gammaBB) = 0;
};

class SlaterExchange : public DFTFunctional {
  double CxVx;  //TF LDA Prefactor (for Vx)  
  double small;    
  double d1over3 ;
  double d4over3 ;
public:
  SlaterExchange(double X = 1.0, double eps = 1e-10);
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
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
  VWNIII(double X = 1.0, double eps = 1e-10);
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
  double Eveps0VWN(double &A_x, double &b_x, double &Q, double &X, double &x0_x, double &X_x0);
  double Eveps1VWN(double &A_x, double &b1, double &b2, double &b3);
  double Eveps2VWN(double A_x, double &b_x, double &c_x, double &X, double &x0_x);
  void popVWNconst();
  void popVWNdens(double rhoA, double rhoB);
};

class VWNV : public VWNIII {
public:
  VWNV(double X = 1.0, double eps = 1e-10);
  double alpha;
  double Xa;
  double Qa;
  double X_x0a;
  double beta;
  double df2_spindensity;
  double delta_eps_etha;
  double db_dr;
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
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
  double small;    
public:
  BEightEight(double X = 1.0, double eps = 1e-10);
  double g0B88 (double x, double &sinhx, double &bx);
  double g1B88 (double x, double &sinhx, double &bx);
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
};

class lyp : public DFTFunctional {
  double d1over3;
  double Cfact;  //TF LDA Prefactor (for Vx)  
  double small;    
  double a     ;    
  double b     ;    
  double c     ;    
  double d     ;    
public:
  struct denspow {
    double rhoT;          // rhoA + rhoB 
    double spindensity;   // (rhoA - rhoB)/rhoT
    double rho2;          // ^(2) 
    double rho1over3;     // ^(1/3) 
    double rho4over3;    // ^(4/3) 
    double rhoA8over3;   // ^(8/3) alpha  
    double rhoB8over3;   // ^(8/3) beta   
    double zeta ;        // (1+ d / rho1over3)^(-1)   Aux variable
    double delta0;        // (c + d*zeta)/ (rho1over3) (A27) /rearrenged
    double delta1;        // ddelta0/d rho   (A33)
    double omega0;        // Eq A26   
    double omega1;        // domega/d rho (A32)   
    double dLYPdgAA;       // Eq A23(delLYP/delgammaAA) ;
    double dLYPdgBB   ;    // Eq A23(delLYP/delgammaBB) 
    double dLYPdgAB   ;    // Eq A24(delLYP/delgammaAB) 
    double d2LYPdrhgAA ;   // Eq A29(del2LYP/delgammaAA*delrhoA)
    double d2LYPdrhgAB ;   // Eq A30(del2LYP/delgammaAB*delrhoA)
    double d2LYPdrhgBB ;   // Eq A31(del2LYP/delgammaBB*delrhoA)
    };
  lyp(double X = 1.0, double eps = 1e-10);
  void popLYPdens(const double &rhoA, const double &rhoB, denspow &denquant);
  void popDensPow(const double &rhoA, const double &rhoB, denspow &denquant);  
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB, bool small);
};
#endif
