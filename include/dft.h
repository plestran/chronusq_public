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
  double small = std::numeric_limits<double>::epsilon(); //numerical accuracy
#ifdef CQ_ENABLE_LIBXC
  xc_func_type func;
#endif

  std::string name;

  DFTFunctional(double X = 1.0, double eps = 1e-10){
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
    inline void operator*=(double x){
      this->eps       *= x;
      this->ddrhoA    *= x;
      this->ddrhoB    *= x;
      this->ddgammaAA *= x;
      this->ddgammaAB *= x;
      this->ddgammaBB *= x;
    }
  };

  virtual DFTInfo eval(const double &rA,const double &rB) = 0;
  virtual DFTInfo eval(const double &rA,const double &rB,const double &gammaAA, const double &gammaAB) = 0;
  virtual DFTInfo eval(const double &rA,const double &rB,const double &gammaAA, const double &gammaAB, const double &gammaBB) = 0;
};

class SlaterExchange : public DFTFunctional {
  double CxVx;  //TF LDA Prefactor (for Vx)  
//A  double small;    
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
//A  double small;    
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
  double b_p  ;
  double b_f  ;
  double c_p  ;
  double c_f  ;
  double x0_p ;
  double x0_f ;
  double x0_a ;
  double X_x0p ; 
  double X_x0f ; 
  double Qp       ; 
  double Qf       ; 
// Density related variable
  struct denspow {
    double rhoT;          // rhoA + rhoB 
    double spindensity;   // (rhoA - rhoB)/rhoT
    double f0_spindensity ;
    double df_spindensity ;
    double r_s      ;
    double r_s_sqrt ;
    double r_s_32   ;
    double Xp       ; 
    double Xf       ; 
    double eps_p    ;
    double eps_f    ;
    double delta_eps_1 ;
    double S1          ;
    double S2          ;
    double M3_A        ;
    double M3_B        ;
    double spindensity_3;
    double spindensity_4;
    double df2_spindensity;
//  VWN5 only
    double Xa;
    double alpha;
    double beta;
    double delta_eps_etha;
    double db_dr;
    double S3          ;
    double M1          ;
    double S4          ;
    double S5          ;
    };

  VWNIII(double X = 1.0, double eps = 1e-10);
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
  double Eveps0VWN(double &A_x, double &b_x, double &Q, double &X, double &x0_x, double &X_x0, denspow &denquant);
  double Eveps1VWN(double &A_x, double &b1, double &b2, double &b3, denspow &denquant);
  double Eveps2VWN(double A_x, double &b_x, double &c_x, double &X, double &x0_x, denspow &denquant);
  void popVWNconst();
  void popVWNdens(const double &rhoA, const double &rhoB, denspow &denquant);
};

class VWNV : public VWNIII {
public:
  VWNV(double X = 1.0, double eps = 1e-10);
//VWN5
// Const
  double b1p      ; 
  double b2p      ; 
  double b3p      ; 
  double b1f      ; 
  double b2f      ; 
  double b3f      ; 
  double A_a  ;
  double b_a  ;
  double c_a  ;
  double Qa;
  double X_x0a;
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
  void popVWNconst();
  void popVWNdens(const double &rhoA, const double &rhoB, denspow &denquant);
};

class BEightEight : public DFTFunctional {
  double  d1over3;
  double  d4over3;
  double CxVx;  //TF LDA Prefactor (for Vx)  
  double beta;    
//A  double small;    
public:
  BEightEight(double X = 1.0, double eps = 1e-10);
  double g0B88 (double x, double &sinhx, double &bx);
  double g1B88 (double x, double &sinhx, double &bx);
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
};

class PBE : public DFTFunctional {
  double  d1over3;
  double  d4over3;
  double CxVx;  //TF LDA Prefactor (for Vx)  
  double beta;    
//A  double small;    
public:
  PBE(double X = 1.0, double eps = 1e-10);
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
};

class lyp : public DFTFunctional {
  double d1over3;
  double Cfact;  //TF LDA Prefactor (for Vx)  
//A  double small;    
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
    double omegainf;      // omega1/omega0 (to avoid problem)   
    double dLYPdgAA;       // Eq A23(delLYP/delgammaAA) ;
    double dLYPdgBB   ;    // Eq A23(delLYP/delgammaBB) 
    double dLYPdgAB   ;    // Eq A24(delLYP/delgammaAB) 
    double d2LYPdrhoAdgAA ;   // Eq A29(del2LYP/delgammaAA*delrhoA)
    double d2LYPdrhoAdgAB ;   // Eq A30(del2LYP/delgammaAB*delrhoA)
    double d2LYPdrhoAdgBB ;   // Eq A31(del2LYP/delgammaBB*delrhoA)
    double d2LYPdrhoBdgBB ;   // Eq A29(del2LYP/delgammaBB*delrhoA)
    double d2LYPdrhoBdgAB ;   // Eq A30(del2LYP/delgammaAB*delrhoB)
    double d2LYPdrhoBdgAA ;   // Eq A31(del2LYP/delgammaAA*delrhoB)
    };
  lyp(double X = 1.0, double eps = 1e-10);
  void popLYPdens(const double &rhoA, const double &rhoB, denspow &denquant);
  void popDensPow(const double &rhoA, const double &rhoB, denspow &denquant);
  double d2LYPdrhoXdgXX(const double &rhoX, const double &rhoY,const double &dLYPdgXX, denspow &denquant);   
  double d2LYPdrhoXdgXY(const double &rhoX, const double &rhoY,const double &dLYPdgXY, denspow &denquant);
  double d2LYPdrhoXdgYY(const double &rhoX, const double &rhoY,const double &dLYPdgYY, denspow &denquant);
  DFTInfo eval(const double &rhoA, const double &rhoB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB);
  DFTInfo eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB, bool small, DFTInfo &info);
};
#endif
