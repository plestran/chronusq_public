#include<global.h>
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
public:
  SlaterExchange();
  DFTInfo eval(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB);
};

class VWN3 : public DFTFunctional {
  double small;    
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
  double over3;
public:
  VWN3();
  DFTInfo eval(double rhoA, double rhoB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB);
  DFTInfo eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB);
  double EvepsVWN(int iop, double A_x, double b_x, double c_x, double x0_x, double rho);
};

