#include<dft.h>

VWNIII::VWNIII(double X, double eps):
DFTFunctional(X,eps){
// General Constants
//  this->small = 1.0e-16; 
  this->over2 = 0.5;
  this->over3 = 1.0/3.0;
  this->over4 = 1.0/4.0;
  this->over6 = 1.0/6.0;
  this->fourover3 = 4.0/3.0;
// Parameter for the fit (according Eq 4.4 Vosko Can. J. Phys. 1980
  this->A1    =   0.0621814; // In the text page 1207 (A^P)
  this->A_p   =  this->A1/2.0; // to Hartree
  this->A_f   =  this->A1/4.0; // to Hartree
//  this->A_a   = -(1.0/(6.0*math.pi*math.pi)) ;// to hartree already
  this->b_p   =  13.0720;   // into text page 1207
  this->c_p   =  42.7198;   // into text page 1207
  this->x0_p  =  -0.409286; // into text page 1207
  this->b_f   =  20.1231;   // into text pagr 1207
  this->c_f   = 101.578;   // into text pagr 1207
  this->x0_f  =  -0.743294;  // into text pagr 1207
//  this->b_a   =   1.13107;   // intext page 1209
//  this->c_a   =  13.0045;    // intext page 1209
//  this->x0_a  =  -0.00475840; // intext page 1209
  this->popVWNconst();

  this->name = "VWN III";
#ifdef CQ_ENABLE_LIBXC
  // This is a dummy initialization to VWN 5
  // FIXME: LibXC claims to have a bunch of different
  // VWN functionals, maybe check to see if one of them
  // is what we call VWN 3?
  xc_func_init(&this->func,XC_LDA_C_VWN,XC_POLARIZED);
#endif
};

void VWNIII::popVWNconst(){
//  this->b1p = (this->b_p*this->x0_p - this->c_p)/(this->c_p*this->x0_p); 
//  this->b2p = (this->x0_p - this->b_p)/(this->c_p*this->x0_p); 
//  this->b3p = (-1.0)/(this->c_p*this->x0_p); 
  this->Qp  = std::pow((4.0*this->c_p - this->b_p*this->b_p),over2); 
  this->X_x0p     = this->x0_p*this->x0_p + this->b_p*this->x0_p + this->c_p; 
//  this->b1f = (this->b_f*this->x0_f - this->c_f)/(this->c_f*this->x0_f); 
//  this->b2f = (this->x0_f - this->b_f)/(this->c_f*this->x0_f); 
//  this->b3f = (-1.0)/(this->c_f*this->x0_f); 
  this->Qf  = std::pow((4.0*this->c_f - this->b_f*this->b_f),over2); 
  this->X_x0f     = this->x0_f*this->x0_f + this->b_f*this->x0_f + this->c_f; 
};

void VWNIII::popVWNdens(const double &rhoA,const double &rhoB, denspow &denquant){
  denquant.rhoT          = rhoA + rhoB;
  denquant.spindensity   = (rhoA - rhoB) / denquant.rhoT;
//  denquant.spindensity_3 = denquant.spindensity*denquant.spindensity*denquant.spindensity;
//  denquant.spindensity_4 = denquant.spindensity*denquant.spindensity_3;

  denquant.f0_spindensity = 0.0;
  if (std::abs(denquant.spindensity) >= this->small)   
     denquant.f0_spindensity += std::pow((1.0+denquant.spindensity),this->fourover3); 
  if (std::abs(denquant.spindensity) >= this->small)   
     denquant.f0_spindensity += std::pow((1.0-(denquant.spindensity)),this->fourover3); 
  denquant.f0_spindensity += -2.0;
  denquant.f0_spindensity /= (-2.0+std::pow((2.0),this->fourover3)); 
//  denquant.f1_spindensity  = std::pow((1.0+denquant.spindensity),this->fourover3); 
//  denquant.f1_spindensity += std::pow((1.0-denquant.spindensity),this->fourover3); 
//  denquant.f1_spindensity /= (2.0) ;
  denquant.df_spindensity = 0.0;
  if (std::abs(denquant.spindensity)   >= this->small) 
     denquant.df_spindensity += std::pow((1.0+denquant.spindensity),this->over3); 
  if (std::abs(denquant.spindensity) >= this->small) 
     denquant.df_spindensity -= std::pow((1.0-denquant.spindensity),this->over3); 
  denquant.df_spindensity *= this->fourover3;
  denquant.df_spindensity /= (-2.0+std::pow((2.0),this->fourover3));  
  denquant.r_s      = std::pow(((3.0)/(4.0*math.pi*denquant.rhoT)),this->over3);
  denquant.r_s_sqrt = std::pow(((3.0)/(4.0*math.pi*denquant.rhoT)),this->over6);
  denquant.r_s_32   = std::pow(((3.0)/(4.0*math.pi*denquant.rhoT)),this->over2);
  denquant.Xp        = denquant.r_s + this->b_p*(denquant.r_s_sqrt) + this->c_p; 
  denquant.Xf        = denquant.r_s + this->b_f*(denquant.r_s_sqrt) + this->c_f; 
};


double VWNIII::Eveps0VWN(double &A_x, double &b_x, double &Q, double &X, 
  double &x0_x, double &X_x0, denspow &denquant){
//    From Reference Vosko en Al., Can. J. Phys., 58, 1200 (1980). VWNIII and VWN5 interpolation formula   
//    IOP 0 -> Eq 4.4 
  double val      = 0.0;
   val = A_x *
       ( std::log(denquant.r_s/X) + 
         2.0*b_x*(std::atan(Q/(2.0*denquant.r_s_sqrt + b_x)))/Q  -
         b_x*x0_x*(1.0/X_x0) * ( 
                               (std::log( (std::pow((denquant.r_s_sqrt-x0_x),2.0))/X )) +
                               (2.0*(b_x + 2.0*x0_x)*(1.0/Q)*(std::atan(Q/(2.0*denquant.r_s_sqrt + b_x))) ) 
                               ) 
       );
   return val;

}

double VWNIII::Eveps1VWN(double &A_x, double &b1, double &b2, double &b3, denspow &denquant){
//    From Reference Vosko en Al., Can. J. Phys., 58, 1200 (1980). VWNIII and VWN5 interpolation formula   
//    IOP 1 Eq. 4.3 (finishing the derivate of eps , rs factor already included)
  double val      = 0.0;
  val = A_x* ( (1.0 + b1*denquant.r_s_sqrt)/(1.0 + b1*denquant.r_s_sqrt + b2*denquant.r_s + b3*denquant.r_s_32));
  return val;

}

double VWNIII::Eveps2VWN(double A_x, double &b_x, double &c_x, double &X, double &x0_x, denspow &denquant){
//    From Reference Vosko en Al., Can. J. Phys., 58, 1200 (1980). VWNIII and VWN5 interpolation formula   
//    IOP 2 Analitic Derv of Eq 4.4 (note this one has to be moltiplied outside by rs to get the final needed term)
  double val      = 0.0;
  double tmp1 = denquant.r_s_sqrt - x0_x;  //dxx0
  double tmp2 = b_x * denquant.r_s_sqrt * x0_x; // bxx0
  double tmp3 = X ; //c +bx+r
  val  = A_x*(c_x*tmp1 - tmp2);
  val /= (denquant.r_s*tmp3*tmp1);
   return val;
}

DFTFunctional::DFTInfo VWNIII::eval(const double &rhoA, const double &rhoB){
}

DFTFunctional::DFTInfo VWNIII::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB){
};

DFTFunctional::DFTInfo VWNIII::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
  DFTFunctional::DFTInfo info;
   if ((rhoA + rhoB) < this->small) { return info;}
  denspow RhoTQuant;
   this->popVWNdens(rhoA, rhoB, RhoTQuant);
   if(std::abs(RhoTQuant.spindensity) > this->small) {
//   Open Shell Case
//   Used Linear Interpolation between parg and ferr 
//   Eq 2.4 and its analytic derivative for VWNIII
     RhoTQuant.eps_p = 
       this->Eveps0VWN(this->A_p,this->b_p,this->Qp,RhoTQuant.Xp,this->x0_p,this->X_x0p,RhoTQuant);  
     RhoTQuant.eps_f = 
       this->Eveps0VWN(this->A_f,this->b_f,this->Qf,RhoTQuant.Xf,this->x0_f,this->X_x0f,RhoTQuant); 
     RhoTQuant.delta_eps_1 = RhoTQuant.eps_f - (RhoTQuant.eps_p);
     info.eps  = RhoTQuant.eps_p + RhoTQuant.delta_eps_1*RhoTQuant.f0_spindensity;
     RhoTQuant.S1 =  -(RhoTQuant.r_s)*this->over3*
       this->Eveps2VWN(this->A_p,this->b_p,this->c_p,RhoTQuant.Xp,this->x0_p,RhoTQuant);
     RhoTQuant.S2 =  -(RhoTQuant.r_s)*this->over3*RhoTQuant.f0_spindensity*
      (this->Eveps2VWN(this->A_f,this->b_f,this->c_f,RhoTQuant.Xf,this->x0_f,RhoTQuant) - 
       this->Eveps2VWN(this->A_p,this->b_p,this->c_p,RhoTQuant.Xp,this->x0_p,RhoTQuant) 
       );
     RhoTQuant.M3_A    =   1.0 - (RhoTQuant.spindensity); 
     RhoTQuant.M3_B    = -(1.0 + RhoTQuant.spindensity);
     info.ddrhoA   = RhoTQuant.S1 + RhoTQuant.S2 + info.eps;
     info.ddrhoB   = RhoTQuant.S1 + RhoTQuant.S2 + info.eps;     
     info.ddrhoA  += RhoTQuant.delta_eps_1*RhoTQuant.M3_A*RhoTQuant.df_spindensity;
     info.ddrhoB  += RhoTQuant.delta_eps_1*RhoTQuant.M3_B*RhoTQuant.df_spindensity;
// Closed Shell
   } else {
     info.eps = 
     this->Eveps0VWN(this->A_p,this->b_p,this->Qp,RhoTQuant.Xp,this->x0_p,this->X_x0p,RhoTQuant);
     info.ddrhoA  = -(this->over3)*
     RhoTQuant.r_s*
     this->Eveps2VWN(this->A_p,this->b_p,this->c_p,RhoTQuant.Xp,this->x0_p,RhoTQuant);
     info.ddrhoA += info.eps ;
     info.ddrhoB  = info.ddrhoA ;
   }
     info.eps *= (rhoA+rhoB);
  return info;
};

