#include<dft.h>

VWNV::VWNV(){
// General Constants
  this->small = 1.0e-10; 
  this->over2 = 0.5;
  this->over3 = 1.0/3.0;
  this->over4 = 1.0/4.0;
  this->over6 = 1.0/6.0;
  this->fourover3 = 4.0/3.0;
// Parameter for the fit (according Eq 4.4 Vosko Can. J. Phys. 1980
  this->A1    =   0.0621814; // In the text page 1207 (A^P)
  this->A_p   =  this->A1/2.0; // to Hartree
  this->A_f   =  this->A1/4.0; // to Hartree
  this->A_a   = -(1.0/(6.0*math.pi*math.pi)) ;// to hartree already
// Updated
  this->b_p   = 3.72744 ; // Caption Table 5
  this->c_p   = 12.9352 ; // Caption Table 5
  this->x0_p  = -0.10498 ; // Caption Table 5
  this->b_f   = 7.06042 ; // Caption Table 5
  this->c_f   = 18.0578;  // Caption Table 5
  this->x0_f  = -0.32500 ; // Caption Table 5
//
  this->b_a   =   1.13107;   // intext page 1209
  this->c_a   =  13.0045;    // intext page 1209
  this->x0_a  =  -0.00475840; // intext page 1209
  this->popVWNconst();
};

void VWNV::popVWNconst(){
  this->b1p = (this->b_p*this->x0_p - this->c_p)/(this->c_p*this->x0_p); 
  this->b2p = (this->x0_p - this->b_p)/(this->c_p*this->x0_p); 
  this->b3p = (-1.0)/(this->c_p*this->x0_p); 
  this->Qp  = std::pow((4.0*this->c_p - this->b_p*this->b_p),over2); 
  this->X_x0p     = this->x0_p*this->x0_p + this->b_p*this->x0_p + this->c_p; 
  this->b1f = (this->b_f*this->x0_f - this->c_f)/(this->c_f*this->x0_f); 
  this->b2f = (this->x0_f - this->b_f)/(this->c_f*this->x0_f); 
  this->b3f = (-1.0)/(this->c_f*this->x0_f); 
  this->Qf  = std::pow((4.0*this->c_f - this->b_f*this->b_f),over2); 
  this->X_x0f     = this->x0_f*this->x0_f + this->b_f*this->x0_f + this->c_f; 
  this->Qa  = std::pow((4.0*this->c_a - this->b_a*this->b_a),over2); 
  this->X_x0a     = this->x0_a*this->x0_a + this->b_a*this->x0_a + this->c_a; 
};

void VWNV::popVWNdens(double rhoA, double rhoB){
  this->rhoT          = rhoA + rhoB;
  this->spindensity   = (rhoA - rhoB) / this->rhoT;
  this->spindensity_4 = std::pow(this->spindensity,4.0);
  this->spindensity_3 = std::pow(this->spindensity,3.0);

  this->f0_spindensity = 0.0;
  if ((1.0+this->spindensity) >= this->small)   f0_spindensity += std::pow((1.0+this->spindensity),this->fourover3); 
  if ((1.0-this->spindensity) >= this->small)   f0_spindensity += std::pow((1.0-(this->spindensity)),this->fourover3); 
  this->f0_spindensity += -2.0;
  this->f0_spindensity /= (-2.0+std::pow((2.0),this->fourover3)); 
  this->f1_spindensity  = std::pow((1.0+this->spindensity),this->fourover3); 
  this->f1_spindensity += std::pow((1.0-this->spindensity),this->fourover3); 
  this->f1_spindensity /= (2.0) ;
  this->df_spindensity = 0.0;
  if ((1.0+this->spindensity)   >= this->small) this->df_spindensity += std::pow((1.0+this->spindensity),this->over3); 
  if ((1.0-(this->spindensity)) >= this->small) this->df_spindensity -= std::pow((1.0-this->spindensity),this->over3); 
  this->df_spindensity *= this->fourover3;
  this->df_spindensity /= (-2.0+std::pow((2.0),this->fourover3));  
  if ((1.0+this->spindensity) >= this->small) this->df2_spindensity +=  std::pow((1.0+this->spindensity),(-2.0/3.0)); 
  if ((1.0-this->spindensity) >= this->small) this->df2_spindensity +=  std::pow((1.0-this->spindensity),(-2.0/3.0)); 
  this->df2_spindensity *= (2.0/9.0);
  this->df2_spindensity /= (-1.0+std::pow((2.0),(1.0/3.0)));
  this->r_s      = std::pow(((3.0)/(4.0*math.pi*this->rhoT)),this->over3);
  this->r_s_sqrt = std::pow(((3.0)/(4.0*math.pi*this->rhoT)),this->over6);
  this->r_s_32   = std::pow(((3.0)/(4.0*math.pi*this->rhoT)),this->over2);
  this->Xp        = this->r_s + this->b_p*(this->r_s_sqrt) + this->c_p; 
  this->Xf        = this->r_s + this->b_f*(this->r_s_sqrt) + this->c_f; 
  this->Xa        = this->r_s + this->b_a*(this->r_s_sqrt) + this->c_a; 
};


double VWNV::Eveps0VWN(double &A_x, double &b_x, double &Q, double &X, 
  double &x0_x, double &X_x0){
//    From Reference Vosko en Al., Can. J. Phys., 58, 1200 (1980). VWNV and VWN5 interpolation formula   
//    IOP 0 -> Eq 4.4 
  double val      = 0.0;
   val = A_x *
       ( std::log(this->r_s/X) + 
         2.0*b_x*(std::atan(Q/(2.0*this->r_s_sqrt + b_x)))/Q  -
         b_x*x0_x*(1.0/X_x0) * ( 
                               (std::log( (std::pow((this->r_s_sqrt-x0_x),2.0))/X )) +
                               (2.0*(b_x + 2.0*x0_x)*(1.0/Q)*(std::atan(Q/(2.0*this->r_s_sqrt + b_x))) ) 
                               ) 
       );
   return val;

}

double VWNV::Eveps1VWN(double &A_x, double &b1, double &b2, double &b3){
//    From Reference Vosko en Al., Can. J. Phys., 58, 1200 (1980). VWNV and VWN5 interpolation formula   
//    IOP 1 Eq. 4.3 (finishing the derivate of eps , rs factor already included)
  double val      = 0.0;
  val = A_x* ( (1.0 + b1*this->r_s_sqrt)/(1.0 + b1*this->r_s_sqrt + b2*this->r_s + b3*this->r_s_32));
  return val;

}

double VWNV::Eveps2VWN(double A_x, double &b_x, double &c_x, double &X, double &x0_x){
//    From Reference Vosko en Al., Can. J. Phys., 58, 1200 (1980). VWNV and VWN5 interpolation formula   
//    IOP 2 Analitic Derv of Eq 4.4 (note this one has to be moltiplied outside by rs to get the final needed term)
  double val      = 0.0;
  this->tmp1 = this->r_s_sqrt - x0_x;  //dxx0
  this->tmp2 = b_x * this->r_s_sqrt * x0_x; // bxx0
  this->tmp3 = X ; //c +bx+r
  val  = A_x*(c_x*this->tmp1 - this->tmp2);
  val /= (this->r_s*this->tmp3*this->tmp1);
   return val;
}

DFTFunctional::DFTInfo VWNV::eval(double rhoA, double rhoB){
  DFTFunctional::DFTInfo info;
   this->popVWNdens(rhoA, rhoB);
   if(std::abs(this->spindensity) > this->small) {
//   Open Shell Case
//   Used Linear Interpolation between parg and ferr 
//   Eq 2.4 and its analytic derivative for VWNV
     this->eps_p =  Eveps0VWN(this->A_p,this->b_p,this->Qp,this->Xp,this->x0_p,this->X_x0p);  
     this->eps_f =  Eveps0VWN(this->A_f,this->b_f,this->Qf,this->Xf,this->x0_f,this->X_x0f); 
//     Used improved Interpolation between parg and ferr 
//     Eq 3.2  and 3.3 and its analytic derivative for VWN5
     this->alpha =  Eveps0VWN(this->A_a,this->b_a,this->Qa,this->Xa,this->x0_a,this->X_x0a);  
     this->delta_eps_1 = this->eps_f - (this->eps_p);
     this->beta  = this->df2_spindensity* this->delta_eps_1 / this->alpha;
     this->beta  += -1.0;
     this->delta_eps_etha = this->alpha;
     this->delta_eps_etha *= (this->f0_spindensity/this->df2_spindensity);
     this->delta_eps_etha *= (1.0 + this->beta*this->spindensity_4);
     info.eps  = this->eps_p +  this->delta_eps_etha;
//   dbeta/dr
     this->db_dr = - this->delta_eps_1 * 
       this->Eveps2VWN(this->A_a,this->b_a,this->c_a,this->Xa,this->x0_a);
     this->db_dr += 
       (Eveps2VWN(this->A_f,this->b_f,this->c_f,this->Xf,this->x0_f) -
        Eveps2VWN(this->A_p,this->b_p,this->c_p,this->Xp,this->x0_p)) 
        * this->alpha;
     this->db_dr *= this->df2_spindensity;
     this->db_dr /= this->alpha*this->alpha;
//   S1 = da/dr * (f(zeta))*(1+zeta^4*beta)/ df2_spindens(0.0)
     this->S1  = this->f0_spindensity;
     this->S1 *= (1.0 + this->beta * this->spindensity_4);
     this->S1 *= 
       this->Eveps2VWN(this->A_a,this->b_a,this->c_a,this->Xa,this->x0_a);
     this->S1 /= this->df2_spindensity;
//   S2 = d eps_p/ dr
     this->S2  =
     this->Eveps2VWN(this->A_p,this->b_p,this->c_p,this->Xp,this->x0_p);
     this->S3  =  this->alpha;
     this->S3 *= (1.0 + this->beta*this->spindensity_4);
     this->S3 *= this->df_spindensity; 
     this->S3 /= this->df2_spindensity;
//   M1 alpha * f(zeta)/this->df2_spindens(0.0)
     this->M1  = this->alpha;
     this->M1 *= this->f0_spindensity; 
     this->M1 /= this->df2_spindensity;
//   S4  zeta^4 * dbeta/dr
     this->S4  = this->spindensity_4 * this->db_dr;
//   S5  4 * beta * zeta^3
     this->S5  = 4.0 * this->beta * this->spindensity_3;
     this->M3_A    =   1.0 - (this->spindensity); 
     this->M3_B    = -(1.0 + this->spindensity);
     info.ddrhoA   = this->r_s*this->over3*
      (this->S1 + this->S2 + this->M1*this->S4);

     info.ddrhoB   =  info.ddrhoA;     
     info.ddrhoA  +=  this->M3_A*(this->M1*this->S5 + this->S3);
     info.ddrhoB  +=  this->M3_B*(this->M1*this->S5 + this->S3);
     info.ddrhoA  += info.eps;
     info.ddrhoB  += info.eps ;
// Closed Shell
   } else {
     info.eps =  Eveps0VWN(this->A_p,this->b_p,this->Qp,this->Xp,this->x0_p,this->X_x0p);
     info.ddrhoA  = -(this->over3)*this->r_s*Eveps2VWN(this->A_p,this->b_p,this->c_p,this->Xp,this->x0_p);
     info.ddrhoA += info.eps ;
   }
  return info;
}

DFTFunctional::DFTInfo VWNV::eval(double rhoA, double rhoB, double gammaAA, double gammaAB){
};

DFTFunctional::DFTInfo VWNV::eval(double rhoA, double rhoB, double gammaAA, double gammaAB, double gammaBB){
};

