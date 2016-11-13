/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#include<dft.h>

VWNV::VWNV(double X, double eps):
VWNIII(X,eps){

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

  this->name = "VWN V";
#ifdef CQ_ENABLE_LIBXC
  xc_func_init(&this->func,XC_LDA_C_VWN,XC_POLARIZED);
#endif
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

void VWNV::popVWNdens(const double &rhoA,const double &rhoB, denspow &denquant){
  denquant.rhoT          = rhoA + rhoB;
  denquant.spindensity   = (rhoA - rhoB) / denquant.rhoT;
  denquant.spindensity_3 = denquant.spindensity*denquant.spindensity*denquant.spindensity;
  denquant.spindensity_4 = denquant.spindensity*denquant.spindensity_3;
  denquant.f0_spindensity = 0.0;

  if (std::abs(denquant.spindensity) >= this->small)   
     denquant.f0_spindensity += std::pow((1.0+denquant.spindensity),this->fourover3); 
  if (std::abs(denquant.spindensity) >= this->small)   
     denquant.f0_spindensity += std::pow((1.0-(denquant.spindensity)),this->fourover3); 
  denquant.f0_spindensity += -2.0;
  denquant.f0_spindensity /= (-2.0+std::pow((2.0),this->fourover3)); 
  denquant.df_spindensity = 0.0;
  if (std::abs(denquant.spindensity)   >= this->small) 
     denquant.df_spindensity += std::pow((1.0+denquant.spindensity),this->over3); 
  if (std::abs(denquant.spindensity) >= this->small) 
     denquant.df_spindensity -= std::pow((1.0-denquant.spindensity),this->over3); 
  denquant.df_spindensity *= this->fourover3;
  denquant.df_spindensity /= (-2.0+std::pow((2.0),this->fourover3));  
  denquant.df2_spindensity = 2.0; 
  denquant.df2_spindensity *= (2.0/9.0);
  denquant.df2_spindensity /= (-1.0+std::pow((2.0),(1.0/3.0)));
  denquant.r_s      = std::pow(((3.0)/(4.0*math.pi*denquant.rhoT)),this->over3);
  denquant.r_s_sqrt = std::pow(((3.0)/(4.0*math.pi*denquant.rhoT)),this->over6);
  denquant.r_s_32   = std::pow(((3.0)/(4.0*math.pi*denquant.rhoT)),this->over2);
  denquant.Xp       = denquant.r_s + this->b_p*(denquant.r_s_sqrt) + this->c_p; 
  denquant.Xf       = denquant.r_s + this->b_f*(denquant.r_s_sqrt) + this->c_f; 
  denquant.Xa       = denquant.r_s + this->b_a*(denquant.r_s_sqrt) + this->c_a; 
};

DFTFunctional::DFTInfo VWNV::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
  DFTFunctional::DFTInfo info;
   if ((rhoA + rhoB) < this->small) { return info;}
   denspow RhoTQuant;
   this->popVWNdens(rhoA, rhoB, RhoTQuant);
   if(std::abs(RhoTQuant.spindensity) > this->small) {
//   Open Shell Case
//   Used Linear Interpolation between parg and ferr 
//   Eq 2.4 and its analytic derivative for VWNV
     RhoTQuant.eps_p =  Eveps0VWN(this->A_p,this->b_p,this->Qp,RhoTQuant.Xp,this->x0_p,this->X_x0p,RhoTQuant);  
     RhoTQuant.eps_f =  Eveps0VWN(this->A_f,this->b_f,this->Qf,RhoTQuant.Xf,this->x0_f,this->X_x0f,RhoTQuant); 
//     Used improved Interpolation between parg and ferr 
//     Eq 3.2  and 3.3 and its analytic derivative for VWN5
     RhoTQuant.alpha =  Eveps0VWN(this->A_a,this->b_a,this->Qa,RhoTQuant.Xa,this->x0_a,this->X_x0a,RhoTQuant);  
     RhoTQuant.delta_eps_1 = RhoTQuant.eps_f - (RhoTQuant.eps_p);
     RhoTQuant.beta  = RhoTQuant.df2_spindensity* RhoTQuant.delta_eps_1 / RhoTQuant.alpha;
     RhoTQuant.beta  += -1.0;
     RhoTQuant.delta_eps_etha = RhoTQuant.alpha;
     RhoTQuant.delta_eps_etha *= (RhoTQuant.f0_spindensity/RhoTQuant.df2_spindensity);
     RhoTQuant.delta_eps_etha *= (1.0 + RhoTQuant.beta*RhoTQuant.spindensity_4);
       info.eps  = RhoTQuant.eps_p +  RhoTQuant.delta_eps_etha;
//   dbeta/dr
     RhoTQuant.db_dr = - RhoTQuant.delta_eps_1 * 
       this->Eveps2VWN(this->A_a,this->b_a,this->c_a,RhoTQuant.Xa,this->x0_a,RhoTQuant);
     RhoTQuant.db_dr += 
       (Eveps2VWN(this->A_f,this->b_f,this->c_f,RhoTQuant.Xf,this->x0_f,RhoTQuant) -
        Eveps2VWN(this->A_p,this->b_p,this->c_p,RhoTQuant.Xp,this->x0_p,RhoTQuant)) 
        * RhoTQuant.alpha;
     RhoTQuant.db_dr *= RhoTQuant.df2_spindensity;
     RhoTQuant.db_dr /= RhoTQuant.alpha*RhoTQuant.alpha;
//   S1: = da/dr * (f(zeta))*(1+zeta^4*beta)/ df2_spindens(0.0)
     RhoTQuant.S1  = RhoTQuant.f0_spindensity;
     RhoTQuant.S1 *= (1.0 + RhoTQuant.beta * RhoTQuant.spindensity_4);
     RhoTQuant.S1 *= 
       this->Eveps2VWN(this->A_a,this->b_a,this->c_a,RhoTQuant.Xa,this->x0_a,RhoTQuant);
     RhoTQuant.S1 /= RhoTQuant.df2_spindensity;
//   S2: d eps_p/ dr
     RhoTQuant.S2  =
       this->Eveps2VWN(this->A_p,this->b_p,this->c_p,RhoTQuant.Xp,this->x0_p,RhoTQuant);
     RhoTQuant.S3  =  RhoTQuant.alpha;
     RhoTQuant.S3 *= (1.0 + RhoTQuant.beta*RhoTQuant.spindensity_4);
     RhoTQuant.S3 *= RhoTQuant.df_spindensity; 
     RhoTQuant.S3 /= RhoTQuant.df2_spindensity;
//   M1: alpha * f(zeta)/this->df2_spindens(0.0)
     RhoTQuant.M1  = RhoTQuant.alpha;
     RhoTQuant.M1 *= RhoTQuant.f0_spindensity; 
     RhoTQuant.M1 /= RhoTQuant.df2_spindensity;
//   S4:  zeta^4 * dbeta/dr
     RhoTQuant.S4  = RhoTQuant.spindensity_4 * RhoTQuant.db_dr;
//   S5:  4 * beta * zeta^3
     RhoTQuant.S5  = 4.0 * RhoTQuant.beta * RhoTQuant.spindensity_3;
     RhoTQuant.M3_A    = 1.0 - (RhoTQuant.spindensity); 
     RhoTQuant.M3_B    = -(1.0 + RhoTQuant.spindensity);
     info.ddrhoA   = -RhoTQuant.r_s*this->over3*
      (RhoTQuant.S1 + RhoTQuant.S2 + RhoTQuant.M1*RhoTQuant.S4);

     info.ddrhoB   =  info.ddrhoA;     
     info.ddrhoA  +=  RhoTQuant.M3_A*(RhoTQuant.M1*RhoTQuant.S5 + RhoTQuant.S3);
     info.ddrhoB  +=  RhoTQuant.M3_B*(RhoTQuant.M1*RhoTQuant.S5 + RhoTQuant.S3);
     info.ddrhoA  += info.eps;
     info.ddrhoB  += info.eps ;
// Closed Shell
   } else {
     info.eps =  Eveps0VWN(this->A_p,this->b_p,this->Qp,RhoTQuant.Xp,this->x0_p,this->X_x0p,RhoTQuant);
     info.ddrhoA  = -(this->over3)*RhoTQuant.r_s*
                    Eveps2VWN(this->A_p,this->b_p,this->c_p,RhoTQuant.Xp,this->x0_p,RhoTQuant);
     info.ddrhoA += info.eps ;
     info.ddrhoB  = info.ddrhoA ;
   }
     info.eps *= RhoTQuant.rhoT; 
  info *= this->scalingFactor; 
  return info;
};

