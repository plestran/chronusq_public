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

PBE::PBE(double X, double eps):
DFTFunctional(X,eps){
// Factor for not spinpolarized. Different from the one used in Slater for spin polarized. Eq2 from Barone J. Chem. Phys. 108, 664 (1998); 
  this->CxVx      = -(3.0/2.0)*(std::pow((3.0/(4.0*math.pi)),(1.0/3.0))); 
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;
  this-> mu       =  0.2195102405736; // PRL 77, 18 pg 3865 kappa page 3867
  this-> kappa    =  0.8040004238; // PRL 77, 18 pg 3865 kappa page 3867
  this-> Ckf      =  std::pow((math.pi*math.pi*6.0),(this->d1over3)); // Precator of Fermi wave number but it has 6 instead of 3
  this->name = "PBE_X";
#ifdef CQ_ENABLE_LIBXC
  xc_func_init(&this->func,XC_GGA_X_PBE,XC_POLARIZED);
#endif
};

double PBE::g0 (double x,double &den){
// Eq 14 PRL 77, 18 pg 3865 kappa page 3866
// Note,I factored in the LDA factor 
  double gx;
  den  = this->mu * x * x;
  den /= this->kappa;
  den += 1.0; 
  gx = -this->kappa;
  gx /= den;
  gx += 1.0 + this->kappa;
  return this->CxVx*gx;
};  //End Form g(x) or F9in Barone)

double PBE::g1 (double x, double &den){
// Derivative rispect to x( s in the paper)  of Eq 14 PRL 77, 18 pg 3865 kappa page 3866
// Note.I factored in the LDA factor 
  double gx;
  gx = 2.0 * this->mu * x ;
  gx /= (den*den);
  return this->CxVx*gx;
};  //End Form g function for B88 Exchange

DFTFunctional::DFTInfo PBE::eval(const double &rhoA, const double &rhoB){
};

DFTFunctional::DFTInfo PBE::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
  DFTFunctional::DFTInfo info;
#ifndef CQ_ENABLE_LIBXC
  double spindensity   = (rhoA - rhoB); 
  spindensity         /= (rhoA + rhoB);
  double rhoA1ov3 = std::pow(rhoA,this->d1over3);
  double rhoB1ov3 ;
  double rhoA4ov3 = rhoA1ov3 * rhoA; 
  double rhoB4ov3 ; 
  if(std::abs(spindensity) > this->small) {
    rhoB1ov3 = std::pow(rhoB,this->d1over3);
    rhoB4ov3 = rhoB1ov3 * rhoB; 
  }
// x = s of paper PRL 77, 18 pg 3865 kappa page 3866 
// s = | MOD(GRAD(rho)) / 2 * Kf * rho
// where Kf (Fermi wavenumber) = Ckf * rho^(1/3) 
  double xA   = std::sqrt(gammaAA) / (2.0 * this->Ckf*rhoA4ov3); 
  double xB ;
  double den;
  double g0 = this->g0(xA,den);
  double g1 = this->g1(xA,den);
  info.eps        = rhoA4ov3*g0;
  info.ddrhoA     = g0 - xA*g1;
  info.ddrhoA    *= this->d4over3*rhoA1ov3;
  info.ddgammaAA  = g1/(4.0*this->Ckf*std::sqrt(gammaAA));
  if(std::abs(spindensity) > this->small) {
    //Open Shell
    xB   = std::sqrt(gammaBB) / (2.0 * this->Ckf*rhoB4ov3); 
    g0 = this->g0(xB,den);
    g1 = this->g1(xB,den);
    info.eps       += rhoB4ov3*g0;
    info.ddrhoB     = g0 - xB*g1;
    info.ddrhoB    *= this->d4over3*rhoB1ov3;
    info.ddgammaBB  = g1/(4.0*this->Ckf*std::sqrt(gammaBB));

  } else {
    //Closed Shell
    info.eps *= 2.0;
    info.ddrhoB     = info.ddrhoA;
    info.ddgammaBB  = info.ddgammaAA;
  }
  info.ddgammaAB  = 0.0;
#else
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

  info.eps *= (rhoA + rhoB); //?

#endif
  info *= this->scalingFactor; 
  return info;
};

DFTFunctional::DFTInfo PBE::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaBB){
};

