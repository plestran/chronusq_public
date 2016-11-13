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
#include <dft.h>

SlaterExchange::SlaterExchange(double X, double eps):
DFTFunctional(X,eps){
// Memo Factor to be added at the end for numerical stability
  this->CxVx  = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));  
//  this->small = 1.0e-16; 
  this-> d1over3  = 1.0/3.0;
  this-> d4over3  = 4.0/3.0;

  this->name = "Slater";
#ifdef CQ_ENABLE_LIBXC
  xc_func_init(&this->func,XC_LDA_X,XC_POLARIZED);
#endif
};

DFTFunctional::DFTInfo SlaterExchange::eval(const double &rhoA, const double &rhoB){

};

DFTFunctional::DFTInfo SlaterExchange::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB){
};

DFTFunctional::DFTInfo SlaterExchange::eval(const double &rhoA, const double &rhoB, const double &gammaAA, const double &gammaAB, const double &gammaBB){
// alpha in A2 is 2/3
// the final expression for the delE/DelphoSigma = - (2)^(1/3) * (3/pi)^(1/3) * rho_sigma
  DFTFunctional::DFTInfo info;
  double rhoT            = rhoA + rhoB;
#ifndef CQ_ENABLE_LIBXC
  double spindensity     = (rhoA - rhoB) / rhoT;
  double eps_spin;
  double tmp1Sppow;
  double tmp2Sppow;
  info.eps        = std::pow(rhoT,this->d1over3)*this->CxVx;      
  if (std::abs(spindensity) > this->small){
// Open Shell
    tmp1Sppow  = std::pow((1.0+spindensity),this->d1over3);
    tmp2Sppow  = std::pow((1.0-spindensity),this->d1over3);
    eps_spin   = tmp1Sppow*(1.0+spindensity) ;      
    eps_spin  += tmp2Sppow*(1.0-spindensity) ;      
    eps_spin   /= 2.0;      
    info.ddrhoA = 
      this->d4over3*info.eps*tmp1Sppow; 
    info.ddrhoB = 
      this->d4over3*info.eps*tmp2Sppow;   
    info.eps *= eps_spin;
  } else {   
// Closed Shell
    info.ddrhoA     = this->d4over3*info.eps;       
    info.ddrhoB     = info.ddrhoA;       
  }
#else
  std::array<double,2> rho = {rhoA,rhoB};
  std::array<double,2> vxc;
  
  xc_lda_exc_vxc(&this->func,1,&rho[0],&info.eps,&vxc[0]);
  info.ddrhoA = vxc[0];
  info.ddrhoB = vxc[1];
#endif
  info.eps  *= rhoT; 
  info *= this->scalingFactor; 

  return info;
};

