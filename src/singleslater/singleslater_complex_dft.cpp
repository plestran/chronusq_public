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
#include <singleslater.h>
namespace ChronusQ {
//----------------------------//
// form the Vxc matrix        //
//----------------------------//

/*
template<>
double SingleSlater<dcomplex>::formBeckeW(cartGP gridPt, int iAtm){
  return 0.0;
};
template<>
double SingleSlater<dcomplex>::normBeckeW(cartGP gridPt){
  return 0.0;
};

template<>
double SingleSlater<dcomplex>::f_spindens(int iop, double spindensity){
  return 0.0;

};

template<>
double SingleSlater<dcomplex>::EvepsVWN(int iop, double A_x, double b_x, double c_x, double x0_x, double rho){
  return 0.0;
};

template<>
double SingleSlater<dcomplex>::df_spindens(double spindensity){
  return 0.0;

};

template<>
double SingleSlater<dcomplex>::df2_spindens(double spindensity){
  return 0.0;

};

template<>
double SingleSlater<dcomplex>::spindens(double rho_A, double rho_B){
  return 0.0;

};

template<>
double SingleSlater<dcomplex>::gB88(int nDer, double x){
  return 0.0;

};
*/

/*
template<>
void SingleSlater<dcomplex>::genSparseBasisMap(){
};

template<>
void SingleSlater<dcomplex>::genSparseRcrosP(){
};
*/

/*
template<>
std::array<double,6> SingleSlater<dcomplex>::formVExSlater (double rho, double spindensity){
  return {0.0,0.0,0.0,0.0,0.0,0.0};
}; 
template<>
double SingleSlater<dcomplex>::omegaLYP(int iop, double c, double d, double rho){  //function used in LYP Correlation (Eq. A26, A32 Ref.LYP2)
  return 0.0;
};//End omegaLYP

template<>
double SingleSlater<dcomplex>::deltaLYP(int iop, double c, double d, double rho){  //function used in LYP Correlation (Eq. A27, A33 Ref.LYP2)
  return 0.0;
};//End deltaLYP

template<>
double SingleSlater<dcomplex>::derLYP(int iop, double a, double b, double c, double d, double rhoA, double rhoB){  //function used in LYP Correlation (A28, A29,A30,A31 Ref.LYP2)
  return 0.0;
};//End derLYP

template<>
double SingleSlater<dcomplex>::der2LYP(int iop, double a, double b, double c, double d, double rhoA, double rhoB, double der1){  
  return 0.0;
};//End der2LYP


template<>
std::array<double,6> SingleSlater<dcomplex>::formVCLYP (double rhoA, double rhoB, 
     double gammaAA, double gammaBB, double gammaAB){
  return {0.0,0.0,0.0,0.0,0.0,0.0};
}; 

template<>
std::array<double,6> SingleSlater<dcomplex>::formVExB88 (double rhoA, double rhoB, 
     double drhoA, double drhoB){
  return {0.0,0.0,0.0,0.0,0.0,0.0};
}; 

template<>
std::array<double,6> SingleSlater<dcomplex>::formVCVWN (double rho, double spindensity){
  return {0.0,0.0,0.0,0.0,0.0,0.0};

}; 

template<>
std::array<double,6> SingleSlater<dcomplex>::formVC (double rho, double spindensity){
  return {0.0,0.0,0.0,0.0,0.0,0.0};

}; 

template<>
std::array<double,6> SingleSlater<dcomplex>::formVCGGA (double rhoA, double rhoB,
     double gammaAA, double gammaBB, double gammaAB){
  return {0.0,0.0,0.0,0.0,0.0,0.0};
}; 

template<>
std::array<double,6> SingleSlater<dcomplex>::formVExGGA (double rhoA, double rhoB,
     double drhoA, double drhoB){
  return {0.0,0.0,0.0,0.0,0.0,0.0};
}; 

template<>
std::array<double,6> SingleSlater<dcomplex>::formVEx(double rho, double spindensity){
  return {0.0,0.0,0.0,0.0,0.0,0.0};

}; 

template<>
void SingleSlater<dcomplex>::formVXC(){

};

template<>
void SingleSlater<dcomplex>::formVXC_store(){

};


template<>
void SingleSlater<dcomplex>::evalVXC(cartGP gridPt, double weight, std::vector<bool> mapRad_,
       double & energyX, double & energyC, RealMatrix * VXA, RealMatrix * VXB, RealMatrix * VCA, RealMatrix * VCB){

};

template<>
void SingleSlater<dcomplex>::evalVXC_store(int iAtm, int ipts, double & energyX, 
       double & energyC, RealMatrix * VXA, RealMatrix * VXB, RealMatrix * VCA, 
       RealMatrix * VCB, RealMatrix *STmp, RealMatrix *dSTmpX, RealMatrix *dSTmpY, RealMatrix *dSTmpZ){
};
*/






} // Namespace ChronusQ
