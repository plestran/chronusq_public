/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
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
#include <grid.h>
#include <cerr.h>

namespace ChronusQ{

  double OneDGrid::f_val(double rad){
  // Test Function to be integrated by One-dimensional grid
  // INT[0,1] r^2 * exp(-r^2);
         return (std::pow(rad,2.0))*(std::exp(-((std::pow(rad,2.0)))));
         } 

  double OneDGrid::f2_val(double elevation,double azimut){
  // Test Function to be integrated by One-dimensional grid over a solid angle
         return (15.0)*(std::pow(sin(elevation),4.0))/(32.0*math.pi);
         } 

  double OneDGrid::integrate(){
  // Integrate a test function for a one dimensional grid radial part
   double sum = 0.0;
     std::cout << "Number of One-grid points= "<< this->nPts_  <<std::endl;
     for(int i = 0; i < this->nPts_; i++){
       if(this->intas2GPt_){
         sum += (this->f2_val(bg::get<0>(this->grid2GPts_[i]),bg::get<1>(this->grid2GPts_[i])))*(this->weights_[i]);
         }else{
         sum += (this->f_val(this->gridPts_[i]))*(this->weights_[i]);
       }
      }
     if(this->intas2GPt_){
       return 4.0*math.pi*sum;
       }else{
       return sum*this->norm_;
     }
  }
void OneDGrid::printGrid(){
    for(int i = 0; i < this->nPts_; i++){
//  Printing for mathematica
    if(this->intas2GPt_){
      cout << "{ 1.0, "<<bg::get<1>(this->grid2GPts_[i])<<", " <<bg::get<0>(this->grid2GPts_[i]) <<"}, "<< endl;
      }else{
      cout << this->gridPts_[i] << endl;
      }    
}

};

  // Class Functions Declaration

  // Function Gauss-Chebyshev 1st kind 
  void GaussChebyshev1stGrid::genGrid(){
     // Gauss-Chebyshev 1st kind grid
     // Int {a,b} f(x) = Int {-1,1} g(x')/sqrt(1-x'^2) ~ Sum [1, NGrid] weights[i] * g(zeta[i])
     // Where g(zeta[i]) = [(b-a)/2] * sqrt(1-zeta[i]^2) f( [(b-a)/2]*zeta[i] + [(b+a)/2] )
     // and   zeta[i] = gridPts_[i] = cos ( [(2*i-1)*math.pi/2*NGrid])
     // weights[i] = math.pi/NGrid
     // Note I have included the "sqrt(1-zeta[i]^2)" of the transformation in the actual weights[i] .
     // Note in c++ i starts from 0, so has been shifted i = i+1 
    for(int i = 0; i < this->nPts_; i++) {
      this->gridPts_[i] = cos(( (2.0*(i+1)-1.0)/(2.0*this->nPts_))*math.pi);
      this->weights_[i] = (math.pi/this->nPts_)*(sqrt(1-(this->gridPts_[i]*this->gridPts_[i]))) ;
//    the weights are only (math.pi/this->nPts_), the second term is including the integrand transformation factor in Eq.25.4.38 Abramowitz Handbook     
//    std::cout << i <<" " <<this->gridPts_[i] << " "<< this->weights_[i]  << std::endl;
      }
    this->transformPts(); 
  }

  void GaussChebyshev1stGrid::transformPts(){
    for(int i = 0; i < this->nPts_; i++)
      this->gridPts_[i] = 
        (this->range_[1] - this->range_[0]) / 2.0 * this->gridPts_[i] +
        (this->range_[1] + this->range_[0]) / 2.0;
  } 

// Function definition for Lebedev
   void LebedevGrid::genGrid(){
     if(this->nPts_ == 6){
       gen6_A1(0,1.0,(1.0/6.0));
      }else if(this->nPts_== 18){
       gen6_A1(0,1.0,0.6666666666666667e-1);
       gen12_A2(6,(1.0/std::sqrt(2.0)),0.7500000000000000e-1);
      }else{
      CErr("Number of points not available in Lebedev quadrature");
      }
  };


    void LebedevGrid::transformPts(){
};

void LebedevGrid::gen6_A1(int num, long double a, long double v){
    
    cartGP tmpCart;
    tmpCart.set<0>(a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(0.0);
    this->weights_[num+0] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+0]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(0.0);
    this->weights_[num+1] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+1]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(a);
    tmpCart.set<2>(0.0);
    this->weights_[num+2] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+2]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(0.0);
    this->weights_[num+3] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+3]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(a);
    this->weights_[num+4] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+4]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(-a);
    this->weights_[num+5] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+5]));
}

void LebedevGrid::gen12_A2(int num, long double a, long double v){
    
    cartGP tmpCart;

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(a);
    tmpCart.set<2>(a);
    this->weights_[num+0] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+0]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(a);
    this->weights_[num+1] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+1]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(a);
    tmpCart.set<2>(-a);
    this->weights_[num+2] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+2]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(-a);
    this->weights_[num+3] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+3]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(a);
    this->weights_[num+4] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+4]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(a);
    this->weights_[num+5] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+5]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(-a);
    this->weights_[num+6] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+6]));
 
    tmpCart.set<0>(-a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>( -a);
    this->weights_[num+7] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+7]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(a);
    tmpCart.set<2>(0.0);
    this->weights_[num+8] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+8]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>( a);
    tmpCart.set<2>(0.0);
    this->weights_[num+9] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+9]));

    tmpCart.set<0>( a);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(0.0);
    this->weights_[num+10] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+10]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(0.0);
    this->weights_[num+11] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+11]));

}

}; // namespace ChronusQ
