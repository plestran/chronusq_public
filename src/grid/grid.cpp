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

namespace ChronusQ{

  class GaussChebyshev1stGrid : public OneDGrid {
    public:
      GaussChebyshev1stGrid(
        int npts = 0, double beg = 0.0, double end = 0.0):
        OneDGrid(npts,beg,end){
          this->gridPts_ = new double[this->nPts_];        ///< Grid Points
          this->weights_ = new double[this->nPts_];        ///< Weights
        };
  // Class Functions
    void genGrid();                                      
    void transformPts();
  }; // class GaussChebyshev1stGrid

  // Class Functions Declaration

  void GaussChebyshev1stGrid::genGrid(){
     // Int {a,b} f(x) = Int {-1,1} g(x')/sqrt(1-x'^2) ~ Sum [1, NGrid] weights[i] * g(zeta[i])
     // Where g(zeta[i]) = [(b-a)/2] * sqrt(1-zeta[i]^2) f( [(b-a)/2]*zeta[i] + [(b+a)/2] )
     // and   zeta[i] = gridPts_[i] = cos ( [(2*i-1)*PI/2*NGrid])
     // weights[i] = PI/NGrid
     // Note I have included the "sqrt(1-zeta[i]^2)" of the transformation in the actual weights[i] .
     // Note in c++ i starts from 0, so has been shifted i = i+1 
    double PI=4.0*atan(1.0);
    for(int i = 0; i < this->nPts_; i++) {
      this->gridPts_[i] = cos(( (2.0*(i+1)-1.0)/(2.0*this->nPts_))*PI);
      this->weights_[i] = (PI/this->nPts_)*(sqrt(1-(this->gridPts_[i]*this->gridPts_[i]))) ;
//    the weights are only (PI/this->nPts_), the second term is including the integrand transformation factor in Eq.25.4.38 Abramowitz Handbook     
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

}; // namespace ChronusQ
