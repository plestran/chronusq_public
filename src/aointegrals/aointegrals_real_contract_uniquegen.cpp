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
#include <aointegrals.h>
namespace ChronusQ{
  /**
   * Forms the unique contributions for a 34 contraction
   * given the permutational symmetry of real ERIs
   *
   * G[X](μ,v) += (μ v | λ σ) X(σ,λ)
   *
   */
  template<>
  void AOIntegrals::Gen34Contract(RealMatrix &G, const RealMatrix &X, int bf1, int bf2, int bf3, int bf4, double v){
    G(bf1,bf2) += X(bf4,bf3)*v;
    G(bf3,bf4) += X(bf2,bf1)*v;
    G(bf2,bf1) += X(bf3,bf4)*v;
    G(bf4,bf3) += X(bf1,bf2)*v;
  } // Gen34Contract
  
  /**
   * Forms the unique contributions for a 23 contraction
   * given the permutational symmetry of real ERIs
   *
   * G[X](μ,v) -= 2*fact*(μ σ | λ v) X(σ,λ)
   *
   */
  template<> 
  void AOIntegrals::Gen23Contract(RealMatrix &G, const RealMatrix &X, int bf1, int bf2, int bf3, int bf4, double v, double fact){
    G(bf1,bf3) -= fact*X(bf2,bf4)*v;
    G(bf2,bf4) -= fact*X(bf1,bf3)*v;
    G(bf1,bf4) -= fact*X(bf2,bf3)*v;
    G(bf2,bf3) -= fact*X(bf1,bf4)*v;
  
    G(bf3,bf1) -= fact*X(bf4,bf2)*v;
    G(bf4,bf2) -= fact*X(bf3,bf1)*v;
    G(bf4,bf1) -= fact*X(bf3,bf2)*v;
    G(bf3,bf2) -= fact*X(bf4,bf1)*v;
  } // Gen23 Contract
  
  /**
   * Forms the unique contributions for a 24 contraction
   * given the permutational symmetry of real ERIs
   *
   * G[X](μ,v) += (μ σ | v λ) X(σ,λ)  (Mulliken)
   *
   * ** or **
   *
   * G[X](μ,v) += <μ v | σ λ> X(σ,λ)  (Dirac)
   *
   */
  template<>
  void AOIntegrals::Gen24Contract(RealMatrix &G, const RealMatrix &X, int bf1, int bf2, int bf3, int bf4, double v){
    G(bf1,bf3) += X(bf2,bf4)*v;
    G(bf2,bf3) += X(bf1,bf4)*v;
    G(bf1,bf4) += X(bf2,bf3)*v;
    G(bf2,bf4) += X(bf1,bf3)*v;
  
    G(bf3,bf1) += X(bf4,bf2)*v; 
    G(bf4,bf1) += X(bf3,bf2)*v;
    G(bf3,bf2) += X(bf4,bf1)*v;
    G(bf4,bf2) += X(bf3,bf1)*v;
  } // Gen24Contract
}; // namespace ChronusQ
