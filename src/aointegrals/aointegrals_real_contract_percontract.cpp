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
   *  Given a shell quartet block of ERIs (row-major), the following contraction is performed
   *
   *  G(μ,v) = (μ v | λ σ) X(σ,λ) - 0.5 (μ σ | λ v) X(σ,λ)
   *
   *  This assumes that G is actually an arbitrary spin block of G and that X is actually
   *  X-Total, i.e. X = X-Alpha + X-Beta (hence the 0.5 on the exchange part)
   */
  template<>
  void AOIntegrals::Restricted34HerContract(RealMatrix &G, const RealMatrix &X, int n1, 
    int n2, int n3, int n4, int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double* buff, 
    double deg){

    for(int i = 0, ijkl = 0 ; i < n1; ++i) {
      int bf1 = bf1_s + i;
      for(int j = 0; j < n2; ++j) {
        int bf2 = bf2_s + j;
        for(int k = 0; k < n3; ++k) {
          int bf3 = bf3_s + k;
          for(int l = 0; l < n4; ++l, ++ijkl) {
            int bf4 = bf4_s + l;
            double v = buff[ijkl]*deg;
  
            // Coulomb
            this->Gen34Contract(G,X,bf1,bf2,bf3,bf4,v);
  
            // Exchange
            this->Gen23Contract(G,X,bf1,bf2,bf3,bf4,v,0.25);
          }
        }
      }
    }
  } // Restricted34HerContract 

  template<>
  void AOIntegrals::UnRestricted34HerContract(RealMatrix &GAlpha, const RealMatrix &XAlpha, 
    RealMatrix &GBeta, const RealMatrix &XBeta, const RealMatrix &XTotal, int n1, int n2, 
    int n3, int n4, int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double* buff, 
    double deg){

    for(int i = 0, ijkl = 0 ; i < n1; ++i) {
      int bf1 = bf1_s + i;
      for(int j = 0; j < n2; ++j) {
        int bf2 = bf2_s + j;
        for(int k = 0; k < n3; ++k) {
          int bf3 = bf3_s + k;
          for(int l = 0; l < n4; ++l, ++ijkl) {
            int bf4 = bf4_s + l;
            double v = buff[ijkl]*deg;
  
            // Coulomb
            this->Gen34Contract(GAlpha,XTotal,bf1,bf2,bf3,bf4,v);
            this->Gen34Contract(GBeta, XTotal,bf1,bf2,bf3,bf4,v);
  
            // Exchange
            this->Gen23Contract(GAlpha,XAlpha,bf1,bf2,bf3,bf4,v,0.5);
            this->Gen23Contract(GBeta, XBeta, bf1,bf2,bf3,bf4,v,0.5);
          }
        }
      }
    }
  } // UnRestricted34HerContract
  
  template<>
  void AOIntegrals::General24CouContract(RealMatrix &G, const RealMatrix &X, int n1, int n2, 
    int n3, int n4, int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double * buff, 
    double deg){

    for(int i = 0, ijkl = 0 ; i < n1; ++i) {
      int bf1 = bf1_s + i;
      for(int j = 0; j < n2; ++j) {
        int bf2 = bf2_s + j;
        for(int k = 0; k < n3; ++k) {
          int bf3 = bf3_s + k;
          for(int l = 0; l < n4; ++l, ++ijkl) {
            int bf4 = bf4_s + l;
            double v = buff[ijkl]*deg;
            this->Gen24Contract(G,X,bf1,bf2,bf3,bf4,v);
          }
        }
      }
    }
  } // General24CouContract

}
