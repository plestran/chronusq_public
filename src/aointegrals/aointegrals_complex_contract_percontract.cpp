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
  void AOIntegrals::Restricted34Contract(ComplexMatrix &G, const ComplexMatrix &X, int n1, 
    int n2, int n3, int n4, int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double* buff, 
    double deg){
    
    for(int i = 0, bf1 = bf1_s, ijkl = 0 ; i < n1; ++i, ++bf1) {
      for(int j = 0, bf2 = bf2_s; j < n2; ++j, ++bf2) {
        for(int k = 0, bf3 = bf3_s; k < n3; ++k, ++bf3) {
          for(int l = 0, bf4 = bf4_s; l < n4; ++l, ++bf4, ++ijkl) {
            double v = buff[ijkl]*deg;
  
            // Coulomb
            this->Gen34Contract(G,X,bf1,bf2,bf3,bf4,v);
  
            // Exchange
            this->Gen23Contract(G,X,bf1,bf2,bf3,bf4,v,0.25);
          }
        }
      }
    }
  } // Restricted34Contract 

  template<>
  void AOIntegrals::UnRestricted34Contract(ComplexMatrix &GAlpha, const ComplexMatrix &XAlpha, 
    ComplexMatrix &GBeta, const ComplexMatrix &XBeta, const ComplexMatrix &XTotal, int n1, int n2, 
    int n3, int n4, int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double* buff, 
    double deg){

    for(int i = 0, bf1 = bf1_s, ijkl = 0 ; i < n1; ++i, ++bf1) {
      for(int j = 0, bf2 = bf2_s; j < n2; ++j, ++bf2) {
        for(int k = 0, bf3 = bf3_s; k < n3; ++k, ++bf3) {
          for(int l = 0, bf4 = bf4_s; l < n4; ++l, ++bf4, ++ijkl) {
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
  } // UnRestricted34Contract

  template<>
  void AOIntegrals::Spinor34Contract(ComplexMatrix &G, const ComplexMatrix &X, int n1, int n2, 
    int n3, int n4, int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double* buff, 
    double deg){

    for(int i = 0, bf1 = bf1_s, ijkl = 0 ; i < n1; ++i, bf1 += 2) {
      for(int j = 0, bf2 = bf2_s; j < n2; ++j, bf2 += 2) {
        for(int k = 0, bf3 = bf3_s; k < n3; ++k, bf3 += 2) {
          for(int l = 0, bf4 = bf4_s; l < n4; ++l, bf4 += 2, ++ijkl) {
            double v = buff[ijkl]*deg;
  
            // Coulomb
          //this->Gen34Contract(G,X,bf1,bf2,bf3,bf4,v);
          //this->Gen34Contract(G,X,bf1+1,bf2+1,bf3+1,bf4+1,v);
            this->GenCouContractSpinor(G,X,bf1,bf2,bf3,bf4,v); 

            // Exchange
/*
            this->Gen23Contract(G,X,bf1,bf2,bf3,bf4,v,0.5);
            this->Gen23Contract(G,X,bf1+1,bf2+1,bf3+1,bf4+1,v,0.5);
//          this->Gen23Contract(G,X,bf1+1,bf2,bf3+1,bf4,v,0.5);
//          this->Gen23Contract(G,X,bf1,bf2+1,bf3,bf4+1,v,0.5);
*/
            this->GenExchContractSpinor(G,X,bf1,bf2,bf3,bf4,v,0.5);
          }
        }
      }
    }
  } // Spinor34Contract 
  
  template<>
  void AOIntegrals::General24CouContract(ComplexMatrix &G, const ComplexMatrix &X, int n1, int n2, 
    int n3, int n4, int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double * buff, 
    double deg){

    for(int i = 0, bf1 = bf1_s, ijkl = 0 ; i < n1; ++i, ++bf1) {
      for(int j = 0, bf2 = bf2_s; j < n2; ++j, ++bf2) {
        for(int k = 0, bf3 = bf3_s; k < n3; ++k, ++bf3) {
          for(int l = 0, bf4 = bf4_s; l < n4; ++l, ++bf4, ++ijkl) {
            double v = buff[ijkl]*deg;
            this->Gen24Contract(G,X,bf1,bf2,bf3,bf4,v);
          }
        }
      }
    }
  } // General24CouContract

  template<>
  void AOIntegrals::Spinor24CouContract(ComplexMatrix &G, const ComplexMatrix &X, int n1, int n2, 
    int n3, int n4, int bf1_s, int bf2_s, int bf3_s, int bf4_s, const double * buff, 
    double deg){

    for(int i = 0, bf1 = bf1_s, ijkl = 0 ; i < n1; ++i, bf1 += 2) {
      for(int j = 0, bf2 = bf2_s; j < n2; ++j, bf2 += 2) {
        for(int k = 0, bf3 = bf3_s; k < n3; ++k, bf3 += 2) {
          for(int l = 0, bf4 = bf4_s; l < n4; ++l, bf4 += 2, ++ijkl) {
            double v = buff[ijkl]*deg;
            this->Spinor24Contract(G,X,bf1,bf2,bf3,bf4,v);
          }
        }
      }
    }
  } // Spinor24CouContract

} // namespace ChronusQ
