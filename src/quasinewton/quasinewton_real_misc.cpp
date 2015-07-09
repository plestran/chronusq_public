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
#include <quasinewton.h>

namespace ChronusQ {
  template<>
  void QuasiNewton<double>::Orth(RealCMMap &A){
    int N   = A.cols();
    int M   = A.rows();
    int LDA = M;
    int INFO;
    double * AMAT = A.data();
    double * TAU = this->LAPACK_SCR;
    this->WORK = TAU + N;

    dgeqrf_(&M,&N,AMAT,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    dorgqr_(&M,&N,&N,AMAT,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    
  }

  template<>
  void QuasiNewton<double>::Orth(RealCMMatrix &A){
    int N   = A.cols();
    int M   = A.rows();
    int LDA = M;
    int INFO;
    double * AMAT = A.data();
    double * TAU = this->LAPACK_SCR;
    this->WORK = TAU + N;

    dgeqrf_(&M,&N,AMAT,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    dorgqr_(&M,&N,&N,AMAT,&LDA,TAU,this->WORK,&this->LWORK,&INFO);
    
  }

  template<>
  void QuasiNewton<double>::metBiOrth(RealCMMap &A, const RealCMMatrix &Met){
    int N = A.cols();
    RealCMMatrix AX = Met*A;
    for(auto i = 0; i < N; i++){
      double inner = A.col(i).dot(AX.col(i));
      int sgn = inner / std::abs(inner);
      inner = sgn*std::sqrt(sgn*inner);
      A.col(i) /= inner;
      for(auto j = i + 1; j < N; j++){
        A.col(j) -= A.col(i) * A.col(i).dot(AX.col(j));
      } 
    }
  }

  template<>
  void QuasiNewton<double>::eigSrt(RealCMMap &V, RealVecMap &E){
    auto N = V.cols();
    while( N != 0){
      auto newn = 0;
      for(auto i = 1; i < N; i++){
        if( E(i-1) > E(i) ){
          auto tmp = E(i);
          E(i) = E(i-1);
          E(i-1) = tmp;
          V.col(i).swap(V.col(i-1));
          newn = i;
        }
      }
      N = newn;
    }
  }
}; // namespace ChronusQ
