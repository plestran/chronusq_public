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
template<typename T>
void QuasiNewton2<T>::metBiOrth(TMap &A, const TMap &Met){
/*
  int N = A.cols();
  TMap AX(this->ASuperMem_,Met.rows(),N);
  AX = Met*A;
  for(auto i = 0; i < N; i++){
    double inner = 0;
    for(auto k = 0; k < N; k++)
      inner += std::norm(AX(k,i));
    int sgn = inner / std::abs(inner);
    inner = sgn*std::sqrt(sgn*inner);
    A.col(i) /= inner;
    for(auto j = i + 1; j < N; j++){
      A.col(j) -= A.col(i) * A.col(i).dot(AX.col(j));
    } 
  }
*/

/*
  int N = A.cols();
  TMap AX(this->ASuperMem_,Met.rows(),N);
  AX = Met*A;
  double inner = A.col(0).dot(AX.col(0));
  int sgn = inner / std::abs(inner);
  inner = sgn * std::sqrt(sgn*inner);
  A.col(0) /= inner;

  for(auto i = 1; i < N; i++){
    for(auto j = 0; j < i; j++) {
      inner
    }
  }
*/

  int N = A.cols();

  auto metInner = [&] (int i, int j) -> double {
    double inner = 0;
    for(auto k = 0; k < N; k++)
    for(auto l = 0; l < N; l++) {
      inner += A(k,i) * Met(k,l) * A(l,j);
    }
    return inner;
  };

  auto metNorm = [&] (int i) -> double {
    double inner = metInner(i,i);
    int sgn = inner / std::abs(inner);
    inner = sgn * std::sqrt(sgn*inner);
    return inner;
  };

  A.col(0) /= metNorm(0);
  for(auto i = 1; i < N; i++) {
    for(auto j = 0; j < i; j++) {
      A.col(i) -= metInner(i,j) * A.col(j);
    }
    A.col(i) /= metNorm(i);
  }
  
}

template <typename T>
void QuasiNewton2<T>::eigSrt(TMap &V, RealVecMap &E){
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

