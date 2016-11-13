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



/*
 * ChronusQ Plugin for Eigen
 */

template<typename OtherDerived>
inline Scalar frobInner(const MatrixBase<OtherDerived>& other) const {
  return (derived().cwiseProduct(other.derived())).sum();
}

template<typename Dervd>
inline double diffNorm( const Dervd &A, const Dervd &B ){
  return (A-B).norm();
}

template<typename Dervd>
inline double diffNormI( const Dervd &A ){
  return (A - Dervd::Identity(A.rows(),A.cols())).norm();
}

template<typename Dervd>
inline double selfInner( const Dervd &A){
  return A.dot(A);
}

inline void printMATLAB(std::ostream &out) const {
  out << "[";
  for(auto i = 0; i < derived().rows(); i++) {
//  out << "[";
    for(auto j = 0; j < derived().cols(); j++){
      if(std::is_same<std::complex<double>,Scalar>::value)
        out << "complex";
      out << derived()(i,j);
      if( j != derived().cols() - 1) out << " , ";
    }
    out << ";" << std::endl;
//  out << "]";
  }
  out << "]";
};
