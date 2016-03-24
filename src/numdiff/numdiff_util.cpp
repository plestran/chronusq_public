/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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
#include <numdiff.h>
namespace ChronusQ{
  template<>
  void NumericalDifferentiation<double>::checkPhase(RealMatrix &A, 
    RealMatrix &B){

    if(A.cols() != B.cols() || A.rows() != B.rows())
      CErr("Phase Checking only works with matricies of the same dimension",
      this->singleSlater_undisplaced_->fileio()->out);

    for(auto j = 0; j < A.cols(); j++){
      int sgn1(1), sgn2(1);
      for(auto i = 0; i < A.rows(); i++){
        if(std::abs(A(i,j)) > 1e-8) {
          if(std::abs(B(i,j)) < 1e-8) continue;
          else { // check B
            sgn1 = A(i,j) / std::abs(A(i,j));
            sgn2 = B(i,j) / std::abs(B(i,j));
            break;
          } // get signs
        } // check A
      } // loop i
      if(sgn1 != sgn2) 
        B.col(j) *= -1.0;
    } // loop j
  };

  template<>
  void NumericalDifferentiation<double>::checkPhase(RealMatrix &M1,
     RealMatrix &M2, RealMatrix &Inner_1_2){
  
    double tol = 1e-1;
    RealMatrix O_1_2(Inner_1_2);
  
    for(auto I = 0; I < Inner_1_2.rows(); I++)
    for(auto J = 0; J < Inner_1_2.cols(); J++){
      if(std::abs(O_1_2(I,J)) < tol) O_1_2(I,J) = 0.0;
      else if(O_1_2(I,J) > 0.0)      O_1_2(I,J) = 1.0;
      else                           O_1_2(I,J) = -1.0;
    }

    for(auto I = 0; I < Inner_1_2.rows(); I++)
    for(auto J = 0; J < Inner_1_2.cols(); J++){
      if(I != J && O_1_2(I,J) != 0.0) O_1_2(I,J) = 0.0;
    }
  
    RealMatrix TMP = M2 * O_1_2;
    M2 = TMP;
  };

  template<>
  void NumericalDifferentiation<double>::checkPhase(SingleSlater<double> &ss1,
    SingleSlater<double> &ss2, RealMatrix &SMO_1_2){
  
    this->checkDegeneracies(ss1);
    this->checkDegeneracies(ss2);

/*
    double tol = 1e-1;
    RealMatrix O_1_2(SMO_1_2);
  
    for(auto mu = 0; mu < SMO_1_2.rows(); mu++)
    for(auto nu = 0; nu < SMO_1_2.rows(); nu++){
      if(std::abs(O_1_2(mu,nu)) < tol) O_1_2(mu,nu) = 0.0;
      else if(O_1_2(mu,nu) > 0.0)      O_1_2(mu,nu) = 1.0;
      else                             O_1_2(mu,nu) = -1.0;
    }
  
//  prettyPrint(cout,O_1_2,"O");
//  prettyPrint(cout,SMO_1_2,"SMO");
    RealMatrix TMP = (*ss2.moA()) * O_1_2;
    (*ss2.moA()) = TMP;
*/
    this->checkPhase((*ss1.moA()),(*ss2.moA()),SMO_1_2);
  };

  template<>
  void NumericalDifferentiation<double>::checkPhase(Response<double> &resp1,
    Response<double> &resp2){
  
    RealMatrix T1,T2;
    if(this->respType_ == RESPONSE_TYPE::CIS){
      T1 = resp1.transDen<SINGLETS>().block(0,0,
          resp1.nMatDim<SINGLETS>(),this->responseNRoots_);
      T2 = resp2.transDen<SINGLETS>().block(0,0,
          resp2.nMatDim<SINGLETS>(),this->responseNRoots_);
    } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
      T1 = resp1.transDen<A_PPTDA_SINGLETS>().block(0,0,
          resp1.nMatDim<A_PPTDA_SINGLETS>(),this->responseNRoots_);
      T2 = resp2.transDen<A_PPTDA_SINGLETS>().block(0,0,
          resp2.nMatDim<A_PPTDA_SINGLETS>(),this->responseNRoots_);
    }

      cout << "  Checking | T - T' | Before Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " 
           << diffNorm(T1,T2) 
           << endl;  
    this->checkPhase(T1,T2);
      cout << "  Checking | T - T' | After Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " 
           << diffNorm(T1,T2) 
           << endl;  

    if(this->respType_ == RESPONSE_TYPE::CIS){
      resp1.transDen<SINGLETS>().block(0,0,
        resp1.nMatDim<SINGLETS>(),this->responseNRoots_) = T1;
      resp2.transDen<SINGLETS>().block(0,0,
        resp2.nMatDim<SINGLETS>(),this->responseNRoots_) = T2;
    } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
      resp1.transDen<A_PPTDA_SINGLETS>().block(0,0,
        resp1.nMatDim<A_PPTDA_SINGLETS>(),this->responseNRoots_) = T1;
      resp2.transDen<A_PPTDA_SINGLETS>().block(0,0,
        resp2.nMatDim<A_PPTDA_SINGLETS>(),this->responseNRoots_) = T2;
    }

  };

}
