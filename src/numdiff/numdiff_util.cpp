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
  
    RealMatrix T1 = resp1.transDen()[0].block(0,0,resp1.nMatDim()[0],
      this->responseNRoots_);
    RealMatrix T2 = resp2.transDen()[0].block(0,0,resp2.nMatDim()[0],
      this->responseNRoots_);

      cout << "  Checking | T - T' | Before Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " 
           << diffNorm(T1,T2) 
           << endl;  
    this->checkPhase(T1,T2);
      cout << "  Checking | T - T' | After Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " 
           << diffNorm(T1,T2) 
           << endl;  

    resp1.transDen()[0].block(0,0,resp1.nMatDim()[0],this->responseNRoots_) = 
      T1;
    resp2.transDen()[0].block(0,0,resp2.nMatDim()[0],this->responseNRoots_) = 
      T2;

  };

}
