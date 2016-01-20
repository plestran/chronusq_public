template<typename T>
void QuasiNewton<T>::metBiOrth(TMap &A, const TMap &Met){
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
}

template<>
void QuasiNewton<double>::eigSrt(RealMap &V, RealVecMap &E){
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

