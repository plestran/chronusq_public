template<typename T>
bool QuasiNewton2<T>::stdHermetianDiag(char JOBV, char UPLO, int N, T *A, 
  double *E) {

  int INFO;
  this->LWORK = 3*N;
  this->WORK = this->memManager_->template malloc<T>(this->LWORK);

  if(std::is_same<double,T>::value) {
    dsyev_(&JOBV,&UPLO,&N,reinterpret_cast<double*>(A),&N,E,
      reinterpret_cast<double*>(this->WORK),&this->LWORK,&INFO);
  } else if(std::is_same<dcomplex,T>::value) {

    this->RWORK_ = 
      this->memManager_->template malloc<double>(std::max(1,3*N-2));

    zheev_(&JOBV,&UPLO,&N,reinterpret_cast<dcomplex*>(A),&N,E,
      reinterpret_cast<dcomplex*>(this->WORK),&this->LWORK,this->RWORK_,
      &INFO);
    this->memManager_->free(this->WORK,std::max(1,3*N-2));

  }

  this->memManager_->free(this->WORK,this->LWORK);
  return (INFO == 0);
}

