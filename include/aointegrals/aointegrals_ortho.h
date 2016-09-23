template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 double>::value,int>::type>
void AOIntegrals::Ortho1Trans(Op& op1, Op& op2){
  double* Scratch = 
    this->memManager_->malloc<double>(op1.size());
  RealMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.noalias() = (*this->ortho1_) * op1;
  op2.noalias() = SCR * this->ortho1_->transpose();
  

  this->memManager_->free(Scratch,op1.size());
}

template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 dcomplex>::value,int>::type>
void AOIntegrals::Ortho1Trans(Op& op1, Op& op2){
  dcomplex* Scratch = 
    this->memManager_->malloc<dcomplex>(op1.size());
  ComplexMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.real() = (*this->ortho1_) * op1.real();
  SCR.imag() = (*this->ortho1_) * op1.imag();
  op2.real() = SCR.real() * this->ortho1_->transpose();
  op2.imag() = SCR.imag() * this->ortho1_->transpose();
  

  this->memManager_->free(Scratch,op1.size());
}

template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 double>::value,int>::type>
void AOIntegrals::Ortho1TransT(Op& op1, Op& op2){
  double* Scratch = 
    this->memManager_->malloc<double>(op1.size());
  RealMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.noalias() = this->ortho1_->transpose() * op1;
  op2.noalias() = SCR * (*this->ortho1_);
  

  this->memManager_->free(Scratch,op1.size());
}

template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 dcomplex>::value,int>::type>
void AOIntegrals::Ortho1TransT(Op& op1, Op& op2){
  dcomplex* Scratch = 
    this->memManager_->malloc<dcomplex>(op1.size());
  ComplexMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.real() = this->ortho1_->transpose() * op1.real();
  SCR.imag() = this->ortho1_->transpose() * op1.imag();
  op2.real() = SCR.real() * (*this->ortho1_);
  op2.imag() = SCR.imag() * (*this->ortho1_);
  

  this->memManager_->free(Scratch,op1.size());
}


template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 double>::value,int>::type>
void AOIntegrals::Ortho2Trans(Op& op1, Op& op2){
  double* Scratch = 
    this->memManager_->malloc<double>(op1.size());
  RealMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.noalias() = this->ortho2_->transpose() * op1;
  op2.noalias() = SCR * (*this->ortho2_);
  

  this->memManager_->free(Scratch,op1.size());
}

template<typename Op,
  typename std::enable_if<
    std::is_same<typename Op::Scalar,
                 dcomplex>::value,int>::type>
void AOIntegrals::Ortho2Trans(Op& op1, Op& op2){
  dcomplex* Scratch = 
    this->memManager_->malloc<dcomplex>(op1.size());
  ComplexMap SCR(Scratch,op1.rows(),op1.rows());

  SCR.real() = this->ortho2_->transpose() * op1.real();
  SCR.imag() = this->ortho2_->transpose() * op1.imag();
  op2.real() = SCR.real() * (*this->ortho2_);
  op2.imag() = SCR.imag() * (*this->ortho2_);
  

  this->memManager_->free(Scratch,op1.size());
}

