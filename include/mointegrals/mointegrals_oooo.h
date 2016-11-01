template <typename T>
void MOIntegrals<T>::formOOOO() {
  if(this->haveMOOOOO_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nVAVA() * this->pscf_.nVAVA();
  int NVA = this->wfn_->nVA();
  int NOA = this->wfn_->nOA();
  int NVB = this->wfn_->nVB();
  int NOB = this->wfn_->nOB();

  if(is2C) {
    form2COOOO();
    return;
  }

  /// Evaluate OOOO (AAAA) MO Integrals

  // Allocate and zero out OOOO (AAAA) storage
  OOOOAAAA_ = this->memManager_->template malloc<T>(NOA*NOA*NOA*NOA);
  std::fill_n(OOOOAAAA_,NOA*NOA*NOA*NOA,0.0);

  // Allocate and zero out Itermediates
  T * I_1_A   = this->memManager_->template malloc<T>(NB*NB*NB*NOA);
  T * I_2_AA  = this->memManager_->template malloc<T>(NB*NB*NOA*NOA);
  T * I_3_AAA = this->memManager_->template malloc<T>(NB*NOA*NOA*NOA);

  std::fill_n(I_1_A  ,NB*NB*NB*NOA,0.0);
  std::fill_n(I_2_AA ,NB*NB*NOA*NOA,0.0);
  std::fill_n(I_3_AAA,NB*NOA*NOA*NOA,0.0);

  // Conjugate for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // First Quarter (AAAA) transformation (mn|ls) -> (in|ls) (A)
  rank4w2Contract(1,this->wfn_->moA()->data(),NB,
    &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
    I_1_A,NOA);

  // Conjugate back for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Second Quarter (AAAA) transformation (in|ls) (A) -> (ij|ls) (AA)
  rank4w2Contract(2,this->wfn_->moA()->data(),NB,  
    I_1_A,NOA,NB,NB,NB,
    I_2_AA,NOA);

  // Conjugate for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Third Quarter (AAAA) transformation (ij|ls) (AA) -> (ij|ks) (AAA)
  rank4w2Contract(3,this->wfn_->moA()->data(),NB,  
    I_2_AA,NOA,NOA,NB,NB,
    I_3_AAA,NOA);

  // Conjugate back for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Fourth Quarter (AAAA) transformation (ij|ks) (AAA) -> (ij|kl) (AAAA)
  rank4w2Contract(4,this->wfn_->moA()->data(),NB,  
    I_3_AAA,NOA,NOA,NOA,NB,
    OOOOAAAA_,NOA);

  this->memManager_->free(I_3_AAA,NB*NOA*NOA*NOA);
  this->memManager_->free(I_1_A  ,NB*NB*NB*NOA);

  if(isOpenShell) {
    /// Evaluate OOOO (AABB) MO Integrals

    // Allocate and zero out OOOO (AABB) storage
    OOOOAABB_ = this->memManager_->template malloc<T>(NOA*NOA*NOB*NOB);
    std::fill_n(OOOOAABB_,NOA*NOA*NOB*NOB,0.0);

    // Allocate and zero out Itermediates
    T * I_3_AAB = this->memManager_->template malloc<T>(NB*NOA*NOA*NOB);

    std::fill_n(I_3_AAB,NB*NOA*NOA*NOB,0.0);

    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // Third Quarter (AABB) transformation (ij|ls) (AA) -> (ij|ks) (AAB)
    rank4w2Contract(3,this->wfn_->moB()->data(),NB,  
      I_2_AA,NOA,NOA,NB,NB,
      I_3_AAB,NOB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
 
    // Fourth Quarter (AABB) transformation (ij|ks) (AAB) -> (ij|kl) (AABB)
    rank4w2Contract(4,this->wfn_->moB()->data(),NB,  
      I_3_AAB,NOA,NOA,NOB,NB,
      OOOOAABB_,NOB);
    
    this->memManager_->free(I_3_AAB,NB*NOA*NOA*NOB);
  }

  this->memManager_->free(I_2_AA ,NB*NB*NOA*NOA);

  if(isOpenShell) {
    /// Evaluate OOOO (BBBB) MO Integrals
 
    // Allocate and zero out OOOO (BBBB) storage
    OOOOBBBB_ = this->memManager_->template malloc<T>(NOB*NOB*NOB*NOB);
    std::fill_n(OOOOBBBB_,NOB*NOB*NOB*NOB,0.0);
 
    // Allocate and zero out Itermediates
    T * I_1_B   = this->memManager_->template malloc<T>(NB*NB*NB*NOB);
    T * I_2_BB  = this->memManager_->template malloc<T>(NB*NB*NOB*NOB);
    T * I_3_BBB = this->memManager_->template malloc<T>(NB*NOB*NOB*NOB);
 
    std::fill_n(I_1_B  ,NB*NB*NB*NOB,0.0);
    std::fill_n(I_2_BB ,NB*NB*NOB*NOB,0.0);
    std::fill_n(I_3_BBB,NB*NOB*NOB*NOB,0.0);

    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
 
    // First Quarter (BBBB) transformation (mn|ls) -> (in|ls) (B)
    rank4w2Contract(1,this->wfn_->moB()->data(),NB,
      &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
      I_1_B,NOB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
 
    // Second Quarter (BBBB) transformation (in|ls) (B) -> (ij|ls) (BB)
    rank4w2Contract(2,this->wfn_->moB()->data(),NB,  
      I_1_B,NOB,NB,NB,NB,
      I_2_BB,NOB);

    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
 
    // Third Quarter (BBBB) transformation (ij|ls) (BB) -> (ij|ks) (BBB)
    rank4w2Contract(3,this->wfn_->moB()->data(),NB,  
      I_2_BB,NOB,NOB,NB,NB,
      I_3_BBB,NOB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
 
    // Fourth Quarter (BBBB) transformation (ij|ks) (BBB) -> (ij|kl) (BBBB)
    rank4w2Contract(4,this->wfn_->moB()->data(),NB,  
      I_3_BBB,NOB,NOB,NOB,NB,
      OOOOBBBB_,NOB);
 
    this->memManager_->free(I_3_BBB,NB*NOB*NOB*NOB);
    this->memManager_->free(I_2_BB ,NB*NB*NOB*NOB);
    this->memManager_->free(I_1_B  ,NB*NB*NB*NOB);
  }

  this->haveMOOOOO_ = true;
} // formOOOO

template <typename T> 
void MOIntegrals<T>::form2COOOO() {

  int NB = this->wfn_->nBasis();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  OOOO_ = this->memManager_->template malloc<T>(NO*NO*NO*NO);
  std::fill_n(OOOO_,NO*NO*NO*NO,0.);


  for(auto l = 0; l < NO; l++)
  for(auto k = 0; k < NO; k++)
  for(auto j = 0; j < NO; j++)
  for(auto i = 0; i < NO; i++)

  for(auto sg = 0; sg < 2*NB; sg+=2)
  for(auto lm = 0; lm < 2*NB; lm+=2)
  for(auto nu = 0; nu < 2*NB; nu+=2)
  for(auto mu = 0; mu < 2*NB; mu+=2){

    OOOO_[i + j*NO + k*NO*NO + l*NO*NO*NO] +=
      (
      std::conj((*wfn_->moA())(mu,i)) * (*wfn_->moA())(nu,j) *
      std::conj((*wfn_->moA())(lm,k)) * (*wfn_->moA())(sg,l) +

      std::conj((*wfn_->moA())(mu+1,i)) * (*wfn_->moA())(nu+1,j) *
      std::conj((*wfn_->moA())(lm,k))   * (*wfn_->moA())(sg,l) +

      std::conj((*wfn_->moA())(mu,i))   * (*wfn_->moA())(nu,j) *
      std::conj((*wfn_->moA())(lm+1,k)) * (*wfn_->moA())(sg+1,l) +

      std::conj((*wfn_->moA())(mu+1,i)) * (*wfn_->moA())(nu+1,j) *
      std::conj((*wfn_->moA())(lm+1,k)) * (*wfn_->moA())(sg+1,l) 
      ) * (*wfn_->aointegrals()->aoERI_)(mu/2,nu/2,lm/2,sg/2);
  }

  this->haveMOOOOO_ = true;
}; // form2COOOO

template <typename T>
void MOIntegrals<T>::formFullOOOO(){
  this->formOOOO();
  if(this->wfn_->nTCS() == 2) return;

  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->OOOO_ = this->memManager_->template malloc<T>(NO*NO*NO*NO);
  std::fill_n(this->OOOO_,NO*NO*NO*NO,0.0);

  for(auto D = 0; D < NO; D+=2)
  for(auto C = 0; C < NO; C+=2)
  for(auto B = 0; B < NO; B+=2)
  for(auto A = 0; A < NO; A+=2){
    int a = A / 2;
    int b = B / 2;
    int c = C / 2;
    int d = D / 2;

    this->OOOO_[A + B*NO + C*NO*NO + D*NO*NO*NO] = 
      this->OOOOAAAA_[a + b*NOA + c*NOA*NOA + d*NOA*NOA*NOA];

    this->OOOO_[(A+1) + (B+1)*NO + C*NO*NO + D*NO*NO*NO] = 
      this->OOOOAABB_[c + d*NOA + a*NOA*NOA + b*NOA*NOA*NOB];

    this->OOOO_[A + B*NO + (C+1)*NO*NO + (D+1)*NO*NO*NO] = 
      this->OOOOAABB_[a + b*NOA + c*NOA*NOA + d*NOA*NOA*NOB];

    this->OOOO_[(A+1) + (B+1)*NO + (C+1)*NO*NO + (D+1)*NO*NO*NO] = 
      this->OOOOBBBB_[a + b*NOB + c*NOB*NOB + d*NOB*NOB*NOB];
  }
} // formFullOOOO

