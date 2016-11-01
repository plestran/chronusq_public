template <typename T>
void MOIntegrals<T>::formVVVV() {
  if(this->haveMOVVVV_) return;

  bool is2C        = this->wfn_->nTCS() == 2;
  bool isOpenShell = !this->wfn_->isClosedShell;
  int NB           = this->wfn_->nBasis();
  int intDim       = this->pscf_.nVAVA() * this->pscf_.nVAVA();
  int NVA = this->wfn_->nVA();
  int NOA = this->wfn_->nOA();
  int NVB = this->wfn_->nVB();
  int NOB = this->wfn_->nOB();

  if(is2C) {
    form2CVVVV();
    return;
  }

  /// Evaluate VVVV (AAAA) MO Integrals

  // Allocate and zero out VVVV (AAAA) storage
  VVVVAAAA_ = this->memManager_->template malloc<T>(NVA*NVA*NVA*NVA);
  std::fill_n(VVVVAAAA_,NVA*NVA*NVA*NVA,0.0);

  // Allocate and zero out Itermediates
  T * I_anls_A   = this->memManager_->template malloc<T>(NB*NB*NB*NVA);
  T * I_abls_AA  = this->memManager_->template malloc<T>(NB*NB*NVA*NVA);
  T * I_abcs_AAA = this->memManager_->template malloc<T>(NB*NVA*NVA*NVA);

  std::fill_n(I_anls_A  ,NB*NB*NB*NVA,0.0);
  std::fill_n(I_abls_AA ,NB*NB*NVA*NVA,0.0);
  std::fill_n(I_abcs_AAA,NB*NVA*NVA*NVA,0.0);

  // Conjugate for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // First Quarter (AAAA) transformation (mn|ls) -> (an|ls) (A)
  rank4w2Contract(1,this->wfn_->moA()->data() + NOA*NB,NB,
    &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
    I_anls_A,NVA);

  // Conjugate back for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Second Quarter (AAAA) transformation (an|ls) (A) -> (ab|ls) (AA)
  rank4w2Contract(2,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_anls_A,NVA,NB,NB,NB,
    I_abls_AA,NVA);

  // Conjugate for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Third Quarter (AAAA) transformation (ab|ls) (AA) -> (ab|cs) (AAA)
  rank4w2Contract(3,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_abls_AA,NVA,NVA,NB,NB,
    I_abcs_AAA,NVA);

  // Conjugate back for contraction if complex
  if(std::is_same<dcomplex,T>::value) {
    (*this->wfn_->moA()) = this->wfn_->moA()->conjugate();
  }

  // Fourth Quarter (AAAA) transformation (ab|cs) (AAA) -> (ab|cd) (AAAA)
  rank4w2Contract(4,this->wfn_->moA()->data() + NOA*NB,NB,  
    I_abcs_AAA,NVA,NVA,NVA,NB,
    VVVVAAAA_,NVA);

  this->memManager_->free(I_abcs_AAA,NB*NVA*NVA*NVA);
  this->memManager_->free(I_anls_A  ,NB*NB*NB*NVA);

  if(isOpenShell) {
    /// Evaluate VVVV (AABB) MO Integrals

    // Allocate and zero out VVVV (AABB) storage
    VVVVAABB_ = this->memManager_->template malloc<T>(NVA*NVA*NVB*NVB);
    std::fill_n(VVVVAABB_,NVA*NVA*NVB*NVB,0.0);

    // Allocate and zero out Itermediates
    T * I_abcs_AAB = this->memManager_->template malloc<T>(NB*NVA*NVA*NVB);

    std::fill_n(I_abcs_AAB,NB*NVA*NVA*NVB,0.0);


    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // Third Quarter (AABB) transformation (ab|ls) (AA) -> (ab|cs) (AAB)
    rank4w2Contract(3,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_abls_AA,NVA,NVA,NB,NB,
      I_abcs_AAB,NVB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
 
    // Fourth Quarter (AABB) transformation (ab|cs) (AAB) -> (ab|cd) (AABB)
    rank4w2Contract(4,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_abcs_AAB,NVA,NVA,NVB,NB,
      VVVVAABB_,NVB);
    
    this->memManager_->free(I_abcs_AAB,NB*NVA*NVA*NVB);
  }

  this->memManager_->free(I_abls_AA ,NB*NB*NVA*NVA);

  if(isOpenShell) {
    /// Evaluate VVVV (BBBB) MO Integrals
 
    // Allocate and zero out VVVV (BBBB) storage
    VVVVBBBB_ = this->memManager_->template malloc<T>(NVB*NVB*NVB*NVB);
    std::fill_n(VVVVBBBB_,NVB*NVB*NVB*NVB,0.0);
 
    // Allocate and zero out Itermediates
    T * I_anls_B   = this->memManager_->template malloc<T>(NB*NB*NB*NVB);
    T * I_abls_BB  = this->memManager_->template malloc<T>(NB*NB*NVB*NVB);
    T * I_abcs_BBB = this->memManager_->template malloc<T>(NB*NVB*NVB*NVB);
 
    std::fill_n(I_anls_B  ,NB*NB*NB*NVB,0.0);
    std::fill_n(I_abls_BB ,NB*NB*NVB*NVB,0.0);
    std::fill_n(I_abcs_BBB,NB*NVB*NVB*NVB,0.0);
 
    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // First Quarter (BBBB) transformation (mn|ls) -> (an|ls) (B)
    rank4w2Contract(1,this->wfn_->moB()->data() + NOB*NB,NB,
      &this->wfn_->aointegrals()->aoERI_->storage()[0],NB,NB,NB,NB,
      I_anls_B,NVB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
 
    // Second Quarter (BBBB) transformation (an|ls) (B) -> (ab|ls) (BB)
    rank4w2Contract(2,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_anls_B,NVB,NB,NB,NB,
      I_abls_BB,NVB);
 
    // Conjugate for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }

    // Third Quarter (BBBB) transformation (ab|ls) (BB) -> (ab|cs) (BBB)
    rank4w2Contract(3,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_abls_BB,NVB,NVB,NB,NB,
      I_abcs_BBB,NVB);

    // Conjugate back for contraction if complex
    if(std::is_same<dcomplex,T>::value) {
      (*this->wfn_->moB()) = this->wfn_->moB()->conjugate();
    }
 
    // Fourth Quarter (BBBB) transformation (ab|cs) (BBB) -> (ab|cd) (BBBB)
    rank4w2Contract(4,this->wfn_->moB()->data() + NOB*NB,NB,  
      I_abcs_BBB,NVB,NVB,NVB,NB,
      VVVVBBBB_,NVB);
 
    this->memManager_->free(I_abcs_BBB,NB*NVB*NVB*NVB);
    this->memManager_->free(I_abls_BB ,NB*NB*NVB*NVB);
    this->memManager_->free(I_anls_B  ,NB*NB*NB*NVB);
  }

  this->haveMOVVVV_ = true;
} // formVVVV

template <typename T> 
void MOIntegrals<T>::form2CVVVV() {

  int NB = this->wfn_->nBasis();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  VVVV_ = this->memManager_->template malloc<T>(NV*NV*NV*NV);
  std::fill_n(VVVV_,NV*NV*NV*NV,0.);


  for(auto d = 0; d < NV; d++)
  for(auto c = 0; c < NV; c++)
  for(auto b = 0; b < NV; b++)
  for(auto a = 0; a < NV; a++)

  for(auto sg = 0; sg < 2*NB; sg+=2)
  for(auto lm = 0; lm < 2*NB; lm+=2)
  for(auto nu = 0; nu < 2*NB; nu+=2)
  for(auto mu = 0; mu < 2*NB; mu+=2){

    VVVV_[a + b*NV + c*NV*NV + d*NV*NV*NV] +=
      (
      std::conj((*wfn_->moA())(mu,a+NO)) * (*wfn_->moA())(nu,b+NO) *
      std::conj((*wfn_->moA())(lm,c+NO)) * (*wfn_->moA())(sg,d+NO) +

      std::conj((*wfn_->moA())(mu+1,a+NO)) * (*wfn_->moA())(nu+1,b+NO) *
      std::conj((*wfn_->moA())(lm,c+NO))   * (*wfn_->moA())(sg,d+NO) +

      std::conj((*wfn_->moA())(mu,a+NO))   * (*wfn_->moA())(nu,b+NO) *
      std::conj((*wfn_->moA())(lm+1,c+NO)) * (*wfn_->moA())(sg+1,d+NO) +

      std::conj((*wfn_->moA())(mu+1,a+NO)) * (*wfn_->moA())(nu+1,b+NO) *
      std::conj((*wfn_->moA())(lm+1,c+NO)) * (*wfn_->moA())(sg+1,d+NO) 
      ) * (*wfn_->aointegrals()->aoERI_)(mu/2,nu/2,lm/2,sg/2);
  }

  this->haveMOVVVV_ = true;
}; // form2CVVVV

template <typename T>
void MOIntegrals<T>::formFullVVVV(){
  this->formVVVV();
  if(this->wfn_->nTCS() == 2) return;

  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();

  this->VVVV_ = this->memManager_->template malloc<T>(NV*NV*NV*NV);
  std::fill_n(this->VVVV_,NV*NV*NV*NV,0.0);

  for(auto D = 0; D < NV; D+=2)
  for(auto C = 0; C < NV; C+=2)
  for(auto B = 0; B < NV; B+=2)
  for(auto A = 0; A < NV; A+=2){
    int a = A / 2;
    int b = B / 2;
    int c = C / 2;
    int d = D / 2;

    this->VVVV_[A + B*NV + C*NV*NV + D*NV*NV*NV] = 
      this->VVVVAAAA_[a + b*NVA + c*NVA*NVA + d*NVA*NVA*NVA];

    this->VVVV_[(A+1) + (B+1)*NV + C*NV*NV + D*NV*NV*NV] = 
      this->VVVVAABB_[c + d*NVA + a*NVA*NVA + b*NVA*NVA*NVB];

    this->VVVV_[A + B*NV + (C+1)*NV*NV + (D+1)*NV*NV*NV] = 
      this->VVVVAABB_[a + b*NVA + c*NVA*NVA + d*NVA*NVA*NVB];

    this->VVVV_[(A+1) + (B+1)*NV + (C+1)*NV*NV + (D+1)*NV*NV*NV] = 
      this->VVVVBBBB_[a + b*NVB + c*NVB*NVB + d*NVB*NVB*NVB];
  }
} // formFullVVVV

