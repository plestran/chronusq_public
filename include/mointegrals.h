#ifndef INCLUDE_MOINTS
#define INCLUDE_MOINTS

#include <global.h>
#include <cerr.h>
#include <wavefunction.h>
#include <pscf.h>

namespace ChronusQ {
template <typename T>
class MOIntegrals {
  typedef Tensor<T,Range4d> TTensor4d; 
  typedef Tensor<T,Range2d> TTensor2d; 
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMatrix;
  typedef Eigen::Map<TMatrix> TMap;
  typedef Eigen::Matrix<T,Dynamic,1,ColMajor> TVec;

  CQMemManager    *memManager_;
  WaveFunction<T> *wfn_; ///< WaveFunction object for MOs
  DummyPostSCF<T>      pscf_; ///< Dummy PSCF Object for dimensions


  T* VVVV_;
  T* VVVO_;
  T* VVOO_;
  T* VOVO_;
  T* VOOO_;
  T* OOOO_;

  T* VVVVAAAA_;
  T* VVVOAAAA_;
  T* VVOOAAAA_;
  T* VOVOAAAA_;
  T* VOOOAAAA_;
  T* OOOOAAAA_;

  T* VVVVAABB_;
  T* VVVOAABB_;
  T* VVOOAABB_;
  T* VOVOAABB_;
  T* VOOOAABB_;
  T* OOOOAABB_;

  T* VVOOBBAA_;
  T* VOOOBBAA_;
  T* VVVOBBAA_;

  T* VVVVBBBB_;
  T* VVVOBBBB_;
  T* VVOOBBBB_;
  T* VOVOBBBB_;
  T* VOOOBBBB_;
  T* OOOOBBBB_;
  // Private functions
  void getLocalMOs();

  bool haveMOVVVV_;
  bool haveMOVVVO_;
  bool haveMOVVOO_;
  bool haveMOVOVO_;
  bool haveMOVOOO_;
  bool haveMOOOOO_;
  bool haveLocMO_;

public:

  MOIntegrals() :
    wfn_(NULL),
    memManager_(NULL) ,
    haveMOVVVV_(false),
    haveMOVVVO_(false),
    haveMOVVOO_(false),
    haveMOVOVO_(false),
    haveMOVOOO_(false),
    haveMOOOOO_(false),
    haveLocMO_ (false)
    { ; };

  
  inline void communicate(WaveFunction<T> &wfn, CQMemManager &memManager){
    this->wfn_        = &wfn;
    this->memManager_ = &memManager;
  }

  inline void initMeta() {
    this->pscf_.communicate(*this->wfn_,*this->memManager_);
    this->pscf_.initMeta();
//  this->getLocalMOs();
  }

  void formVVVV();
  void formVVVO();
  void formVVOO();
  void formVOVO();
  void formVOOO();
  void formOOOO();

  void form2CVVVV();
  void form2CVVVO();
  void form2CVVOO();
  void form2CVOVO();
  void form2CVOOO();
  void form2COOOO();

  void formFullVVVV();
  void formFullVVVO();
  void formFullVVOO();
  void formFullVOVO();
  void formFullVOOO();
  void formFullOOOO();

  void testMOInts();

  T* VVVV() const { return this->VVVV_; };
  T* VVVO() const { return this->VVVO_; };
  T* VVOO() const { return this->VVOO_; };
  T* VOVO() const { return this->VOVO_; };
  T* VOOO() const { return this->VOOO_; };
  T* OOOO() const { return this->OOOO_; };

  T* VVVVAAAA() const { return this->VVVVAAAA_; };
  T* VVVOAAAA() const { return this->VVVOAAAA_; };
  T* VVOOAAAA() const { return this->VVOOAAAA_; };
  T* VOVOAAAA() const { return this->VOVOAAAA_; };
  T* VOOOAAAA() const { return this->VOOOAAAA_; };
  T* OOOOAAAA() const { return this->OOOOAAAA_; };
                                             
  T* VVVVAABB() const { return this->VVVVAABB_; };
  T* VVVOAABB() const { return this->VVVOAABB_; };
  T* VVOOAABB() const { return this->VVOOAABB_; };
  T* VOVOAABB() const { return this->VOVOAABB_; };
  T* VOOOAABB() const { return this->VOOOAABB_; };
  T* OOOOAABB() const { return this->OOOOAABB_; };
                                             
  T* VVOOBBAA() const { return this->VVOOBBAA_; };
  T* VOOOBBAA() const { return this->VOOOBBAA_; };
  T* VVVOBBAA() const { return this->VVVOBBAA_; };
                                             
  T* VVVVBBBB() const { return this->VVVVBBBB_; };
  T* VVVOBBBB() const { return this->VVVOBBBB_; };
  T* VVOOBBBB() const { return this->VVOOBBBB_; };
  T* VOVOBBBB() const { return this->VOVOBBBB_; };
  T* VOOOBBBB() const { return this->VOOOBBBB_; };
  T* OOOOBBBB() const { return this->OOOOBBBB_; };
};

template <typename T, typename U, typename V>
void rank4w2Contract(int INDX, T* A, int LDA, 
  U* B, int LDB1, int LDB2, int LDB3, int LDB4,
  V* C, int LDX) {

  if(INDX == 1) {
    for(auto l = 0; l < LDB4; l++)
    for(auto k = 0; k < LDB3; k++)  
    for(auto j = 0; j < LDB2; j++)
    for(auto x = 0; x < LDX ; x++)
    for(auto i = 0; i < LDB1; i++)
      C[x + j*LDX + k*LDX*LDB2 + l*LDX*LDB2*LDB3] +=
        A[i + x*LDA] *
        B[i + j*LDB1 + k*LDB1*LDB2 + l*LDB1*LDB2*LDB3];
  } else if(INDX == 2) {
    for(auto l = 0; l < LDB4; l++)
    for(auto k = 0; k < LDB3; k++)  
    for(auto x = 0; x < LDX ; x++)
    for(auto j = 0; j < LDB2; j++)
    for(auto i = 0; i < LDB1; i++)
      C[i + x*LDB1 + k*LDX*LDB1 + l*LDX*LDB1*LDB3] +=
        A[j + x*LDA] *
        B[i + j*LDB1 + k*LDB1*LDB2 + l*LDB1*LDB2*LDB3];
  } else if(INDX == 3) {
    for(auto l = 0; l < LDB4; l++)
    for(auto x = 0; x < LDX ; x++)
    for(auto k = 0; k < LDB3; k++)  
    for(auto j = 0; j < LDB2; j++)
    for(auto i = 0; i < LDB1; i++)
      C[i + j*LDB1 + x*LDB1*LDB2 + l*LDX*LDB1*LDB2] +=
        A[k + x*LDA] *
        B[i + j*LDB1 + k*LDB1*LDB2 + l*LDB1*LDB2*LDB3];
  } else if(INDX == 4) {
    for(auto x = 0; x < LDX ; x++)
    for(auto l = 0; l < LDB4; l++)
    for(auto k = 0; k < LDB3; k++)  
    for(auto j = 0; j < LDB2; j++)
    for(auto i = 0; i < LDB1; i++)
      C[i + j*LDB1 + k*LDB1*LDB2 + x*LDB1*LDB2*LDB3] +=
        A[l + x*LDA] *
        B[i + j*LDB1 + k*LDB1*LDB2 + l*LDB1*LDB2*LDB3];
  }
};
};
#endif
