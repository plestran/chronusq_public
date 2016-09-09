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

  CQMemManager    *memManager_;
  WaveFunction<T> *wfn_; ///< WaveFunction object for MOs
  DummyPostSCF<T>      pscf_; ///< Dummy PSCF Object for dimensions

/*
  // Spin-Orbital MO Integral Storage
  std::unique_ptr<TTensor4d> VVVV_;
  std::unique_ptr<TTensor4d> VVVO_;
  std::unique_ptr<TTensor4d> VVOO_;
  std::unique_ptr<TTensor4d> VOVO_;
  std::unique_ptr<TTensor4d> VOOO_;
  std::unique_ptr<TTensor4d> OOOO_;

  // Spacial-Orbital MO Integral Storage
  std::unique_ptr<TTensor4d> VAVAVAVA_;
  std::unique_ptr<TTensor4d> VAVAVAOA_;
  std::unique_ptr<TTensor4d> VAVAOAOA_;
  std::unique_ptr<TTensor4d> VAOAVAOA_;
  std::unique_ptr<TTensor4d> VAOAOAOA_;
  std::unique_ptr<TTensor4d> OAOAOAOA_;

  std::unique_ptr<TTensor4d> VAVAVBVB_;
  std::unique_ptr<TTensor4d> VAVAVBOB_;
  std::unique_ptr<TTensor4d> VAVAOBOB_;
  std::unique_ptr<TTensor4d> VAOAVBOB_;
  std::unique_ptr<TTensor4d> VAOAOBOB_;
  std::unique_ptr<TTensor4d> OAOAOBOB_;

  std::unique_ptr<TTensor4d> VBVBVBVB_;
  std::unique_ptr<TTensor4d> VBVBVBOB_;
  std::unique_ptr<TTensor4d> VBVBOBOB_;
  std::unique_ptr<TTensor4d> VBOBVBOB_;
  std::unique_ptr<TTensor4d> VBOBOBOB_;
  std::unique_ptr<TTensor4d> OBOBOBOB_;


  // Local copies of MO coefficients
  std::unique_ptr<TTensor2d> locMOAOcc_;
  std::unique_ptr<TTensor2d> locMOBOcc_;
  std::unique_ptr<TTensor2d> locMOAVir_;
  std::unique_ptr<TTensor2d> locMOBVir_;

  // Local copies of (Complex) MO coefficients
  std::unique_ptr<TTensor2d> reLocMOAOcc_;
  std::unique_ptr<TTensor2d> reLocMOBOcc_;
  std::unique_ptr<TTensor2d> reLocMOAVir_;
  std::unique_ptr<TTensor2d> reLocMOBVir_;
  std::unique_ptr<TTensor2d> imLocMOAOcc_;
  std::unique_ptr<TTensor2d> imLocMOBOcc_;
  std::unique_ptr<TTensor2d> imLocMOAVir_;
  std::unique_ptr<TTensor2d> imLocMOBVir_;
*/

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
    this->getLocalMOs();
  }

  void formVVVV();
  void formVVVO();
  void formVVOO();
  void formVOVO();
  void formVOOO();
  void formOOOO();

  void formFullVVVV();
  void formFullVVVO();
  void formFullVVOO();
  void formFullVOVO();
  void formFullVOOO();
  void formFullOOOO();

  void testMOInts();
};
};
#endif
