#ifndef INCLUDED_PSCF
#define INCLUDED_PSCF
#include <quantum.h>
#include <wavefunction.h>

template <typename T>
class PostSCF : public Quantum<T>{
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMatrix;
  typedef Eigen::Map<TMatrix> TMap;
  typedef Eigen::Matrix<T,Dynamic,1,ColMajor> TVec;

protected:
  WaveFunction<T> *reference_;
  AOIntegrals     *aointegrals_;
  FileIO          *fileio_;

  // Single Particle Dims
  int nO_;
  int nV_;
  int nOA_;
  int nVA_;
  int nOB_;
  int nVB_;
  int nBasis_;
  
  // Composite Particle Dims
  int nOO_;
  int nOV_;
  int nVV_;
  int nOAOA_;
  int nOAOB_;
  int nOBOB_;
  int nOAVA_;
  int nOAVB_;
  int nOBVA_;
  int nOBVB_;
  int nVAVA_;
  int nVAVB_;
  int nVBVB_;

  int nOO_SLT_;
  int nVV_SLT_;
  int nOAOA_SLT_;
  int nOBOB_SLT_;
  int nVAVA_SLT_;
  int nVBVB_SLT_;

  int nOO_LT_;
  int nVV_LT_;
  int nOAOA_LT_;
  int nOBOB_LT_;
  int nVAVA_LT_;
  int nVBVB_LT_;

public:

  PostSCF() : Quantum<T>(),
    reference_   (NULL),
    aointegrals_ (NULL),
    fileio_      (NULL){ ; };

  PostSCF(PostSCF &other):
    Quantum<T>(dynamic_cast<Quantum<T>&>(other)),
    reference_   (other.reference_  ),
    aointegrals_ (other.aointegrals_),
    fileio_      (other.fileio_     ){ 
      this->initMeta(); 
  };

  inline void communicate(WaveFunction &wfn) {
    this->reference_ = &wfn;
  };

  inline void initMeta() {
    this->aointegrals_ = this->reference_->aointegrals();
    this->fileio_      = this->reference_->fileio();

    this->nOA_ = this->reference_->nOA();
    this->nOB_ = this->reference_->nOB();
    this->nVA_ = this->reference_->nVA();
    this->nVB_ = this->reference_->nVB();
    this->nO_  = this->reference_->nO();
    this->nV_  = this->reference_->nV();

    this->nBasis_ = this->reference_->nBasis();

    this->nOO_   = this->nO_  * this->nO_;
    this->nOAOA_ = this->nOA_ * this->nOA_;
    this->nOAOB_ = this->nOA_ * this->nOB_;
    this->nOBOB_ = this->nOB_ * this->nOB_;
    this->nVV_   = this->nV_  * this->nV_;
    this->nVAVA_ = this->nVA_ * this->nVA_;
    this->nVAVB_ = this->nVA_ * this->nVB_;
    this->nVBVB_ = this->nVB_ * this->nVB_;
    this->nOV_   = this->nO_  * this->nV_;
    this->nOAVA_ = this->nOA_ * this->nVA_;
    this->nOAVB_ = this->nOA_ * this->nVB_;
    this->nOBVA_ = this->nOB_ * this->nVA_;
    this->nOBVB_ = this->nOB_ * this->nVB_;

    this->nOO_SLT_   = this->nO_  * (this->nO_  - 1)/2;
    this->nOAOA_SLT_ = this->nOA_ * (this->nOA_ - 1)/2;
    this->nOBOB_SLT_ = this->nOB_ * (this->nOB_ - 1)/2;
    this->nVV_SLT_   = this->nV_  * (this->nV_  - 1)/2;
    this->nVAVA_SLT_ = this->nVA_ * (this->nVA_ - 1)/2;
    this->nVBVB_SLT_ = this->nVB_ * (this->nVB_ - 1)/2;

    this->nOO_LT_   = this->nO_  * (this->nO_  - 1)/2;
    this->nOAOA_LT_ = this->nOA_ * (this->nOA_ - 1)/2;
    this->nOBOB_LT_ = this->nOB_ * (this->nOB_ - 1)/2;
    this->nVV_LT_   = this->nV_  * (this->nV_  - 1)/2;
    this->nVAVA_LT_ = this->nVA_ * (this->nVA_ - 1)/2;
    this->nVBVB_LT_ = this->nVB_ * (this->nVB_ - 1)/2;
  };

}

#endif
