#ifndef INCLUDED_PSCF
#define INCLUDED_PSCF

template <typename T>
class PostSCF : {
  // Single Particle Dims
  int nO_;
  int nV_;
  int nOA_;
  int nVA_;
  int nOB_;
  int nVB_;
  
  // Composite Particle Dims
  int nOO_;
  int nOV_;
  int nVV_;
  int nOAOA_;
  int nOAOB_;
  int nOBOB_;
  int nOAVA_;
  int nOAVB_;
  int nOBVB_;
  int nVAVA_;
  int nVAVB_;
  int nVBVB_;

  int nOO_SLT_;
  int nOV_SLT_;
  int nVV_SLT_;
  int nOAOA_SLT_;
  int nOAOB_SLT_;
  int nOBOB_SLT_;
  int nOAVA_SLT_;
  int nOAVB_SLT_;
  int nOBVB_SLT_;
  int nVAVA_SLT_;
  int nVAVB_SLT_;
  int nVBVB_SLT_;

}

#endif
