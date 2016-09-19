#include <qn.h>
namespace ChronusQ {
template<>
void QuasiNewton2<double>::invertSuperMetric(const int NTrial) {
  int N = 2*NTrial;
  int LWORK = N*N;
  int INFO;

  int * iPiv = this->memManager_->template malloc<int>(N+1);
  double * WORK = this->memManager_->template malloc<double>(LWORK);

  dgetrf_(&N,&N,this->SSuperMem_,&N,iPiv,&INFO);
  dgetri_(&N,this->SSuperMem_,&N,iPiv,WORK,&LWORK,&INFO);

  this->memManager_->free(iPiv,N+1);
  this->memManager_->free(WORK,LWORK);
};

};
