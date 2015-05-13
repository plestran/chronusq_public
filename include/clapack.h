/**********************************************************************************************************
 *	C++ declaration of LAPACK drivers
 **********************************************************************************************************/
#include <complex>

extern "C" {
void dgemm_(const char*,const char*,const int*,const int*,const int*,double*,double*,const int*,double*,const int*,double*,double*,const int*);
void dgeev_(char*,char*,int*,double*,int*,double*,double*,double*,int*,double*,int*,double*,int*,int*);
void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);
void dlaswp_(int*,double*,int*,int*,int*,int*,int*);
void zgemm_(const char*,const char*,const int*,const int*,const int*,dcomplex*,dcomplex*,const int*,dcomplex*,const int*,dcomplex*,dcomplex*,const int*);
void zgeev_(char*,char*,int*,dcomplex*,int*,dcomplex*,dcomplex*,int*,dcomplex*,int*,dcomplex*,int*,double*,int*);
void zgesv_(int*,int*,dcomplex*,int*,int*,dcomplex*,int*,int*);

void dsymm_(char*,char*,const int*,const int*,double*,double*,const int*,double*,const int*,double*,double*,const int*);
void zsymm_(char*,char*,const int*,const int*,dcomplex*,dcomplex*,const int*,dcomplex*,const int*,dcomplex*,dcomplex*,const int*);
void zhemm_(char*,char*,const int*,const int*,dcomplex*,dcomplex*,const int*,dcomplex*,const int*,dcomplex*,dcomplex*,const int*);
void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
void zheev_(char*,char*,int*,dcomplex*,int*,double*,dcomplex*,int*,double*,int*);
void pkconv_(int*,char*,int*,double*,int*,double*);
void eigsrt_(int*,double*,double*,double*,double*,int*,int*);
void dgegv_(char*,char*,int*,double*,int*,double*,int*,double*,double*,double*,double*,int*,double*,int*,double*,int*,int*);
void dsygv_(int*,char*,char*,int*,double*,int*,double*,int*,double*,double*,int*,int*);
}
