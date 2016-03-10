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
void zhegv_(int*,char*,char*,int*,dcomplex*,int*,dcomplex*,int*,double*,dcomplex*,int*,double*,int*);
void dgeqrf_(int *, int *, double *, int *, double *, double *, int *, int *);
void dorgqr_(int *, int *, int *, double *, int *, double *, double *, int *, int *);
void zgeqrf_(int *, int *, dcomplex *, int *, dcomplex *, dcomplex *, int *, int *);
void zungqr_(int *, int *, int *, dcomplex *, int *, dcomplex *, dcomplex *, int *, int *);
void dgesvd_(char*,char*,int*,int*,double*,int*,double*,double*,int*,double*,int*,double*,int*,int*);
void dgetrf_(int*,int*,double*,int*,int*,int*);
void dgetri_(int*,double*,int*,int*,double*,int*,int*);
void zgetrf_(int*,int*,dcomplex*,int*,int*,int*);
void zgetri_(int*,dcomplex*,int*,int*,dcomplex*,int*,int*);
void dpotrf_(char*,int*,double*,int*,int*); 
void dpotri_(char*,int*,double*,int*,int*); 
void dgecon_(char*,int*,double*,int*,double*,double*,double*,int*,int*);
void dgels_(char*,int*,int*,int*,double*,int*,double*,int*,double*,int*,int*);
void dgelss_(int*,int*,int*,double*,int*,double*,int*,double*,double*,int*,
  double*,int*,int*);
}

