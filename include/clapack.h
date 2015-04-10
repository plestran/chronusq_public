/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explictly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
/**********************************************************************************************************
 *	C++ declaration of LAPACK drivers
 **********************************************************************************************************/
#include <complex>
//using namespace std;
typedef std::complex<double> dcomplex;

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
