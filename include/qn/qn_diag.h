/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
template<typename T>
bool QuasiNewton2<T>::stdHermetianDiag(char JOBV, char UPLO, int N, T *A, 
  double *E) {

  int INFO;
  int LWORK = 3*N;
  T* WORK = this->memManager_->template malloc<T>(LWORK);

  if(std::is_same<double,T>::value) {
    dsyev_(&JOBV,&UPLO,&N,reinterpret_cast<double*>(A),&N,E,
      reinterpret_cast<double*>(WORK),&LWORK,&INFO);
  } else if(std::is_same<dcomplex,T>::value) {

    double *RWORK = 
      this->memManager_->template malloc<double>(std::max(1,3*N-2));

    zheev_(&JOBV,&UPLO,&N,reinterpret_cast<dcomplex*>(A),&N,E,
      reinterpret_cast<dcomplex*>(WORK),&LWORK,RWORK,&INFO);
    this->memManager_->free(RWORK,std::max(1,3*N-2));

  }

  this->memManager_->free(WORK,LWORK);
  return (INFO == 0);
}

template<typename T>
bool QuasiNewton2<T>::stdNonHermetianDiag(char JOBVR, char JOBVL, int N, T *A, 
  T *E, T *VR, T *VL) {

  int INFO;
  int LWORK = 4*N;
  T *WORK = this->memManager_->template malloc<T>(LWORK);

  if(std::is_same<double,T>::value) {
    double* ER = reinterpret_cast<double*>(E);
    double* EI = reinterpret_cast<double*>(E) + N;

    dgeev_(&JOBVL,&JOBVR,&N,reinterpret_cast<double*>(A),&N,ER,EI,
      reinterpret_cast<double*>(VL),&N,reinterpret_cast<double*>(VR),&N,
      reinterpret_cast<double*>(WORK),&LWORK,&INFO);

  } else if(std::is_same<dcomplex,T>::value) {

    double *RWORK = 
      this->memManager_->template malloc<double>(std::max(1,2*N));

    zgeev_(&JOBVR,&JOBVL,&N,reinterpret_cast<dcomplex*>(A),&N,
      reinterpret_cast<dcomplex*>(E),
      reinterpret_cast<dcomplex*>(VL),&N,reinterpret_cast<dcomplex*>(VR),&N,
      reinterpret_cast<dcomplex*>(WORK),&LWORK,RWORK,&INFO);

    this->memManager_->free(RWORK,std::max(1,2*N));

  }

  this->memManager_->free(WORK,LWORK);
  return (INFO == 0);
}
