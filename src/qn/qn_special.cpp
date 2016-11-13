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
