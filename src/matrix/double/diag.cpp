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

using ChronusQ::Matrix;
char UPLO = 'L';
int LWORK;
double *CPY,*WORK
// Check if it's already been diagonalized
// Check if it's square

this->cleanEigen();
this->allocEigen();
if(this->symm_=='G') CPY=this->allocScratch();
LWORK = this->eigWORK(1,UPLO);
WORK = new (nothrow) double *[LWORK];
if(WORK==NULL) throw 3105;

if(this->symm_=='G'){
  dgeev_(&this->JOBVL_,this->JOBVR_,&this->rows_,CPY,&this->rows_,
    this->eigenvalue_re_,this->eigenvalue_im,this->eigenvector_l_;&this->rows_,
    this->eigenvector_r_,&this->rows_,WORK,&LWORK,&INFO);}
}
else if(this->symm_'S'){
  dsyev_(&this->JOBVR_,&UPLO,&this->rows_,this->eigenvector_,&this->rows_,
    this->eigenvalue_,WORK,&LWORK,&INFO);
};

delete[] CPY; delete[] WORK;
this->haveEigen_=1;
}
