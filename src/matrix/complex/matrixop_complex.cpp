/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
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
#include "global.h"
#include "matrix.h"
using Chronus::Matrix;

namespace ChronusQ {
/*
int dotProd(Matrix *a, Matrix *b, Matrix *c) {
  PetscMPIInt rank,size;
  PetscInt    namelen;
  char        processor_name[MPI_MAX_PROCESSOR_NAME];
  Mat A,B,C;
  PetscInt i, j;
  PetscErrorCode ierr;

  PetscInitialize(&argc,&argv,(char *)0,help);
  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,a->nRow()*a->nCol(),a->nRow()*a->nCol());CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);

  for (i=0; i<a->nRow(); i++) 
    for (j=0; j<a->nCol(); j++) {
      ierr = MatSetValues(A,1,&i,1,&j,&(*a)(i,j),INSERT_VALUES); CHKERRQ(ierr);
    };

  MatCreate(PETSC_COMM_WORLD, &B);
  MatSetSizes(A,b->nRow(),b->nCol(),b->nRow(),b->nCol());

  MatCreate(PETSC_COMM_WORLD, &C);
  MatSetSizes(A,c->nRow(),c->nCol(),c->nRow(),c->nCol());

  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Get_processor_name(processor_name,&namelen);
  cout<<"Process "<<rank<<" of "<<size<<" is on "<<processor_name<<endl;

  PetscFinalize();
  return 0;
};
*/

Matrix* iden(const int n, const int m){
  Matrix* m = new Matrix(n,m,"Identity","STD");
  for(int i = 0; i < m; i++){ m->data_[i+n*i] = 1.0;};
  return m;
};
Matrix* tBasis(const Matrix *C, const Matrix *X){
  Matrix *prod = new Matrix(X->rows_,X->cols_);
  double *scr = new double[->]
};
} // namespace ChronusQ
