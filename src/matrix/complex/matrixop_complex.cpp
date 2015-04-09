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
