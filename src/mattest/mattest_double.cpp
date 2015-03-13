#include "matrix.h"

void double_test(int N, int M){
  // Allocate Matricies
  Matrix<double> *A = new Matrix<double>(N,N,"A","STD");
  Matrix<double> *B = new Matrix<double>(N,M,"B","STD");
  Matrix<double> *C = new Matrix<double>(N,1,"C","STD");
  Matrix<double> *D = new Matrix<double>(M,1,"D","STD");

  Matrix<double> *AB = new Matrix<double>(N,M,"AB = A*B","STD");
  Matrix<double> *AC = new Matrix<double>(N,1,"AC = A*C","STD");
  Matrix<double> *BD = new Matrix<double>(N,1,"BD = B*D","STD");

  Matrix<double> *A_Vec = new Matrix<double>(N,N,"A Eigenvector","STD");
  Matrix<double> *A_Eig = new Matrix<double>(N,N,"A Eigenvalues","STD");

  Matrix<double> *AT  = new Matrix<double>(N,N,"Back-transformed -> A","STD");
  Matrix<double> *AET = new Matrix<double>(N,N,"Back-transformed -> A_Eig","STD");

  Matrix<double> *SqrtA = new Matrix<double>(N,N,"[A]^2","STD");
  Matrix<double>  *ExpA = new Matrix<double>(N,N,"Exp[A]","STD");

  Matrix<double> *VcT = new Matrix<double>(N,N,"VcT");
  Matrix<double> *tmp = new Matrix<double>(N,N,"tmp");

  // Clear Contents
  A->setSymm('S');
  A->clearAll();
  B->clearAll();
  C->clearAll();
  D->clearAll();
  A_Eig->clearAll();
  A_Vec->clearAll();
  VcT->clearAll();
  tmp->clearAll();

  // Initialize Test Matricies
  for(int i = 0; i < N; i++)
  for(int j = 0; j < N; j++){
    (*A)(i,j) = i + j;
  }
  A->scale(0.01);

  for(int i = 0; i < N; i++)
  for(int j = 0; j < M; j++){
    (*B)(i,j) = i * j;
  }

  
  for(int i = 0; i < N; i++){
    (*C)(i,0) = i + (double)1.0/(i+1);
  }

  for(int i = 0; i < M; i++){
    (*D)(i,0) = i + (double)1.0/(i*i+1);
  }

  // Perform Operations
  (*AB) = (*A)*(*B);
  (*AC) = (*A)*(*C);
  (*BD) = (*B)*(*D);
  
  
//A->setSymm('G');
  try{ A->diag();}
  catch(int msg){ cout << msg << endl;}
  
  BD->printAll();
//try{B->diag();}
//catch(int msg) {CErr(msg);}; 
  
  try{
  (*A_Vec) =     A->eigenvector();
    A_Eig->setDag(A->eigenvalue());  
  } catch (int msg){ cout << msg << endl;};
  try{ (*AT) = A_Eig->transNT(*A_Vec);}
  catch(int msg){ cout << msg << endl;}
  (*AET) = A->transTN(*A_Vec);
  A->printAll();
  AT->printAll();
  A_Eig->printAll();
  AET->printAll();
  (*SqrtA) = (*A)^2.0;
  (*ExpA) = exp(*A);
  SqrtA->printAll();
  ExpA->printAll();
  A_Vec->printAll();
  (*VcT) = A_Vec;
  VcT->transpose();
  (*tmp) = (*VcT)*(*A_Vec);
  tmp->printAll();
  
  delete A;
  delete B;
  delete C;
  delete A_Vec;
  delete A_Eig;
  delete SqrtA;
  delete ExpA;
}
