#include "global.h"
#include "fileio.h"
#include "molecule.h"
#include "atoms.h"
#include "matrix.h"
#include "basisset.h"
#include "aointegrals.h"
#include "workers.h"
#include "cerr.h"

void double_test(int,int);
void complex_test(int,int);

int main(int argc,char *argv[]) {
/*
  clock_t    start,finish;
  cout.precision(8);
  cout.fill(' ');
  cout.setf(ios::right,ios::adjustfield);
  cout.setf(ios::fixed,ios::floatfield);

  //Initialize MPI
  GlobalMPI globalMPI;
//  OOMPI_COMM_WORLD.Init(argc, argv);
//  globalMPI.myid = OOMPI_COMM_WORLD.Rank();
//  globalMPI.size = OOMPI_COMM_WORLD.Size();
//  MPI_Get_processor_name(globalMPI.nodeName,&(globalMPI.nodeNameLen));

//  if(globalMPI.myid==0) {
    start = clock();
    controller(argc, argv, &globalMPI);
//  } else {
//    worker(&globalMPI);
//  };

//  if(globalMPI.myid==0) {
    finish = clock();
    //cout<<"The total job running time is:"<<(finish-start)/CLOCKS_PER_SEC<<" seconds."<<endl;
//  }
//  OOMPI_COMM_WORLD.Finalize();  
*/
  cout << setw(25) << left << "Performing Tests for:"<< left <<setw(15) <<"N = "<< setw(4) << left << 5 <<setw(15)<<left <<" M = "<<setw(4)<<left<<7  << "  ... ";
  double_test(5,7);
  cout << "Done" << endl;
  cout << setw(25) << left << "Performing Tests for:"<< left <<setw(15) <<"N = "<< setw(4) << left << 10 <<setw(15)<<left <<" M = "<<setw(4)<<left<<12  << "  ... ";
  double_test(10,12);
  cout << "Done" << endl;
  cout << setw(25) << left << "Performing Tests for:"<< left <<setw(15) <<"N = "<< setw(4) << left << 50 <<setw(15)<<left <<" M = "<<setw(4)<<left<<105  << "  ... ";
  double_test(50,105);
  cout << "Done" << endl;
  cout << setw(25) << left << "Performing Tests for:"<< left <<setw(15) <<"N = "<< setw(4) << left << 2 <<setw(15)<<left <<" M = "<<setw(4)<<left<<105  << "  ... ";
  double_test(2,105);
  cout << "Done" << endl;
  cout << "HERE IN MAIN" << endl;
  complex_test(5,5);
/*
  int n = 2;
  int m;
  m = 3;
  srand(time(NULL));
  Matrix<double> *A = new Matrix<double>(n*n,n*n,"A","STD");
  Matrix<double> *B = new Matrix<double>(n,n,"B","STD");
  Matrix<double> *C = new Matrix<double>(n*n,1), *D;
  A->clearAll();
  (*A)(0,0) = 6.0;
  (*A)(1,1) = 6.0;
  (*A)(0,1) = 4.0;
  (*A)(1,0) = 4.0;
  (*A)(2,2) = 6.0;
  (*A)(0,2) = 4.0;
  (*A)(2,0) = 4.0;
  (*A)(3,3) = 6.0;
  (*A)(0,3) = 4.0;
  (*A)(3,0) = 4.0;
  A->scaleDag(2.0);
  A->scale(0.002);
  for(int i=0;i<n;i++){
    for(int j=i;j<n;j++){
      (*A)(i,j) = (double)rand()/RAND_MAX;
      (*A)(j,i) = (*A)(i,j); 
    }
  }
//A->setSymm('S');
//B->setName("B mat");
  (*B)(0,0) = 2.5;
  (*B)(1,0) = 1.0;
  (*B)(1,1) = 2.5;
  (*B)(0,1) = 1.0;
  B->vectorize();
//A->printAll();
//C = (*A)*(*((*A)*(*A)));
//B->printAll();
  (*C) = (*A)*(*B);
  C->setName("C");
//C->printAll();
//C->~Matrix<double>();
//C = (*B)*(*A);
//C->setName("C");
//C->printAll();
//D = (*C)*double(5.3);
//D->printAll();
  Matrix<double> *I = new Matrix<double>(n*n,n*n);
  for(int i = 0; i < n*n; i++){ (*I)(i,i) = 1.0;};
  clock_t start,finish;
  A->setSymm('S');
//I->setSymm('S');
  start = clock();
  I->printAll();
  A->printAll();
//A->diag(I);
  try{ A->diag();}
  catch(int msg){ CErr(msg);};
  A->eSort();
  finish = clock();
//cout<<"The total job running time is:"<<setprecision(8)<<(double)(finish-start)/CLOCKS_PER_SEC<<" seconds."<<endl;
//try{ A->diag(I);}
//catch(int msg){ CErr(msg);};
  Matrix<double> *E = new Matrix<double>(n*n,n*n);
  Matrix<double> *F = new Matrix<double>(n*n,n*n);
  F->clearAll();
  E->clearAll();
  (*E)=(A->eigenvector());
//(*E1)=(A->eigenvector_l());
  F->setDag( A->eigenvalue());
  E->setName("A Eigenvector");
  F->setName("A Eigenvalues");
  E->printAll();
  F->printAll();
  Matrix<double> *G = new Matrix<double>(n*n,n*n);
  Matrix<double> *H = new Matrix<double>(n*n,n*n);
  Matrix<double> *K = new Matrix<double>(n*n,n*n);
  (*H) = A->transTN((*E));
  H->setName("Trans Eig");
  H->printAll();
  (*K) = exp(*A);
  K->setName("Exp");
  K->printAll();
*/
  return 0;
};

