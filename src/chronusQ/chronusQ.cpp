#include "global.h"
#include "fileio.h"
#include "molecule.h"
#include "atoms.h"
#include "matrix.h"
#include "basisset.h"
#include "aointegrals.h"
#include "workers.h"
using namespace ChronusQ;

int main(int argc,char *argv[]) {
  clock_t start,finish;
  cout.precision(6);
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
    atlas(argc, argv, &globalMPI);
//  } else {
//    worker(&globalMPI);
//  };

//  if(globalMPI.myid==0) {
    finish = clock();
    //cout<<"The total job running time is:"<<(finish-start)/CLOCKS_PER_SEC<<" seconds."<<endl;
//  }
//  OOMPI_COMM_WORLD.Finalize();  

  return 0;
};


