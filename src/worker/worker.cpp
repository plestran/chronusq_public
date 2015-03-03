#include "workers.h"

int worker(GlobalMPI *globalMPI) {
//  cout<<"Worker "<<globalMPI->myid<<" is on "<<globalMPI->nodeName<<endl;

  Molecule  *molecule=new Molecule();

//  molecule->mpiRecv(0);

  return 1;
};


