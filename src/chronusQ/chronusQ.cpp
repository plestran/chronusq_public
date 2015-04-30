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
#include <global.h>
#include <fileio.h>
#include <molecule.h>
#include <atoms.h>
#include <basisset.h>
#include <aointegrals.h>
#include <workers.h>
using namespace ChronusQ;

int main(int argc,char *argv[]) {
  std::chrono::high_resolution_clock::time_point start,finish;
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
    start = std::chrono::high_resolution_clock::now();
    std::string str(argv[1]);
    atlas(argc, str, &globalMPI);
//  } else {
//    worker(&globalMPI);
//  };

//  if(globalMPI.myid==0) {
    finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> JobT = finish - start;
    cout<<"The total job running time is:  "<< JobT.count() <<" seconds."<<endl;
//  }
//  OOMPI_COMM_WORLD.Finalize();  

  return 0;
};


