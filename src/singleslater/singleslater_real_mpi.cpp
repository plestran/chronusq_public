/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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
#include <singleslater.h>
namespace ChronusQ {
/*************************/
/* MPI Related Routines  */
/*************************/
template<>
void SingleSlater<double>::mpiSend(int toID,int tag) {
  //OOMPI_COMM_WORLD[toID].Send(this->nAtoms_,tag);
  //OOMPI_COMM_WORLD[toID].Send(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiSend(toID,tag);
};
template<>
void SingleSlater<double>::mpiRecv(int fromID,int tag) {
  //OOMPI_COMM_WORLD[fromID].Recv(this->nAtoms_,tag);
  //this->index_=new int[this->nAtoms_];
  //this->cart_ =new RealMatrix(3, this->nAtoms_, "Molecule");
  //OOMPI_COMM_WORLD[fromID].Recv(this->index_,this->nAtoms_,tag);
  //this->cart_->mpiRecv(fromID,tag);
};
} // namespace ChronusQ

