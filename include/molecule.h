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
#ifndef INCLUDED_MOLECULES
#define INCLUDED_MOLECULES
#include "global.h"
#include "fileio.h"
#include "atoms.h"
#include "matrix.h"
#include "controls.h"

/****************************/
/* Error Messages 2000-2999 */
/****************************/

namespace ChronusQ {
class Molecule {
  int      nAtoms_;      // number of atoms in the system
  int      charge_;      // total charge
  int      spin_;        // spin multiplicity
  int      nTotalE_;     // total number of electrons
  int      size_;        // size of the object in terms of sizeof(char)
  int     *index_;       // index of atom in the atoms[] array
  double   energyNuclei_;// nuclear repulsion energy
  Matrix<double>  *cart_;        // cartesian coordinates

public:

  // constructor
  Molecule(int nAtoms=0,FileIO *fileio=NULL){ if(nAtoms>0) iniMolecule(nAtoms,fileio);};
  ~Molecule(){
    delete[] index_;
    delete   cart_;
  };
  void iniMolecule(int,FileIO*);

  // print
  void printInfo(FileIO*,Controls*);

  // access to private data
  inline int index(int i) { return this->index_[i];};
  inline int nAtoms() {return this->nAtoms_;};
  inline Matrix<double> *cart() {return this->cart_;}
  inline int charge() {return this->charge_;}
  inline int spin() {return this->spin_;}
  inline int nTotalE() {return this->nTotalE_;};
  inline void readCharge(int charge) {this->charge_=charge;};
  inline void readSpin(int spin) {this->spin_=spin;};
  inline int size() { return this->size_;};
  inline double energyNuclei() { return this->energyNuclei_;};

  // read from input file
  void readMolecule(FileIO*);

  // read|write scratch|binary files
  void ioRead(FileIO*);
  void ioWrite(FileIO*);

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagMolecule);
  void mpiRecv(int,int tag=tagMolecule);
};
} // namespace ChronusQ
#endif
