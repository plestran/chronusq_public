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
#include <global.h>
#include <cerr.h>
#include <fileio.h>
#include <atoms.h>
#include <controls.h>

/****************************/
/* Error Messages 2000-2999 */
/****************************/

namespace ChronusQ {
class Molecule {
  int                          nAtoms_;      // number of atoms in the system
  int                          charge_;      // total charge
  int                          multip_;      // spin multiplicity
  int                          nTotalE_;     // total number of electrons
  int                         *index_;       // center index in the atoms[] array
  double                       energyNuclei_;// nuclear repulsion energy
  std::unique_ptr<RealMatrix>  cart_;        // cartesian coordinates
  std::unique_ptr<VectorXd>  COM_;         // center of mass coordinate or center of nuclear charges 
  std::unique_ptr<RealMatrix>  momentOfInertia_; // Moment of inertia
  std::unique_ptr<RealMatrix>  rIJ_;             // Interatomic distance matrix
public:

  // constructor
  Molecule(){ 
    this->nTotalE_ = 0;
  };

  Molecule(Atoms atm, std::ostream &out){
    this->iniMolecule(1,out);
    auto n = HashAtom(atm.symbol,atm.massNumber);
    if(n!=-1) index_[0] = n;
    else
      CErr("Error: invalid atomic symbol or mass number!",out);
    nTotalE_ = atm.atomicNumber;
    (*cart_)(0,0) = 0.0;
    (*cart_)(1,0) = 0.0;
    (*cart_)(2,0) = 0.0;
    energyNuclei_ = 0.0;
    this->computeRij();
    this->toCOM(0);
    this->computeI();
  }
  ~Molecule(){
    delete[] index_;
  };
  void iniMolecule(int,std::ostream &out=cout);
//APS Compute center of mass (or center of nuclear charges) of a molecule
  void toCOM(int Iop0);    
//APE
  void computeI();
  void computeRij();
  // print
  void printInfo(std::ostream &);

  // Python API
  void Wrapper_printInfo(FileIO &);

  // access to private data
  inline int index(int i) { return this->index_[i];};
  inline int nAtoms() {return this->nAtoms_;};
  inline RealMatrix* cart() {return this->cart_.get();}
  inline RealMatrix* rIJ() {return this->rIJ_.get();}
  inline int charge() {return this->charge_;}
  inline int multip() {return this->multip_;}
  inline int nTotalE() {return this->nTotalE_;};
  inline void readCharge(int charge) {this->charge_=charge; this->nTotalE_ -= charge;};
  inline void readMultip(int multip) {this->multip_=multip;};
  inline double energyNuclei() { return this->energyNuclei_;};

  // read from input file
  void readMolecule(FileIO *, std::istream &);

  // read|write scratch|binary files
  void ioRead(FileIO *);
  void ioWrite(FileIO *);

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagMolecule);
  void mpiRecv(int,int tag=tagMolecule);
};
} // namespace ChronusQ
#endif
