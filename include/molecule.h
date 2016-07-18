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
  std::vector<libint2::Shell>  finiteWidthNuclei_;

  // Misc options
  int                          printLevel_;  // Level of print
public:

  // constructor
  Molecule(){ this->loadDefaults();};

  Molecule(Atoms atm, std::ostream &out) : Molecule(){
    this->nAtoms_ = 1;
    this->alloc(out);
    auto n = HashAtom(atm.symbol,atm.massNumber);
    if(n != -1) index_[0] = n;
    else
      CErr("Error: invalid atomic symbol or mass number!",out);
    this->nTotalE_ = atm.atomicNumber;
    this->computeRij();
    this->toCOM(0);
    this->computeI();
  }

  ~Molecule(){
    delete[] index_;
  };

 
  inline void loadDefaults(){
    this->nAtoms_          = 0;
    this->charge_          = 0; 
    this->multip_          = 0;
    this->nTotalE_         = 0;
    this->index_           = NULL;
    this->energyNuclei_    = 0.0;
    this->cart_            = nullptr;
    this->COM_             = nullptr;
    this->momentOfInertia_ = nullptr;
    this->rIJ_             = nullptr;

    this->printLevel_      = 1;
  }
  
  void alloc(std::ostream &out=cout);
  void toCOM(int Iop0);    
  void computeI();
  void computeRij();
  void computeNucRep();
  void printInfo(std::ostream & out=cout );
  void generateFiniteWidthNuclei();
  inline void convBohr(){ (*this->cart_) /= phys.bohr; }

  // Python API
  void Wrapper_printInfo(FileIO &);
  void Wrapper_alloc(FileIO &);
  boost::python::list Wrapper_cart();

  // Reference Access
  inline int& index(int i) { return this->index_[i];};
  inline int  atomicZ(int i) { return elements[index_[i]].atomicNumber;};

  // Getters
  inline int nAtoms() {return this->nAtoms_;};
  inline int charge() {return this->charge_;}
  inline int multip() {return this->multip_;}
  inline int nTotalE() {return this->nTotalE_;};
  inline int printLevel(){return this->printLevel_;};

  inline double energyNuclei() { return this->energyNuclei_;};

  inline RealMatrix* cart() {return this->cart_.get();}
  inline RealMatrix* rIJ() {return this->rIJ_.get();}
  inline VectorXd*   COM() {return this->COM_.get();}
  inline std::vector<std::array<double,3>> cartArray() {
    std::vector<std::array<double,3>> atomicCenters;
 
    for(auto iAtm = 0; iAtm < this->nAtoms(); iAtm++){
      atomicCenters.push_back(
          {(*this->cart())(0,iAtm),
           (*this->cart())(1,iAtm),
           (*this->cart())(2,iAtm)}
      );
    };

    return atomicCenters;
  };

  inline libint2::Shell& nucShell(int i){ return this->finiteWidthNuclei_[i];};

  // Setters
  inline void setCharge(int i) {this->charge_ = i; this->nTotalE_ -= i;};
  inline void setMultip(int i) {this->multip_ = i;};
  inline void setNAtoms(int i) {this->nAtoms_ = i;};
  inline void setIndex(int i, int n){ this->index_[i] = n; };
  inline void setCart(int i, double x, double y, double z){
    (*this->cart_)(0,i) = x;
    (*this->cart_)(1,i) = y;
    (*this->cart_)(2,i) = z;
  }
  inline void setNTotalE(int i){   this->nTotalE_    = i;};
  inline void setPrintLevel(int i){this->printLevel_ = i;};



  // read from input file
  void readMolecule(FileIO *, std::istream &);

  // read|write scratch|binary files
  void ioRead(FileIO *);
  void ioWrite(FileIO *);

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int);
  void mpiRecv(int,int);
};
} // namespace ChronusQ
#endif
