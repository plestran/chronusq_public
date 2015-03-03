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
#endif
