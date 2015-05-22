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
#ifndef INCLUDED_SINGLESLATER
#define INCLUDED_SINGLESLATER
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <controls.h>
#include <aointegrals.h>

/****************************/
/* Error Messages 5000-5999 */
/****************************/

namespace ChronusQ {
class SingleSlater {
  int      nBasis_;
  int      nTT_;
  int      nAE_;
  int      nBE_;
  int      RHF_;
  int      nOccA_;
  int      nOccB_;
  int      nVirA_;
  int      nVirB_;
  int      spin_;
  int    **R2Index_;
  std::unique_ptr<RealMatrix>  densityA_;
  std::unique_ptr<RealMatrix>  densityB_;
  std::unique_ptr<RealMatrix>  fockA_;
  std::unique_ptr<RealMatrix>  fockB_;
  std::unique_ptr<RealMatrix>  coulombA_;
  std::unique_ptr<RealMatrix>  coulombB_;
  std::unique_ptr<RealMatrix>  exchangeA_;
  std::unique_ptr<RealMatrix>  exchangeB_;
  std::unique_ptr<RealMatrix>  moA_;
  std::unique_ptr<RealMatrix>  moB_;
  std::unique_ptr<RealMatrix>  PTA_;
  std::unique_ptr<RealMatrix>  PTB_;
  std::unique_ptr<RealMatrix>  dipole_;
  std::unique_ptr<RealMatrix>  quadpole_;
  std::unique_ptr<RealTensor3d>  octpole_;
  BasisSet *    basisset_;
  Molecule *    molecule_;
  FileIO *      fileio_;
  Controls *    controls_;
  AOIntegrals * aointegrals_;

public:
 
  bool	haveMO;
  bool	haveDensity; 
  bool	haveCoulomb;
  bool	haveExchange;
  bool  havePT;

  double   energyOneE;
  double   energyTwoE;
  double   energyNuclei;
  double   totalEnergy;

  // constructor & destructor
  SingleSlater(){;};
  ~SingleSlater() {;};
  // pseudo-constructor
  void iniSingleSlater(Molecule *,BasisSet *,
                       AOIntegrals *,FileIO *,
                       Controls *);

  //set private data
  inline void setNBasis(int nBasis) { this->nBasis_ = nBasis;};
  inline void setNAE(int nAE)    { this->nAE_ = nAE;};
  inline void setNBE(int nBE)    { this->nBE_ = nBE;};
  inline void setRHF(int RHF)    { this->RHF_ = RHF;};

  // access to private data
  inline int nBasis() { return this->nBasis_;};
  inline int nAE()    { return this->nAE_;};
  inline int nBE()    { return this->nBE_;};
  inline int nOccA()  { return this->nOccA_;};
  inline int nOccB()  { return this->nOccB_;}
  inline int nVirA()  { return this->nVirB_;};
  inline int nVirB()  { return this->nVirB_;};
  inline int RHF()    { return this->RHF_; };
  inline int spin()   { return this->spin_; };
  inline RealMatrix* densityA() { return this->densityA_.get();};
  inline RealMatrix* densityB() { return this->densityB_.get();};
  inline RealMatrix* fockA()    { return this->fockA_.get();};
  inline RealMatrix* fockB()    { return this->fockB_.get();};
  inline RealMatrix* coulombA() { return this->coulombA_.get();};
  inline RealMatrix* coulombB() { return this->coulombB_.get();};
  inline RealMatrix* exchangeA(){ return this->exchangeA_.get();};
  inline RealMatrix* exchangeB(){ return this->exchangeB_.get();};
  inline RealMatrix* moA()      { return this->moA_.get();};
  inline RealMatrix* moB()      { return this->moB_.get();};

  void formGuess();	        // form the intial guess of MO's
  void formDensity();		// form the density matrix
  void formFock();	        // form the Fock matrix
  void formCoulomb();		// form the Coulomb matrix
  void formExchange();		// form the exchange matrix
  void formPT();
  void readGuessIO();       	// read the initial guess of MO's from the input stream
  void readGuessGauMatEl(GauMatEl&); // read the intial guess of MO's from Gaussian raw matrix element file
  void readGuessGauFChk(std::string &);	// read the initial guess of MO's from the Gaussian formatted checkpoint file
  void computeEnergy();         // compute the total electronic energy
  void computeMultipole();      // compute multipole properties
  void SCF();  
  void printEnergy(); 
  void printMultipole();
  void printInfo();
  void printDensityinf();
  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagSingleSlater);
  void mpiRecv(int,int tag=tagSingleSlater);
};
} // namespace ChronusQ
#endif
