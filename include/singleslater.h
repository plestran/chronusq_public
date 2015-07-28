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
template<typename T>
class SingleSlater {
  typedef Eigen::Matrix<T,Dynamic,Dynamic,RowMajor> TMatrix;
  int      nBasis_;
  int      nShell_;
  int      nTT_;
  int      nAE_;
  int      nBE_;
  int      Ref_;
  int      nOccA_;
  int      nOccB_;
  int      nVirA_;
  int      nVirB_;
  int      multip_;
  int    **R2Index_;

  int      nTCS_; ///< Integer to scale the dimension of matricies for TCS's

  // Internal Storage
  std::unique_ptr<TMatrix>  densityA_;   ///< Alpha or Full (TCS) Density Matrix
  std::unique_ptr<TMatrix>  densityB_;   ///< Beta Density Matrix
  std::unique_ptr<TMatrix>  fockA_;      ///< Alpha or Full (TCS) Fock Matrix
  std::unique_ptr<TMatrix>  fockB_;      ///< Beta Fock Matrix
  std::unique_ptr<TMatrix>  coulombA_;   ///< Alpha or Full (TCS) Coulomb Matrix
  std::unique_ptr<TMatrix>  coulombB_;   ///< Beta Coulomb Matrix
  std::unique_ptr<TMatrix>  exchangeA_;  ///< Alpha or Full (TCS) Exchange Matrix
  std::unique_ptr<TMatrix>  exchangeB_;  ///< Beta Exchange Matrix
  std::unique_ptr<TMatrix>  moA_;        ///< Alpha or Full (TCS) MO Coefficient Matrix
  std::unique_ptr<TMatrix>  moB_;        ///< Beta MO Coefficient Matrix
  std::unique_ptr<TMatrix>  epsA_;       ///< Alpha or Full (TCS) Fock Eigenenergies
  std::unique_ptr<TMatrix>  epsB_;       ///< Beta Fock Eigenenergie
  std::unique_ptr<TMatrix>  PTA_;        ///< Alpha or Full (TCS) Perturbation Tensor
  std::unique_ptr<TMatrix>  PTB_;        ///< Beta Perturbation Tensor
  std::unique_ptr<RealMatrix>  dipole_;  ///< Electric Dipole Moment
  std::unique_ptr<RealMatrix>  quadpole_; ///< Electric Quadrupole Moment
  std::unique_ptr<RealMatrix>  tracelessQuadpole_; ///< Traceless Electric Quadrupole Moment
  std::unique_ptr<RealTensor3d>  octpole_; ///< Electric Octupole Moment
  BasisSet *    basisset_;               ///< Basis Set
  Molecule *    molecule_;               ///< Molecular specificiations
  FileIO *      fileio_;                 ///< Access to output file
  Controls *    controls_;               ///< General ChronusQ flow parameters
  AOIntegrals * aointegrals_;            ///< Molecular Integrals over GTOs (AO basis)

  int lenX_;
  int lenXp_;
  int lenF_;
  int lenP_;
  int lenB_;
  int lenCoeff_;
  int LWORK_;
  int lenLambda_;
  int lenDelF_;
  int lenOccNum_;
  int lenScr_;

  T *SCF_SCR;
  T *XMem_;
  T *FpAlphaMem_;
  T *FpBetaMem_;
  T *POldAlphaMem_;
  T *POldBetaMem_;
  T *ErrorAlphaMem_;
  T *ErrorBetaMem_;
  T *FADIIS_;
  T *FBDIIS_;
  T *WORK_;
  T *XpMem_;
  T *lambdaMem_;
  T *delFMem_;
  T *PNOMem_;
  T *occNumMem_;
  

  void initSCFMem();
  void initMemLen();
  void initSCFPtr();
  void formX();
  void formNO();
  void diagFock();
  void evalConver();

public:
 
  enum{
    RHF,
    UHF,
    CUHF,
    TCS
  }; ///< Supported references
  bool	haveMO;      ///< Have MO coefficients?
  bool	haveDensity; ///< Computed Density? (Not sure if this is used anymore)
  bool	haveCoulomb; ///< Computed Coulomb Matrix?
  bool	haveExchange;///< Computed Exchange Matrix?
  bool  havePT;      ///< Computed Perturbation Tensor?
  bool  isClosedShell;
  bool  isConverged;
  bool doCUHF;

  double   energyOneE; ///< One-bodied operator tensors traced with Density
  double   energyTwoE; ///< Two-bodied operator tensors traced with Density
  double   energyNuclei; ///< N-N Repulsion Energy
  double   totalEnergy; ///< Sum of all energetic contributions

  // constructor & destructor
  SingleSlater(){;};
  ~SingleSlater() {
  //if(this->SCF_SCR != NULL) delete [] this->SCF_SCR;
  };

  template<typename U>
  SingleSlater(SingleSlater<U> *); ///< Copy Constructor
  // pseudo-constructor
  void iniSingleSlater(Molecule *,BasisSet *,AOIntegrals *,FileIO *,Controls *);
  //set private data
  inline void setNBasis(int nBasis) { this->nBasis_ = nBasis;};
  inline void setNAE(int nAE)    { this->nAE_ = nAE;};
  inline void setNBE(int nBE)    { this->nBE_ = nBE;};
  inline void setRef(int Ref)    { this->Ref_ = Ref;};

  // access to private data
  inline int nBasis() { return this->nBasis_;};
  inline int nTT()     { return this->nTT_;};
  inline int nShell() { return this->nShell_;};
  inline int nAE()    { return this->nAE_;};
  inline int nBE()    { return this->nBE_;};
  inline int nOccA()  { return this->nOccA_;};
  inline int nOccB()  { return this->nOccB_;}
  inline int nVirA()  { return this->nVirA_;};
  inline int nVirB()  { return this->nVirB_;};
  inline int Ref()    { return this->Ref_; };
  inline int multip()  { return this->multip_;};
  inline int nOVA()    { return nOccA_*nVirA_;};
  inline int nOVB()    { return nOccB_*nVirB_;};
  inline TMatrix* densityA() { return this->densityA_.get();};
  inline TMatrix* densityB() { return this->densityB_.get();};
  inline TMatrix* fockA()    { return this->fockA_.get();};
  inline TMatrix* fockB()    { return this->fockB_.get();};
  inline TMatrix* coulombA() { return this->coulombA_.get();};
  inline TMatrix* coulombB() { return this->coulombB_.get();};
  inline TMatrix* exchangeA(){ return this->exchangeA_.get();};
  inline TMatrix* exchangeB(){ return this->exchangeB_.get();};
  inline TMatrix* moA()      { return this->moA_.get();};
  inline TMatrix* moB()      { return this->moB_.get();};
  inline TMatrix* epsA()     { return this->epsA_.get();};
  inline TMatrix* epsB()     { return this->epsB_.get();};
  inline TMatrix* PTA()      { return this->PTA_.get();};
  inline TMatrix* PTB()      { return this->PTB_.get();};
  inline RealMatrix* dipole(){ return this->dipole_.get();};
  inline RealMatrix* quadpole(){ return this->quadpole_.get();};
  inline RealMatrix* tracelessQuadpole(){ return this->tracelessQuadpole_.get();};
  inline RealTensor3d* octpole(){ return this->octpole_.get();};
  inline BasisSet *    basisset(){return this->basisset_;};
  inline Molecule *    molecule(){return this->molecule_;};
  inline FileIO *      fileio(){return this->fileio_;};
  inline Controls *    controls(){return this->controls_;};
  inline AOIntegrals * aointegrals(){return this->aointegrals_;};

  void formGuess();	        // form the intial guess of MO's
  void formDensity();		// form the density matrix
  void formFock();	        // form the Fock matrix
  void formCoulomb();		// form the Coulomb matrix
  void formExchange();		// form the exchange matrix
  void formPT();
  void matchord();              // match Guassian order of guess
  void readGuessIO();       	// read the initial guess of MO's from the input stream
  void readGuessGauMatEl(GauMatEl&); // read the intial guess of MO's from Gaussian raw matrix element file
  void readGuessGauFChk(std::string &);	// read the initial guess of MO's from the Gaussian formatted checkpoint file
  void computeEnergy();         // compute the total electronic energy
  void computeMultipole();      // compute multipole properties
  void SCF();  
  void SCF2();
  void CDIIS(int,T*,T*,T*,T*);
  void printEnergy(); 
  void printMultipole();
  void printInfo();
  void printDensityInfo(double,double,double);
  void printDensityInfo(double,double);

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagSingleSlater);
  void mpiRecv(int,int tag=tagSingleSlater);
};

#include <singleslater_alloc.h>
#include <singleslater_print.h>
#include <singleslater_fock.h>
#include <singleslater_misc.h>


} // namespace ChronusQ
#endif
