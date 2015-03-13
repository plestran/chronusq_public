#ifndef INCLUDED_SINGLESLATER
#define INCLUDED_SINGLESLATER
#include "global.h"
#include "matrix.h"
#include "molecule.h"
#include "controls.h"
#include "aointegrals.h"

/****************************/
/* Error Messages 5000-5999 */
/****************************/

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
  Matrix<double>  *densityA_;
  Matrix<double>  *densityB_;
  Matrix<double>  *fockA_;
  Matrix<double>  *fockB_;
  Matrix<double>  *coulombA_;
  Matrix<double>  *coulombB_;
  Matrix<double>  *exchangeA_;
  Matrix<double>  *exchangeB_;
  Matrix<double>  *moA_;
  Matrix<double>  *moB_;
  BasisSet     	*basisset_;
  Molecule    	*molecule_;
  FileIO       	*fileio_;
  Controls     	*controls_;
  AOIntegrals   *aointegrals_;

public:
 
  bool	haveMO;
  bool	haveDensity; 
  bool	haveCoulomb;
  bool	haveExchange;

  double   energyOneE;
  double   energyTwoE;
  double   energyNuclei;
  double   totalEnergy;

  // constructor & destructor
  SingleSlater(){;};
  ~SingleSlater() {
    if(   densityA_!=NULL) delete densityA_;
    if(      fockA_!=NULL) delete fockA_;
    if(   coulombA_!=NULL) delete coulombA_;
    if(  exchangeA_!=NULL) delete exchangeA_;
    if(        moA_!=NULL) delete moA_;
    if(!RHF_) {
      if(      moB_!=NULL) delete moB_;
      if(exchangeB_!=NULL) delete exchangeB_;
      if( coulombB_!=NULL) delete coulombB_;
      if(    fockB_!=NULL) delete fockB_;
      if( densityB_!=NULL) delete densityB_;
    };
  };
  // pseudo-constructor
  void iniSingleSlater(Molecule*,BasisSet*,AOIntegrals*,FileIO*,Controls*);

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
  inline Matrix<double> *densityA() { return this->densityA_;};
  inline Matrix<double> *densityB() { return this->densityB_;};
  inline Matrix<double> *fockA()    { return this->fockA_;};
  inline Matrix<double> *fockB()    { return this->fockB_;};
  inline Matrix<double> *coulombA() { return this->coulombA_;};
  inline Matrix<double> *coulombB() { return this->coulombB_;};
  inline Matrix<double> *exchangeA(){ return this->exchangeA_;};
  inline Matrix<double> *exchangeB(){ return this->exchangeB_;};
  inline Matrix<double> *moA()      { return this->moA_;};
  inline Matrix<double> *moB()      { return this->moB_;};

  void formGuess();	        // form the intial guess of MO's
  void formDensity();		// form the density matrix
  void formFock();	        // form the Fock matrix
  void formCoulomb();		// form the Coulomb matrix
  void formExchange();		// form the exchange matrix
  void readGuessIO();       	// read the initial guess of MO's from the input stream
  void readGuessGauFChk(char*);	// read the initial guess of MO's from the Gaussian formatted checkpoint file
  void computeEnergy();         // compute the total electronic energy
  void SCF();  
  void printEnergy(); 
  void printInfo();

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagSingleSlater);
  void mpiRecv(int,int tag=tagSingleSlater);
};
#endif
