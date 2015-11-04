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
#include <grid.h>
/****************************/
/* Error Messages 5000-5999 */
/****************************/

namespace ChronusQ {
template<typename T>
class SingleSlater {
  typedef Eigen::Matrix<T,Dynamic,Dynamic,RowMajor> TMatrix;

/*
  struct MetaData_ {
    int nBasis;      ///< Number of basis functions (inherited from BasisSet)
    int nShell;      ///< Number of basis shells (inherited from BasisSet)
    int nTT;         ///< LT NBasis dimenstion (N(N+1)/2)
    int nAE;         ///< Total number of alpha electrons
    int nBE;         ///< Total number of beta electrons
    int nOccA;       ///< Number of occupied alpha electrons
    int nOccB;       ///< Number of occupied beta electrons
    int nVirA;       ///< Number of virtual (unoccupied) alpha electrons
    int nVirB;       ///< Number of virtual (unoccupied) beta electrons
    int multip;      ///< Spin multiplicity (inherited from Molecule)
    int nTCS;        ///< Integer to scale the dimension of matricies for TCS's
    int maxSCFIter;  ///< Maximum number of SCF Cycles
    int maxMultipole;///< Maximum multipole order (inherited from AOIntegrals)
    int printLevel;  ///< Print level (Verbosity)

    bool isClosedShell; ///< Boolean to run through closed shell machinery
    bool isHF;          ///< Boolean of whether or not it's a HF reference
    bool isDFT;         ///< Boolean of whether or not it's a DFT reference
    bool isPrimary;

    double denTol; ///< SCF tolerence on the density
    double eneTol; ///< SCF tolerence on the energy

    int Ref;         ///< Specifies the Reference via enum

    std::string SCFType;      ///< String containing SCF Type (R/C) (R/U/G/CU)
    std::string SCFTypeShort; ///< String containing SCF Type (R/C) (R/U/G/CU)
    std::string algebraicField;      ///< String Real/Complex/(Quaternion)
    std::string algebraicFieldShort; ///< String Real/Complex/(Quaternion)

    std::array<double,3> elecField; ///< Static electric field


    int        ** R2Index;     ///< Not sure? artifact of in-house integrals
    BasisSet    * basisset;    ///< Basis Set
    Molecule    * molecule;    ///< Molecular specificiations
    FileIO      * fileio;      ///< Access to output file
    Controls    * controls;    ///< General ChronusQ flow parameters
    AOIntegrals * aointegrals; ///< Molecular Integrals over GTOs (AO basis)
  };
*/

  int      nBasis_;
  int      nShell_;
  int      nTT_;
  int      nAE_;
  int      nBE_;
  int      Ref_;
  int      CorrKernel_;
  int      ExchKernel_;
  int      DFTKernel_;
  int      nOccA_;
  int      nOccB_;
  int      nVirA_;
  int      nVirB_;
  int      multip_;
  int    **R2Index_;

  int      nTCS_; ///< Integer to scale the dimension of matricies for TCS's
  int      guess_;

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
  std::unique_ptr<RealMatrix>  epsA_;       ///< Alpha or Full (TCS) Fock Eigenenergies
  std::unique_ptr<RealMatrix>  epsB_;       ///< Beta Fock Eigenenergie
  std::unique_ptr<TMatrix>  PTA_;        ///< Alpha or Full (TCS) Perturbation Tensor
  std::unique_ptr<TMatrix>  PTB_;        ///< Beta Perturbation Tensor
  std::unique_ptr<TMatrix>  vXA_;        ///< Alpha or Full (TCS) VX
  std::unique_ptr<TMatrix>  vXB_;        ///< Beta VXC
  std::unique_ptr<TMatrix>  vCorA_;        ///< Alpha or Full Vcorr
  std::unique_ptr<TMatrix>  vCorB_;        ///< Beta Vcorr
  std::unique_ptr<RealMatrix>  dipole_;  ///< Electric Dipole Moment
  std::unique_ptr<RealMatrix>  quadpole_; ///< Electric Quadrupole Moment
  std::unique_ptr<RealMatrix>  tracelessQuadpole_; ///< Traceless Electric Quadrupole Moment
  std::unique_ptr<RealTensor3d>  octpole_; ///< Electric Octupole Moment
  BasisSet *    basisset_;               ///< Basis Set
  Molecule *    molecule_;               ///< Molecular specificiations
  FileIO *      fileio_;                 ///< Access to output file
  Controls *    controls_;               ///< General ChronusQ flow parameters
  AOIntegrals * aointegrals_;            ///< Molecular Integrals over GTOs (AO basis)
  TwoDGrid    * twodgrid_   ;            ///< 3D grid (1Rad times 1 Ang) 

  std::string SCFType_;                  ///< String containing SCF Type (R/C) (R/U/G/CU)
  std::string SCFTypeShort_;             ///< String containing SCF Type (R/C) (R/U/G/CU)
  std::string algebraicField_;           ///< String Real/Complex/(Quaternion)
  std::string algebraicFieldShort_;      ///< String Real/Complex/(Quaternion)
  std::array<double,3> elecField_;
  std::vector<double> mullPop_; ///< mulliken partial charge
  double Sx_, Sy_, Sz_, Ssq_;

  // Lengths of scratch partitions (NOT MEANT TO BE COPIED)
  int lenX_;
  int lenXp_;
  int lenF_;
  int lenP_;
  int lenB_;
  int lenCoeff_;
  int LWORK_;
  int LRWORK_;
  int lenLambda_;
  int lenDelF_;
  int lenOccNum_;
  int lenScr_;
  int lenRealScr_;

  // Pointers of scratch partitions (NOT MEANT TO BE COPIED)
  double *REAL_SCF_SCR;
  double *occNumMem_;
  double *RWORK_;
  double *SCpyMem_;
  double *SEVlMem_;
  double *SEVcMem_;
  double *LowdinWORK_;

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
  

  // Various functions the perform SCF and SCR allocation
  void initSCFMem();       ///< Initialize scratch memory for SCF
  void allocAlphaScr();    ///< Allocate scratch for Alpha related quantities
  void allocBetaScr();     ///< Allocate scratch for Beta related quantities
  void allocCUHFScr();     ///< Allocate scratch for CUHF realted quantities
  void allocLAPACKScr();   ///< Allocate LAPACK scratch space
  void allocLowdin();      ///< Allocate space for Lowin intermediates
  void cleanupSCFMem();    ///< Cleanup scratch memoty for SCF
  void cleanupAlphaScr();  ///< Cleanup scratch for Alpha related quantities
  void cleanupBetaScr();   ///< Cleanup scratch for Beta related quantities
  void cleanupCUHFScr();   ///< Cleanup scratch for CUHF realted quantities
  void cleanupLAPACKScr(); ///< Cleanup LAPACK scratch space
  void cleanupLowdin();    ///< Cleanup space for Lowin intermediates
  void complexMem();       ///< Add scratch space for Complex intermediates (?)
  void initMemLen();       ///< Populate lengths of scratch partitions
  void initSCFPtr();       ///< NULL-out pointers to scratch partitions
  void formX();            ///< Form orthonormal basis transformation matrix
  void formNO();           ///< Form Natural Orbitals
  void diagFock();         ///< Diagonalize Fock Matrix
  void mixOrbitalsSCF();   ///< Mix the orbitals for Complex / TCS SCF
  void evalConver(int);    ///< Evaluate convergence criteria for SCF

  double denTol_;
  double eneTol_;
  int maxSCFIter_;
  int maxMultipole_;

  int printLevel_;

  void allocOp();
  void allocAlphaOp();
  void allocBetaOp();
  void allocDFT();
  void allocAlphaDFT();
  void allocBetaDFT();
  void allocMultipole();

  inline void checkWorkers(){
    if(this->fileio_  == NULL) 
      CErr("Fatal: Must initialize SingleSlater with FileIO Object");
    if(this->basisset_ == NULL) 
      CErr("Fatal: Must initialize SingleSlater with BasisSet Object",
           this->fileio_->out);
    if(this->molecule_ == NULL) 
      CErr("Fatal: Must initialize SingleSlater with Molecule Object",
           this->fileio_->out);
    if(this->controls_ == NULL) 
      CErr("Fatal: Must initialize SingleSlater with Controls Object",
           this->fileio_->out);
    if(this->aointegrals_== NULL)
      CErr("Fatal: Must initialize SingleSlater with AOIntegrals Object",
           this->fileio_->out);
    
  }

  inline void checkMeta(){
    this->checkWorkers();
    if(this->nBasis_ == 0 || this->nShell_ == 0)
      CErr(
        "Fatal: SingleSlater Object Initialized with NBasis = 0 or NShell = 0",
        this->fileio_->out);

    if((this->molecule_->nTotalE() % 2 == 0 && this->multip_ % 2 == 0) ||
       (this->molecule_->nTotalE() % 2 != 0 && this->multip_ % 2 != 0))
      CErr(std::string("Fatal: The specified multiplicity is impossible within")
           +std::string(" the number of electrons given"),this->fileio_->out);

    // This is often called before the method is set
//  if(this->Ref_ == _INVALID) 
//    CErr("Fatal: SingleSlater reference not set!",this->fileio_->out);
  }


public:
 
  enum REFERENCE {
    _INVALID,
    RHF,
    UHF,
    CUHF,
    TCS,
    RKS,
    UKS,
    CUKS,
    GKS
  }; ///< Supported references

  enum GUESS {
    SAD,
    CORE,
    READ
  }; // Supported Guess Types

  enum DFT {
    NODFT,
    USERDEFINED,
    LSDA
  };

  enum CORR {
    NOCORR,
    VWN3,
    VWN5
  };

  enum EXCH {
    NOEXCH,
    EXACT,
    SLATER
  };

  bool	haveMO;      ///< Have MO coefficients?
  bool	haveDensity; ///< Computed Density? (Not sure if this is used anymore)
  bool	haveCoulomb; ///< Computed Coulomb Matrix?
  bool	haveExchange;///< Computed Exchange Matrix?
  bool  havePT;      ///< Computed Perturbation Tensor?
  bool  isClosedShell;
  bool  isConverged;
  bool  isHF;
  bool  isDFT;
  bool  isPrimary;
  bool  doDIIS;

  double   energyOneE; ///< One-bodied operator tensors traced with Density
  double   energyTwoE; ///< Two-bodied operator tensors traced with Density
  double   energyNuclei; ///< N-N Repulsion Energy
  double   totalEnergy; ///< Sum of all energetic contributions

  double   totalEx;     ///< LDA Exchange
  double   totalEcorr;  ///< Total VWN Energy
  double   eps_corr;    ///< VWN Correlation Energy Density
  double   mu_corr;     ///<  VWN Correlation Potential
  double   mu_corr_B;   ///<  VWN Correlation Potential (beta)

  int      nSCFIter;

  // constructor & destructor
  SingleSlater(){
    // Zero out integers to be set
    this->nBasis_  = 0;
    this->nShell_  = 0;
    this->nTT_     = 0;
    this->nAE_     = 0;
    this->nBE_     = 0;
    this->nOccA_   = 0;
    this->nOccB_   = 0;
    this->nVirA_   = 0;
    this->nVirB_   = 0;
    this->multip_  = 0;
    this->nSCFIter = 0;

    // Initialize Smart Pointers
    this->densityA_          = nullptr;   
    this->densityB_          = nullptr;   
    this->fockA_             = nullptr;      
    this->fockB_             = nullptr;      
    this->coulombA_          = nullptr;   
    this->coulombB_          = nullptr;   
    this->exchangeA_         = nullptr;  
    this->exchangeB_         = nullptr;  
    this->moA_               = nullptr;        
    this->moB_               = nullptr;        
    this->epsA_              = nullptr;    
    this->epsB_              = nullptr;    
    this->PTA_               = nullptr;        
    this->PTB_               = nullptr;        
    this->vXA_              = nullptr;       
    this->vXB_              = nullptr;       
    this->vCorA_              = nullptr;       
    this->vCorB_              = nullptr;       
    this->dipole_            = nullptr;  
    this->quadpole_          = nullptr;
    this->tracelessQuadpole_ = nullptr; 
    this->octpole_           = nullptr;

    // Initialize Raw Pointers
    this->R2Index_     = NULL;
    this->basisset_    = NULL;               
    this->molecule_    = NULL;               
    this->fileio_      = NULL;                 
    this->controls_    = NULL;               
    this->aointegrals_ = NULL;            

    // Initialize Booleans
    this->isConverged   = false;
    this->haveCoulomb   = false;
    this->haveExchange  = false;
    this->haveDensity   = false;
    this->haveMO        = false;
    this->havePT        = false;
    this->isClosedShell = false;


    // Standard Values
    this->Ref_         = _INVALID;
    this->CorrKernel_  = NOCORR;
    this->ExchKernel_  = NOEXCH;
    this->DFTKernel_   = NODFT;
    this->denTol_      = 1e-10;
    this->eneTol_      = 1e-12;
    this->maxSCFIter_  = 256;
    //this->maxSCFIter_  = 128;
    this->nTCS_        = 1;
    this->maxMultipole_ = 3;
    this->elecField_   = {0.0,0.0,0.0};
    this->printLevel_  = 1;
    this->isPrimary    = true;
    this->doDIIS       = true;
    this->isHF         = true;
    this->isDFT         = false;
    this->guess_       = SAD;

  };
  ~SingleSlater() {
  //if(this->SCF_SCR != NULL) delete [] this->SCF_SCR;
  };

  template<typename U>
  SingleSlater(SingleSlater<U> *); ///< Copy Constructor

  // pseudo-constructor
  void iniSingleSlater(Molecule *,BasisSet *,AOIntegrals *,FileIO *,Controls *);

  // Link up to all of the other worker classes
  inline void communicate(Molecule &mol, BasisSet&basis, AOIntegrals &aoints, 
    FileIO &fileio, Controls &controls){

    this->molecule_    = &mol;
    this->basisset_    = &basis;
    this->fileio_      = &fileio;
    this->controls_    = &controls;
    this->aointegrals_ = &aoints;
  }

  // Initialize Meta data from other worker classes
  inline void initMeta(){
    
    this->checkWorkers();

    this->nBasis_      = this->basisset_->nBasis();
    this->nTT_         = this->nBasis_ * (this->nBasis_ + 1) / 2;
    this->multip_      = this->molecule_->multip();
    this->nShell_      = this->basisset_->nShell();
    this->maxMultipole_= this->aointegrals_->maxMultipole();
    this->energyNuclei = this->molecule_->energyNuclei();

    int nTotalE  = this->molecule_->nTotalE();
    int nSingleE = this->multip_ - 1;
    this->nOccB_ = (nTotalE - nSingleE)/2;
    this->nVirB_ = this->nBasis_ - this->nOccB_;
    this->nOccA_ = this->nOccB_ + nSingleE;
    this->nVirA_ = this->nBasis_ - this->nOccA_;
    this->nAE_   = this->nOccA_;
    this->nBE_   = this->nOccB_;
  }
  void alloc();

  //set private data
  inline void setNBasis(int nBasis) { this->nBasis_ = nBasis;};
  inline void setNAE(int nAE)    { this->nAE_ = nAE;};
  inline void setNBE(int nBE)    { this->nBE_ = nBE;};
  inline void setRef(int Ref)    { this->Ref_ = Ref;};
//inline void setField(double x, double y, double z){
//  this->elecField_[0] = x;
//  this->elecField_[1] = y;
//  this->elecField_[2] = z;
//}
  inline void setField(std::array<double,3> field){
    this->elecField_ = field;
  }
  inline void Wrapper_setField(double x, double y, double z){
    this->setField({{x,y,z}});
  }
  inline void setNTCS(int i){ this->nTCS_ = i;};
  inline void setMaxMultipole(int i){ this->maxMultipole_ = i;};
  inline void setPrintLevel(int i){ this->printLevel_ = i;};
  inline void setSCFDenTol(double x){ this->denTol_ = x;};
  inline void setSCFEneTol(double x){ this->eneTol_ = x;};
  inline void setSCFMaxIter(int i){ this->maxSCFIter_ = i;};
  inline void setGuess(int i){ this->guess_ = i;};
  inline void isNotPrimary(){this->isPrimary = false;};
  inline void setCorrKernel(int i){this->CorrKernel_ = i;};
  inline void setExchKernel(int i){this->ExchKernel_ = i;};
  inline void setDFTKernel( int i){this->DFTKernel_  = i;};

  // access to private data
  inline int nBasis() { return this->nBasis_;};
  inline int nTCS()   { return this->nTCS_;};
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
  inline int CorrKernel(){return this->CorrKernel_;};
  inline int ExchKernel(){return this->ExchKernel_;};
  inline int DFTKernel(){ return this->DFTKernel_ ;};
  inline int printLevel(){ return this->printLevel_;};
  inline std::vector<double> mullPop(){ return this->mullPop_;};
  inline std::array<double,3> elecField(){ return this->elecField_;  };
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
  inline TMatrix* vXA()      { return this->vXA_.get();};
  inline TMatrix* vXB()      { return this->vXB_.get();};
  inline TMatrix* vCorA()      { return this->vCorA_.get();};
  inline TMatrix* vCorB()      { return this->vCorB_.get();};
  inline RealMatrix* epsA()     { return this->epsA_.get();};
  inline RealMatrix* epsB()     { return this->epsB_.get();};
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
  inline TwoDGrid *    twodgrid(){return this->twodgrid_;};
  
  inline std::string SCFType(){return this->SCFType_;};
  inline int         guess(){return this->guess_;};

  void formGuess();	        // form the intial guess
  void SADGuess();
  void COREGuess();
  void READGuess();
  void placeAtmDen(std::vector<int>, SingleSlater<double> &);           // Place the atomic densities into total densities for guess
  void scaleDen();              // Scale the unrestricted densities for correct # electrons
  void formDensity();		// form the density matrix
  void formFock();	        // form the Fock matrix
  void formCoulomb();		// form the Coulomb matrix
  void formExchange();		// form the exchange matrix
  void formPT();
  void formVXC();               // Form DFT VXC Term
  void computeSExpect();
  void formCor (double rho, double spindensity); // Form DFT correlarion potential 
  double EvepsVWN(int iop,double a_x, double b_x, double c_x, double x0_x, double rho ); // Form DFT correlarion potential 
  void formEx(double rho); // Form DFT exchange
  double f_spindens(int iop, double spindens);  // define f(spindendity)
  double df_spindens(double spindens);  // define df(spindendity)/dspindensity
  double df2_spindens(double spindens);  // define df2(spindendity)/dspindensity2
  double spindens(double rho_A, double rho_B);  // define spindendity
  double formBeckeW(cartGP gridPt, int iAtm);            // Evaluate Becke Weights
  double normBeckeW(cartGP gridPt);            // Evaluate Becke Weights
  void   buildVxc(cartGP gridPt, double weight);            // function to build the Vxc therm
  void matchord();              // match Guassian order of guess
  void readGuessIO();       	// read the initial guess of MO's from the input stream
  void readGuessGauMatEl(GauMatEl&); // read the intial guess of MO's from Gaussian raw matrix element file
  void readGuessGauFChk(std::string &);	// read the initial guess of MO's from the Gaussian formatted checkpoint file
  void computeEnergy();         // compute the total electronic energy
  void computeMultipole();      // compute multipole properties
  void SCF();  
  void CDIIS();
  void CpyFock(int);
  void GenDComm(int);
  void mullikenPop();
  void printEnergy(); 
  void printMultipole();
  void printInfo();
  void printDensityInfo(double,double,double);
  void printDensityInfo(double,double);
  void printSCFHeader(ostream &output=cout);
  void printSCFIter(int,double,double,double);
  void printDensity();
  void printPT();
  void printFock();
  void getAlgebraicField();
  void writeSCFFiles();
  void checkReadReference();
  
  inline void genMethString(){
    if(this->Ref_ == _INVALID) 
      CErr("Fatal: SingleSlater reference not set!",this->fileio_->out);

    this->getAlgebraicField(); 
    this->SCFType_      = this->algebraicField_      + " ";
    this->SCFTypeShort_ = this->algebraicFieldShort_ + "-";
    
    std::string generalReference;
    std::string generalRefShort;
    if(this->isHF){
      generalReference = "Hartree-Fock";
      generalRefShort  = "HF";
    } else if(this->isDFT) {
      generalReference = "Kohn-Sham";
      if(this->DFTKernel_ == USERDEFINED)
        generalRefShort  = "KS";
      else if(this->DFTKernel_ == LSDA) {
        generalReference += " (LSDA)";
        generalRefShort  = "LSDA";
      }
    }

    if(this->Ref_ == RHF) {
      this->SCFType_      += "Restricted " + generalReference; 
      this->SCFTypeShort_ += "R" + generalRefShort;
    } else if(this->Ref_ == UHF) {
      this->SCFType_      += "Unrestricted " + generalReference; 
      this->SCFTypeShort_ += "U" + generalRefShort;
    } else if(this->Ref_ == CUHF) {
      this->SCFType_      += "Constrained Unrestricted " + generalReference; 
      this->SCFTypeShort_ += "CU" + generalRefShort;
    } else if(this->Ref_ == TCS) {
      this->SCFType_      += "Generalized " + generalReference; 
      this->SCFTypeShort_ += "G" + generalRefShort;
    }
  }

  // Python API
  void Wrapper_iniSingleSlater(Molecule&,BasisSet&,AOIntegrals&,FileIO&,
    Controls&); 
  boost::python::list Wrapper_dipole();
  boost::python::list Wrapper_quadrupole();
  boost::python::list Wrapper_octupole();

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagSingleSlater);
  void mpiRecv(int,int tag=tagSingleSlater);
};

#include <singleslater_alloc.h>
#include <singleslater_guess.h>
#include <singleslater_print.h>
#include <singleslater_fock.h>
#include <singleslater_misc.h>
#include <singleslater_scf.h>
#include <singleslater_dft.h>


} // namespace ChronusQ
#endif
