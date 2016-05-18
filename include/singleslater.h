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
#ifndef INCLUDED_SINGLESLATER
#define INCLUDED_SINGLESLATER
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <aointegrals.h>
#include <grid.h>
#include <quantum.h>
#include <dft.h>

/****************************/
/* Error Messages 5000-5999 */
/****************************/

namespace ChronusQ {
template<typename T>
class SingleSlater : public Quantum<T> {
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMatrix;
  typedef Eigen::Map<TMatrix> TMap;

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

  int      guess_;
  int      nDIISExtrap_;
  int      iDIISStart_;

  // DFT Parameters
  int weightScheme_;
  int dftGrid_;
  int nRadDFTGridPts_;
  int nAngDFTGridPts_;
  double nElectrons_;

  std::vector<std::unique_ptr<DFTFunctional>> dftFunctionals_;


  std::unique_ptr<TMap>  NBSqScratch_;
  std::unique_ptr<TMap>  NBSqScratch2_;

  // Internal Storage
  std::unique_ptr<TMatrix>  coulombA_;   ///< deprecated 
  std::unique_ptr<TMatrix>  coulombB_;   ///< deprecated 
  std::unique_ptr<TMatrix>  exchangeA_;  ///< deprecated 
  std::unique_ptr<TMatrix>  exchangeB_;  ///< deprecated 

  // Fock Matrix
  std::unique_ptr<TMap>  fockA_;      ///< Alpha or Full (TCS) Fock Matrix
  std::unique_ptr<TMap>  fockB_;      ///< Beta Fock Matrix
  std::unique_ptr<TMap>  fockScalar_;
  std::unique_ptr<TMap>  fockMx_;
  std::unique_ptr<TMap>  fockMy_;
  std::unique_ptr<TMap>  fockMz_;

  // Orthonormal Fock
  std::unique_ptr<TMap>  fockOrthoA_;
  std::unique_ptr<TMap>  fockOrthoB_;
  std::unique_ptr<TMap>  fockOrthoScalar_;
  std::unique_ptr<TMap>  fockOrthoMx_;
  std::unique_ptr<TMap>  fockOrthoMy_;
  std::unique_ptr<TMap>  fockOrthoMz_;


  // Fock Eigensystem
  std::unique_ptr<TMap>  moA_;        ///< Alpha or Full MO Coefficients
  std::unique_ptr<TMap>  moB_;        ///< Beta MO Coefficient Matrix
  std::unique_ptr<RealMap>  epsA_;       ///< Alpha or Full Eigenenergies
  std::unique_ptr<RealMap>  epsB_;       ///< Beta Fock Eigenenergie

  std::unique_ptr<TMap>  PTA_;        ///< Alpha or Full Perturbation Tensor
  std::unique_ptr<TMap>  PTB_;        ///< Beta Perturbation Tensor
  std::unique_ptr<TMap>  PTScalar_;
  std::unique_ptr<TMap>  PTMx_;
  std::unique_ptr<TMap>  PTMy_;
  std::unique_ptr<TMap>  PTMz_;

  // Orthonormal Density
  std::unique_ptr<TMap>  onePDMOrthoA_;        
  std::unique_ptr<TMap>  onePDMOrthoB_;        
  std::unique_ptr<TMap>  onePDMOrthoScalar_;
  std::unique_ptr<TMap>  onePDMOrthoMx_;
  std::unique_ptr<TMap>  onePDMOrthoMy_;
  std::unique_ptr<TMap>  onePDMOrthoMz_;

  std::unique_ptr<TMap>  vXA_;        ///< Alpha or Full (TCS) VX
  std::unique_ptr<TMap>  vXB_;        ///< Beta VXC
  std::unique_ptr<TMap>  vCorA_;        ///< Alpha or Full Vcorr
  std::unique_ptr<TMap>  vCorB_;        ///< Beta Vcorr

  std::vector<RealSparseMatrix>  sparsedmudX_; ///< basis derivative (x)
  std::vector<RealSparseMatrix>  sparsedmudY_; ///< basis derivative (y)
  std::vector<RealSparseMatrix>  sparsedmudZ_; ///< basis derivative (z)
  std::vector<RealSparseMatrix>  sparsemudotX_; ///< basis * (x) 
  std::vector<RealSparseMatrix>  sparsemudotY_; ///< basis * (y) 
  std::vector<RealSparseMatrix>  sparsemudotZ_; ///< basis * (z) 
  std::vector<RealSparseMatrix> sparseMap_;     // BasisFunction Map 
  std::vector<RealSparseMatrix> sparseWeights_; // Weights Map
  std::vector<RealSparseMatrix> sparseDoRho_; // Evaluate density Map

  BasisSet *    basisset_;         ///< Basis Set
  Molecule *    molecule_;         ///< Molecular specificiations
  FileIO *      fileio_;           ///< Access to output file
  AOIntegrals * aointegrals_;      ///< Molecular Integrals over GTOs (AO basis)
  TwoDGrid    * twodgrid_   ;      ///< 3D grid (1Rad times 1 Ang) 

  std::string SCFType_;      ///< String containing SCF Type (R/C) (R/U/G/CU)
  std::string SCFTypeShort_; ///< String containing SCF Type (R/C) (R/U/G/CU)
  std::string algebraicField_;     ///< String Real/Complex/(Quaternion)
  std::string algebraicFieldShort_;///< String Real/Complex/(Quaternion)
  std::array<double,3> elecField_;
  std::vector<double> mullPop_; ///< mulliken partial charge
//  double Sx_, Sy_, Sz_, Ssq_;

  // Lengths of scratch partitions (NOT MEANT TO BE COPIED)
  int lenF_;
  int lenP_;
  int lenB_;
  int lenCoeff_;
  int LWORK_;
  int LRWORK_;
  int lenLambda_;
  int lenDelF_;
  int lenOccNum_;

  // Pointers of scratch partitions (NOT MEANT TO BE COPIED)
  double *occNumMem_;
  double *RWORK_;

  T *FpAlphaMem_;
  T *FpBetaMem_;
  T *POldAlphaMem_;
  T *POldBetaMem_;
  T *ErrorAlphaMem_;
  T *ErrorBetaMem_;
  T *FADIIS_;
  T *FBDIIS_;
  T *WORK_;
  T *lambdaMem_;
  T *delFMem_;
  T *PNOMem_;

  T *NBSQScr1_;
  T *NBSQScr2_;
  T *ScalarScr1_;
  T *ScalarScr2_;
  T *MzScr1_;
  T *MzScr2_;
  T *MyScr1_;
  T *MyScr2_;
  T *MxScr1_;
  T *MxScr2_;

  

  // Various functions the perform SCF and SCR allocation
  void initSCFMem();       ///< Initialize scratch memory for SCF
  void allocAlphaScr();    ///< Allocate scratch for Alpha related quantities
  void allocBetaScr();     ///< Allocate scratch for Beta related quantities
  void allocCUHFScr();     ///< Allocate scratch for CUHF realted quantities
  void allocLAPACKScr();   ///< Allocate LAPACK scratch space
  void cleanupSCFMem();    ///< Cleanup scratch memoty for SCF
  void cleanupAlphaScr();  ///< Cleanup scratch for Alpha related quantities
  void cleanupBetaScr();   ///< Cleanup scratch for Beta related quantities
  void cleanupCUHFScr();   ///< Cleanup scratch for CUHF realted quantities
  void cleanupLAPACKScr(); ///< Cleanup LAPACK scratch space
  void initMemLen();       ///< Populate lengths of scratch partitions
  void initSCFPtr();       ///< NULL-out pointers to scratch partitions
  void formNO();           ///< Form Natural Orbitals
  void diagFock();         ///< Diagonalize Fock Matrix
  void mixOrbitalsSCF();   ///< Mix the orbitals for Complex / TCS SCF
  void evalConver(int);    ///< Evaluate convergence criteria for SCF

  void initSCFMem2();       ///< Initialize scratch memory for SCF (2)
  void diagFock2();         ///< Diagonalize Fock Matrix
  void orthoFock();
  void fockCUHF();
  void orthoDen();
  void cleanupSCFMem2();
  void copyDen();
  void genDComm2(int);
  void backTransformMOs();

  double denTol_;
  double eneTol_;
  int maxSCFIter_;
  bool  fixPhase_;

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
    if(this->aointegrals_== NULL)
      CErr("Fatal: Must initialize SingleSlater with AOIntegrals Object",
           this->fileio_->out);
    if(this->memManager_ == NULL)
      CErr("Fatal: Must initialize SingleSlater with CQMemManager Object",
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
    VWN5,
    LYP
  };

  enum EXCH {
    NOEXCH,
    EXACT,
    SLATER,
    B88
  };

  enum DFT_GRID {
    EULERMACL,
    GAUSSCHEB
  };

  enum DFT_WEIGHT_SCHEME{
    BECKE,
    FRISCH
  };
 
  bool	haveMO;      ///< Have MO coefficients?
  bool	haveDensity; ///< Computed Density? (Not sure if this is used anymore)
  bool	haveCoulomb; ///< Computed Coulomb Matrix?
  bool	haveExchange;///< Computed Exchange Matrix?
  bool	screenVxc   ;///< Do the screening for Vxc?
  bool  havePT;      ///< Computed Perturbation Tensor?
  bool  isConverged;
  bool  isHF;
  bool  isDFT;
  bool  isPrimary;
  bool  doDIIS;
  bool  isGGA;

  double   energyOneE; ///< One-bodied operator tensors traced with Density
  double   energyTwoE; ///< Two-bodied operator tensors traced with Density
  double   energyNuclei; ///< N-N Repulsion Energy
  double   totalEnergy; ///< Sum of all energetic contributions

  double   totalEx;     ///< LDA Exchange
  double   totalEcorr;  ///< Total VWN Energy
  double   eps_corr;    ///< VWN Correlation Energy Density
  double   mu_corr;     ///<  VWN Correlation Potential
  double   mu_corr_B;   ///<  VWN Correlation Potential (beta)
  double   epsScreen;   ///<  Screening value for both basis and Bweight
  double   epsConv;     ///<  Threshold value for converging cutoff radius given epsScreen
  int      maxiter;     ///<  Maximum number of iteration to find cutoff radius
  int      ngpts;       ///<  Total number of grid pts
//// T
//  std::chrono::duration<double> duration_1;
//  std::chrono::duration<double> duration_2;
//  std::chrono::duration<double> duration_3;
//  std::chrono::duration<double> duration_4;
//  std::chrono::duration<double> duration_5;
//  std::chrono::duration<double> duration_6;
//  std::chrono::duration<double> duration_7;
//  std::chrono::duration<double> duration_8;
  int      nSCFIter;



  // constructor & destructor
  SingleSlater() : Quantum<T>(){
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
    this->ngpts    = 0;

    // Initialize Smart Pointers
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

    // Initialize Raw Pointers
    this->R2Index_     = NULL;
    this->basisset_    = NULL;               
    this->molecule_    = NULL;               
    this->fileio_      = NULL;                 
    this->aointegrals_ = NULL;            

    // Initialize Booleans
    this->isConverged   = false;
    this->haveCoulomb   = false;
    this->haveExchange  = false;
    this->haveDensity   = false;
    this->haveMO        = false;
    this->havePT        = false;


    // Standard Values
    this->Ref_          = _INVALID;
    this->CorrKernel_   = NOCORR;
    this->ExchKernel_   = NOEXCH;
    this->DFTKernel_    = NODFT;
    this->denTol_       = 1e-8;
    this->eneTol_       = 1e-10;
    this->maxSCFIter_   = 256;
    this->nDIISExtrap_  = 7;
    this->iDIISStart_   = 4;

    this->elecField_   = {0.0,0.0,0.0};
    this->printLevel_  = 1;
    this->isPrimary    = true;
    this->doDIIS       = true;
    this->isHF         = true;
    this->isDFT         = false;
    this->fixPhase_     = true;
    this->guess_       = SAD;

    this->weightScheme_ = BECKE;
    this->dftGrid_      = GAUSSCHEB;
    this->screenVxc     = true;
    this->epsScreen     = 1.0e-10;
    this->nRadDFTGridPts_ = 100;
    this->nAngDFTGridPts_ = 302;
    this->isGGA = false;

    this->dftFunctionals_.emplace_back(new SlaterExchange());

    // FIXME: maybe hardcode these?
    this->epsConv       = 1.0e-7;
    this->maxiter       = 50;

  };
  ~SingleSlater() { ; };

  SingleSlater(SingleSlater *other) : 
    Quantum<T>(dynamic_cast<Quantum<T>&>(*other)){

    this->nBasis_ = other->nBasis();
    this->nTT_    = other->nTT();
    this->nAE_    = other->nAE();
    this->nBE_    = other->nBE(); 
    this->nOccA_  = other->nOccA();
    this->nOccB_  = other->nOccB();
    this->nVirA_  = other->nVirA();
    this->nVirB_  = other->nVirB();
    this->multip_   = other->multip();
    this->energyNuclei = other->energyNuclei;
    this->Ref_    = other->Ref();
    this->haveDensity = true;
    this->haveMO	    = true;
    this->havePT      = true;
    this->printLevel_ = other->printLevel();
    this->doDIIS = other->doDIIS;
    this->isHF   = other->isHF;
    this->isDFT  = other->isDFT;
    this->guess_ = other->guess();
    this->elecField_   = (other->elecField());
    this->basisset_    = other->basisset();    
    this->molecule_    = other->molecule();
    this->fileio_      = other->fileio();
    this->aointegrals_ = other->aointegrals();

    auto NB = this->nBasis_*this->nTCS_;
    auto NBSq = NB*NB;
    this->allocOp();

    (*this->fockA_) = *other->fockA();
    (*this->moA_  ) = *other->moA();
    (*this->PTA_  ) = *other->PTA();

    if(!this->isClosedShell || this->nTCS_ == 2){
      (*this->fockB_) = *other->fockB();
      (*this->moB_  ) = *other->moB();
      (*this->PTB_  ) = *other->PTB();
    }
    
  };

  template<typename U>
  SingleSlater(U *);

  // Link up to all of the other worker classes
  inline void communicate(Molecule &mol, BasisSet&basis, AOIntegrals &aoints, 
    FileIO &fileio, CQMemManager &memManager){

    this->molecule_    = &mol;
    this->basisset_    = &basis;
    this->fileio_      = &fileio;
    this->aointegrals_ = &aoints;
    this->memManager_  = &memManager;
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
  inline void setNBasis(int nBasis)       { this->nBasis_ = nBasis;    };
  inline void setNAE(int nAE)             { this->nAE_ = nAE;          };
  inline void setNBE(int nBE)             { this->nBE_ = nBE;          };
  inline void setRef(int Ref)             { this->Ref_ = Ref;          };
  inline void setPrintLevel(int i)        { this->printLevel_ = i;     };
  inline void setSCFDenTol(double x)      { this->denTol_ = x;         };
  inline void setSCFEneTol(double x)      { this->eneTol_ = x;         };
  inline void setSCFMaxIter(int i)        { this->maxSCFIter_ = i;     };
  inline void setGuess(int i)             { this->guess_ = i;          };
  inline void isNotPrimary()              { this->isPrimary = false;   };
  inline void setCorrKernel(int i)        { this->CorrKernel_ = i;     };
  inline void setExchKernel(int i)        { this->ExchKernel_ = i;     };
  inline void setDFTKernel( int i)        { this->DFTKernel_  = i;     };
  inline void setDFTWeightScheme(int i)   { this->weightScheme_ = i;   };
  inline void setDFTGrid(int i)           { this->dftGrid_ = i;        };
  inline void setDFTNGridPts(int i, int j){ this->nRadDFTGridPts_ = i;
                                            this->nAngDFTGridPts_ = j; };
  inline void setDFTNRad(int i)           { this->nRadDFTGridPts_ = i;};
  inline void setDFTNAng(int i)           { this->nAngDFTGridPts_ = i; };
  inline void setDFTScreenTol(double x)   { this->epsScreen = x;       };
  inline void turnOffDFTScreening()       { this->screenVxc = false;   }; 

  inline void setField(std::array<double,3> field){ 
    this->elecField_ = field;  
  };
  inline void Wrapper_setField(double x, double y, double z){ 
    this->setField({{x,y,z}}); 
  };
  


  // access to private data
  inline int nBasis()    { return this->nBasis_;                  };
  inline int nTT()       { return this->nTT_;                     };
  inline int nShell()    { return this->nShell_;                  };
  inline int nAE()       { return this->nAE_;                     };
  inline int nBE()       { return this->nBE_;                     };
  inline int nOccA()     { return this->nOccA_;                   };
  inline int nOccB()     { return this->nOccB_;                   };     
  inline int nVirA()     { return this->nVirA_;                   };
  inline int nVirB()     { return this->nVirB_;                   };
  inline int Ref()       { return this->Ref_;                     };      
  inline int multip()    { return this->multip_;                  };
  inline int nOVA()      { return nOccA_*nVirA_;                  };
  inline int nOVB()      { return nOccB_*nVirB_;                  };
  inline int CorrKernel(){ return this->CorrKernel_;              };
  inline int ExchKernel(){ return this->ExchKernel_;              };
  inline int DFTKernel() { return this->DFTKernel_ ;              };
  inline int printLevel(){ return this->printLevel_;              };
  inline std::vector<double> mullPop()   { return this->mullPop_; };
  inline std::array<double,3> elecField(){ return this->elecField_; };

  inline TMap* fockA()                { return this->fockA_.get();    };
  inline TMap* fockB()                { return this->fockB_.get();    };
  inline TMap* fockScalar()           { return this->fockScalar_.get();};
  inline TMap* fockMz()           { return this->fockMz_.get();};
  inline TMap* fockMy()           { return this->fockMy_.get();};
  inline TMap* fockMx()           { return this->fockMx_.get();};
  inline TMap* coulombA()             { return this->coulombA_.get(); };
  inline TMap* coulombB()             { return this->coulombB_.get(); };
  inline TMap* exchangeA()            { return this->exchangeA_.get();};
  inline TMap* exchangeB()            { return this->exchangeB_.get();};
  inline TMap* moA()                  { return this->moA_.get();      };
  inline TMap* moB()                  { return this->moB_.get();      };
  inline TMap* vXA()                  { return this->vXA_.get();      };
  inline TMap* vXB()                  { return this->vXB_.get();      };
  inline TMap* vCorA()                { return this->vCorA_.get();    };
  inline TMap* vCorB()                { return this->vCorB_.get();    };
  inline RealMap* epsA()              { return this->epsA_.get();     };
  inline RealMap* epsB()              { return this->epsB_.get();     };
  inline TMap* PTA()                  { return this->PTA_.get();      };
  inline TMap* PTB()                  { return this->PTB_.get();      };

  inline BasisSet     * basisset()       { return this->basisset_;       };
  inline Molecule     * molecule()       { return this->molecule_;       };
  inline FileIO       * fileio()         { return this->fileio_;         };
  inline AOIntegrals  * aointegrals()    { return this->aointegrals_;    };
  inline TwoDGrid     * twodgrid()       { return this->twodgrid_;       };
  inline std::string SCFType()           { return this->SCFType_;        };
  inline int         guess()             { return this->guess_;          };

  inline void checkDFTType(){
    if(this->CorrKernel_ == LYP || this->ExchKernel_ == B88)
      this->isGGA = true;
  };
  void formGuess();	        // form the intial guess
  void SADGuess();
  void COREGuess();
  void READGuess();
  void placeAtmDen(std::vector<int>, SingleSlater<double> &); // Place the atomic densities into total densities for guess
  void scaleDen();              // Scale the unrestricted densities for correct # electrons
  void formDensity();		// form the density matrix
  void formFock();	        // form the Fock matrix
  void formCoulomb();		// form the Coulomb matrix
  void formExchange();		// form the exchange matrix
  void formPT();
  void formVXC();               // Form DFT VXC Term
  void formVXC_store();               // Form DFT VXC Term
  void genSparseBasisMap();     // Generate Basis Set Mapping
  void genSparseRcrosP();      //  Generate R cros P int
  void evalVXC_store(int, int, double &, double &,RealMatrix *, RealMatrix *,
         RealMatrix *, RealMatrix *, RealMatrix *,RealMatrix *,RealMatrix *,
         RealMatrix *);
  void evalVXC(cartGP, double, std::vector<bool>, double &, double &,
      RealMatrix *, RealMatrix *, RealMatrix *, RealMatrix*); // evaluate DFT VXC Matrix Term( at a given pts)
  std::array<double,6 > formVC (double, double);    // Form DFT correlarion density,potential (A and B)
  std::array<double,6 > formVCGGA (double, double, double, double, double);    // Form DFT GGA correlarion density,potential (A and B)
  std::array<double,6 > formVCVWN (double, double); // Form DFT VWN correlation (VWN3 and VWN5)
  std::array<double,6 > formVCLYP (double, double, double, double, double); // Form DFT LYP correlation 
  std::array<double,6 > formVEx (double, double); // Form DFT exchange density, potential (A and B)  
  std::array<double,6 > formVExGGA (double, double, double, double);    // Form DFT GGA Exchange density,potential (A and B)
  std::array<double,6 > formVExSlater (double, double); // Form DFT Slater exchange
  std::array<double,6 > formVExB88 (double, double, double, double); // Form DFT Becke88 exchange
  double EvepsVWN(int iop,double a_x, double b_x, double c_x, double x0_x, 
      double rho ); // Form DFT correlarion potential 
  double gB88(int, double);                   //funtion used in B88 Exchange
  double omegaLYP(int, double, double, double);                    //function used in LYP Correlation 
  double derLYP(int, double, double, double, double, double, double);//function used in LYP Correlation 
  double der2LYP(int, double, double, double, double, double, double, double);//function used in LYP Correlation 
  double deltaLYP(int, double, double, double);  //function used in LYP Correlation
  double f_spindens(int iop, double spindens);  // define f(spindendity)
  double df_spindens(double spindens);  // define df(spindendity)/dspindensity
  double df2_spindens(double spindens);  // define df2(spindendity)/dspindensity2
  double spindens(double rho_A, double rho_B);   // define spindendity
  double formBeckeW(cartGP gridPt, int iAtm);    // Evaluate Becke Weights
  double normBeckeW(cartGP gridPt);             // Normalize Becke Weights
//  void   buildVxc(cartGP gridPt, double weight, std::vector<bool> mapRad_);            // function to build the Vxc therm
  void matchord();              // match Guassian order of guess
  void readGuessIO();       	// read the initial guess of MO's from the input stream
  void readGuessGauMatEl(GauMatEl&); // read the intial guess of MO's from Gaussian raw matrix element file
  void readGuessGauFChk(std::string &);	// read the initial guess of MO's from the Gaussian formatted checkpoint file
  void computeEnergy();         // compute the total electronic energy
  void computeMultipole();      // compute multipole properties
//void computeSExpect();        // compute <S> <S^2>
  void computeSSq();
  inline void computeProperties(){
    this->computeMultipole();
    this->computeSExpect(*this->aointegrals_->overlap_);
  };
  void SCF();  
  void CDIIS();
  void CpyFock(int);
  void GenDComm(int);
  void mullikenPop();
  void fixPhase();

  void SCF2();


  void printEnergy(); 
  void printMultipole();
  void printSExpect();
  inline void printProperties() {
    this->printSExpect();
    this->printMultipole();
  };
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
  boost::python::list Wrapper_dipole();
  boost::python::list Wrapper_quadrupole();
  boost::python::list Wrapper_octupole();

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag);
  void mpiRecv(int,int tag);
};

#include <singleslater/singleslater_alloc.h>
#include <singleslater/singleslater_guess.h>
#include <singleslater/singleslater_print.h>
#include <singleslater/singleslater_fock.h>
#include <singleslater/singleslater_misc.h>
#include <singleslater/singleslater_scf.h>
#include <singleslater/singleslater_diis.h>
#include <singleslater/singleslater_properties.h>
//#include <singleslater_dft.h>


} // namespace ChronusQ
#endif
