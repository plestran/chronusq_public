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
#include <grid.h>
#include <grid2.h>
//#include <quantum.h>
#include <wavefunction.h>
#include <dft.h>

/****************************/
/* Error Messages 5000-5999 */
/****************************/

namespace ChronusQ {

struct SCFConvergence {
  double EDelta;
  double PSRMS;
  double PMRMS;
};

template<typename T>
class KernelIntegrand {
  typedef Eigen::Matrix<T,Dynamic,Dynamic> TMatrix;
  public:
  TMatrix VXCScalar;
  TMatrix VXCMz;
  TMatrix VXCMy;
  TMatrix VXCMx;
  double Energy;

  KernelIntegrand(size_t N) : VXCScalar(N,N), VXCMz(N,N), Energy(0.0){ 
    VXCScalar.setZero();
    VXCMz.setZero();
  };
};

enum DIIS_ALGORITHM {
  NO_DIIS_SET,
  NO_DIIS,
  CDIIS,
  EDIIS,
  CEDIIS
};

template<typename T>
class SingleSlater : public WaveFunction<T> {
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMatrix;
  typedef Eigen::Map<TMatrix> TMap;
  typedef Eigen::Matrix<T,Dynamic,1,ColMajor> TVec;

/*
  // SingleSlater class dependencies
  BasisSet *    basisset_;         ///< Basis Set
  Molecule *    molecule_;         ///< Molecular specificiations
  FileIO *      fileio_;           ///< Access to output file
  AOIntegrals * aointegrals_;      ///< Molecular Integrals over GTOs (AO basis)
*/

  // Misc Metadata
  int printLevel_;

/*
  // Privite Metadata regarding basis dimensions
  int      nBasis_;
  int      nShell_;
  int      nTT_;

  // Private Metadata regarding molecular specification
  int      nAE_;
  int      nBE_;
  int      nOccA_;
  int      nOccB_;
  int      nVirA_;
  int      nVirB_;
  int      multip_;
*/

  // Private Metadata regarding the reference
  int      Ref_;
  std::string SCFType_;      ///< String containing SCF Type (R/C) (R/U/G/CU)
  std::string SCFTypeShort_; ///< String containing SCF Type (R/C) (R/U/G/CU)
  std::string algebraicField_;     ///< String Real/Complex/(Quaternion)
  std::string algebraicFieldShort_;///< String Real/Complex/(Quaternion)

  // Private Metadata regarding SCF control
  // Fock Formation
  bool doIncFock_;
  int  nIncFock_;

  // SCF convergence parameters
  double denTol_;
  double eneTol_;
  int maxSCFIter_;

  // DIIS related parameters
  int            nDIISExtrap_;
  int            iDIISStart_;
  DIIS_ALGORITHM diisAlg_;


  // General extrapolation parameters
  int nKeep_;

  // Level Shifting Parameters
  int      iStartLevelShift_;
  int      nLevelShift_;
  double   levelShiftParam_;

  // Misc SCF control
  bool  fixPhase_;

  int      guess_;

  int    **R2Index_;

  // Private metadata regarding DFT
  int DFTKernel_;
  int weightScheme_;
  int dftGrid_;
  int nRadDFTGridPts_;
  int nAngDFTGridPts_;
  double nElectrons_;


  std::unique_ptr<TMap>  NBSqScratch_;
  std::unique_ptr<TMap>  NBSqScratch2_;
  std::unique_ptr<TMap>  NBSqScratch3_;
  std::unique_ptr<TMap>  NBSqScratch4_;
  std::unique_ptr<TMap>  NBTSqScratch_;

  // Internal Storage
//std::unique_ptr<TMap>  coulombA_;   ///< deprecated 
//std::unique_ptr<TMap>  coulombB_;   ///< deprecated 
//std::unique_ptr<TMap>  exchangeA_;  ///< deprecated 
//std::unique_ptr<TMap>  exchangeB_;  ///< deprecated 

  // Fock Matrix
//std::unique_ptr<TMap>  fockA_;      ///< Alpha or Full (TCS) Fock Matrix
//std::unique_ptr<TMap>  fockB_;      ///< Beta Fock Matrix
  std::unique_ptr<TMap>  fockScalar_;
  std::unique_ptr<TMap>  fockMx_;
  std::unique_ptr<TMap>  fockMy_;
  std::unique_ptr<TMap>  fockMz_;
  std::vector<TMap*>     fock_;

  // Orthonormal Fock
//std::unique_ptr<TMap>  fockOrthoA_;
//std::unique_ptr<TMap>  fockOrthoB_;
  std::unique_ptr<TMap>  fockOrthoScalar_;
  std::unique_ptr<TMap>  fockOrthoMx_;
  std::unique_ptr<TMap>  fockOrthoMy_;
  std::unique_ptr<TMap>  fockOrthoMz_;
  std::vector<TMap*>     fockOrtho_;


  // Fock Eigensystem
/*
  std::unique_ptr<TMap>  moA_;        ///< Alpha or Full MO Coefficients
  std::unique_ptr<TMap>  moB_;        ///< Beta MO Coefficient Matrix
*/
/*
  std::unique_ptr<RealMap>  epsA_;       ///< Alpha or Full Eigenenergies
  std::unique_ptr<RealMap>  epsB_;       ///< Beta Fock Eigenenergie
*/

//std::unique_ptr<TMap>  PTA_;        ///< Alpha or Full Perturbation Tensor
//std::unique_ptr<TMap>  PTB_;        ///< Beta Perturbation Tensor
  std::unique_ptr<TMap>  PTScalar_;
  std::unique_ptr<TMap>  PTMx_;
  std::unique_ptr<TMap>  PTMy_;
  std::unique_ptr<TMap>  PTMz_;
  std::vector<TMap*>     PT_;

  // Orthonormal Density
//std::unique_ptr<TMap>  onePDMOrthoA_;        
//std::unique_ptr<TMap>  onePDMOrthoB_;        
  std::unique_ptr<TMap>  onePDMOrthoScalar_;
  std::unique_ptr<TMap>  onePDMOrthoMx_;
  std::unique_ptr<TMap>  onePDMOrthoMy_;
  std::unique_ptr<TMap>  onePDMOrthoMz_;
  std::vector<TMap*>     onePDMOrtho_;

//std::unique_ptr<TMap>  vXA_;        ///< Alpha or Full (TCS) VX
//std::unique_ptr<TMap>  vXB_;        ///< Beta VXC
//std::unique_ptr<TMap>  vCorA_;        ///< Alpha or Full Vcorr
//std::unique_ptr<TMap>  vCorB_;        ///< Beta Vcorr
  std::unique_ptr<TMap>  vXCScalar_;
  std::unique_ptr<TMap>  vXCMz_;
  std::unique_ptr<TMap>  vXCMy_;
  std::unique_ptr<TMap>  vXCMx_;
  std::vector<TMap*>     vXC_;



  std::array<double,3> elecField_;
  std::vector<double> mullPop_; ///< mulliken partial charge


  // Pointers of scratch partitions (NOT MEANT TO BE COPIED)
  double *occNumMem_;

  T *lambdaMem_;
  T *delFMem_;
  T *PNOMem_;


  // Storage Files for Extrapolation
  H5::DataSet *FScalarDIIS_;
  H5::DataSet *FMzDIIS_;
  H5::DataSet *FMyDIIS_;
  H5::DataSet *FMxDIIS_;
  std::vector<H5::DataSet*> FDIIS_;

  H5::DataSet *DScalarDIIS_;
  H5::DataSet *DMzDIIS_;
  H5::DataSet *DMyDIIS_;
  H5::DataSet *DMxDIIS_;
  std::vector<H5::DataSet*> DDIIS_;

  H5::DataSet *EScalarDIIS_;
  H5::DataSet *EMzDIIS_;
  H5::DataSet *EMyDIIS_;
  H5::DataSet *EMxDIIS_;
  std::vector<H5::DataSet*> CommDIIS_;

  H5::DataSet *PTScalarDIIS_;
  H5::DataSet *PTMzDIIS_;
  H5::DataSet *PTMyDIIS_;
  H5::DataSet *PTMxDIIS_;
  std::vector<H5::DataSet*> PTDIIS_;

  H5::DataSet *ADMPGradScalar_;
  H5::DataSet *ADMPGradMz_;
  H5::DataSet *ADMPGradMy_;
  H5::DataSet *ADMPGradMx_;
  std::vector<H5::DataSet*> ADMPGrad_;

  H5::DataSet *DMSErrScalar_;
  H5::DataSet *DMSErrMz_;
  H5::DataSet *DMSErrMy_;
  H5::DataSet *DMSErrMx_;
  std::vector<H5::DataSet*> DMSErr_;


  // Storage Files for most recent Density (for SCF comparison)
  H5::DataSet *DScalarOld_;
  H5::DataSet *DMzOld_;
  H5::DataSet *DMyOld_;
  H5::DataSet *DMxOld_;

  // Storage Files for most recent Fock (for SCF comparison)
  H5::DataSet *PTScalarOld_;
  H5::DataSet *PTMzOld_;
  H5::DataSet *PTMyOld_;
  H5::DataSet *PTMxOld_;

  // Storage for FP
  H5::DataSet *FPScalar_;
  H5::DataSet *FPMz_;
  H5::DataSet *FPMy_;
  H5::DataSet *FPMx_;

  // Storage Files for DeltaD
  H5::DataSet *DeltaDScalar_;
  H5::DataSet *DeltaDMz_;
  H5::DataSet *DeltaDMy_;
  H5::DataSet *DeltaDMx_;
 
  // New DIIS Functions
  void genDIISCom(int);
  void CDIIS4(int);

  // Various functions the perform SCF and SCR allocation
  void initSCFPtr();       ///< NULL-out pointers to scratch partitions
  void formNO();           ///< Form Natural Orbitals
  void mixOrbitalsSCF();   ///< Mix the orbitals for Complex / TCS SCF

  void fockCUHF();
  void copyDen();
  void backTransformMOs();

  void doImagTimeProp(double); ///< Propagate the wavefunction in imaginary time 


  void allocOp();
  void allocScr();
//void allocAlphaOp();
//void allocBetaOp();
  void allocDFT();
//void allocAlphaDFT();
//void allocBetaDFT();
  void allocMultipole();

/*
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
*/


public:
 
  enum REFERENCE {
    _INVALID,
    RHF,
    UHF,
    CUHF,
    TCS,
    X2C
//  RKS,
//  UKS,
//  CUKS,
//  GKS
  }; ///< Supported references

  enum GUESS {
    ONLY,
    SAD,
    CORE,
    READ,
    RANDOM
  }; // Supported Guess Types

  enum DFT {
    NODFT,
    USERDEFINED,
    LSDA,
    SVWN3,
    SVWN5,
    BLYP
  };


 
//bool	haveMO;      ///< Have MO coefficients?
//bool	haveDensity; ///< Computed Density? (Not sure if this is used anymore)
//bool	haveCoulomb; ///< Computed Coulomb Matrix?
//bool	haveExchange;///< Computed Exchange Matrix?
//bool  havePT;      ///< Computed Perturbation Tensor?
  bool  isConverged;
  bool  isHF;
  bool  isDFT;
  bool  isPrimary;
  bool  doDIIS;

  bool  doITP; ///< Do Imaginary Time Propagation (ITP)?
  float dt;    ///< Timestep for Imaginary Time Propagation

  bool	screenVxc   ;///< Do the screening for Vxc?
  bool  isGGA;

  double   energyOneE; ///< One-bodied operator tensors traced with Density
  double   energyTwoE; ///< Two-bodied operator tensors traced with Density
//double   energyNuclei; ///< N-N Repulsion Energy
  double   energyExc;
//double   totalEnergy; ///< Sum of all energetic contributions

//double   totalEx;     ///< LDA Exchange
//double   totalEcorr;  ///< Total VWN Energy
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

  std::vector<std::unique_ptr<DFTFunctional>> dftFunctionals_;


  // constructor & destructor
  SingleSlater() : WaveFunction<T>(),
/*
    nBasis_ (0),
    nShell_ (0),
    nTT_    (0),
    nAE_    (0),
    nBE_    (0),
    nOccA_  (0),
    nOccB_  (0),
    nVirA_  (0),
    nVirB_  (0),
    multip_ (0),
*/
    nSCFIter(0),
    ngpts   (0),
/*
    moA_    (nullptr),        
    moB_    (nullptr),        
    epsA_   (nullptr),    
    epsB_   (nullptr),    
*/
    vXCScalar_ (nullptr),       
    vXCMz_ (nullptr),       
    vXCMy_ (nullptr),       
    vXCMx_ (nullptr),       
    R2Index_     (NULL),
/*
    basisset_    (NULL),               
    molecule_    (NULL),               
    fileio_      (NULL),                 
    aointegrals_ (NULL),            
*/
    isConverged  (false),
    isGGA        (false),
    energyOneE   (0.0),
    energyTwoE   (0.0),
//  energyNuclei (0.0),
//  energyExc    (0.0),
//  totalEnergy  (0.0){
    energyExc    (0.0){


    // Standard Values
    this->Ref_              = _INVALID;
    this->DFTKernel_        = NODFT;
    this->denTol_           = 1e-8;
    this->eneTol_           = 1e-10;
    this->maxSCFIter_       = 256;
    this->iStartLevelShift_ = 0;
    this->nLevelShift_      = 4;
    this->levelShiftParam_  = 2.42;

    this->elecField_   = {0.0,0.0,0.0};
    this->printLevel_  = 1;
    this->isPrimary    = true;
    this->doITP        = false;
    this->dt           = 0.1;
    this->isHF         = true;
    this->isDFT        = false;
    this->fixPhase_    = true;
    this->guess_       = SAD;

    // SCF Fock Formation
    this->doIncFock_   = false;
    this->nIncFock_    = 20;

    // Extrapolation
    this->doDIIS       = true;
    this->doDMS        = false;
    this->nDIISExtrap_ = 6;
    this->iDIISStart_  = 0;
    this->diisAlg_     = DIIS_ALGORITHM::NO_DIIS_SET;

    // DFT
    this->weightScheme_ = ATOMIC_PARTITION::BECKE;
    this->dftGrid_      = GRID_TYPE::EULERMAC;
    this->screenVxc     = true;
    this->epsScreen     = 1.0e-10;
    this->nRadDFTGridPts_ = 100;
    this->nAngDFTGridPts_ = 302;


    // FIXME: maybe hardcode these?
    this->epsConv       = 1.0e-7;
    this->maxiter       = 50;

  };

  ~SingleSlater() { ; };

  SingleSlater(SingleSlater *other) : 
    WaveFunction<T>(dynamic_cast<WaveFunction<T>&>(*other)),

/*
    nBasis_ ( other->nBasis() ),
    nTT_    ( other->nTT() ),
    nAE_    ( other->nAE() ),
    nBE_    ( other->nBE() ), 
    nOccA_  ( other->nOccA() ),
    nOccB_  ( other->nOccB() ),
    nVirA_  ( other->nVirA() ),
    nVirB_  ( other->nVirB() ),
    multip_ ( other->multip() ),
*/
//  energyNuclei ( other->energyNuclei ),
    Ref_    ( other->Ref() ),
    printLevel_  ( other->printLevel() ),
    doDIIS  ( other->doDIIS ),
    isHF    ( other->isHF ),
    isDFT   ( other->isDFT ),
    guess_  ( other->guess() ),
    elecField_   ( other->elecField() )
/*
    basisset_    ( other->basisset() ),    
    molecule_    ( other->molecule() ),
    fileio_      ( other->fileio() ),
    aointegrals_ ( other->aointegrals() ) {
*/
    {
    this->alloc();

/*
    (*this->fockA_) = *other->fockA();
    (*this->moA_  ) = *other->moA();
    (*this->PTA_  ) = *other->PTA();

    if(!this->isClosedShell || this->nTCS_ == 2){
      (*this->fockB_) = *other->fockB();
      (*this->moB_  ) = *other->moB();
      (*this->PTB_  ) = *other->PTB();
    }
*/
    for(auto iF = 0; iF < this->fock_.size(); iF++){
      *this->fock_[iF] = *other->fock_[iF];
      *this->PT_[iF]   = *other->PT_[iF];

      // Copy over the orthonormal density for RT calculations
      *this->onePDMOrtho_[iF] = *other->onePDMOrtho_[iF];
    }
/*
    *this->moA_ = *other->moA_;
    if(this->nTCS_ == 2 and !this->isClosedShell)
      *this->moB_ = *other->moB_;
*/
    
  };

  template<typename U>
  SingleSlater(U *);

/*
  // Link up to all of the other worker classes
  inline void communicate(Molecule &mol, BasisSet&basis, AOIntegrals &aoints, 
    FileIO &fileio, CQMemManager &memManager){

    this->molecule_    = &mol;
    this->basisset_    = &basis;
    this->fileio_      = &fileio;
    this->aointegrals_ = &aoints;
    this->memManager_  = &memManager;
  }
*/

  // Initialize Meta data from other worker classes
  inline void initMeta(){
    WaveFunction<T>::initMeta();
  //this->checkWorkers();

/*
    this->nBasis_      = this->basisset_->nBasis();
    this->nTT_         = this->nBasis_ * (this->nBasis_ + 1) / 2;
    this->multip_      = this->molecule_->multip();
    this->nShell_      = this->basisset_->nShell();
*/
    this->maxMultipole_= this->aointegrals_->maxMultipole();
//  this->energyNuclei = this->molecule_->energyNuclei();

/*
    int nTotalE  = this->molecule_->nTotalE();
    int nSingleE = this->multip_ - 1;
    this->nOccB_ = (nTotalE - nSingleE)/2;
    this->nVirB_ = this->nBasis_ - this->nOccB_;
    this->nOccA_ = this->nOccB_ + nSingleE;
    this->nVirA_ = this->nBasis_ - this->nOccA_;
    this->nAE_   = this->nOccA_;
    this->nBE_   = this->nOccB_;
*/

    if(this->doDIIS && this->diisAlg_ == DIIS_ALGORITHM::NO_DIIS_SET)
      this->diisAlg_ = DIIS_ALGORITHM::CDIIS;

    if(this->isDFT){
      this->epsScreen /= this->molecule_->nAtoms() * 
        this->nRadDFTGridPts_ * this->nAngDFTGridPts_;
      this->basisset_->radcut(this->epsScreen,this->maxiter,this->epsConv);
    }
  }

  void alloc();
  void dealloc();

  //set private data
/*
  inline void setNBasis(int nBasis)       { this->nBasis_ = nBasis;    };
  inline void setNAE(int nAE)             { this->nAE_ = nAE;          };
  inline void setNBE(int nBE)             { this->nBE_ = nBE;          };
*/
  inline void setRef(int Ref)             { this->Ref_ = Ref;          };
  inline void setPrintLevel(int i)        { this->printLevel_ = i;     };
  inline void setSCFDenTol(double x)      { this->denTol_ = x;         };
  inline void setSCFEneTol(double x)      { this->eneTol_ = x;         };
  inline void setSCFMaxIter(int i)        { this->maxSCFIter_ = i;     };
  inline void setGuess(int i)             { this->guess_ = i;          };
  inline void isNotPrimary()              { this->isPrimary = false;   };
  inline void setDFTKernel( int i)        { this->DFTKernel_  = i;     };
  inline void setDFTWeightScheme(int i)   { this->weightScheme_ = i;   };
  inline void setDFTGrid(int i)           { this->dftGrid_ = i;        };
  inline void setDFTNGridPts(int i, int j){ this->nRadDFTGridPts_ = i;
                                            this->nAngDFTGridPts_ = j; };
  inline void setDFTNRad(int i)           { this->nRadDFTGridPts_ = i;};
  inline void setDFTNAng(int i)           { this->nAngDFTGridPts_ = i; };
  inline void setDFTScreenTol(double x)   { this->epsScreen = x;       };
  inline void turnOffDFTScreening()       { this->screenVxc = false;   }; 

  inline void setITPdt(double x)          { this->dt = x;              }; 

  inline void setField(std::array<double,3> field){ 
    this->elecField_ = field;  
  };
  inline void Wrapper_setField(double x, double y, double z){ 
    this->setField({{x,y,z}}); 
  };
  


  // access to private data
/*
  inline int nBasis()    { return this->nBasis_;                  };
  inline int nTT()       { return this->nTT_;                     };
  inline int nShell()    { return this->nShell_;                  };
  inline int nAE()       { return this->nAE_;                     };
  inline int nBE()       { return this->nBE_;                     };
  inline int nOccA()     { return this->nOccA_;                   };
  inline int nOccB()     { return this->nOccB_;                   };     
  inline int nVirA()     { return this->nVirA_;                   };
  inline int nVirB()     { return this->nVirB_;                   };
*/
  inline int Ref()       { return this->Ref_;                     };      
//inline int multip()    { return this->multip_;                  };
/*
  inline int nOVA()      { return nOccA_*nVirA_;                  };
  inline int nOVB()      { return nOccB_*nVirB_;                  };
*/
  inline int DFTKernel() { return this->DFTKernel_ ;              };
  inline int printLevel(){ return this->printLevel_;              };
  inline std::vector<double> mullPop()   { return this->mullPop_; };
  inline std::array<double,3> elecField(){ return this->elecField_; };

//inline TMap* fockA()                { return this->fockA_.get();    };
//inline TMap* fockB()                { return this->fockB_.get();    };
  inline TMap* fockScalar()           { return this->fockScalar_.get();};
  inline TMap* fockMz()           { return this->fockMz_.get();};
  inline TMap* fockMy()           { return this->fockMy_.get();};
  inline TMap* fockMx()           { return this->fockMx_.get();};
  inline std::vector<TMap*>& fock(){ return this->fock_;};
  inline std::vector<TMap*>& fockOrtho(){ return this->fockOrtho_;};
//inline TMap* coulombA()             { return this->coulombA_.get(); };
//inline TMap* coulombB()             { return this->coulombB_.get(); };
//inline TMap* exchangeA()            { return this->exchangeA_.get();};
//inline TMap* exchangeB()            { return this->exchangeB_.get();};
/*
  inline TMap* moA()                  { return this->moA_.get();      };
  inline TMap* moB()                  { return this->moB_.get();      };
*/
//inline TMap* vXA()                  { return this->vXA_.get();      };
//inline TMap* vXB()                  { return this->vXB_.get();      };
//inline TMap* vCorA()                { return this->vCorA_.get();    };
//inline TMap* vCorB()                { return this->vCorB_.get();    };
  inline TMap* vXCScalar()            { return this->vXCScalar_.get();};
  inline TMap* vXCMz()                { return this->vXCMz_.get();    };
  inline TMap* vXCMy()                { return this->vXCMy_.get();    };
  inline TMap* vXCMx()                { return this->vXCMx_.get();    };
/*
  inline RealMap* epsA()              { return this->epsA_.get();     };
  inline RealMap* epsB()              { return this->epsB_.get();     };
*/
//inline TMap* PTA()                  { return this->PTA_.get();      };
//inline TMap* PTB()                  { return this->PTB_.get();      };
  inline TMap* PTScalar()           { return this->PTScalar_.get();};
  inline TMap* PTMz()           { return this->PTMz_.get();};
  inline TMap* PTMy()           { return this->PTMy_.get();};
  inline TMap* PTMx()           { return this->PTMx_.get();};
  inline std::vector<TMap*>& PT(){ return this->PT_;};

  inline TMap* onePDMOrthoScalar() { return this->onePDMOrthoScalar_.get();};
  inline TMap* onePDMOrthoMz()     { return this->onePDMOrthoMz_.get();};
  inline TMap* onePDMOrthoMy()     { return this->onePDMOrthoMy_.get();};
  inline TMap* onePDMOrthoMx()     { return this->onePDMOrthoMx_.get();};
  inline std::vector<TMap*>& onePDMOrtho(){ return this->onePDMOrtho_;};
/*
  inline BasisSet     * basisset()       { return this->basisset_;       };
  inline Molecule     * molecule()       { return this->molecule_;       };
  inline FileIO       * fileio()         { return this->fileio_;         };
  inline AOIntegrals  * aointegrals()    { return this->aointegrals_;    };
*/
//inline TwoDGrid     * twodgrid()       { return this->twodgrid_;       };
  inline std::string SCFType()           { return this->SCFType_;        };
  inline int         guess()             { return this->guess_;          };

  void formGuess();	        // form the intial guess
  void SADGuess();
  void COREGuess();
  void READGuess();
  void RandomGuess();
  void placeAtmDen(std::vector<int>, SingleSlater<double> &); // Place the atomic densities into total densities for guess
  void scaleDen();              // Scale the unrestricted densities for correct # electrons
  void formDensity();		// form the density matrix
  void formFock(bool increment=false);	        // form the Fock matrix
  void formCoulomb();		// form the Coulomb matrix
  void formExchange();		// form the exchange matrix
  void formPT();
//void formVXC();               // Form DFT VXC Term
//void formVXC_store();               // Form DFT VXC Term
  void formVXC_new();               // Form DFT VXC Term
//void genSparseBasisMap();     // Generate Basis Set Mapping
//void genSparseRcrosP();      //  Generate R cros P int
//void evalVXC_store(int, int, double &, double &,RealMatrix *, RealMatrix *,
//       RealMatrix *, RealMatrix *, RealMatrix *,RealMatrix *,RealMatrix *,
//       RealMatrix *);
//void evalVXC(cartGP, double, std::vector<bool>, double &, double &,
//    RealMatrix *, RealMatrix *, RealMatrix *, RealMatrix*); // evaluate DFT VXC Matrix Term( at a given pts)
//std::array<double,6 > formVC (double, double);    // Form DFT correlarion density,potential (A and B)
//std::array<double,6 > formVCGGA (double, double, double, double, double);    // Form DFT GGA correlarion density,potential (A and B)
//std::array<double,6 > formVCVWN (double, double); // Form DFT VWN correlation (VWN3 and VWN5)
//std::array<double,6 > formVCLYP (double, double, double, double, double); // Form DFT LYP correlation 
//std::array<double,6 > formVEx (double, double); // Form DFT exchange density, potential (A and B)  
//std::array<double,6 > formVExGGA (double, double, double, double);    // Form DFT GGA Exchange density,potential (A and B)
//std::array<double,6 > formVExSlater (double, double); // Form DFT Slater exchange
//std::array<double,6 > formVExB88 (double, double, double, double); // Form DFT Becke88 exchange
//double EvepsVWN(int iop,double a_x, double b_x, double c_x, double x0_x, 
//    double rho ); // Form DFT correlarion potential 
//double gB88(int, double);                   //funtion used in B88 Exchange
//double omegaLYP(int, double, double, double);                    //function used in LYP Correlation 
//double derLYP(int, double, double, double, double, double, double);//function used in LYP Correlation 
//double der2LYP(int, double, double, double, double, double, double, double);//function used in LYP Correlation 
//double deltaLYP(int, double, double, double);  //function used in LYP Correlation
//double f_spindens(int iop, double spindens);  // define f(spindendity)
//double df_spindens(double spindens);  // define df(spindendity)/dspindensity
//double df2_spindens(double spindens);  // define df2(spindendity)/dspindensity2
//double spindens(double rho_A, double rho_B);   // define spindendity
//double formBeckeW(cartGP gridPt, int iAtm);    // Evaluate Becke Weights
//double normBeckeW(cartGP gridPt);             // Normalize Becke Weights
//  void   buildVxc(cartGP gridPt, double weight, std::vector<bool> mapRad_);            // function to build the Vxc therm
  void matchord();              // match Guassian order of guess
  void readGuessIO();       	// read the initial guess of MO's from the input stream
  void readGuessGauMatEl(GauMatEl&); // read the intial guess of MO's from Gaussian raw matrix element file
  void readGuessGauFChk(std::string &);	// read the initial guess of MO's from the Gaussian formatted checkpoint file
  void computeEnergy();         // compute the total electronic energy
  void computeMultipole();      // compute multipole properties
  void computeSSq();
  inline void computeProperties(){
    this->computeMultipole();
    this->computeSExpect(*this->aointegrals_->overlap_);
  };
  void CpyFock(int);
  void GenDComm(int);
  void mullikenPop();
  void fixPhase();

  void levelShift();
  void levelShift2();


  // DFT Setup Routines
  inline void addSlater(){
    this->dftFunctionals_.emplace_back(new SlaterExchange());
  };
  inline void addB88(){
    this->dftFunctionals_.emplace_back(new BEightEight());
    this->isGGA = true;
  };
  inline void addLYP(){
    this->dftFunctionals_.emplace_back(new lyp()); 
    this->isGGA = true;
  };
  inline void addVWN5(){
    this->dftFunctionals_.emplace_back(new VWNV()); 
  };
  inline void addVWN3(){
    this->dftFunctionals_.emplace_back(new VWNIII()); 
  };


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
    } else if(this->nTCS_ == 2) {
      this->SCFType_      += "Generalized " + generalReference; 
      this->SCFTypeShort_ += "G" + generalRefShort;
    } else if(this->Ref_ == X2C) {
      this->SCFType_      += "Exact Two-Component " + generalReference; 
      this->SCFTypeShort_ += "X2C-"+generalRefShort;
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



  void SCF3();
  void initSCFMem3();
  void initDIISFiles();
  void orthoFock3();
  void unOrthoDen3();
  SCFConvergence evalConver3();
  void backTransformMOs3();
  void cleanupSCFMem3();
  void cpyDenDIIS(int);
  void cpyFockDIIS(int);
  void readDIIS(H5::DataSet*,int,T*);
  void writeDIIS(H5::DataSet*,int,T*);
  void diagFock2();         ///< Diagonalize Fock Matrix
  void cpyAOtoOrthoDen();

  void gatherOrthoFock();
  void gatherFock();
  void gatherOrthoDen();

  void populateMO4Diag();

  void mixOrbitals2C();
  void mixOrbitalsComplex();

  void formADMPGrad(int);
  void formDMSErr(int);
  void DMSExtrap(int);
  void initDMSFiles();
  void initADMPFiles();
  void McWeeny(int);
  bool doDMS;

  void formDeltaD();
  void copyDeltaDtoD();
  void copyDOldtoD();
  void copyPT();
  void incPT();

  void formFP();
  
};

} // namespace ChronusQ
#endif
