/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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
#ifndef  INCLUDED_AOINTEGRAL
#define  INCLUDED_AOINTEGRAL
//#include <gsl/gsl_sf_erf.h>
#include <global.h>
#include <cerr.h>
#include <basisset.h>
#include <molecule.h>
#include <fileio.h>
#include <tools.h>
#include <grid.h>
#include <memory.h>

#define MaxFmTPt 3201
#define MaxTotalL 18

namespace ChronusQ {
/**
 *  \brief Information pertaining to classical nuclei 
 *
 *  Mainly for use with in-house integral code
 */
#define MAXNAOS 21

//xslis
struct ShellPair{
  typedef double real_t;
  ShellCQ  *iShell;
  ShellCQ  *jShell;
  int	  lTotal;		    // total angular momenta of the shell pair
  int     nPGTOPair;
  int     ibf_s;
  int     jbf_s;
  int     isphbf_s;             //spherical basis function start
  int     jsphbf_s;
  int     icarbf_s;
  int     jcarbf_s;
  int     isphsize;             // number of spherical gaussian
  int     jsphsize;
  int     icarsize;
  int     jcarsize;
  std::array<real_t,3> A;	// x,y,z coordinate of center A
  std::array<real_t,3> B;	// x,y,z coordinate of center B
  std::array<real_t,3> AB;	// x,y,z distance between centers xA-xB, yA-yB, zA-zB
  std::vector<real_t> KAB;
  std::vector<real_t> UAB;
  std::vector<real_t> zetaa;
  std::vector<real_t> zetab;
  std::vector<real_t> Zeta;	// the total of exponents (alpha+beta) 
  std::vector<real_t> Xi;	// the alpha*beta/(alpha+beta) 
  std::vector<real_t> invZeta;	// the inverse of the total of exponents 1/(alpha+beta) 
  std::vector<real_t> halfInvZeta;	// the inverse of the total of exponents 0.5/(alpha+beta) 
  std::vector<real_t> ss;
  std::vector<real_t> norm;
  std::vector<real_t> beta;
  std::vector<std::array<real_t,3>> P;
  std::vector<std::array<real_t,3>> PA;
  std::vector<std::array<real_t,3>> PB;
  std::vector<std::array<real_t,3>> PZeta;	// P*(alpha+beta)

  std::vector<real_t> ssV;
  std::vector<real_t> ssT;
  std::vector<std::array<real_t,3>> ssL;
};
//xslie

struct MolecularConstants{
  int nAtom; ///< number of nuclei
  int atomZ[MAXATOMS]; ///< Classical charges of nuclei
  double cart[3][MAXATOMS]; ///< Cartesian co-ordinates for classical nuclei
};

/**
 *  \brief Struct to contain information about basis set shell-pair constants
 *
 *  Mainly for use with in-houe integral code
 */
struct PairConstants{
  /* one-center info */
  int nPGTOs[2];   ///< number of PGTOs
  int aoIndex[2];  ///< location of PGTO in ao[]
  int L[2];        ///< total angular momentum of each AO
  int nBasis[2];   ///< number of Basis Function
  double deltaAB[3];
  double intSmall;

  /* two-center info */
  double ssPairTotal;
  double ssPair[MAXCONTRACTION][MAXCONTRACTION];
  bool   ssNonzero[MAXCONTRACTION][MAXCONTRACTION];
  double deltaPA[3][MAXCONTRACTION][MAXCONTRACTION];
  double deltaPB[3][MAXCONTRACTION][MAXCONTRACTION];
  double deltaPZ[3][MAXCONTRACTION][MAXCONTRACTION][MAXATOMS];

  /* parameters required in hRR and vRR */
  double Sa0Par[MAXCONTRACTION][MAXCONTRACTION];
  double TabPar1[MAXCONTRACTION][MAXCONTRACTION];
  double TabPar2[MAXCONTRACTION][MAXCONTRACTION];
  double TabPar3[MAXCONTRACTION][MAXCONTRACTION];
  double Ta0Par3[MAXCONTRACTION][MAXCONTRACTION];
  double FmU[10][MAXCONTRACTION][MAXCONTRACTION][MAXATOMS];   // [ss]^m = FmU * vK
};

/**
 *  \brief Struct to contain information about basis set shell-quartet constants
 *
 *  Mainly for use with in-houe integral code
 */
struct QuartetConstants{
  /* four-center info */
  double deltaWP[3][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];
  double deltaWQ[3][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];
  double normQuartet[MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];

  /* parameters required in vRR */
  double a0c0Par1[MAXCONTRACTION][MAXCONTRACTION];
  double a000Par1[MAXCONTRACTION][MAXCONTRACTION];
  double a0c0Par2[MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];
  double a000Par2[MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];
  double a0c0Par3[MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];

  /* [ss|ss]^m = Fmt * vKK */
  double FmT[MAXANGULARMOMENTUM][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];
};

/****************************/
/* Error Messages 7000-7999 */
/****************************/
/**
 *  \brief Class to contains molecular integrals over basis functions.
 *
 *  A general class that contained both molecular integrals over (gaussian) basis functions
 *  as well as associated meta data that pertains to their computation. This class acts as
 *  the driver for both the in-house and Libint (Valeev) gaussian integral engines. It aims
 *  to provide a high-level interface for the computation of molecular integrals in the
 *  atomic orbital basis for ease of implementation of new quantum mechanical methods
 */
class AOIntegrals{
  int       nBasis_; ///< Number of basis functions \f$N_{b}\f$
  int       nTT_; ///< Reduced number of basis functions (lower triangle) \f$ N_b (N_b+1) / 2\f$
  int       maxMultipole_;
  int       maxNumInt_;
  int       printLevel_;
  int       nSphBasis_;
  int       nCartBasis_;

  double thresholdSchwartz_;
  double thresholdS_;
  double thresholdAB_;

  BasisSet *    	basisSet_; ///< Pointer to primary basis set
  BasisSet *     DFbasisSet_; ///< Pointer to density fitting basis set
  Molecule *   	molecule_; ///< Pointer to molecule specification
  FileIO *      	fileio_; ///< Pointer to FileIO
  TwoDGrid *            twodgrid_; ///< 3D grid (1Rad times 1 Ang) 
  CQMemManager *        memManager_;

  int       **R2Index_;
  std::unique_ptr<PairConstants>        pairConstants_; ///< Smart pointer to struct containing shell-pair meta-data
  std::unique_ptr<MolecularConstants>   molecularConstants_; ///< Smart pointer to struct containing molecular struture meta-data
  std::unique_ptr<QuartetConstants>     quartetConstants_; ///< Smart pointer to struct containing shell-quartet meta-data

  void OneEDriver(libint2::Operator); ///< General wrapper for one-electron integrals using Libint integral engine
  void breakUpMultipole();

  inline void checkWorkers(){
    if(this->memManager_  == NULL) 
      CErr("Fatal: Must initialize AOIntegrals with CQMemManager Object");
    if(this->fileio_  == NULL) 
      CErr("Fatal: Must initialize AOIntegrals with FileIO Object");
    if(this->basisSet_ == NULL) 
      CErr("Fatal: Must initialize AOIntegrals with BasisSet Object",
           this->fileio_->out);
    if(this->molecule_ == NULL) 
      CErr("Fatal: Must initialize AOIntegrals with Molecule Object",
           this->fileio_->out);
  }

  inline void checkMeta(){
    this->checkWorkers();
    if(this->nBasis_ == 0)
      CErr(
        "Fatal: AOIntegrals Object Initialized with NBasis = 0 or NShell = 0",
        this->fileio_->out);
  }

public:
  // these should be protected
  std::unique_ptr<RealMatrix>    pVp_; ///pVdotp integral
  // Two-Body Integrals
  std::unique_ptr<RealMap>    twoEC_; ///< Two-body Coulomb integrals  (deprecated)
  std::unique_ptr<RealMap>    twoEX_; ///< Two-body Exchange integrals (deprecated)
  std::unique_ptr<RealMap>    schwartz_; ///< Schwartz bounds for ERI screening
  std::unique_ptr<RealTensor4d>  aoERI_; ///< Rank-4 ERI tensor over primary basis functions \f$ (\mu \nu \vert \lambda\delta )\f$

  // One-Body Integrals
  std::unique_ptr<RealMap>    coreH_; ///< Core Hamiltonian \f$ h = T + V \f$
  std::unique_ptr<RealMap>    oneEmx_;
  std::unique_ptr<RealMap>    oneEmy_;
  std::unique_ptr<RealMap>    oneEmz_;
  std::unique_ptr<RealMap>    overlap_; ///< Overlap matrix \f$ S_{\mu\nu} = \langle \mu \vert \nu \rangle \f$
  std::unique_ptr<RealMap>    kinetic_; ///< Kinetic energy tensor \f$ T_{\mu\nu} = \langle \mu \vert \Delta \vert \nu \rangle \f$
  std::unique_ptr<RealMap>    potential_; ///< Potential energy tensor \f$ V_{\mu\nu} = \sum_A \left\langle \mu \vert r_{1A}^{-1}\vert \nu\right\rangle\f$
  std::unique_ptr<RealMap>    ortho1_;
  std::unique_ptr<RealMap>    ortho2_;

  // Resolution of the Identity
  std::unique_ptr<RealTensor3d>  aoRII_; ///< Rank-3 DFI tensor over density-fitting basis functions \f$ ( \mu\nu \vert X ) \f$
  std::unique_ptr<RealTensor2d>  aoRIS_; ///< Rank-2 Metric overlap tensor over density-fitting basis functions \f$\left( X \vert r_{12}^{-1} \vert Y \right)\f$

  // Multipole Integrals
  std::unique_ptr<RealTensor3d>  elecDipole_; ///< Electric dipole matrix \f$\vec{\mu}_{\nu\sigma}=\langle\nu\vert\vec{r}\vert\sigma\rangle\f$
  std::unique_ptr<RealTensor3d>  elecQuadpole_;///< Electric quadrupole matrix \f$Q_{\mu\nu}^{ij}=\langle\mu\vert r_i r_j \vert\nu\rangle\f$
  std::unique_ptr<RealTensor3d>  elecOctpole_;///< Electric octupole matrix \f$O_{\mu\nu}^{ijk}=\langle\mu\vert r_i r_j r_k \vert\nu\rangle\f$
  std::unique_ptr<RealTensor3d>  angular_;///Angular momentum integral, also known as R cross P
  std::unique_ptr<RealTensor3d>  RcrossDel_; ///< R cross Del matrix \f$\vec{\mu}_{\nu\sigma}=\langle\nu\vert\vec{r} \Del \vert\sigma\rangle\f$
  std::unique_ptr<RealTensor3d>  SOCoupling_;///Spin orbital integral


  // Storage for separation of multipoles
  std::vector<ConstRealMap> elecDipoleSep_;
  std::vector<ConstRealMap> elecQuadpoleSep_;
  std::vector<ConstRealMap> elecOctpoleSep_;
  std::vector<ConstRealMap> orbitalDipoleSep_; ///< r x del


  bool	haveAOTwoE; ///< Whether or not the two-bodied molecular integrals have been evaluated (for in-core integrals)
  bool	haveAOOneE; ///< Whether or not the one-body molecular integrals have been evaluated
  bool  haveSchwartz; ///< Whether or not the Schwartz bound tensor has been evaluated for the primary basis set
  bool  haveRIS; ///< Whether or not the DFI tensor has been evaluated for the density-fiting basis set
  bool  haveRII; ///< Whether or not the Metric overlap tensor has been evaluated for the density-fitting basis set
  bool  haveTRII;
  bool  isPrimary; ///< Whether or not this is the primary calculation (i.e. not guess)
  bool  useFiniteWidthNuclei; ///< Whether or not to use finite width nuclei in the calculation of the one-body potential
  bool  doX2C; //Whether to do the X2C transformation
  int   twoEFudge; // Apply fudge factor before to PVP integrals (2), or after to spin hamiltonian (1), or not at all (0).


  // Timing Stats
  std::chrono::duration<double> OneED; ///< Timer for core Hamiltonian assembly
  std::chrono::duration<double> SED; ///< Timer for overlap matrix evaluation
  std::chrono::duration<double> TED; ///< Timer for kinetic energy tensor evaluation
  std::chrono::duration<double> VED; ///< Timer for potential energy tensor evaluation
  std::chrono::duration<double> SchwartzD; ///< Timer for Schwartz bound evaluation
  std::chrono::duration<double> OrthoD; ///< Timer for orthonormal transformation formation
  std::chrono::duration<double> X2CD; ///< Timer for X2C transformation  

  enum INTEGRAL_ALGORITHM {
    DIRECT,
    INCORE,
    DENFIT
  };

  int integralAlgorithm;

  AOIntegrals(){
    this->nBasis_ = 0;
    this->nTT_    = 0;
    this->nSphBasis_ = 0;
    this->nCartBasis_ = 0;

    this->R2Index_  = NULL;

    this->molecule_   = NULL; 
    this->basisSet_   = NULL; 
    this->fileio_     = NULL; 
    this->DFbasisSet_ = NULL;
    this->memManager_ = NULL;

    this->pairConstants_      = nullptr;
    this->molecularConstants_ = nullptr;
    this->quartetConstants_   = nullptr;

    this->twoEC_        = nullptr;
    this->twoEX_        = nullptr;
    this->coreH_        = nullptr;
    this->oneEmx_       = nullptr;
    this->oneEmy_       = nullptr;
    this->oneEmz_       = nullptr;
    this->overlap_      = nullptr;
    this->kinetic_      = nullptr;
    this->potential_    = nullptr;
    this->angular_      = nullptr;
    this->schwartz_     = nullptr;
    this->aoERI_        = nullptr;
    this->aoRII_        = nullptr;
    this->aoRIS_        = nullptr;
    this->elecDipole_   = nullptr;
    this->elecQuadpole_ = nullptr;
    this->elecOctpole_  = nullptr;
    this->RcrossDel_   = nullptr;
    this->SOCoupling_   = nullptr;
    this->pVp_          = nullptr;

    this->haveAOTwoE   = false;
    this->haveAOOneE   = false;
    this->haveSchwartz = false;
    this->haveRIS      = false;
    this->haveRII      = false;
    this->haveTRII     = false;
    this->doX2C	       = false;
    this->twoEFudge    = 1;

    // Standard Values
    this->maxMultipole_     = 3;
    this->integralAlgorithm = DIRECT;
    this->isPrimary         = true;
    this->useFiniteWidthNuclei = false;
    this->printLevel_       = 1;
    this->thresholdS_ =        1.0e-10;
    this->thresholdAB_ =       1.0e-6;
    this->thresholdSchwartz_ = 1.0e-14;
    this->orthoType = LOWDIN;
  };

  ~AOIntegrals(){
/*
    // Free up memory from memory manager
    cout << "HERE" << endl;
    prettyPrintSmart(cout,*this->coreH_,"H");
    this->memManager_->free(this->coreH_->data(),
      this->nBasis_*this->nBasis_);
    cout << "HERE" << endl;
    this->memManager_->free(this->overlap_->data(),
      this->nBasis_*this->nBasis_);
    cout << "HERE" << endl;
    this->memManager_->free(this->kinetic_->data(),
      this->nBasis_*this->nBasis_);
    cout << "HERE" << endl;
    this->memManager_->free(this->potential_->data(),
      this->nBasis_*this->nBasis_);
    cout << "HERE" << endl;
    this->memManager_->free(this->schwartz_->data(),
      this->nBasis_*this->nBasis_);
    cout << "HERE" << endl;
    this->memManager_->free(this->ortho1_->data(),this->nBasis_*this->nBasis_);
    cout << "HERE" << endl;
    this->memManager_->free(this->ortho2_->data(),this->nBasis_*this->nBasis_);
    cout << "HERE" << endl;

    if(this->doX2C) {
      this->memManager_->free(this->oneEmx_->data(),
        this->nBasis_*this->nBasis_);
      this->memManager_->free(this->oneEmy_->data(),
        this->nBasis_*this->nBasis_);
      this->memManager_->free(this->oneEmz_->data(),
        this->nBasis_*this->nBasis_);
    }
    cout << "HERE" << endl;
*/
  };
  

  inline void communicate(Molecule &mol, BasisSet &basis, FileIO &fileio,
      CQMemManager &memManager){
    this->molecule_   = &mol;
    this->basisSet_   = &basis;
    this->fileio_     = &fileio;
    this->memManager_ = &memManager;
    this->DFbasisSet_ = NULL;
  }

  inline void initMeta(){
    this->checkWorkers();
    this->nBasis_ = this->basisSet_->nBasis();
    this->nTT_    = this->nBasis_ * (this->nBasis_+1) / 2;
  }

  void alloc();
  void allocOp();
  void allocMultipole();
  void allocOrth();



  // IO
  void writeOneE();

  // Getters
  inline int maxMultipole(){ return this->maxMultipole_;}
  inline CQMemManager* memManager(){ return this->memManager_;};

  // Setters

  inline void setAlgorithm(const std::string &alg){
    if(!alg.compare("DIRECT"))
      this->integralAlgorithm = DIRECT;
    else if(!alg.compare("INCORE"))
      this->integralAlgorithm = INCORE;
  }

  inline void setMaxMultipole(int i){ this->maxMultipole_ = i;}
//inline void setAlgorithm(int i)   { this->integralAlgorithm = i;}
  inline void setPrintLevel(int i){ this->printLevel_ = i;}

  inline double &twoEC(int i, int j, int k, int l){
    return (*twoEC_)(this->R2Index_[i][j],this->R2Index_[k][l]);
  }; ///< Access element i,j,k,l of Coulomb tensor
  inline double &twoEC(int ij, int kl){
    return (*twoEC_)(ij,kl);
  }; ///< Access mapped element ij,kl of Coulomb tensor
  inline double &twoEX(int i, int j, int k, int l){
    return (*twoEX_)(this->R2Index_[i][j],this->R2Index_[k][l]);
  };///< Access element i,j,k,l of Exchange tensor
  inline double &twoEX(int ij, int kl){
    return (*twoEX_)(ij,kl);
  }; ///< Access mapped element i,j,k,l of Exchange tensor
  
  void computeFmTTaylor(double*,double,int,int);
  void generateFmTTable();
  void iniQuartetConstants(ShellPair*,ShellPair*);
  void iniPairConstants(ShellPair*);
  void iniMolecularConstants();

  void printTimings();
  void computeAOTwoE(); // build two-electron AO integral matrices
  void computeAOOneE(); // build one-electron AO integral matrices
  void finiteWidthPotential();
  void formP2Transformation();
  void formPVP(std::vector<std::reference_wrapper<RealMap>>&);


  enum ORTHOTYPE {
    LOWDIN,
    CHOLESKY
  };

  ORTHOTYPE orthoType;
  inline void computeOrtho(){
    if     (this->orthoType == LOWDIN)   this->computeLowdin();
    else if(this->orthoType == CHOLESKY) this->computeCholesky();
  };
  void computeLowdin();
  void computeCholesky();


  void DKH0(); // compute DKH0 relativistic correction to kinetic energy
  void printOneE();
#ifdef USE_LIBINT
  void computeSchwartz();
  void computeAORII();
  void computeAORIS();
  void transformAORII();


  enum ERI_CONTRACTION_TYPE {
    COULOMB,
    EXCHANGE,
    PAIR
  };

  template<typename Op, typename T> 
  void newTwoEContractDirect(
      const std::vector<std::reference_wrapper<Op>>&,
      std::vector<std::reference_wrapper<Op>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<T>&
  );

  template<typename Op, typename T> 
  void newTwoEContractIncore(
      const std::vector<std::reference_wrapper<Op>>&,
      std::vector<std::reference_wrapper<Op>>&,
      std::vector<ERI_CONTRACTION_TYPE>&,
      std::vector<T>&);
  
  template<typename Op, typename T>
  inline void newTwoEContract(
      const std::vector<std::reference_wrapper<Op>> &X,
      std::vector<std::reference_wrapper<Op>> &AX,
      std::vector<ERI_CONTRACTION_TYPE> &contractionList,
      std::vector<T> &scalingFactors) {
    if(this->integralAlgorithm == DIRECT)
      this->newTwoEContractDirect(X,AX,contractionList,scalingFactors);
    else if(this->integralAlgorithm == INCORE)
      this->newTwoEContractIncore(X,AX,contractionList,scalingFactors);
    else
      CErr("Only DIRECT and INCORE (N4) Integrals are avaliable",
          this->fileio_->out);
  };

  template<typename Op>
  inline void subElecDipole(Op &op, std::array<double,3> &field){
    std::array<double,3> null{{0,0,0}};
    if(field == null) return;
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++)
      op.real() -= field[iXYZ] * this->elecDipoleSep_[iXYZ];
  };





  void compareRI();
#endif
//----------------------------------------//
// member functions in integrals_onee.cpp //
//----------------------------------------//
//xslis
  std::array<std::array<double,MaxTotalL+1>,MaxFmTPt> FmTTable_;
  std::vector<ChronusQ::ShellPair> shellPairs_;
  int  nShellPair_;
  int  nShellQuartet_;
  void createShellPairs();

  void computeOverlapS();
  double hRRSab(ChronusQ::ShellPair*,int,int*,int,int*);
  double vRRSa0(ChronusQ::ShellPair*,int,int*,int);
  void computeKineticT();
  double vRRTab(ChronusQ::ShellPair*,int,int*,int,int*);
  double RRTab(ChronusQ::ShellPair*,int,int*,int,int*);
  double vRRiPPTab(ChronusQ::ShellPair*,int,int*,int,int*,int);
  double hRRiPPSab(ChronusQ::ShellPair*,int,int*,int,int*,int);
  void computePotentialV();
  double hRRVab(ChronusQ::ShellPair*,MolecularConstants*,int,int*,int,int*);
  double vRRVa0(ChronusQ::ShellPair*,double*,double*,int,int,int*,int,int);
//xslie
//SS
  void computeAngularL(std::vector<ChronusQ::ShellPair>&,RealTensor3d&);
  double Labmu(ChronusQ::ShellPair*,RealMatrix*,RealMatrix*,RealTensor3d*,int,int*,int,int*,int,int);
  void computeSL();
  double hRRiPPVab(ChronusQ::ShellPair*,int,int*,int,int*,double*,int,int,int);
  double Slabmu(ChronusQ::ShellPair*,RealMatrix*,RealMatrix*,RealTensor3d*,double*,int,int*,int,int*,int,int,int,int); 
  double vRRV0b(ChronusQ::ShellPair*,double*,double*,int,int,int*,int,int);
  void computepVdotp();
  double pVpab(ChronusQ::ShellPair*,double*,int,int*,int,int*,int,int,int); 
//  void cart2sphtrans(ChronusQ::ShellPair*, RealMatrix*, RealMatrix& , double);
  std::vector<double> cart2SphTrans(ChronusQ::ShellPair*, double*);
  std::complex <double> car2sphcoeff(int , int , int*);

//SS
//----------------------------------------//
// member functions in integrals_twoe.cpp //
//----------------------------------------//
  double twoehRRabcd(int*,ShellPair*,ShellPair*,int,int*,int,int*,int,int*,int,int*);
  double twoehRRa0cd(int*,ShellPair*,ShellPair*,int,int*,int,int*,int,int*);
  double twoepp00(int*,ShellPair*,ShellPair*,int,int*,int,int*,int,int*,int,int*);
  double twoeppp0(int*,ShellPair*,ShellPair*,int,int*,int,int*,int,int*,int,int*);
  double twoepppp(int*,ShellPair*,ShellPair*,int,int*,int,int*,int,int*,int,int*);
  double twoevRRa0c0(ShellPair*,ShellPair*,int,int,int*,int,int*,int*,int*,int*,int*);
  double twoevRRa000(ShellPair*,ShellPair*,int,int,int*,int*,int*,int*,int*);
  double twoeSSSS0(int*,ShellPair*,ShellPair*);

  // Misc Utility Functions that deal with AOIntegrals
  static RealMatrix genSpx(BasisSet&, BasisSet&);

//template<typename Op>
//void Ortho1Trans(Op&,Op&);
//template<typename Op> Ortho1Trans<Op,double>(Op&,Op&);
//template<typename Op> Ortho1Trans<Op,dcomplex>(Op&,Op&);

  template<typename Op,
    typename std::enable_if<
      std::is_same<typename Op::Scalar,
                   double>::value,int>::type = 0>
  void Ortho1Trans(Op&,Op&);

  template<typename Op,
    typename std::enable_if<
      std::is_same<typename Op::Scalar,
                   dcomplex>::value,int>::type = 0>
  void Ortho1Trans(Op&,Op&);

  template<typename Op,
    typename std::enable_if<
      std::is_same<typename Op::Scalar,
                   double>::value,int>::type = 0>
  void Ortho1TransT(Op&,Op&);

  template<typename Op,
    typename std::enable_if<
      std::is_same<typename Op::Scalar,
                   dcomplex>::value,int>::type = 0>
  void Ortho1TransT(Op&,Op&);

  template<typename Op,
    typename std::enable_if<
      std::is_same<typename Op::Scalar,
                   double>::value,int>::type = 0>
  void Ortho2Trans(Op&,Op&);

  template<typename Op,
    typename std::enable_if<
      std::is_same<typename Op::Scalar,
                   dcomplex>::value,int>::type = 0>
  void Ortho2Trans(Op&,Op&);
};


} // namespace ChronusQ

#endif
