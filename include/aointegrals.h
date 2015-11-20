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
#ifndef  INCLUDED_AOINTEGRAL
#define  INCLUDED_AOINTEGRAL
//#include <gsl/gsl_sf_erf.h>
#include <global.h>
#include <cerr.h>
#include <basisset.h>
#include <molecule.h>
#include <fileio.h>
#include <controls.h>
#include <tools.h>

#define MaxFmTPt 3201
#define MaxTotalL 18

namespace ChronusQ {
/**
 *  \brief Information pertaining to classical nuclei 
 *
 *  Mainly for use with in-house integral code
 */
#define MAXNAOS 21

struct ShellPair{
  typedef double real_t;
  ShellCQ  *iShell;
  ShellCQ  *jShell;
  int	  lTotal;		    // total angular momenta of the shell pair
  int     nPGTOPair;
  std::array<real_t,3> A;	// x,y,z coordinate of center A
  std::array<real_t,3> B;	// x,y,z coordinate of center B
  std::array<real_t,3> AB;	// x,y,z distance between centers xA-xB, yA-yB, zA-zB
  std::vector<real_t> KAB;
  std::vector<real_t> UAB;
  std::vector<real_t> Zeta;	// the total of exponents (alpha+beta) 
  std::vector<real_t> invZeta;	// the inverse of the total of exponents 1/(alpha+beta) 
  std::vector<real_t> halfInvZeta;	// the inverse of the total of exponents 0.5/(alpha+beta) 
  std::vector<real_t> ss;
  std::vector<real_t> norm;
  std::vector<std::array<real_t,3>> P;
  std::vector<std::array<real_t,3>> PA;
  std::vector<std::array<real_t,3>> PB;
  std::vector<std::array<real_t,3>> PZeta;	// P*(alpha+beta)
};

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
  int       nTCS_;
  int       nTT_; ///< Reduced number of basis functions (lower triangle) \f$ N_b (N_b+1) / 2\f$
  int       maxMultipole_;
  int       **R2Index_;
  double	**FmTTable_;

  BasisSet *    	basisSet_; ///< Pointer to primary basis set
  BasisSet *     DFbasisSet_; ///< Pointer to density fitting basis set
  Molecule *   	molecule_; ///< Pointer to molecule specification
  FileIO *      	fileio_; ///< Pointer to FileIO
  Controls *    	controls_; ///< Pointer to job control

  std::unique_ptr<PairConstants>        pairConstants_; ///< Smart pointer to struct containing shell-pair meta-data
  std::unique_ptr<MolecularConstants>   molecularConstants_; ///< Smart pointer to struct containing molecular struture meta-data
  std::unique_ptr<QuartetConstants>     quartetConstants_; ///< Smart pointer to struct containing shell-quartet meta-data

//dbwys
#ifdef USE_LIBINT
  void OneEDriver(libint2::OneBodyEngine::integral_type); ///< General wrapper for one-electron integrals using Libint integral engine
#endif
//dbwye


  inline void checkWorkers(){
    if(this->fileio_  == NULL) 
      CErr("Fatal: Must initialize AOIntegrals with FileIO Object");
    if(this->basisSet_ == NULL) 
      CErr("Fatal: Must initialize AOIntegrals with BasisSet Object",
           this->fileio_->out);
    if(this->molecule_ == NULL) 
      CErr("Fatal: Must initialize AOIntegrals with Molecule Object",
           this->fileio_->out);
    if(this->controls_ == NULL) 
      CErr("Fatal: Must initialize AOIntegrals with Controls Object",
           this->fileio_->out);
  }

  inline void checkMeta(){
    this->checkWorkers();
    if(this->nBasis_ == 0)
      CErr(
        "Fatal: SingleSlater Object Initialized with NBasis = 0 or NShell = 0",
        this->fileio_->out);
  }

public:
  // these should be protected
  std::unique_ptr<RealMatrix>    twoEC_; ///< Two-body Coulomb integrals 
  std::unique_ptr<RealMatrix>    twoEX_; ///< Two-body Exchange integrals
  std::unique_ptr<RealMatrix>    oneE_; ///< Core Hamiltonian \f$ h = T + V \f$
  std::unique_ptr<RealMatrix>    overlap_; ///< Overlap matrix \f$ S_{\mu\nu} = \langle \mu \vert \nu \rangle \f$
  std::unique_ptr<RealMatrix>    kinetic_; ///< Kinetic energy tensor \f$ T_{\mu\nu} = \langle \mu \vert \Delta \vert \nu \rangle \f$
  std::unique_ptr<RealMatrix>    kineticP_; ///< Kinetic energy tensor in momentum space \f$ T_{\mu\nu} = \langle \mu \vert \Delta \vert \nu \rangle \f$
  std::unique_ptr<RealMatrix>    potential_; ///< Potential (nuclear attraction) energy tensor \f$ V_{\mu\nu} = \sum_A \left\langle \mu \vert r_{1A}^{-1}\vert \nu\right\rangle\f$
  std::unique_ptr<RealMatrix>    schwartz_; ///< Schwartz bounds for ERI screening
  std::unique_ptr<RealTensor4d>  aoERI_; ///< Rank-4 ERI tensor over primary basis functions \f$ (\mu \nu \vert \lambda\delta )\f$
  std::unique_ptr<RealTensor3d>  aoRII_; ///< Rank-3 DFI tensor over density-fitting basis functions \f$ ( \mu\nu \vert X ) \f$
  std::unique_ptr<RealTensor2d>  aoRIS_; ///< Rank-2 Metric overlap tensor over density-fitting basis functions \f$\left( X \vert r_{12}^{-1} \vert Y \right)\f$
  std::unique_ptr<RealTensor3d>  elecDipole_; ///< Electric dipole matrix \f$\vec{\mu}_{\nu\sigma}=\langle\nu\vert\vec{r}\vert\sigma\rangle\f$
  std::unique_ptr<RealTensor3d>  elecQuadpole_;///< Electric quadrupole matrix \f$Q_{\mu\nu}^{ij}=\langle\mu\vert r_i r_j \vert\nu\rangle\f$
  std::unique_ptr<RealTensor3d>  elecOctpole_;///< Electric octupole matrix \f$O_{\mu\nu}^{ijk}=\langle\mu\vert r_i r_j r_k \vert\nu\rangle\f$

  bool		haveAOTwoE; ///< Whether or not the two-bodied molecular integrals have been evaluated (for in-core integrals)
  bool		haveAOOneE; ///< Whether or not the one-body molecular integrals have been evaluated
  bool          haveSchwartz; ///< Whether or not the Schwartz bound tensor has been evaluated for the primary basis set
  bool          haveRIS; ///< Whether or not the DFI tensor has been evaluated for the density-fiting basis set
  bool          haveRII; ///< Whether or not the Metric overlap tensor has been evaluated for the density-fitting basis set
  bool          haveTRII;
//bool          allocERI;
//bool          doDF;
  bool          isPrimary;


  // Timing Stats
  std::chrono::duration<double> OneED; ///< High-precision timing for core Hamiltonian assembly
  std::chrono::duration<double> SED; ///< High-precision timing for overlap matrix evaluation
  std::chrono::duration<double> TED; ///< High-precision timing for kinetic energy tensor evaluation
  std::chrono::duration<double> VED; ///< High-precision timing for potential energy tensor evaluation
  std::chrono::duration<double> CoulD; ///< High-precision timing for Coulomb tensor evaluation
  std::chrono::duration<double> ExchD; ///< High-precision timing for Exchange tensor evaluation
  std::chrono::duration<double> PTD; ///< High-precision timing for Perturbation tensor evaluation, \f$G[P]\f$
  std::chrono::duration<double> SchwartzD; ///< High-precision timing for Schwartz bound evaluation
  std::chrono::duration<double> DenShBlkD; ///< High-precision timing for Density shell-block norm evaluation

  enum INTEGRAL_ALGORITHM {
    DIRECT,
    INCORE,
    DENFIT
  };

  int integralAlgorithm;

  AOIntegrals(){
    this->nBasis_ = 0;
    this->nTT_    = 0;

    this->R2Index_  = NULL;
    this->FmTTable_ = NULL;

    this->molecule_   = NULL; 
    this->basisSet_   = NULL; 
    this->fileio_     = NULL; 
    this->controls_   = NULL; 
    this->DFbasisSet_ = NULL;

    this->pairConstants_      = nullptr;
    this->molecularConstants_ = nullptr;
    this->quartetConstants_   = nullptr;

    this->twoEC_        = nullptr;
    this->twoEX_        = nullptr;
    this->oneE_         = nullptr;
    this->overlap_      = nullptr;
    this->kinetic_      = nullptr;
    this->kineticP_     = nullptr;
    this->potential_    = nullptr;
    this->schwartz_     = nullptr;
    this->aoERI_        = nullptr;
    this->aoRII_        = nullptr;
    this->aoRIS_        = nullptr;
    this->elecDipole_   = nullptr;
    this->elecQuadpole_ = nullptr;
    this->elecOctpole_  = nullptr;

    this->haveAOTwoE   = false;
    this->haveAOOneE   = false;
    this->haveSchwartz = false;
    this->haveRIS      = false;
    this->haveRII      = false;
    this->haveTRII     = false;
//  this->allocERI     = false;
//  this->doDF         = false;

    // Standard Values
    this->nTCS_             = 1;
    this->maxMultipole_     = 3;
    this->integralAlgorithm = DIRECT;
    this->isPrimary         = true;
  };
  ~AOIntegrals(){;};
  
  void iniAOIntegrals(Molecule *,BasisSet *,
                      FileIO *,Controls *,
                      BasisSet * DFbasisSet=NULL); ///< Initialization function

  inline void communicate(Molecule &mol, BasisSet &basis, FileIO &fileio, 
                          Controls &controls){
    this->molecule_   = &mol;
    this->basisSet_   = &basis;
    this->fileio_     = &fileio;
    this->controls_   = &controls;
    this->DFbasisSet_ = NULL;
  }

  inline void initMeta(){
    this->checkWorkers();
    this->nBasis_ = this->basisSet_->nBasis();
    auto NTCSxNBASIS = this->nTCS_*this->nBasis_;
    this->nTT_    = NTCSxNBASIS * (NTCSxNBASIS+1) / 2;
  }

  void alloc();
  void allocOp();
  void allocMultipole();


  // IO
  void writeOneE();

  // Getters
  inline int nTCS(){ return this->nTCS_;}
  inline int maxMultipole(){ return this->maxMultipole_;}

  // Setters
  inline void setNTCS(int i)        { this->nTCS_         = i;}
  inline void setMaxMultipole(int i){ this->maxMultipole_ = i;}
  inline void setAlgorithm(int i)   { this->integralAlgorithm = i;}

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
//--------------------------------------------//
// member functions in integrals_builders.cpp //
//--------------------------------------------//
  void computeAOTwoE(); // build two-electron AO integral matrices
  void computeAOOneE(); // build one-electron AO integral matrices
  void DKH0(); // compute DKH0 relativistic correction to kinetic energy
  void printOneE();
#ifdef USE_LIBINT
  void computeSchwartz();
  void computeAORII();
  void computeAORIS();
  void transformAORII();
  template<typename T> void twoEContractDirect(bool,bool,bool,bool,bool,const T&,T&,const T&,T&);
  template<typename T> void twoEContractN4(bool,bool,bool,bool,bool,const T &,T &,const T &, T &);
  template<typename T> void twoEContractDF(bool,bool,bool,const T &,T &,const T &, T &);
  template<typename T>
    void multTwoEContractDirect(int, bool,bool,bool,bool,bool,const std::vector<T> &,std::vector<T> &,
                                const std::vector<T> &,std::vector<T> &);
  template<typename T> 
    void multTwoEContractN4(int, bool,bool,bool,const std::vector<T> &,std::vector<T> &,
                            const std::vector<T> &,std::vector<T> &);
  template<typename T> 
    void multTwoEContractDF(int, bool,bool,bool,const std::vector<T> &,std::vector<T> &,
                            const std::vector<T> &,std::vector<T> &);
  template<typename TMat,typename T> 
    void Restricted34Contract(bool,TMat&, const TMat &, int,int,int,int,int,int,int,int,
                                 const T*,T);
  template<typename TMat,typename T> 
    void UnRestricted34Contract(bool,TMat&, const TMat &, TMat&, const TMat &, const TMat &, 
                                   int,int,int,int,int,int,int,int,const T*,T);
  template<typename TMat,typename T>
    void Spinor34Contract(bool,TMat&,const TMat&,int,int,int,int,int,int,int,int,const T*,T);
  template<typename TMat,typename T> 
    void General24CouContract(TMat&, const TMat &, int,int,int,int,int,int,int,int,
                              const T*,T);
  template<typename TMat,typename T> 
    void Spinor24CouContract(TMat&, const TMat &, int,int,int,int,int,int,int,int,
                              const T*,T);
  template<typename TMat,typename T> void Gen34Contract(TMat&,const TMat&,int,int,int,int,T);
  template<typename TMat,typename T> void Gen23Contract(TMat&,const TMat&,int,int,int,int,T,
                                                        double);
  template<typename TMat,typename T> void Gen24Contract(TMat&,const TMat&,int,int,int,int,T);
  template<typename TMat,typename T> 
    void Spinor24Contract(TMat&,const TMat&,int,int,int,int,T);
  template<typename TMat,typename T> void GenCouContractSpinor(TMat&,const TMat&,int,int,int,
                                                               int,T);
  template<typename TMat,typename T> void GenExchContractSpinor(TMat&,const TMat&,int,int,
                                                                int,int,T,double);

  void compareRI();
#endif
//----------------------------------------//
// member functions in integrals_onee.cpp //
//----------------------------------------//
//ShellPair         *shellPairs_;
  std::vector<ShellPair> shellPairs_;
  int  nShellPair_;
  int  nShellQuartet_;
  void createShellPairs();

  void computeOverlapS(); // Depreciated
  double hRRSab(ShellPair*,int,int*,int,int*);
  double vRRSa0(ShellPair*,int,int*,int);
  void computeKineticT(); // Depreciated
  double vRRTab(ShellPair*,int,int*,int,int*,int*,int*);
  double vRRTa0(ShellPair*,int,int*,int*,int*);
  void computePotentialV(); // Depreciated
  double oneehRRTSab(ShellPair*,int,int*,int,int*,int*,int*);
  double oneehRRSab(ShellPair*,int,int*,int,int*);
  double oneehRRVab(ShellPair*,int,int*,int,int*);
  double oneevRRSa0(ShellPair*,int,int*,int*,int*);
  double oneevRRVa0(ShellPair*,int*,int,int,int*,int*,int*);
  double oneevRRTab(ShellPair*,int,int*,int,int*,int*,int*);
  double oneevRRTa0(ShellPair*,int,int*,int*,int*);
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

  // Python API
  void Wrapper_iniAOIntegrals(Molecule&,BasisSet&,FileIO&,Controls&); 
};
} // namespace ChronusQ

#endif
