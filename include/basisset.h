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
#ifndef INCLUDED_BASISSET
#define INCLUDED_BASISSET
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <fileio.h>
#include <controls.h>
#include <tools.h>

/****************************/
/* Error Messages 4000-4999 */
/****************************/


namespace ChronusQ{
struct Shell{
  std::vector<double> expo;
  int     L;                     // angular momentum  -  0 , 1 , 2,  3 , 4
  bool    forceCart_; ///< Whether or not to force cartesian basis functions 
  std::vector<double> coef;
  std::vector<double> cart;
  std::vector<double> norm;
  int     nPGTOs;                // level of contraction of the shell
  int     aoIndex;               // the starting index of the ao's in the shell in the ao[] array 
  int	  nAOs;			 // number of AOs in this shell
  int     center;                // index of the atom to which the shell belongs
};

struct ShellCQ {
  typedef double real_t;

  /// contracted Gaussian = angular momentum + sph/cart flag + contraction coefficients
  int l;
  bool pure;
  std::vector<real_t> coeff;
  std::vector<real_t> alpha; //!< exponents
  std::array<real_t, 3> O;   //!< origin
  std::vector<real_t> max_ln_coeff; //!< maximum ln of (absolute) contraction coefficient for each primitive

  size_t cartesian_size() const { return (l + 1) * (l + 2) / 2; }
  size_t size() const { return pure ? (2 * l + 1) : cartesian_size(); }

  ShellCQ(libint2::Shell& other){
    this->alpha = other.alpha;
    this->l = other.contr[0].l;
    this->coeff = other.contr[0].coeff;
    this->pure = other.contr[0].pure;
    this->O = other.O;
    this->max_ln_coeff = other.max_ln_coeff;
  };


  ShellCQ& move(const std::array<real_t, 3> new_origin) {
    O = new_origin;
    return *this;
  }

  /// embeds normalization constants into contraction coefficients. Do this before computing integrals.
  /// \note Must be done only once.
/*
  void renorm() {
    using libint2::math::df_Kminus1;
    const auto sqrt_Pi_cubed = real_t{5.56832799683170784528481798212};
    const auto np = nprim();
    for(auto& c: contr) {
      assert(c.l <= 15); // due to df_Kminus1[] a 64-bit integer type; kinda ridiculous restriction anyway
      for(auto p=0; p!=np; ++p) {
        assert(alpha[p] >= 0.0);
        if (alpha[p] != 0.) {
          const auto two_alpha = 2 * alpha[p];
          const auto two_alpha_to_am32 = pow(two_alpha,c.l+1) * sqrt(two_alpha);
          const auto norm = sqrt(pow(2,c.l) * two_alpha_to_am32/(sqrt_Pi_cubed * df_Kminus1[2*c.l] ));

          c.coeff[p] *= norm;
        }
      }
    }

    // update max log coefficients
    max_ln_coeff.resize(np);
    for(auto p=0; p!=np; ++p) {
      real_t max_ln_c = - std::numeric_limits<real_t>::max();
      for(auto& c: contr) {
        max_ln_c = std::max(max_ln_c, log(std::abs(c.coeff[p])));
      }
      max_ln_coeff[p] = max_ln_c;
    }
  }
*/
  size_t nprim() const { return alpha.size(); }

  bool operator==(const ShellCQ& other) const {
    return &other == this || (O == other.O && alpha == other.alpha && coeff == other.coeff);
  }
  bool operator!=(const ShellCQ& other) const {
    return not this->operator==(other);
  }

  static char am_symbol(size_t l) {
    static char lsymb[] = "spdfghikmnoqrtuvwxyz";
    assert(l<=19);
    return lsymb[l];
  }
  static unsigned short am_symbol_to_l(char am_symbol) {
    const char AM_SYMBOL = ::toupper(am_symbol);
    switch (AM_SYMBOL) {
      case 'S': return 0;
      case 'P': return 1;
      case 'D': return 2;
      case 'F': return 3;
      case 'G': return 4;
      case 'H': return 5;
      case 'I': return 6;
      case 'K': return 7;
      case 'M': return 8;
      case 'N': return 9;
      case 'O': return 10;
      case 'Q': return 11;
      case 'R': return 12;
      case 'T': return 13;
      case 'U': return 14;
      case 'V': return 15;
      case 'W': return 16;
      case 'X': return 17;
      case 'Y': return 18;
      case 'Z': return 19;
      default: throw "invalid angular momentum label";
    }
  }

  /// @return "unit" Shell, with exponent=0. and coefficient=1., located at the origin
/*
  static ShellCQ unit() {
    ChronusQ::ShellCQ unitshell{
                {0.0}, // exponent
                {{0, false, {1.0}}},
                {{0.0, 0.0, 0.0}} // placed at origin
            };
    unitshell.renorm();
    return unitshell;
  }
*/
};

class BasisSet{
  struct ReferenceShell{
    int atomicNumber;
    int index;
    std::vector<libint2::Shell> shells;
    std::vector<std::vector<double>> unNormalizedCons;
  }; ///< struct to hold information about the basis reference
  int  nBasis_      ; ///< Number of (Gaußian) contracted basis functions
  int  nPrimitive_  ; ///< Number of uncontracted Gaußian primitives
  int  nShell_      ; ///< Number of basis shells
  int  nShellPair_  ; ///< Number of unique basis shell pairs
  int  maxPrim_     ; ///< Maximum number of Gaußian primitives for a single basis function
  int  maxL_        ; ///< Maximum angular momentum for a single basis function
  bool forceCart_   ; ///< Whether or not to force cartesian basis functions 
  double      * radCutSh_ ; ///< CutOff Radius for each Shell
  double      * expPairSh_ ; ///< SS Exp for each Shel Pair
  std::vector<double> basisEvalScr_;
  std::vector<double> basisEvalScr2_;
  std::vector<RealMatrix>    Car2Sph_;///< Matrix transformation Cart -> Sph

  std::vector<int>               nLShell_  ; ///< Maps L value to # of shells of that L
  std::vector<int>               mapSh2Bf_ ; ///< Maps shell number to first basis funtion
  std::vector<int>               mapSh2Cen_; ///< Maps shell number to atomic center
  std::vector<std::array<int,2>> mapCen2Bf_; ///< Maps atomic center to first basis function
  std::vector<ReferenceShell>    refShells_; ///< Set of reference shells for given basis
  std::vector<libint2::Shell>    shells_   ; ///< Local basis storage (in shells)
  std::vector<std::vector<double>> unNormCons_ ; // Unnormalized contraction coefficients
  std::unique_ptr<RealMatrix>    mapPrim2Bf_;///< Matrix transformation Prim -> Bf

  std::string basisPath_; ///< Path to the basis file

  std::unique_ptr<ifstream> basisFile_; ///< The file containing the basis defintion (G94)

  FileIO *fileio_; ///< FileIO Object for output file

  int printLevel_; ///< Level of print for basis set object


public:
//AOCartesian *ao;
//ShellPair   *shellPairs;
//Shell       *shells_old;
//int         *sortedShells;


  enum BASISSETS {
    PopleSTO3G,
    PopleSTO6G,
    Pople321G,
    Pople431G,
    Pople631G,
    Pople631ppGs,
    Pople6311pGs,
    Pople6311pGss,
    Pople6311pG2dp,
    ccpVDZ,
    ccpVTZ,
    def2SVP,
    def2SVPD,
    def2TZVP
  };

  std::map<BASISSETS,std::string> basisMap;
  std::map<std::string,BASISSETS> basisKey;



  bool haveMapSh2Bf  ; ///< (?) The map from shells to basis functions has been made
  bool haveMapSh2Cen ; ///< (?) The map from shells to atomic centers has been made
  bool haveMapCen2Bf ; ///< (?) The map from atomic centers to basis functions has been made

  std::vector<ChronusQ::ShellCQ> shellsCQ ; ///< Local basis storage (in CQ shell format)

  std::unique_ptr<RealMatrix> shBlkNormAlpha; ///< Shell Block (Inf) norm for Alpha matrix
  std::unique_ptr<RealMatrix> shBlkNormBeta;  ///< Shell Block (Inf) norm for Beta  matrix

//TIMING
  std::chrono::duration<double> duration_1;
  std::chrono::duration<double> duration_2;
  std::chrono::duration<double> duration_3;
  std::chrono::duration<double> duration_4;
  std::chrono::duration<double> duration_5;
//TIMING
  /**
   *  Default Constructor
   */ 
  BasisSet(){
    this->nBasis_          = 0      ;
    this->nPrimitive_      = 0      ;
    this->nShell_          = 0      ;
    this->nShellPair_      = 0      ;
    this->maxPrim_         = 0      ;
    this->maxL_            = 0      ;
    this->forceCart_       = false  ;
    this->haveMapSh2Bf     = false  ;
    this->haveMapSh2Cen    = false  ;
    this->basisFile_       = nullptr;
    this->fileio_          = NULL   ;
    this->radCutSh_        = NULL   ; 
    this->expPairSh_        = NULL   ; 
    this->printLevel_      = 1      ;
    this->makeBasisMap();


    this->makeCar2Sph(LIBINT2_MAX_AM);
  };

  /**
   *  Copy Constructor
   */ 
  BasisSet(const BasisSet &basis) : BasisSet(){
    this->forceCart_       = basis.forceCart_       ; 
    this->fileio_          = basis.fileio_      ; 
    this->basisPath_       = basis.basisPath_   ;
    this->refShells_       = basis.refShells_   ;
  }
  /**
   *  Destructor
   */ 
  ~BasisSet(){ // FIXME need to move these over to unique_ptr
/*
    delete[] ao;
    delete[] shellPairs;
    delete[] shells_old;
    delete[] sortedShells;
*/
    delete[] this->radCutSh_;
    delete[] this->expPairSh_;
  };

  inline void communicate(FileIO &fileio){ this->fileio_ = &fileio;};

  // Getters
  inline int     nBasis() {return this->nBasis_;       }; ///< Return # of basis functions
  inline int nPrimitive() {return this->nPrimitive_;   }; ///< Return # of primitive GTOs
  inline int     nShell() {return this->shells_.size();}; ///< Return # of basis shells
  inline int nShellPair() {return this->nShellPair_;   }; ///< Return # of unique shell pairs
  inline int       maxL() {return this->maxL_;         }; ///< Return max angular momentum
  inline int    maxPrim() {return this->maxPrim_;      }; ///< Return max # primitive GTOs
  inline int printLevel() {return this->printLevel_;   }; ///< Return printLevel
  inline double * radCutSh() {return this->radCutSh_;  }; ///< Return radCutSh
  inline double * expPairSh() {return this->expPairSh_;  }; ///< Return expPairSh
  inline  double getradCutSh(int iShell) {return this->radCutSh_[iShell];   }; ///< Return radCutSh
  
  template <typename T> double * basisEval(int,std::array<double,3>,T*);
  template <typename T> double * basisEval(libint2::Shell&,T*);
  template <typename T> double * basisDEval(int,libint2::Shell&,T*);
  double * CarToSpDEval(int iop, int L, double *cart);
  template <typename T> double * basisProdEval(libint2::Shell,libint2::Shell,T*);
  double * basisonFlyProdEval(libint2::Shell s1, int s1size, libint2::Shell s2, int s2size,double rx, double ry, double rz);
//std::vector<bool> MapGridBasis(cartGP& pt);  ///< Create a Mapping of basis over grid points
  void MapGridBasis(std::vector<bool>&,cartGP& pt);  ///< Create a Mapping of basis over grid points
  void     radcut(double thr, int maxiter, double epsConv);   //Populate all shell cut off radius
  void     popExpPairSh();   // Populate expPairSh
  double   fSpAv (int iop,int l, double alpha, double r);   //Evaluate Spheric Average of a Shell
  double   fRmax (int l, double alpha, double thr, double epsConv, int maxiter);   //Evaluate Spheric Average of a Shell
  inline libint2::Shell&     shells(int i) {return this->shells_[i];    };
  inline int                nLShell(int L) {return this->nLShell_[L];   };
  inline int &             mapSh2Bf(int i) {return this->mapSh2Bf_[i];  };
  inline int &             mapSh2Cen(int i) {return this->mapSh2Cen_[i];};
  inline std::array<int,2> mapCen2Bf(int i) {return this->mapCen2Bf_[i];};
  inline RealMatrix *      mapPrim2Bf() {return this->mapPrim2Bf_.get();};
  inline std::string       basisPath(){return this->basisPath_;};

  inline void resetMapSh2Bf() {this->mapSh2Bf_.clear(); this->haveMapSh2Bf  = false;};
  inline void resetMapSh2Cen(){this->mapSh2Cen_.clear();this->haveMapSh2Cen = false;};
  inline void resetMapCen2Bf(){this->mapCen2Bf_.clear();this->haveMapCen2Bf = false;};
  inline void resetLocalShells(){this->shells_.clear(); this->shellsCQ.clear();};
  inline void resetAll(){
    this->nBasis_     = 0;
    this->nPrimitive_ = 0;
    this->nShell_     = 0;
    this->nShellPair_ = 0;
    this->maxPrim_    = 0;
    this->maxL_       = 0;
    this->resetMapSh2Bf();
    this->resetMapSh2Cen();
    this->resetMapCen2Bf();
    this->resetLocalShells();
  }






  

  inline void setBasisPath( std::string str){ this->basisPath_  = str;};
  inline void setPrintLevel(int i          ){ this->printLevel_ = i  ;};
  inline void forceCart(){ this->forceCart_ = !this->forceCart_;};
  inline bool getforceCart(){ return this->forceCart_;};


  void printInfo();   ///< Print all info
  void printMeta();   ///< Print Meta data (nBasis, etc)
  void printHeader(); ///< Print header
  void printBasis();  ///< Print basis definition

  void basisSetRead(FileIO *,Molecule *,Controls *); ///< Parse and construct basis from file
  void findBasisFile(std::string);         ///< Attempt to find basis set file
  void parseGlobal();                      ///< Parse basis set file, generate reference
  void constructLocal(Molecule *);         ///< Construct local basis defintion
  void computeMeta();                      ///< Compute meta data (nBasis, etc)
  void makeMapSh2Bf();                  ///< generate mapSh2Bf
  void makeMapSh2Cen(Molecule *);          ///< generate mapSh2Cen
  void makeMapCen2Bf(Molecule *);          ///< generate mapCen2Bf
  void makeMapPrim2Bf();
  void makeCar2Sph(int L);
  void makeBasisMap();  ///< Generate map from basis enum to pasis path
  void renormShells();                     ///< Renormalize Libint2::Shell set
  std::vector<libint2::Shell> uncontractBasis(); ///< Unconctract the basis
  template<typename TMat> void computeShBlkNorm(const TMat&, RealMatrix &);

  void constructExtrn(Molecule *, BasisSet *); ///< Generate new basis from refernce shells
  void genUCvomLocal(BasisSet *);

  std::pair<double,double> cart2sphCoeff(unsigned,unsigned,unsigned,unsigned,unsigned);

  inline void makeMaps(Molecule* mol){
    this->makeMapSh2Bf();
    this->makeMapSh2Cen(mol);
    this->makeMapCen2Bf(mol);
  }

  // Python API
  void Wrapper_constructLocal(Molecule&);
  void Wrapper_makeMaps(Molecule&);

}; // class BasisSet

template<typename TMat>
void BasisSet::computeShBlkNorm(const TMat &Alpha, RealMatrix &ShBlkNorm) {
  if(!this->haveMapSh2Bf) this->makeMapSh2Bf(); 

  for(auto s1 = 0; s1 < this->nShell_; s1++) {
    int bf1 = this->mapSh2Bf_[s1];
    int n1  = this->shells_[s1].size();
    for(auto s2 = 0; s2 < this->nShell_; s2++) {
      int bf2 = this->mapSh2Bf_[s2];
      int n2  = this->shells_[s2].size();
     
      ShBlkNorm(s1,s2) = Alpha.block(bf1,bf2,n1,n2).template lpNorm<Infinity>();
    }
  }
};
}; // namespace ChronusQ
#endif
