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
  char    name[2];               // name of the shell - "S","P","D","F","G"
  bool	  SP;			 // is this part of an SP shell?
  int     L;                     // angular momentum  -  0 , 1 , 2,  3 , 4
  int     center;                // index of the atom to which the shell belongs
  int     nPGTOs;                // level of contraction of the shell
  int     aoIndex;               // the starting index of the ao's in the shell in the ao[] array 
  int	  nAOs;			 // number of AOs in this shell
  double  coef[MAXCONTRACTION];  // contraction coefficients
  double  expo[MAXCONTRACTION];  // exponents of primitive GTOs
  double  norm[MAXCONTRACTION];  // normalization constants of primitive GTOs (divided by AOCartesian.divConst)
};

struct AOCartesian{
  int     lx;                    // x angular momentum
  int     ly;                    // y angular momentum
  int     lz;                    // z angular momentum
  int	  l[3];			 // x,y,z angular momentum
  int     shIndex;               // index of the shell in the shells[] array to which current ao belongs
  double  divConst;              // normalization constant = Shell.norm/divConst
};

struct BasisPair{
  int     aoIndex[2];
  int     nPGTOs[2];		// number of primitive Gaussian (degree of contract) on each shell
  int     L[2];                 // angular momenta of the first and second shells
  int     centerIndex[2];	// indices of the atom where the first and second shells are centered
};

struct ShellPair{
  int     L[2];                 // angular momenta of the first and second shells
  int	  LTotal;		// total angular momenta of the shell pair
  int     aoIndex[2];		// starting indeices of AOs in the first and second shells
  int     shIndex[2];           // indices of the first and second shell in the shells[] array to which current ao belongs
  int     center[2];		// indices of the atom where the first and second shells are centered
  int     nPGTOs[2];		// number of primitive Gaussian (degree of contract) on each shell
  int     nBasis[2];		// number of basis function on each shell
  int     nAOPair;		// number of unique AO pairs
  int     aoPairIndex[250][2];  // the indices of the ao pair in the ao[] array

  double  divConst[250];	// division constant; 
  double  centerA[3];		// x,y,z coordinate of center A
  double  centerB[3];		// x,y,z coordinate of center B
  double  deltaAB[3];		// x,y,z distance between centers xA-xB, yA-yB, zA-zB
  double  centerP[3][MAXCONTRACTION][MAXCONTRACTION];	// x,y,z coordinate of the combined center of the shellpair (alpha*A+beta*B)/(alpha+beta)
  double  centerPZeta[3][MAXCONTRACTION][MAXCONTRACTION];	// centerP*(alpha+beta)
  double  deltaPA[3][MAXCONTRACTION][MAXCONTRACTION];	// x,y,z distance P-A
  double  deltaPB[3][MAXCONTRACTION][MAXCONTRACTION];	// x,y,z distance P-B
  double  deltaPApPB[3][3][MAXCONTRACTION][MAXCONTRACTION];// PAi+PBj
  double  deltaPAtPB[3][3][MAXCONTRACTION][MAXCONTRACTION];// PAi*PBj
  double  zeta[MAXCONTRACTION][MAXCONTRACTION];		// the total of exponents (alpha+beta) 
  double  inversezeta[MAXCONTRACTION][MAXCONTRACTION];	// the inverse of the total of exponents 0.5/(alpha+beta) 
  double  invzeta[MAXCONTRACTION][MAXCONTRACTION];	// the inverse of the total of exponents 0.5/(alpha+beta) 
  double  KAB[MAXCONTRACTION][MAXCONTRACTION];		// KAB used to compute the [ss|ss] integral
  double  UAB[MAXCONTRACTION][MAXCONTRACTION];		// KAB used to compute the [ss|ss] integral
  double  norm[MAXCONTRACTION][MAXCONTRACTION];		// pairwise normalization constant
};

class BasisSet{
  struct ReferenceShell{
    int atomicNumber;
    int index;
    std::vector<libint2::Shell> shells;
  }; ///< struct to hold information about the basis reference
  int  nBasis_      ; ///< Number of (Gaußian) contracted basis functions
  int  nPrimitive_  ; ///< Number of uncontracted Gaußian primitives
  int  nShell_      ; ///< Number of basis shells
  int  nShellPair_  ; ///< Number of unique basis shell pairs
  int  maxPrim_     ; ///< Maximum number of Gaußian primitives for a single basis function
  int  maxL_        ; ///< Maximum angular momentum for a single basis function
  bool doSph_       ; ///< Whether or not to make the cartesian -> spherical transformation
  double      * radCutSh_ ; ///< CutOff Radius for each Shell

  std::vector<int>               nLShell_  ; ///< Maps L value to # of shells of that L
  std::vector<int>               mapSh2Bf_ ; ///< Maps shell number to first basis funtion
  std::vector<int>               mapSh2Cen_; ///< Maps shell number to atomic center
  std::vector<std::array<int,2>> mapCen2Bf_; ///< Maps atomic center to first basis function
  std::vector<ReferenceShell>    refShells_; ///< Set of reference shells for given basis
  std::vector<libint2::Shell>    shells_   ; ///< Local basis storage (in shells)

  std::string basisPath_; ///< Path to the basis file

  std::unique_ptr<ifstream> basisFile_; ///< The file containing the basis defintion (G94)

  FileIO *fileio_; ///< FileIO Object for output file

  int printLevel_; ///< Level of print for basis set object


public:
  AOCartesian *ao;
  ShellPair   *shellPairs;
  Shell       *shells_old;
  int         *sortedShells;



  bool haveMapSh2Bf  ; ///< (?) The map from shells to basis functions has been made
  bool haveMapSh2Cen ; ///< (?) The map from shells to atomic centers has been made
  bool haveMapCen2Bf ; ///< (?) The map from atomic centers to basis functions has been made

  std::unique_ptr<RealMatrix> shBlkNormAlpha; ///< Shell Block (Inf) norm for Alpha matrix
  std::unique_ptr<RealMatrix> shBlkNormBeta;  ///< Shell Block (Inf) norm for Beta  matrix

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
    this->doSph_           = false  ;
    this->haveMapSh2Bf     = false  ;
    this->haveMapSh2Cen    = false  ;
    this->basisFile_       = nullptr;
    this->fileio_          = NULL   ;
    this->radCutSh_        = NULL   ; 
    this->printLevel_      = 1      ;
  };

  /**
   *  Copy Constructor
   */ 
  BasisSet(const BasisSet &basis) : BasisSet(){
    this->doSph_           = basis.doSph_       ; 
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
  inline  double getradCutSh(int iShell) {return this->radCutSh_[iShell];   }; ///< Return radCutSh
  
  template <typename T> double * basisEval(int,std::array<double,3>,T*);
  template <typename T> double * basisEval(libint2::Shell&,T*);
  template <typename T> double * basisProdEval(libint2::Shell,libint2::Shell,T*);
  std::vector<bool> MapGridBasis(cartGP pt);  ///< Create a Mapping of basis over grid points
  void     radcut(double thr, int maxiter, double epsConv);   //return all shell cut off radius
  double   fSpAv (int iop,int l, double alpha, double r);   //Evaluate Spheric Average of a Shell
  double   fRmax (int l, double alpha, double thr, double epsConv, int maxiter);   //Evaluate Spheric Average of a Shell
  inline libint2::Shell      shells(int i) {return this->shells_[i];    };
  inline int                nLShell(int L) {return this->nLShell_[L];   };
  inline int               mapSh2Bf(int i) {return this->mapSh2Bf_[i];  };
  inline int               mapSh2Cen(int i) {return this->mapSh2Cen_[i];};
  inline std::array<int,2> mapCen2Bf(int i) {return this->mapCen2Bf_[i];};

  inline void resetMapSh2Bf() {this->mapSh2Bf_.clear(); this->haveMapSh2Bf  = false;};
  inline void resetMapSh2Cen(){this->mapSh2Cen_.clear();this->haveMapSh2Cen = false;};
  inline void resetMapCen2Bf(){this->mapCen2Bf_.clear();this->haveMapCen2Bf = false;};
  

  inline void setBasisPath( std::string str){ this->basisPath_  = str;};
  inline void setPrintLevel(int i          ){ this->printLevel_ = i  ;};


  void printInfo();   ///< Print all info
  void printMeta();   ///< Print Meta data (nBasis, etc)
  void printHeader(); ///< Print header
  void printBasis();  ///< Print basis definition

  void basisSetRead(FileIO *,Molecule *,Controls *); ///< Parse and construct basis from file
  void findBasisFile(std::string);         ///< Attempt to find basis set file
  void parseGlobal();                      ///< Parse basis set file, generate reference
  void constructLocal(Molecule *);         ///< Construct local basis defintion
  void computeMeta();                      ///< Compute meta data (nBasis, etc)
  void makeMapSh2Bf(int);                  ///< generate mapSh2Bf
  void makeMapSh2Cen(Molecule *);          ///< generate mapSh2Cen
  void makeMapCen2Bf(int,Molecule *);          ///< generate mapCen2Bf
  void renormShells();                     ///< Renormalize Libint2::Shell set
  std::vector<libint2::Shell> uncontractBasis(); ///< Unconctract the basis
  template<typename TMat> void computeShBlkNorm(bool,int,const TMat*, const TMat*);

  void constructExtrn(Molecule *, BasisSet *); ///< Generate new basis from refernce shells
  void genUCvomLocal(BasisSet *);

  inline void makeMaps(int nTCS, Molecule* mol){
    this->makeMapSh2Bf(nTCS);
    this->makeMapSh2Cen(mol);
    this->makeMapCen2Bf(nTCS,mol);
  }

  // Python API
  void Wrapper_constructLocal(Molecule&);
  void Wrapper_makeMaps(int,Molecule&);

}; // class BasisSet
}; // namespace ChronusQ
#endif
