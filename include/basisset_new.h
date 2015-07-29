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
#ifndef INCLUDED_BASISSETNEW
#define INCLUDED_BASISSETNEW
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
/*
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
*/

class BasisSetNew{
  struct ReferenceShell{
    int atomicNumber;
    int index;
    std::vector<libint2::Shell> shells;
  }; 
  int  nBasis_           = 0    ;
  int  nPrimitive_       = 0    ;
  int  nShell_           = 0    ;
  int  nShellPair_       = 0    ;
  int  maxPrim_          = 0    ;
  int  maxL_             = 0    ;
  bool haveMapSh2Bf_     = false;
  bool haveMapSh2Cen_    = false;
  bool doSph_            = false;

  std::vector<int>            nLShell_;
  std::vector<int>            mapSh2Bf_;
  std::vector<int>            mapSh2Cen_;
  std::vector<libint2::Shell> shells_;
  std::vector<ReferenceShell> refShells_;

  std::string basisPath_;

  std::unique_ptr<ifstream> basisFile_;


  FileIO *fileio_;


public:
/*
  AOCartesian *ao;
  ShellPair   *shellPairs;
  Shell       *shells;
  int         *sortedShells;
*/

  std::unique_ptr<RealMatrix> shBlkNormAlpha; 
  std::unique_ptr<RealMatrix> shBlkNormBeta; 

  BasisSetNew(int nBasis=0, int nShell=0){;};
  ~BasisSetNew(){ // FIXME need to move these over to unique_ptr
/*
    delete[] ao;
    delete[] shellPairs;
    delete[] shells;
    delete[] sortedShells;
*/
    this->basisFile_->close();
  };

  inline int     nBasis()  {return this->nBasis_;};
  inline int nPrimitive()  {return this->nPrimitive_;};
  inline int     nShell()  {return this->shells_.size();};
  inline int nShellPair()  {return this->nShellPair_;};
  inline int nLShell(int L){return this->nLShell_[L];};

  inline void setBasisPath(std::string str){ this->basisPath_ = str;};


  void printInfo();
  void printMeta();
  void printHeader();
  void printBasis();

  void basisSetRead(FileIO *, Molecule *);
  void findBasisFile(std::string);
  void parseGlobal();
  void constructLocal(Molecule *);
  void computeMeta();
  void makeMapSh2Bf();
  void makeMapSh2Cen(Molecule *);
  template<typename TMat> void computeShBlkNorm(bool, const TMat*, const TMat*);

  BasisSetNew constructExtrn(Molecule *);


}; // class BasisSet
}; // namespace ChronusQ
#endif
