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
#include <basisset.h>
#include <molecule.h>
#include <fileio.h>
#include <controls.h>
#include <tools.h>

#define MaxFmTPt 3201
#define MaxTotalL 18

namespace ChronusQ {
struct MolecularConstants{
  int nAtom; //number of Atoms
  int atomZ[MAXATOMS];
  double cart[3][MAXATOMS];
};

struct PairConstants{
  /* one-center info */
  int nPGTOs[2];   //number of PGTOs
  int aoIndex[2];  //location of PGTO in ao[]
  int L[2];        //total angular momentum of each AO
  int nBasis[2];   //number of Basis Function
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
class AOIntegrals{
  int       nBasis_;
  int       nTT_;
  int       **R2Index_;
  double	**FmTTable_;

  std::shared_ptr<BasisSet>    	basisSet_;
  std::shared_ptr<Molecule>   	molecule_;
  std::shared_ptr<FileIO>      	fileio_;
  std::shared_ptr<Controls>    	controls_;

  PairConstants 	    *pairConstants_;
  MolecularConstants	*molecularConstants_;
  QuartetConstants	    *quartetConstants_;

//dbwys
#ifdef USE_LIBINT
  void OneEDriver(libint2::OneBodyEngine::integral_type);
#endif
//dbwye

public:
  // these should be protected
  std::shared_ptr<RealMatrix>     twoEC_;
  std::shared_ptr<RealMatrix>     twoEX_;
  std::shared_ptr<RealMatrix>     oneE_;
  std::shared_ptr<RealMatrix>     overlap_;
  std::shared_ptr<RealMatrix>     kinetic_;
  std::shared_ptr<RealMatrix>     potential_;
  std::shared_ptr<RealMatrix>     schwartz_;
  std::shared_ptr<Tensor<double>> aoERI_;

  bool		haveAOTwoE;
  bool		haveAOOneE;
  bool          haveSchwartz;


  // Timing Stats
  std::chrono::duration<double> OneED;
  std::chrono::duration<double> SED;
  std::chrono::duration<double> TED;
  std::chrono::duration<double> VED;
  std::chrono::duration<double> CoulD;
  std::chrono::duration<double> ExchD;
  std::chrono::duration<double> PTD;
  std::chrono::duration<double> SchwartzD;
  std::chrono::duration<double> DenShBlkD;
 
  AOIntegrals(){;};
  ~AOIntegrals(){
    twoEC_.reset();
    twoEX_.reset();
    oneE_.reset();
    overlap_.reset();
    kinetic_.reset();
    potential_.reset();
    schwartz_.reset();
    aoERI_.reset();
  };
  
  // initialization function
  void iniAOIntegrals(std::shared_ptr<Molecule>,std::shared_ptr<BasisSet>,
                      std::shared_ptr<FileIO>,std::shared_ptr<Controls>);

  inline double &twoEC(int i, int j, int k, int l){
    return (*twoEC_)(this->R2Index_[i][j],this->R2Index_[k][l]);
  };
  inline double &twoEC(int ij, int kl){
    return (*twoEC_)(ij,kl);
  };
  inline double &twoEX(int i, int j, int k, int l){
    return (*twoEX_)(this->R2Index_[i][j],this->R2Index_[k][l]);
  };
  inline double &twoEX(int ij, int kl){
    return (*twoEX_)(ij,kl);
  };
  
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
#ifdef USE_LIBINT
  void computeSchwartz();
#endif
//----------------------------------------//
// member functions in integrals_onee.cpp //
//----------------------------------------//
  void computeOverlapS();
  void computeKineticT();
  void computePotentialV();
  double oneehRRTSab(int,int*,int,int*,int*,int*);
  double oneehRRSab(int,int*,int,int*);
  double oneehRRVab(int,int*,int,int*);
  double oneevRRSa0(int,int*,int*,int*);
  double oneevRRVa0(int*,int,int,int*,int*,int*);
  double oneevRRTab(int,int*,int,int*,int*,int*);
  double oneevRRTa0(int,int*,int*,int*);
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
};
} // namespace ChronusQ

#endif
