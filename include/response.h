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
#ifndef INCLUDED_RESPONSE
#define INCLUDED_RESPONSE

#include <global.h>
#include <cerr.h>
#include <qn.h>
#include <singleslater.h>
#include <mointegrals.h>



namespace ChronusQ {

// Enumerate possible Response Jobs
enum RESPONSE_TYPE {
  NOMETHOD,
  CIS,
  RPA,
  PPRPA,
  PPATDA,
  PPCTDA,
  STAB
};
template<typename T>
class Response : public QNCallable<T> {

  /** Useful TypeDefs **/
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1>                TVec;
  typedef Eigen::Map<TVec>                          TVecMap;
  typedef Eigen::Map<TMat>                          TMap;
  typedef Tensor<T,Range3d>                         TTensor3d;
  typedef Tensor<T,Range4d>                         TTensor4d;

  /** Connection to other CQ classes **/
  FileIO          * fileio_;
  SingleSlater<T> * singleSlater_; 
  MOIntegrals<T>  * mointegrals_;

  /** Meta inherited from SingleSlater **/
  int nBasis_; ///< Number of Basis Functions
  int nTCS_;   ///< Number of components
  int Ref_;    ///< Reference ID

  /** Job Control for response calcluation **/
  bool useIncoreInts;   ///< Use Incore Integrals
  bool doFull;          ///< Do full matrix problem (in-core)
  bool debugIter;       ///< Diagonalize full incore matrix iteratively
  bool doTDA;           ///< Invoke TDA
  RESPONSE_TYPE iMeth_; ///< Response type 
  // nSek, nGuess and nSingleDim inherited from QNCallable
  

  /** Derived Dimensions Post-SCF **/
  // Single Dimensions
  int       nOA_; ///< Number of occupied orbitals (Alpha)
  int       nVA_; ///< Number of virtual orbitals (Alpha)
  int       nOB_; ///< Number of occupied orbitals (Beta)
  int       nVB_; ///< Number of virtual orbitals (Beta)
  int       nO_;  ///< Total number of occupied orbitals NOA + NOB
  int       nV_;  ///< Total number of virtual orbitals NVA + NVB

  // Quadratic Dimenstions
    
  // Occupied-Occupied
  int       nOAOA_;      ///< NOA * NOA
  int       nOBOB_;      ///< NOB * NOB
  int       nOAOB_;      ///< NOA * NOB
  int       nOO_;        ///< NO * NO
  int       nOAOA_SLT_;  ///< NOA * (NOA - 1) / 2
  int       nOBOB_SLT_;  ///< NOB * (NOB - 1) / 2
  int       nOO_SLT_;    ///< NO * (NO - 1) / 2
  int       nOAOA_LT_;   ///< NOA * (NOA + 1) / 2
  int       nOO_LT_;     ///< NO * (NO + 1) / 2

  // Occupied-Virtual
  int       nOAVA_;      ///< NOA * NVA
  int       nOBVB_;      ///< NOB * NVB
  int       nOAVB_;      ///< NOA * NVB
  int       nOBVA_;      ///< NOB * NVA
  int       nOV_;        ///< NO * NV

  // Virtual-Virtual
  int       nVAVA_;      ///< NVA * NVA
  int       nVBVB_;      ///< NVB * NVB
  int       nVAVB_;      ///< NVA * NVB
  int       nVV_;        ///< NV * NV
  int       nVAVA_SLT_;  ///< NVA * (NVA - 1) / 2
  int       nVBVB_SLT_;  ///< NVB * (NVB - 1) / 2
  int       nVV_SLT_;    ///< NV * (NV - 1) / 2
  int       nVAVA_LT_;   ///< NVA * (NVA + 1) / 2
  int       nVV_LT_;     ///< NV * (NV + 1) / 2


  /** Misc values **/
  double rMu_; ///< Level shift, currently only used for PP-methods

  /** Internal Storage of desired quantities **/
  // Misc Required for QN
  std::unique_ptr<VectorXd> rmDiag_; ///< Diagonal elements of response matrix

  // Solution Quantities
  std::vector<TMat>             transDen_;    ///< Transition Density (MO)
  std::unique_ptr<VectorXd>     frequencies_; ///< Transition Frequencies

  // Properties
  std::vector<double>               oscStrength_; ///< Oscillator Strengths
  std::vector<std::array<double,3>> transDipole_; ///< Transition Dipole
  

public:
  /** Constructors **/
  /**
   *  Default Constructor loads default values
   *  Inherits QNCallable default constructor
   */ 
  Response() : QNCallable<T>(){
    // Intialize pointers to NULL
    this->fileio_       = NULL;
    this->mointegrals_  = NULL;
    this->singleSlater_ = NULL;
    this->rmDiag_       = nullptr;
    this->frequencies_  = nullptr;

    // Zero out meta data to be initialized by SingleSlater
    this->nBasis_     = 0;
    this->nTCS_       = 0;
    this->Ref_        = 0;

    // Zero out PSCF dimensions to be built later
    this->nOA_        = 0;
    this->nVA_        = 0;
    this->nOB_        = 0;
    this->nVB_        = 0;
    this->nOAVA_      = 0;      
    this->nOBVB_      = 0;      
    this->nOAVB_      = 0;      
    this->nOBVA_      = 0;      
    this->nVAVA_SLT_  = 0;  
    this->nVBVB_SLT_  = 0;  
    this->nVAVA_LT_   = 0;   
    this->nVAVA_      = 0;      
    this->nVBVB_      = 0;      
    this->nOAOA_SLT_  = 0;  
    this->nOBOB_SLT_  = 0;  
    this->nOAOA_LT_   = 0;   
    this->nOAOA_      = 0;      
    this->nOBOB_      = 0;      
    this->nVAVB_      = 0;      
    this->nOAOB_      = 0;      
    this->nO_         = 0;         
    this->nV_         = 0;         
    this->nOV_        = 0;        
    this->nVV_SLT_    = 0;    
    this->nVV_LT_     = 0;     
    this->nVV_        = 0;        
    this->nOO_SLT_    = 0;    
    this->nOO_LT_     = 0;     
    this->nOO_        = 0;        

    // Zero out Misc values
    this->rMu_ = 0.0;

    // Standard (default) values
    this->iMeth_        = NOMETHOD;
    this->useIncoreInts = false;
    this->doFull        = false;
    this->debugIter     = false;
    this->doTDA         = false;
  };

  // Dummy Destructor
  ~Response(){;};

  /** Communication Routines **/

  /**
   *  Initialize pointers to other CQ workers
   */ 
  inline void communicate(SingleSlater<T> &ss, MOIntegrals<T> &moints, 
    FileIO &fileio) {

    this->singleSlater_ = &ss;
    this->mointegrals_  = &moints;
    this->fileio_       = &fileio;

  }

  /**
   *  Checks that Response object has been initialized with proper CQ workers
   */ 
  inline void checkWorkers(){
    if(this->fileio_  == NULL) 
      CErr("Fatal: Must initialize SDResponse with FileIO Object");
    if(this->singleSlater_ == NULL) 
      CErr("Fatal: Must initialize SDResponse with SingleSlater Object",
           this->fileio_->out);
    if(this->mointegrals_ == NULL) 
      CErr("Fatal: Must initialize SDResponse with MOIntegrals Object",
           this->fileio_->out);
  }


  /** Meta data routines **/
  inline void initMeta() {
    this->checkWorkers();

    this->nBasis_         = this->singleSlater_->nBasis();
    this->nTCS_           = this->singleSlater_->nTCS();
    this->Ref_            = this->singleSlater_->Ref();
    this->nOA_            = this->singleSlater_->nOccA();
    this->nOB_            = this->singleSlater_->nOccB();
    this->nVA_            = this->singleSlater_->nVirA();
    this->nVB_            = this->singleSlater_->nVirB();

    this->nOAVA_          = this->nOA_*this->nVA_;
    this->nOBVB_          = this->nOB_*this->nVB_;
    this->nOAVB_          = this->nOA_*this->nVB_;
    this->nOBVA_          = this->nOB_*this->nVA_;
    this->nVAVA_SLT_      = this->nVA_*(this->nVA_-1)/2;
    this->nVBVB_SLT_      = this->nVB_*(this->nVB_-1)/2;
    this->nVAVA_LT_       = this->nVA_*(this->nVA_+1)/2;
    this->nVAVA_          = this->nVA_*this->nVA_;
    this->nVBVB_          = this->nVB_*this->nVB_;
    this->nOAOA_SLT_      = this->nOA_*(this->nOA_-1)/2;
    this->nOBOB_SLT_      = this->nOB_*(this->nOB_-1)/2;
    this->nOAOA_LT_       = this->nOA_*(this->nOA_+1)/2;
    this->nOAOA_          = this->nOA_*this->nOA_;
    this->nOBOB_          = this->nOB_*this->nOB_;
    this->nVAVB_          = this->nVA_*this->nVB_;
    this->nOAOB_          = this->nOA_*this->nOB_;
    this->nO_             = this->nOA_ + this->nOB_;
    this->nV_             = this->nVA_ + this->nVB_;
    this->nOV_            = this->nO_  * this->nV_;
    this->nVV_SLT_        = this->nV_*(this->nV_-1)/2;
    this->nVV_LT_         = this->nV_*(this->nV_+1)/2;
    this->nVV_            = this->nV_*this->nV_;
    this->nOO_SLT_        = this->nO_*(this->nO_-1)/2;
    this->nOO_LT_         = this->nO_*(this->nO_+1)/2;
    this->nOO_            = this->nO_*this->nO_;
  } //initMeta

  // Setters
  inline void setMeth(RESPONSE_TYPE n) { this->iMeth_ = n; };
  
  // IO Related
  // In-Core Related
  void incoreFOPPA();
  void incoreSOPPA();
  void incoreTOPPA();
  void incorePPRPA();

  // QN Related
  void IterativeResponse();
  void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &){;};
  void linearTransFOPPA(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &);
  void linearTransSOPPA(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &);
  void linearTransTOPPA(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &);
  void linearTransPPRPA(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &);
  void formGuess(){;};
  void formGuessFOPPA();
  void formGuessSOPPA();
  void formGuessTOPPA();
  void formGuessPPRPA();
  void formDiag(){;};
  void formDiagFOPPA();
  void formDiagSOPPA();
  void formDiagTOPPA();
  void formDiagPPRPA();


  // Misc
  void checkValid();      ///< Check validity of calculation with meta data
  void checkValidFOPPA(); ///< Checks specific to FOPPA
  void checkValidSOPPA(); ///< Checks specific to SOPPA
  void checkValidTOPPA(); ///< Checks specific to TOPPA
  void checkValidPPRPA(); ///< Checks specific to PPRPA
  
}; // class Response
}; // namespace ChronusQ

#endif
