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
#ifndef INCLUDED_RESPONSE
#define INCLUDED_RESPONSE

#include <global.h>
#include <cerr.h>
#include <qn.h>
#include <singleslater.h>
#include <mointegrals.h>



namespace ChronusQ {

// Enumerate possible Response Functions
enum RESPONSE_TYPE {
  NOMETHOD,
  CIS,
  RPA,
  PPRPA,
  PPTDA,
  STAB
};

// Enumerate Response Function Classes
enum RESPONSE_CLASS {
  NOCLASS,
  FOPPA,
  SOPPA,
  TOPPA,
  PPPA
};

// Enumerate Response Function Partitioning Schemes
enum RESPONSE_MATRIX_PARTITIONING {
  SPIN_SEPARATED,
  SPIN_ADAPTED    ///< RHF ONLY
};

enum RESPONSE_PARTITION {
  FULL,
  FULL_A_PPTDA,
  FULL_C_PPTDA,
  SINGLETS,
  TRIPLETS,
  AA_PPRPA,
  AB_PPRPA,
  BB_PPRPA,
  AAA_PPTDA,
  AAB_PPTDA,
  ABB_PPTDA,
  CAA_PPTDA,
  CAB_PPTDA,
  CBB_PPTDA,
  A_PPTDA_SINGLETS,
  C_PPTDA_SINGLETS,
  PPRPA_SINGLETS,
  A_PPTDA_TRIPLETS,
  C_PPTDA_TRIPLETS,
  PPRPA_TRIPLETS
};

enum RESPONSE_JOB_TYPE {
  NOJOB,
  EIGEN,
  DYNAMIC
};
template<typename T>
class Response : public QNCallable<T>, public Quantum<T> {

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
  int Ref_;    ///< Reference ID

  /** Job Control for response calcluation **/
  bool useIncoreInts_;   ///< Use Incore Integrals
  bool doFull_;          ///< Do full matrix problem (in-core)
  bool debugIter_;       ///< Diagonalize full incore matrix iteratively
  bool doTDA_;           ///< Invoke TDA
  RESPONSE_TYPE iMeth_;  ///< Response type 
  RESPONSE_CLASS iClass_; ///< Response Function Class
  RESPONSE_JOB_TYPE iJob_; ///< Response Job Type
  // nSek, nGuess and nSingleDim inherited from QNCallable

  std::vector<int> nMatDim_;///< Dimensions of the different Response Matricies
  std::vector<RESPONSE_PARTITION> iMatIter_; ///< Response Matrix Partition
  RESPONSE_MATRIX_PARTITIONING iPart_;///< Type of Response matrix Partitioning
  std::vector<H5::DataSet*> guessFiles_;

  RESPONSE_PARTITION currentMat_;
  std::map<RESPONSE_PARTITION,size_t> partitionIndexMap_;

  bool doSinglets_;      ///< (?) Find NSek Singlet Roots (depends on SA)
  bool doTriplets_;      ///< (?) Find NSek Triplet Roots (depends on SA)

  // Job control specific to PPRPA
  // ** Only make sense with SP **
  bool doAllAlpha_;      ///< (?) Do All Alpha Block ** PPRPA Only **
  bool doMixedAB_;       ///< (?) Do Mixed Alpha-Beta Block ** PPRPA ONLY **
  bool doAllBeta_;       ///< (?) Do All Beta Block ** PPRPA ONLY **
  

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
  std::vector<VectorXd> rmDiag_; ///< Diagonal elements of response matrix

  // Solution Quantities
  std::vector<TMat>             transDen_;    ///< Transition Density (MO)
  std::vector<VectorXd>         frequencies_; ///< Transition Frequencies

  // Properties
  std::vector<double>               oscStrength_; ///< Oscillator Strengths
  std::vector<std::array<double,3>> transDipole_; ///< Transition Dipole
  

public:
  /** Constructors **/
  /**
   *  Default Constructor loads default values
   *  Inherits QNCallable default constructor
   */ 
  Response() : QNCallable<T>(), Quantum<T>(){
    // Intialize pointers to NULL
    this->fileio_       = NULL;
    this->mointegrals_  = NULL;
    this->singleSlater_ = NULL;

    // Zero out meta data to be initialized by SingleSlater
    this->nBasis_     = 0;
//  this->nTCS_       = 0;
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
    this->iMeth_         = NOMETHOD;
    this->iClass_        = NOCLASS;
    this->iJob_          = EIGEN;
    this->useIncoreInts_ = false;
    this->doFull_        = false;
    this->debugIter_     = false;
    this->doTDA_         = false;
    this->iPart_         = SPIN_SEPARATED;
    this->doSinglets_    = true;
    this->doTriplets_    = true;
    this->doAllAlpha_    = true;
    this->doMixedAB_     = true;
    this->doAllBeta_     = false;
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
  void initMeta();
  void initMetaFOPPA();
  void initMetaSOPPA();
  void initMetaTOPPA();
  void initMetaPPRPA();
  void initPSCFDims();
  void initRMu();

  // Setters
  inline void setMeth(RESPONSE_TYPE n) { this->iMeth_  = n;            };
  inline void doFull()                 { this->doFull_ = true;         };
  inline void doTDA()                  { this->doTDA_  = true;         };
  inline void doSA()                   { this->iPart_  = SPIN_ADAPTED; };

  // Getters
  inline SingleSlater<T> * singleSlater(){return this->singleSlater_;};
  inline RESPONSE_TYPE     Meth(){return this->iMeth_;};

//inline std::vector<TMat> transDen(){return this->transDen_;};
//inline std::vector<VectorXd> frequencies(){return this->frequencies_;};
//inline std::vector<int> nMatDim(){return this->nMatDim_;};
  
  template<RESPONSE_PARTITION U>
  inline TMat& transDen(){
    auto search = this->partitionIndexMap_.find(U);
    if(search == this->partitionIndexMap_.end())
      CErr("Requested Transition Density Not Available",
          this->fileio_->out);
    return this->transDen_[search->second];
  };

  template<RESPONSE_PARTITION U>
  inline VectorXd& frequencies(){
    auto search = this->partitionIndexMap_.find(U);
    if(search == this->partitionIndexMap_.end())
      CErr("Requested Freqencies Not Available",
          this->fileio_->out);
    return this->frequencies_[search->second];
  };

  template<RESPONSE_PARTITION U>
  inline int nMatDim(){
    auto search = this->partitionIndexMap_.find(U);
    if(search == this->partitionIndexMap_.end())
      CErr("Requested Transition Density Not Available",
          this->fileio_->out);
    return this->nMatDim_[search->second];
  };
  
  // IO Related
  void printInfo();

  // Run a Response Calculation
  inline void doResponse() {
    // Initialize MetaData
    this->initMeta();

    this->alloc();
    // Print Response Module Info
    this->printInfo();

    if(this->doFull_) this->full();
    else this->IterativeResponse();
  };

  // In-Core Related
  inline void full() {
    if(this->singleSlater_->aointegrals()->integralAlgorithm != 
       AOIntegrals::INCORE)
      CErr("Full Response Problems Require In-Core Integrals",
        this->fileio_->out);
 
    if(this->iClass_ == FOPPA) this->fullFOPPA();
    else if(this->iClass_ == SOPPA) this->fullSOPPA();
    else if(this->iClass_ == TOPPA) this->fullTOPPA();
    else if(this->iClass_ == PPPA)  this->fullPPRPA();
  };
  void fullFOPPA();
  void fullSOPPA(){;}; // NYI
  void fullTOPPA(){;}; // NYI
  void fullPPRPA();

  // QN Related
  inline void IterativeResponse(){
    this->formDiag();
    this->formGuess();
    this->nGuess_     = this->nSek_; // Quick Hack to get things working
    if(this->iMeth_ == RESPONSE_TYPE::RPA){
      this->needsLeft_ = true;
    }

    std::function<H5::DataSet*(const H5::PredType&,std::string&,
      std::vector<hsize_t>&)> fileFactory = 
        std::bind(&FileIO::createScratchPartition,this->fileio_,
        std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);

    for(auto iMat = 0; iMat != iMatIter_.size(); iMat++){
      this->currentMat_ = this->iMatIter_[iMat];
      this->nSingleDim_ = this->nMatDim_[iMat];  
      this->solutionVecR_ = &this->transDen_[iMat];
      this->omega_        = &this->frequencies_[iMat];
      this->guessFile_    = this->guessFiles_[iMat];
      this->diag_         = &this->rmDiag_[iMat];

      QuasiNewton2<T> qn(this,fileFactory);
      if(this->iMeth_ == RESPONSE_TYPE::RPA){
        qn.setMatrixType( QNMatrixType::HERMETIAN_GEP           );
        qn.setAlgorithm(  QNSpecialAlgorithm::SYMMETRIZED_TRIAL );
      }
      qn.run();
      cout << endl << endl;
    }
  };

  // Required for QNCallable
  inline void linearTrans(TMap &VR,TMap &VL,TMap &SR,TMap &SL,
                          TMap &RR,TMap &RL){
    
    if(this->iClass_ == FOPPA) this->linearTransFOPPA(VR,VL,SR,SL,RR,RL);
    else if(this->iClass_ == SOPPA) this->linearTransSOPPA(VR,VL,SR,SL,RR,RL);
    else if(this->iClass_ == TOPPA) this->linearTransTOPPA(VR,VL,SR,SL,RR,RL);
    else if(this->iClass_ == PPPA)  this->linearTransPPRPA(VR,VL,SR,SL,RR,RL);
  };
  void linearTransFOPPA(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &);
  void linearTransSOPPA(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &){;};
  void linearTransTOPPA(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &){;};
  void linearTransPPRPA(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &){;};

  // Required for QNCallable
  void formGuess();

  inline void formDiag(){
    if(this->iClass_ == FOPPA) this->formDiagFOPPA();
    else if(this->iClass_ == SOPPA) this->formDiagSOPPA();
    else if(this->iClass_ == TOPPA) this->formDiagTOPPA();
    else if(this->iClass_ == PPPA)  this->formDiagPPRPA();
  };
  void formDiagFOPPA();
  void formDiagSOPPA(){;};
  void formDiagTOPPA(){;};
  void formDiagPPRPA();


  // Allocation Related
  void alloc();

  // Misc
  void checkValid();      ///< Check validity of calculation with meta data
  void checkValidFOPPA(); ///< Checks specific to FOPPA
  void checkValidSOPPA(); ///< Checks specific to SOPPA
  void checkValidTOPPA(); ///< Checks specific to TOPPA
  void checkValidPPRPA(); ///< Checks specific to PPRPA

  // Transformations of Transition Densities
  inline void formAOTransDen(TVecMap &TVec, TMat &TAOA, TMat &TAOB){
    if(this->iClass_ == FOPPA) this->formAOTransDenFOPPA(TVec,TAOA,TAOB);
    else if(this->iClass_ == PPPA) this->formAOTransDenPPRPA(TVec,TAOA,TAOB);
  };
  inline void formMOTransDen(TVecMap &TVec, TMat &TMOA, TMat &TMOB){
    if(this->iClass_ == FOPPA) this->formMOTransDenFOPPA(TVec,TMOA,TMOB);
    else if(this->iClass_ == PPPA) this->formMOTransDenPPRPA(TVec,TMOA,TMOB);
  };
  void formAOTransDenFOPPA(TVecMap&, TMat&, TMat&);
  void formMOTransDenFOPPA(TVecMap&, TMat&,TMat&);
  void formAOTransDenPPRPA(TVecMap&, TMat&, TMat&);
  void formMOTransDenPPRPA(TVecMap&, TMat&,TMat&);
  void formDensity(){;};
  void placeOccVir(TVecMap&, TMat&, TMat&);
  void placeVirOcc(TVecMap&, TMat&, TMat&);
  void placeOccOcc(TVecMap&, TMat&, TMat&);
  void placeVirVir(TVecMap&, TMat&, TMat&);
  void retrieveOccVir(TVecMap&, TMat&, TMat&);
  void retrieveVirOcc(TVecMap&, TMat&, TMat&);
  void retrieveOccOcc(TVecMap&, TMat&, TMat&);
  void retrieveVirVir(TVecMap&, TMat&, TMat&);

}; // class Response
#include <response_alloc.h>
#include <response_meta.h>
#include <response_io.h>
}; // namespace ChronusQ

#endif
