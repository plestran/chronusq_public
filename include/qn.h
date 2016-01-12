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
#ifndef INCLUDED_QUASINEWTON2 
#define INCLUDED_QUATINEWTON2
#include <global.h>
#include <cerr.h>

namespace ChronusQ { 
  template<typename T>
  class QNCallable {
  protected:
    // Useful Eigen Typedefs
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;
    typedef Eigen::Map<TMat> TMap;
    typedef Eigen::Map<TVec> TVecMap;

    // Dimension Variables
    size_t nSingleDim_;  ///< Dimension of the problem
    size_t nSek_;        ///< Number of desired solution vectors
    size_t nGuess_;      ///< Number of initial guess vectors

    // Booleans for QN setup
    bool generateGuess_;  ///< (?) Generate identity guess
    bool allocSolution_;  ///< (?) Allocate space for solution
    bool needsLeft_;      ///< (?) Left vectors need to be taken into account
    bool isFullMatrix_;   ///< (?) Full matrix passes (DEBUG OPTION)
    bool directDiag_;     ///< (?) Arithmetic with direct access of diagonals
   
    // In-Core storage
    TMat        * fullMatrix_;   ///< Full Matrix (DEBUG PURPOSES)
    TMat        * fullMetric_;   ///< Full Metric (DEBUG PURPOSES)
    TMat        * solutionVecR_; ///< (right) solution vectors
    TMat        * solutionVecL_; ///< (left)  solution vectors
    VectorXd    * omega_;        ///< Frequencies for solution vectors
    VectorXd    * diag_;         ///< Diagonal elements of problem

    // Disk storage
    H5::H5File  * scratchFile_; ///< Scratch file
    H5::DataSet * guessFile_;   ///< File to store guess vectors

  public:
    /**
     *  Defualt Constructor for QNCallable
     *
     *  Default constructor loads defaults
     */ 
    QNCallable() {
      this->nSingleDim_   = 0;
      this->nSek_         = 0;
      this->nGuess_       = 0;
      this->fullMatrix_   = NULL;
      this->solutionVecR_ = NULL;
      this->solutionVecL_ = NULL;
      this->omega_        = NULL;
      this->diag_         = NULL;
      this->scratchFile_  = NULL;
      this->guessFile_    = NULL;

      this->generateGuess_ = false;
      this->allocSolution_ = false;
      this->needsLeft_     = false;
      this->isFullMatrix_  = false;
      this->directDiag_    = false;
    };

    QNCallable(TMat *A, int NSek, int NGuess) : QNCallable(){
      this->fullMatrix_ = A;
      this->nSek_       = NSek;
      this->nGuess_     = NGuess;

      this->nSingleDim_    = A->rows();
      this->generateGuess_ = true;
      this->allocSolution_ = true;
      this->isFullMatrix_  = true;
    };

    QNCallable(TMat *A, int NSek) : QNCallable(A,NSek,0) {
      this->nGuess_ = 2*NSek;
    };

    QNCallable(TMat *A, int NSek, int NGuess, H5::H5File *scr) :
      QNCallable(A,NSek,NGuess) {
      this->scratchFile_ = scr;
    }

    inline void initQN() {
      if(this->generateGuess_) this->generateGuess();
      if(this->allocSolution_) this->allocSolution();
    }

    /** Function Declarations **/

    // Setters
    inline void setNSek(int n)  { this->nSek_   = n;};
    inline void setNGuess(int n){ this->nGuess_ = n;};

    // Getters
    inline int nGuess()     { return this->nGuess_    ; };
    inline int nSek()       { return this->nSek_      ; };
    inline int nSingleDim() { return this->nSingleDim_; };
    
    inline bool needsLeft() { return this->needsLeft_;  };

    inline VectorXd * omega()        { return this->omega_       ; };
    inline TMat     * solutionVecR() { return this->solutionVecR_; };
    inline TMat     * solutionVecL() { return this->solutionVecL_; };
    inline VectorXd * diag()         { return this->diag_        ; };

    inline H5::H5File  * scratchFile(){ return this->scratchFile_; };
    inline H5::DataSet * guessFile()  { return this->guessFile_;   };  

    virtual void linearTrans(TMap &,TMap &,TMap &,TMap &,TMap &,TMap &) = 0;
    virtual void formGuess() = 0;
    virtual void formDiag()  = 0;

    
  }; // class QNCallable

  // Enums for Job Control
  enum QNProblemType {
    DIAGONALIZATION,
    LINEAR_SOLVE
  };

  enum QNMatrixType {
    HERMETIAN,
    HERMETIAN_GEP,
    NON_HERMETIAN
  };

  enum QNGuessType {
    RESIDUAL_DAVIDSON,
    RESIDUAL_OLSEN
  };

  enum QNSpecialAlgorithm {
    NOT_SPECIAL,
    SYMMETRIZED_TRIAL
  };

  template <typename T>
  class QuasiNewton2 {

    // Useful Eigen Typedefs
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;
    typedef Eigen::Map<TMat> TMap;
    typedef Eigen::Map<TVec> TVecMap;

    QNCallable<T> * qnObj_;

    // Dimension Variables
    size_t maxSubSpace_; ///< Maximum dimension of iterative subspace

    // Iteration related variables
    int maxMicroIter_;
    int maxMacroIter_;

    // Various tolerance values
    double residualTol_;

    // Job Control
    QNProblemType      problemType_;
    QNMatrixType       matrixType_; 
    QNGuessType        guessType_;
    QNSpecialAlgorithm specialAlgorithm_;

    // MetaData about QN Calculation
    int  nMicroIter_;
    int  nMacroIter_;
    int  nTotalIter_;
    bool isConverged_;

    // Output
    std::ostream * out_;
    
    // Factory to create new files
    std::function<
      H5::DataSet*(const H5::PredType&,std::string&,std::vector<hsize_t>&)
    > genScrFile_;


    /** Scratch Partitions **/

    // Constant data-type storage
    int        LRWORK; ///< Length of RWORK
    double   * RWORK_; ///< Real workspace for ZGEEV (Complex Only!)
    double   * ERMem_; ///< Storage of Re[omega] for subspace (Real Only!)
    double   * EIMem_; ///< Storage of Im[omega] for subspace (Real Only!)
    dcomplex * ECMem_; ///< Storage of omega for subspace (Complex Only!)

    // Trial vector storage
    T* TRMem_; ///< Right trial vectors
    T* TLMem_; ///< Left trial vectors

    // Matrix/Metric-Vector Product storage
    T* SigmaRMem_; ///< Matrix-Vector Product on Right trial vectors
    T* SigmaLMem_; ///< Matrix-Vector Product on Left trial vectors
    T* RhoRMem_;   ///< Metric-Vector Product on Right trial vectors
    T* RhoLMem_;   ///< Metric-Vector Product on Left trial vectors

    // Full projection storage
    T* XTSigmaRMem_; ///< Full projection of Right Sigma onto subspace
    T* XTSigmaLMem_; ///< Full projection of Left Sigma onto subspace
    T* XTRhoRMem_;   ///< Full projection of Right Rho onto subspace
    T* XTRhoLMem_;   ///< Full projection of Left Rho onto subspace

    // Residuals
    T* ResRMem_; ///< Right Residuals
    T* ResLMem_; ///< Left Residuals

    // Reconstructed Solution Vectors
    T* URMem_; ///< Reconstructed Right solution vectors
    T* ULMem_; ///< Reconstructed Left solution vectors

    // LAPACK Scratch Space
    T* WORK;   ///< LAPACK scratch space
    int LWORK; ///< Length of WORK

    // Special Algorithm Memory
       
    // Symmetrized Trial Vectors
    T* ASuperMem_; ///< Matrix in basis of gerade and ungerade vectors
    T* SSuperMem_; ///< Metrix in basis of gerade and ungedade vectors
      
      
    void allocScr();
    void allocScrSpecial();
    void cleanupScr();
    void cleanupScrSpecial();


    // Trial vector storage (Disk)
    H5::DataSet * TRFile_; ///< Right trial vectors
    H5::DataSet * TLFile_; ///< Left trial vectors

    // Matrix/Metric-Vector Product storage (Disk)
    H5::DataSet * SigmaRFile_; ///< Matrix-Vector Product on Right trial vectors
    H5::DataSet * SigmaLFile_; ///< Matrix-Vector Product on Left trial vectors
    H5::DataSet * RhoRFile_;   ///< Metric-Vector Product on Right trial vectors
    H5::DataSet * RhoLFile_;   ///< Metric-Vector Product on Left trial vectors

/*  Should always be able to store reduced dimension in Core
    // Full projection storage (Disk)
    H5::DataSet * XTSigmaRFile_; ///< Full projection of Right Sigma onto subspace
    H5::DataSet * XTSigmaLFile_; ///< Full projection of Left Sigma onto subspace
    H5::DataSet * XTRhoRFile_;   ///< Full projection of Right Rho onto subspace
    H5::DataSet * XTRhoLFile_;   ///< Full projection of Left Rho onto subspace
*/

    // Residuals (Disk)
    H5::DataSet * ResRFile_; ///< Right Residuals
    H5::DataSet * ResLFile_; ///< Left Residuals


    // Reconstructed Solution Vectors (Disk)
    H5::DataSet * URFile_; ///< Reconstructed Right solution vectors
    H5::DataSet * ULFile_; ///< Reconstructed Left solution vectors

    void iniScratchFiles();
    void writeTrialVectors(const int);
    void readTrialVectors(const int);

  public:
    // Set private data
    void setMatrixType(QNMatrixType      type){this->matrixType_       = type;};
    void setProblemType(QNProblemType    type){this->problemType_      = type;};
    void setGuessType(QNGuessType        type){this->guessType_        = type;};
    void setAlgorithm(QNSpecialAlgorithm type){this->specialAlgorithm_ = type;};

    // Run the QN Calculation
    void run();
    void runMicro();

    // Procedural Functions
    void readGuess();
    void checkOrthogonality(int&);
    void formLinearTrans(const int, const int);
    void fullProjection(const int);
    void reducedDimDiag(const int);
    void reconstructSolution(const int);
    void generateResiduals(const int);
    std::vector<bool> checkConvergence(const int, int&); 
    void formNewGuess(std::vector<bool>&, int&, int&, int&, int&);


    // Orthogonality Check Routines
    void checkLinearDependence(int &);
    void orthogonalize(int);


    // Diagonalization Routines
    void stdHermetianDiag(const int);
    void stdNonHermetianDiag(const int);

    // Guess Generation Routines
    void davResidualGuess(    T,const TMap&,TMap&,const TMap&,TMap&);
    void davStdResidualGuess( T,const TMap&,TMap&);
    void davSymmResidualGuess(T,const TMap&,TMap&,const TMap&,TMap&);

    // Special Algorithm Routines
    void symmetrizeTrial();
    void buildSuperMatricies(const int);

    #include <qn_constructors.h>

  }; // class QuasiNewton2
  #include <qn_memory.h>
  #include <qn_procedural.h>
  #include <qn_special.h>
}; // namespace ChronusQ

#endif

