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
#ifndef INCLUDED_QUASINEWTON
#define INCLUDED_QUASINEWTON
#include <global.h>
#include <cerr.h>
#include <sdresponse.h>
namespace ChronusQ {
/**
 * A class to setup and run various iterative Quasi-Newton type 
 * calculations such as:
 *   Davidson Diagonalization
 *   Quasi-Newton Linear Equation Solution
 */
template <typename T>
  class QuasiNewton{
    /** Member Variables **/
    // Useful Eigen Typedefs
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;
    typedef Eigen::Map<TMat> TMap;
    typedef Eigen::Map<TVec> TVecMap;

    // Boolean logic member variables
    bool isHermitian_;     // Hermitian Scheme
    bool doDiag_;          // Quasi-Newton Diagonalization (Davidson)
    bool doGEP_;           // Generalized Eigenproblem
    bool doLin_;           // Quasi-Newton Linear Equation Solve
    bool doResGuess_;      // Generate new vectors from residual
    bool symmetrizedTrial_;// Symmetrized Trial Vectors (Kauczor et al.)
    bool debug_;           // Enables various debug options
    bool cleanupMem_;      // True if memory cleanup is required
    bool isConverged_;     // True if the iterative calculation converged
    bool genGuess_;        // Generate Standard (identity) Guess
    bool doRestart_;

    // Integer size related member variables
    int N_;                // Dimension of the problem
    int nSek_;             // Number of desired roots for diagonalization
    int nGuess_;           // Numer of given (or to be generated) guess vectors
    int maxSubSpace_;      // Maximum dimension of the iterative subspace

    // Integer iteration related variables
    int maxIter_;          // Maximum number of miro iterations
    int MaxIter_;          // Maximum number of macro iterations
    int nMicroIter_;
    int nIter_;

    // Double precision variables related to convergence tolerence
    double resTol_;        // Residual norm tolerence for convergence

    // Templated pointers to Eigen storage
    TMat * A_;                // Pointer to full matrix to be diagonalized
    RealMatrix * diag_;             // Pointer to diagonal storage
    TMat * solutionVector_;   // Solution vectors at current iteration 
    VectorXd * solutionValues_;   // Solution values (eigen) at current iteration 
    std::unique_ptr<TMat> guessR_;   // Guess vectors (always local copy)
    std::unique_ptr<TMat> guessL_;   // Guess vectors (always local copy)

    SDResponse<T> * sdr_; // Pointer to SDResponse object
    /** Scratch Variables **/
    // Length of memory partitions
    int LenScr        ; 
    int LenRealScr    ;
    int LenSigma      ; 
    int LenRho        ;
    int LenXTSigma    ;
    int LenXTRho      ;
    int LenSuper      ;
    int LenU          ;
    int LenRes        ;
    int LenTVec       ;
    int LEN_LAPACK_SCR; 
    int LWORK         ;
    
    double * REAL_SCR;
    double * RWORK;
    double * RealEMem;
    // Templated pointers for scratch and paritions
    T * SCR        ; 
    T * SigmaRMem  ; 
    T * SigmaLMem  ; 
    T * XTSigmaRMem; 
    T * XTSigmaLMem; 
    T * RhoRMem    ; 
    T * RhoLMem    ; 
    T * XTRhoRMem  ; 
    T * XTRhoLMem  ; 
    T * ASuperMem  ;
    T * SSuperMem  ;
    T * SCPYMem    ;
    T * NHrProdMem ;
    T * BiOrthMem  ;
    T * URMem      ; 
    T * ULMem      ; 
    T * ResRMem    ; 
    T * ResLMem    ; 
    T * TVecRMem   ;         
    T * TVecLMem   ;         
    T * LAPACK_SCR ;
    T * ERMem      ;
    T * EIMem      ;
    T * WORK       ;
    /** Inline Functions **/
    inline void loadDefaults(){
      // Defaults for # of interations
      this->maxIter_          = 128;
      this->MaxIter_          = 20;
      this->nMicroIter_       = 0;
      this->nIter_            = 0;
      
      // Default tolerence 
      this->resTol_           = 5.0e-6;
   
      // Initialize values to be set later
      this->initBoolean();
      this->initSize();
      this->initPtr();
      this->initScrLen();
    };
    inline void initBoolean(){
      this->isHermitian_      = true;  // Default to Hermitian Scheme
      this->doDiag_           = true;  // Defualt to diagonalization
      this->doGEP_            = false; // Default standard eigenproblem if doing a diagonalization
      this->doLin_            = false; // Mildly redundent, as it will currently always be the opposite of doDiag
      this->doResGuess_       = true;  // Default to residual based guess
      this->symmetrizedTrial_ = false; // Default to unmodified trial vectors
      this->debug_            = false; // Default to terse output
      this->cleanupMem_       = false; // Assume that the space for results was allocated elsewhere
      this->isConverged_      = false; // Obviously the calculation doesn't begin converged...
      this->genGuess_         = true;  // Defualt generate identity guess
      this->doRestart_        = false; // class global variable to determine if a restart is needed
    }
    inline void initSize(){
      // Zero out all of the size dependent quantities
      this->N_                = 0;
      this->maxSubSpace_      = 0;
      this->nSek_             = 0;
      this->nGuess_           = 0; 
    }
    inline void initPtr(){
      // Initialize all pointers to some varient of NULL
      this->A_                = NULL;
      this->diag_             = NULL;
      this->solutionVector_   = NULL;
      this->solutionValues_   = NULL;
      this->guessR_           = nullptr;
      this->guessL_           = nullptr;
      this->sdr_              = NULL;
     
      this->SCR         = NULL; 
      this->REAL_SCR    = NULL; 
      this->SigmaRMem   = NULL; 
      this->SigmaLMem   = NULL; 
      this->XTSigmaRMem = NULL; 
      this->XTSigmaLMem = NULL; 
      this->RhoRMem     = NULL; 
      this->RhoLMem     = NULL; 
      this->XTRhoRMem   = NULL; 
      this->XTRhoLMem   = NULL; 
      this->ASuperMem   = NULL;
      this->SSuperMem   = NULL;
      this->SCPYMem     = NULL;
      this->NHrProdMem  = NULL;
      this->BiOrthMem   = NULL;
      this->URMem       = NULL; 
      this->ULMem       = NULL; 
      this->ResRMem     = NULL; 
      this->ResLMem     = NULL; 
      this->TVecRMem    = NULL;         
      this->TVecLMem    = NULL;         
      this->LAPACK_SCR  = NULL;
      this->ERMem       = NULL;
      this->EIMem       = NULL;
      this->WORK        = NULL;
      this->RWORK       = NULL;
      this->RealEMem    = NULL;
    };
  inline void initScrLen(){
/** DETERMINE LENGTH OF SCRATCH SPACE **
 * 
 *
 * (1) Linear transform of A onto right / gerade trial vectors 
 *
 * (2) Reduced dimension of A onto the right / gerade trial vector subspace
 *     (also stores the reduced dimension eigenvectors after diagonalization
 *      via DSYEV)
 *
 * (3) Approximate right / gerade eigenvectors 
 *
 * (4) Residuals of the right / gerade eigenvectors 
 *
 * (5) Right / gerade trial vectors 
 *
 * (6) Temp storage for a single precondictioned eigenvector 
 *
 * (7) Linear transform of the metric S onto the right / gerade trial vectors
 *
 * (8) Reduced dimension of S into the right / gerade trial vector subspace
 *
 * (9) Linear transform of A onto left / ungerade trial vectors 
 *
 * (10) Linear transform of the metric S onto the left /un gerade trial vectors
 *
 * (11) Reduced dimension of A onto the left / ungerade trial vector subspace
 *
 * (12) Reduced dimension of S into the left / ungerade trial vector subspace
 *
 * (13) Approximate left / ungerade eigenvectors 
 *
 * (14) Residuals of the left / ungerade eigenvectors 
 *
 * (15) Left / ungerade trial vectors 
 *
 * (16) Supermatrix of reduced dimension A onto the two subspaces
 *
 * (17) Supermatrix of reduced dimension S onto the two subspaces
 *
 * (23) Space for a copy of the S supermatrix to use for re-orthogonalization
 *
 * (24) Space for the non-hermetian product of S^-1 * A in the reduced dimention
 *
 * (25) Space for SX product for BiOrth wrt the metric
 */
    // Lenth of memory partitions
    this->LenSigma   = this->N_ * this->maxSubSpace_;
    this->LenRho     = this->N_ * this->maxSubSpace_;
    this->LenXTSigma = this->maxSubSpace_ * this->maxSubSpace_;
    this->LenXTRho   = this->maxSubSpace_ * this->maxSubSpace_;
    this->LenSuper   = 4 * this->maxSubSpace_ * this->maxSubSpace_;
    this->LenU       = this->N_ * this->maxSubSpace_;
    this->LenRes     = this->N_ * this->maxSubSpace_;
    this->LenTVec    = this->N_ * this->maxSubSpace_;
    this->LWORK      = 0;
    this->LEN_LAPACK_SCR = 0;
  
    // Init LenScr
    this->LenScr     = 0;
  
    // Memory length to hold intermediates
    this->LenScr += this->LenSigma;   // 1
    this->LenScr += this->LenXTSigma; // 2
    this->LenScr += this->LenU;       // 3
    this->LenScr += this->LenRes;     // 4 
    this->LenScr += this->LenTVec;    // 5
  
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      this->LenScr += this->LenRho;     // 7
      this->LenScr += this->LenXTRho;   // 8 
      this->LenScr += this->LenSigma;   // 9
      this->LenScr += this->LenRho;     // 10
      this->LenScr += this->LenXTSigma; // 11
      this->LenScr += this->LenXTRho;   // 12
      this->LenScr += this->LenU;       // 13
      this->LenScr += this->LenRes;     // 14
      this->LenScr += this->LenTVec;    // 15
      this->LenScr += this->LenSuper;   // 16
      this->LenScr += this->LenSuper;   // 17
      this->LenScr += this->LenSuper;   // 23
      this->LenScr += this->LenSuper;   // 25
    }

    if(!this->isHermitian_ && this->symmetrizedTrial_){
      this->LenScr += this->LenSuper; // 24
    }

    // LAPACK Storage Space Length
    this->initLAPACKScrLen();  
  }
  
    inline int stdSubSpace(){
      // Standard value for the maximum dimension of the
      // iterative subspace min(6*NSek,N/2)
      //return std::min(20*this->nSek_,this->N_/2);
      return std::min(250,this->N_/2);
//    return this->N_;
    };

    inline int stdNGuess(){
      // Standard value for the number of inital guess
      // vectors that are to be generated 2*NSek
      return 2*this->nSek_;
    }
    inline void allocGuess(){ 
      // Allocate space for local copy of the guess vectors
      this->guessR_ = std::unique_ptr<TMat>(new TMat(this->N_,this->nGuess_));
      if(this->doRestart_)
        this->guessL_ = std::unique_ptr<TMat>(new TMat(this->N_,this->nGuess_));
      if(this->genGuess_) this->identGuess();
    };
    inline void identGuess(){
      // Generate the identity (standard) guess
      (*this->guessR_) = TMat::Identity(this->N_,this->nGuess_);
      this->genGuess_ = false; // So we dont generate the guess on restart
    };
  inline void allocScr(){
    // Allocate scratch space
    this->SCR      = new T      [this->LenScr]; 
    this->REAL_SCR = new double [this->LenRealScr];
    std::memset((void*)this->SCR,     0.0,this->LenScr     * sizeof(T));
    std::memset((void*)this->REAL_SCR,0.0,this->LenRealScr * sizeof(double));
  
    // Partition scratch space
    this->SigmaRMem     = this->SCR;
    this->XTSigmaRMem   = this->SigmaRMem   + this->LenSigma;
    this->URMem         = this->XTSigmaRMem + this->LenXTSigma;
    this->ResRMem       = this->URMem       + this->LenU; 
    this->TVecRMem      = this->ResRMem     + this->LenRes;
    this->LAPACK_SCR    = this->TVecRMem    + this->LenTVec;
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      this->RhoRMem       = this->LAPACK_SCR  + this->LEN_LAPACK_SCR;
      this->XTRhoRMem     = this->RhoRMem     + this->LenRho;
      this->SigmaLMem     = this->XTRhoRMem   + this->LenXTRho;
      this->XTSigmaLMem   = this->SigmaLMem   + this->LenSigma;
      this->RhoLMem       = this->XTSigmaLMem + this->LenXTSigma;
      this->XTRhoLMem     = this->RhoLMem     + this->LenRho;
      this->ULMem         = this->XTRhoLMem   + this->LenXTRho;
      this->ResLMem       = this->ULMem       + this->LenU;
      this->TVecLMem      = this->ResLMem     + this->LenRes;
      this->ASuperMem     = this->TVecLMem    + this->LenTVec;
      this->SSuperMem     = this->ASuperMem   + this->LenSuper;
      this->SCPYMem       = this->SSuperMem   + this->LenSuper;
      this->BiOrthMem     = this->SCPYMem     + this->LenSuper;
    }
    if(!this->isHermitian_ && this->symmetrizedTrial_){
      this->NHrProdMem    = this->BiOrthMem   + this->LenSuper;
    }

    this->RWORK = this->REAL_SCR;
  }
  inline void cleanupScr(){
    delete [] this->SCR     ;
    delete [] this->REAL_SCR;
  }
  inline void checkValid(ostream &output=cout){
    if(this->A_ != NULL){
      if(this->A_->rows() != this->A_->cols())
        CErr("Quasi-Newton only supported for square problems!");
    }
    if(this->nGuess_ > this->maxSubSpace_)
      CErr("Number of initial guess vectors exceeds maximum dimension of iterative subspace");
  } // checkValid
    inline void allocSolution(){
      // Allocate space for solution
      this->solutionValues_ = new VectorXd(this->nSek_,1);
      this->solutionVector_ = new TMat(this->n_,this->nSek_);
      this->cleanupMem_ = true;
    };
    /** Quasi-Newton Procedural Functions **/
    void runMicro(ostream &output=cout); 
    void linearTrans(const int, const int);
    void fullProjection(const int);
    void buildSuperMat(const int);
    void redDiag(int,ostream &output=cout);
    void diagMem(int);
    void stdHerDiag(int, ostream &output=cout);
    void symmHerDiag(int, ostream &output=cout);
    void symmNonHerDiag(int, ostream &output=cout);
    void reconstructSolution(const int);
    void genRes(const int);
    std::vector<bool> checkConv(const int, int &,ostream &output=cout);
    void formNewGuess(std::vector<bool> &,int&,int,int&,int&);
    void formResidualGuess(double,const TMap &, TMap &, const TMap &, TMap &);
    void genStdHerResGuess(double, const TMap &, TMap &);
    void genSymmResGuess(double,const TMap &, TMap &, const TMap &, TMap &);
    void setupRestart();
    void Orth(TMap &);
    void Orth(TMat &);
    void metBiOrth(TMap &, const TMat &);
    void eigSrt(TMap &, TVecMap &);
    void initLAPACKScrLen();

  public:
    /** Destructor
     *  
     *
     *  Only useful if memory for the solution was allocated as opposed
     *  to being passed in
     */ 
    ~QuasiNewton(){
      if(this->cleanupMem_){
        delete this->solutionVector_;
        delete this->solutionValues_;
      }
    }
    /** Default Constructor
     *  Default Constructor
     *
     *  Loads the default parameter values (see loadDefaults)
     */ 
    QuasiNewton(){ this->loadDefaults();};
    /** Constructor for Quasi-Newton based on a SDResponse Object
     *  Constructor for Quasi-Newton based on a SDResponse Object
     *
     *  All of the parameters are set based on the SDResponse object
     *  that is passed in. Currently, QuasiNewton knows about the
     *  following from SDResponse:
     *   -  Configuration Interaction Singles            (CIS)
     *   -  Random Phase Approximation                   (RPA)
     *   -  Particle-Particle Random Phase Approximation (pp-RPA)
     *   
     *  Currently, the following developments are under-way and
     *  QuasiNewton will know how to handle them soon enough
     *   -  Standard Frequency Dependent Linear Response (SFDLR)
     *   -  Damped Frequency Dependent Linear Response   (DFDLR)
     */ 
    QuasiNewton(SDResponse<T> * SDR) : QuasiNewton(){
      this->sdr_              = SDR; 
      this->N_                = SDR->nSingleDim();
      this->nSek_             = SDR->nSek();
      this->nGuess_           = SDR->nGuess();
      this->solutionValues_   = SDR->omega();
      this->solutionVector_   = SDR->transDen();
      this->diag_             = SDR->rmDiag();
      this->symmetrizedTrial_ = (SDR->iMeth() == SDResponse<T>::RPA);
      this->doGEP_            = (SDR->iMeth() == SDResponse<T>::RPA);
      this->genGuess_         = (this->nGuess_ == 0);
      this->maxSubSpace_      = this->stdSubSpace();
      this->isHermitian_      = !(SDR->iMeth() == SDResponse<T>::RPA);
      this->initScrLen();

      if(this->genGuess_) this->nGuess_ = this->stdNGuess();
      this->checkValid(SDR->fileio()->out);
      this->allocGuess();
//    if(!this->genGuess_) *this->guessR_ = *SDR->davGuess();
      if(!this->genGuess_) *this->guessR_ = TMat::Identity(this->N_,this->nGuess_);;
      this->allocScr();
    };

    QuasiNewton(bool isH, bool sT, int n, TMat *A,RealMatrix* diag,TMat *Vc, VectorXd *Eig) : QuasiNewton() {
      this->N_ = A->rows();
      this->A_ = A;
      this->nSek_ = n;
      this->solutionValues_ = Eig;
      this->solutionVector_ = Vc;
      this->diag_           = diag;
      this->symmetrizedTrial_ = sT;
      this->doGEP_            = false;
      this->genGuess_         = true;
      this->maxSubSpace_      = this->stdSubSpace();
      this->isHermitian_      = isH;
      this->initScrLen();

      if(this->genGuess_) this->nGuess_ = this->stdNGuess();
      this->checkValid(cout);
      this->allocGuess();
      this->allocScr();
    }
    /** Public inline functions **/
    inline TVec* eigenValues(){return this->solutionValues_;};
    inline TMat* eigenVector(){return this->solutionVector_;};
    inline void printInfo(ostream &output=cout) {
      output << bannerTop << endl;
      output << "Davidson Diagonalization Settings:" << endl << endl;
 
      output << std::setw(50) << std::left << "  Dimension of the Matrix:" << this->N_ << endl;
      output << std::setw(50) << std::left << "  Number of Desired Roots:" << this->nSek_ << endl;
      output << std::setw(50) << std::left << "  Number of Initial Guess Vectors:" << this->nGuess_;
      if(this->nGuess_ == 2*this->nSek_) output << "    (default 2*NSek)";
      output << endl;
      output << std::setw(50) << std::left << "  Maximum Dimension of Iterative Subspace:"
           << this->maxSubSpace_ << endl;
      output << std::setw(50) << std::left << "  Maximum Number of Micro Iterations:"
           << this->maxIter_ << endl;
      output << std::setw(50) << std::left << "  Maximum Number of Macro Iterations:"
           << this->MaxIter_ << endl;
      output << std::setw(50) << std::left << "  Using an Hermitian algorithm?:";
      if(this->isHermitian_) output << "Yes";
      else output << "No";
      output << endl;
      output << std::setw(50) << std::left << "  Full Matrix Passed to for AX:";
      if(this->A_ != nullptr) output << "Yes";
      else output << "No";
      output << endl;
 
      output << endl << bannerEnd << endl;
    }
    void run(ostream &output=cout);
  
    inline int nIter(){ return this->nIter_;};
    
  }; // class QuasiNewton
  /** Run the Quasi-Newton Calculation (templated) **/
  template <typename T>
  void QuasiNewton<T>::run(ostream &output){
    time_t currentTime;
    std::chrono::high_resolution_clock::time_point start,finish;
    std::chrono::duration<double> elapsed;
    output << bannerTop << endl << endl;
    time(&currentTime);
    output << "Quasi-Newton Calculation Started: " << ctime(&currentTime) << endl;
    this->printInfo(output);
    start = std::chrono::high_resolution_clock::now();
    for(auto iter = 0; iter < this->MaxIter_; iter++){
      this->nMicroIter_ = 0;
      this->runMicro(output);
      this->nIter_ += this->nMicroIter_;
      if(this->isConverged_) break;
    }
    this->cleanupScr();
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    time(&currentTime);
    if(!this->isConverged_) CErr("Quasi-Newton Failed to Converge Within Given Criteria!",output);
    output << "Quasi-Newton Calculation Finished: " << ctime(&currentTime);
    output << "Time Elapsed: " << elapsed.count() << " sec" << endl;
    output << bannerEnd << endl << endl;
  }
  /** Linear Transformation **/
  /*  Compute the linear transformation of the matrix (σ) [and possibly the 
   *  metric (ρ)] onto the basis vectors (b)
   *
   *  For Hermitian matricies in general (viz. Davidson J. Comput. Phys. 17 (1975))
   *
   *  σ = E| b >
   *
   *  For RPA (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 74,75)):
   *
   *  σ_g = E| b_g >   σ_u = E| b_u >
   *  ρ_g = E| b_u >   ρ_u = E| b_g >
   *
   *  ** ρ not referenced for CIS although it is passed **
   *
   */
  template<typename T>
  void QuasiNewton<T>::linearTrans(const int NOld,const int NNew){
    TMap NewSR  (this->SigmaRMem+NOld*this->N_,this->N_,NNew);
    TMap NewVecR(this->TVecRMem+ NOld*this->N_,this->N_,NNew);
    
    // Initialize these Eigen Maps so they remain in scope
    TMap NewRhoR(this->SCR,0,0);
    TMap NewRhoL(this->SCR,0,0);
    TMap NewSL  (this->SCR,0,0);
    TMap NewVecL(this->SCR,0,0);
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      new (&NewRhoR) TMap(this->RhoRMem + NOld*this->N_,this->N_,NNew);
      new (&NewRhoL) TMap(this->RhoLMem + NOld*this->N_,this->N_,NNew);
      new (&NewSL  ) TMap(this->SigmaLMem+NOld*this->N_,this->N_,NNew);
      new (&NewVecL) TMap(this->TVecLMem+ NOld*this->N_,this->N_,NNew);
    }
    if(this->sdr_ != NULL){
      if(this->sdr_->iMeth() == SDResponse<T>::CIS || 
         this->sdr_->iMeth() == SDResponse<T>::RPA ||
         this->sdr_->iMeth() == SDResponse<T>::STAB){
        // Linear transformation onto right / gerade
        this->sdr_->formRM3(NewVecR,NewSR,NewRhoL); 
        // Linear trasnformation onto left / ungerade
        if(this->sdr_->iMeth() == SDResponse<T>::RPA){
          this->sdr_->formRM3(NewVecL,NewSL,NewRhoR);
        }
      } else if(this->sdr_->iMeth() == SDResponse<T>::PPRPA  || 
                this->sdr_->iMeth() == SDResponse<T>::PPATDA ||
                this->sdr_->iMeth() == SDResponse<T>::PPCTDA) {
        this->sdr_->formRM4(NewVecR,NewSR,NewRhoL); 
        if(this->sdr_->iMeth() == SDResponse<T>::PPRPA)   
          // Linear trasnformation onto left / ungerade
          this->sdr_->formRM4(NewVecL,NewSL,NewRhoR);
      }
    } else {
      NewSR = (*this->A_) * NewVecR;
      if(!this->isHermitian_ || this->symmetrizedTrial_){
        NewSL = (*this->A_) * NewVecL;
        NewRhoR = NewVecL;
        NewRhoL = NewVecR;
        NewRhoR.block(this->N_/2,0,this->N_/2,NNew) *= -1.0;
        NewRhoL.block(this->N_/2,0,this->N_/2,NNew) *= -1.0;
        //cout << "SR" << endl << std::setprecision(6) << NewSR << endl;
        //cout << "RhoR" << endl << std::setprecision(6) << NewRhoR << endl;
        //cout << "SL" << endl << std::setprecision(6) << NewSL << endl;
        //cout << "RhoL" << endl << std::setprecision(6) << NewRhoL << endl;
      }
    }
  //CErr();
  } // linearTrans
  /** Full projection onto reduced space  **/
  /*  Full projection of the matrix (and the metric) onto the reduced subspace
   *
   *  For Hermitian matricies in general (viz. Davidson J. Comput. Phys. 17 (1975))
   *
   *  E(R) = < b | σ >
   *
   *  For RPA (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 89)):
   *
   *  E(R)_gg = < b_g | σ_g >
   *  E(R)_uu = < b_u | σ_u >
   *  S(R)_gu = < b_g | ρ_g >
   *  S(R)_uu = < b_u | ρ_u >
   *
   */ 
  template<typename T>
  void QuasiNewton<T>::fullProjection(const int NTrial){
    TMap SigmaR   (this->SigmaRMem,  this->N_,NTrial);
    TMap XTSigmaR (this->XTSigmaRMem,NTrial,  NTrial);
    TMap TrialVecR(this->TVecRMem,   this->N_,NTrial);

    // Initialize these Eigen Maps so they remain in scope
    TMap RhoR     (this->RhoRMem,    0,0);
    TMap XTRhoR   (this->XTRhoRMem,  0,0);
    TMap SigmaL   (this->SigmaLMem,  0,0);
    TMap XTSigmaL (this->XTSigmaLMem,0,0);
    TMap RhoL     (this->RhoLMem,    0,0);
    TMap XTRhoL   (this->XTRhoLMem,  0,0);
    TMap TrialVecL(this->TVecLMem,   0,0);
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      new (&RhoR)      TMap(this->RhoRMem,    this->N_,NTrial);
      new (&XTRhoR)    TMap(this->XTRhoRMem,  NTrial,  NTrial);
      new (&SigmaL)    TMap(this->SigmaLMem,  this->N_,NTrial);
      new (&XTSigmaL)  TMap(this->XTSigmaLMem,NTrial,  NTrial);
      new (&RhoL)      TMap(this->RhoLMem,    this->N_,NTrial);
      new (&XTRhoL)    TMap(this->XTRhoLMem,  NTrial,  NTrial);
      new (&TrialVecL) TMap(this->TVecLMem,   this->N_,NTrial);
    }
    XTSigmaR = TrialVecR.adjoint()*SigmaR; // E(R) or E(R)_gg
    if(!isHermitian_ || symmetrizedTrial_){
      XTRhoR   = TrialVecR.adjoint()*RhoR;   // S(R)_gu
      XTSigmaL = TrialVecL.adjoint()*SigmaL; // E(R)_uu
      XTRhoL   = TrialVecL.adjoint()*RhoL;   // S(R)_ug
    }
  } // fullProjection
  /** Build Supermatricies in reduced dimension **/
  /*
   * Set up the reduced dimensional supermatricies
   * viz. Kauczor et al. JCTC p. 1610  (Eq 88)
   *
   * E(R) = [ E(R)_gg   0  ]    S(R) = [   0  S(R)_gu ]
   *        [   0  E(R)_uu ]           [ S(R)_ug   0  ]
   *
   */  
  template<typename T>
  void QuasiNewton<T>::buildSuperMat(const int NTrial){
    // Note this routine only makes sense when two
    // sets of guess vectors are employed
    TMap XTSigmaR(this->XTSigmaRMem,NTrial,  NTrial);
    TMap XTRhoR  (this->XTRhoRMem,  NTrial,  NTrial);
    TMap XTSigmaL(this->XTSigmaLMem,NTrial,  NTrial);
    TMap XTRhoL  (this->XTRhoLMem,  NTrial,  NTrial);
    TMap ASuper  (this->ASuperMem, 2*NTrial,2*NTrial);
    TMap SSuper  (this->SSuperMem, 2*NTrial,2*NTrial);

    ASuper.setZero();
    SSuper.setZero();
    ASuper.block(0,     0,     NTrial,NTrial) = XTSigmaR;
    ASuper.block(NTrial,NTrial,NTrial,NTrial) = XTSigmaL;
    SSuper.block(0,     NTrial,NTrial,NTrial) = XTRhoR;
    SSuper.block(NTrial,0,     NTrial,NTrial) = XTRhoL;
  } // buildSuperMat

  /** Construct the residual vector **/
  /*
   *  For Hermitian matricies in general (viz. Davidson J. Comput. Phys. 17 (1975))
   *
   *  R = A| X > - | X > * ω = | σ_i > * X(R)_i - | X > ω
   *
   *  For RPA (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 92,93)):
   *
   *  R_g = | {σ_g}_i > * {X(R)_g}_i - | {ρ_g}_i > * {X(R)_u}_i * ω
   *  R_u = | {σ_u}_i > * {X(R)_g}_i - | {ρ_u}_i > * {X(R)_g}_i * ω
   */ 
  template<typename T>
  void QuasiNewton<T>::genRes(const int NTrial){
    TMap SigmaR   (this->SigmaRMem,  this->N_,NTrial);
    TMap XTSigmaR (this->XTSigmaRMem,NTrial,  NTrial);
    TMap ResR     (this->ResRMem,    this->N_,NTrial);
    TMap UR       (this->URMem,      this->N_,NTrial);
    TVecMap ER      (this->ERMem,      NTrial         );

    // Initialize these Eigen Maps so they remain in scope
    TMap RhoR     (this->RhoRMem,    0,0);
    TMap SigmaL   (this->SigmaLMem,  0,0);
    TMap XTSigmaL (this->XTSigmaLMem,0,0);
    TMap RhoL     (this->RhoLMem,    0,0);
    TMap ResL     (this->ResLMem,    0,0);
    TMap UL       (this->ULMem,      0,0);
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      new (&RhoR)      TMap(this->RhoRMem,    this->N_,NTrial);
      new (&SigmaL)    TMap(this->SigmaLMem,  this->N_,NTrial);
      new (&XTSigmaL)  TMap(this->XTSigmaLMem,NTrial,  NTrial);
      new (&RhoL)      TMap(this->RhoLMem,    this->N_,NTrial);
      new (&ResL)      TMap(this->ResLMem,    this->N_,NTrial);
      new (&UL)        TMap(this->ULMem,      this->N_,NTrial);
    }
    if(this->isHermitian_ && !this->symmetrizedTrial_) 
      ResR = SigmaR*XTSigmaR - UR*ER.asDiagonal();
    if(symmetrizedTrial_) {
      ResR = SigmaR*XTSigmaR - RhoR*XTSigmaL*ER.asDiagonal();
      ResL = SigmaL*XTSigmaL - RhoL*XTSigmaR*ER.asDiagonal();
    }
  } //genRes
  /** Check for convergence **/ 
  template<typename T>
  std::vector<bool> QuasiNewton<T>::checkConv(const int NTrial, int & NNotConv,ostream &output) {
    TMap ResR     (this->ResRMem,    this->N_,NTrial);
    TVecMap ER      (this->ERMem,      NTrial         );

    // Initialize these Eigen Maps so they remain in scope
    TMap ResL     (this->ResLMem,    0,0);
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      new (&ResL)      TMap(this->ResLMem,    this->N_,NTrial);
    }
    // Vector to store convergence info
    std::vector<bool> resConv;
    NNotConv = 0;

    // Loop over NSek residual vectors. Decide from which residuals
    // will be made perturbed guess vectors
    for(auto k = 0; k < this->nSek_; k++) {
      double NORM = ResR.col(k).norm();
      if(!this->isHermitian_ || this->symmetrizedTrial_) 
        NORM = std::max(NORM,ResL.col(k).norm());
      if(NORM < this->resTol_) resConv.push_back(true);
      else {
        resConv.push_back(false); NNotConv++;
      }
    }
    output << std::fixed << std::setprecision(12);
    output << "  Checking Quasi-Newton Convergence:" << endl;
    output << "    " << std::setw(8)  << " " << std::setw(32) << std::left << "    Roots at Current Iteration:";
    output << std::setw(32) << std::left << "    (Max) Norm of Residual(s):" << endl;
    for(auto k = 0 ; k < this->nSek_; k++){
      double NORM = ResR.col(k).norm();
      if(!this->isHermitian_ || this->symmetrizedTrial_) 
        NORM = std::max(NORM,ResL.col(k).norm());

      output << "    " << std::setw(12) << "State " + std::to_string(k+1) + ":";
      output << std::setw(32) << std::left << std::fixed << (*this->solutionValues_)(k,0);
      output << std::setw(32) << std::left << std::scientific << NORM;
      if(resConv[k]) output << "     Root has converged" << endl;
      else output << "     Root has not converged" << endl;
    }
    return resConv;
  } // checkConv
} // namespace ChronusQ
#endif
