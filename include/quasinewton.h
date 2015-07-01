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
#ifndef INCLUDED_QUASINEWTON
#define INCLUDED_QUATINEWTON
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
    // Useful typedefs for Eigen templates
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
    typedef Eigen::Matrix<T,Dynamic,1> TVec;
    typedef Eigen::Map<TMat> TCMMap;
    typedef Eigen::Map<TVec> TVecMap;

    bool isHermetian_;     // Hermetian Scheme
    bool doDiag_;          // Quasi-Newton Diagonalization (Davidson)
    bool doGEP_;           // Generalized Eigenproblem
    bool doLin_;           // Quasi-Newton Linear Equation Solve
    bool doResGuess_;      // Generate new vectors from residual
    bool symmetrizedTrial_;// Symmetrized Trial Vectors (Kauczor et al. JCTC 7 (2010))
    bool debug_;           // Enables various debug options
    bool cleanupMem_;      // True if memory cleanup is required (if memory was allocated locally)
    bool isConverged_;     // True if the iterative calculation converged
    bool genGuess_;        // Generate Standard (identity) Guess

    int N_;                // Dimension of the problem
    int maxSubSpace_;      // Maximum dimension of the iterative subspace
    int maxIter_;          // Maximum number of miro iterations
    int MaxIter_;          // Maximum number of macro iterations
    int nSek_;             // Number of desired roots in the case of diagonalization (Davidson)
    int nGuess_;           // Numer of given (or to be generated) guess vectors
//  int method_;           // Post-SCF method being performed

    double resTol_;        // Residual norm tolerence for convergence

    TMat * A_;                // Pointer to full matrix to be diagonalized
    TMat * solutionVector_;   // Solution vectors at current iteration (could be local or non-local)
    TVec * solutionValues_;   // Solution values (eigenvalues) at current iteration (could be local or non-local)
    std::unique_ptr<TMat> guess_;   // Guess vectors (always local copy, even if generated elsewhere)

    SDResponse * sdr_;     // Pointer to SDResponse object

    void runMicro(ostream &output=cout);   // Run a set of micro iterations
    void checkValid(ostream &output=cout); // Check if varibles make sense (idiot check)
    void loadDefaults();                   // Initialize the control parameters to defualt values
    #include <quasinewtonscratch.h>
    void resizeMaps(const int,const int,const int);
    void linearTrans();
    void fullProjection();
    void buildSuperMat(const int);
    void redDiag(int,ostream &output=cout);
    void reconstructSolution();
    void genRes();
    void genResGuess();
    std::vector<bool> checkConv(int &,ostream &output=cout);
    void formNewGuess(std::vector<bool> &,int,int,int,int);

    // Allocate space for local copy of the guess vectors
    inline void allocGuess(){ 
      this->guess_ = std::unique_ptr<TMat>(new TMat(this->n_,this->nGuess_));
      if(this->genGuess_) this->identGuess();
    };

    // Allocate space for solution
    inline void allocSolution(){
      this->solutionVector_ = new TMat(this->nSek_,1);
      this->solutionValues_ = new TMat(this->n_,this->nSek_);
      this->cleanupMem_ = true;
    };

    // Generate the identity (standard) guess
    inline void identGuess(){
      (*this->guess_) = TMat::Identity(this->n_,this->nGuess_);
    };

    // Standard value for the maximum dimension of the
    // iterative subspace min(6*NSek,N/2)
    inline int stdSubSpace(){
      return std::min(6*this->nSek_,this->N_/2);
    };

    // Standard value for the number of inital guess
    // vectors that are to be generated 2*NSek
    inline int stdNGuess(){
      return 2*this->nSek_;
    }

  public:
    /**
     *  Destructor
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

    /**
     *  Default Constructor
     *
     *  Loads the default parameter values (see loadDefaults)
     */ 
    QuasiNewton(){ this->loadDefaults();};

    /**
     *  Constructor for debug comparison
     *  
     *  Full matrix is passed in so one can test any number of things:
     *    - Matrix-Vector product (sigma)
     *    - Quality of guess
     *    - Iteration comparison
     *
     *  This has no place in production code and should only be used for
     *  debugging. Storing the full matrix in memory completely defeats
     *  the purpose of iterative Quasi-Newton schemes.
     *
     *  Debug print is turned on by default
     */
    QuasiNewton(TMat* A, int nSek, ostream &output=cout) : QuasiNewton() {
      this->A_           = A;
      this->N_           = A.rows();
      this->nSek_        = nSek;
      this->nGuess_      = this->stdNGuess();
      this->maxSubSpace_ = this->stdSubSpace();
      this->debug_       = true;

      this->checkValid();
      this->allocGuess();
      this->allocSolution();
      this->allocScr();
    };

    /**
     *  Constructor for Quasi-Newton based on a SDResponse Object
     *
     *  All of the parameters are set based on the SDResponse object
     *  that is passed in. Currently, QuasiNewton knows about the
     *  following from SDResponse:
     *   -  Configuration Interaction Singles (CIS)
     *   -  Random Phase Approximation        (RPA)
     *   
     *  Currently, the following developments are under-way and
     *  QuasiNewton will know how to handle them soon enough
     *   -  Standard Frequency Dependent Linear Response (SFDLR)
     *   -  Damped Frequency Dependent Linear Response   (DFDLR)
     */ 
    QuasiNewton(SDResponse * SDR) : QuasiNewton(){
      this->sdr_              = SDR; 
      this->N_                = SDR->nSingleDim();
      this->nSek_             = SDR->nSek();
      this->nGuess_           = SDR->nGuess();
      this->solutionValues_   = SDR->omega();
      this->solutionVector_   = SDR->transDen();
      this->symmetrizedTrial_ = (SDR->iMeth() == SDResponse::RPA);
      this->doGEP_            = (SDR->iMeth() == SDResponse::RPA);
      this->genGuess_         = (this->nGuess_ == 0);

      if(this->genGuess_) this->nGuess_ = this->stdNGuess();
      this->checkValid();
      this->allocGuess();
      if(!this->genGuess_) *this->guess_ = *SDR->davGuess();
      this->allocScr();
    };

   inline TVec* eigenValues(){return this->solutionValues_;};
   inline TMat* eigenVector(){return this->solutionVector_;};
   inline void run(ostream &output=cout){
     time_t currentTime;
     std::chrono::high_resolution_clock::time_point start,finish;
     std::chrono::duration<double> elapsed;
     output << bannerTop << endl << endl;
     time(&currentTime);
     output << "Quasi-Newton Calculation Started: " << ctime(&currentTime) << endl;
     this->printInfo(output);
     start = std::chrono::high_resolution_clock::now();
     this->runMicro(output);
     finish = std::chrono::high_resolution_clock::now();
     elapsed = finish - start;
     time(&currentTime);
     if(!this->isConverged_) CErr("Quasi-Newton Failed to Converge Within Given Criteria!",output);
     output << "Quasi-Newton Calculation Finished: " << ctime(&currentTime);
     output << "Time Elapsed: " << elapsed.count() << " sec" << endl;
     output << bannerEnd << endl << endl;
   }
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
      output << std::setw(50) << std::left << "  Using an Hermetian algorithm?:";
      if(this->hermetian_) output << "Yes";
      else output << "No";
      output << endl;
      output << std::setw(50) << std::left << "  Using LAPACK to diagonalize subspace?:";
      if(this->useLAPACK_) output << "Yes";
      else output << "No";
      output << endl;
      output << std::setw(50) << std::left << "  Full Matrix Passed to for AX:";
      if(this->mat_ != nullptr) output << "Yes";
      else output << "No";
      output << endl;
 
      output << endl << bannerEnd << endl;
    }
  
  }; // class QuasiNewton
  template<typename T>
  void QuasiNewton<T>::loadDefaults(){
    this->initScrLen();
    this->initScrPtr();
    this->isHermetian_      = true;  // Default to Hermetian Scheme
    this->doDiag_           = true;  // Defualt to diagonalization
    this->doGEP_            = false; // Default standard eigenproblem if doing a diagonalization
    this->doLin_            = false; // Mildly redundent, as it will currently always be the opposite of doDiag
    this->doResGuess_       = true;  // Default to residual based guess
    this->symmetrizedTrial_ = false; // Default to unmodified trial vectors
    this->debug_            = false; // Default to terse output
    this->cleanupMem_       = false; // Assume that the space for results was allocated elsewhere
    this->isConverged_      = false; // Obviously the calculation doesn't begin converged...
    this->genGuess_         = true;  // Defualt generate identity guess

    // Zero out all of the size dependent quantities
    this->N_                = 0;
    this->maxSubSpace_      = 0;
    this->nSek_             = 0;
    this->nGuess_           = 0; 

    // Defaults for # of interations
    this->maxIter_          = 128;
    this->MaxIter_          = 20;

    // Initialize all pointers to some varient of NULL
    this->A_                = NULL;
    this->solutionVector_   = NULL;
    this->solutionValues_   = NULL;
    this->guess_            = nullptr;
    this->sdr_              = NULL;

    this->resTol_           = 5.0e-6;
  } // loadDefaults

  template<typename T>
  void QuasiNewton<T>::checkValid(ostream &output){
    if(this->A_ != NULL){
      if(this->A_->rows() != this->A_->cols())
        CErr("Quasi-Newton only supported for square problems!");
    }
    if(this->nGuess_ >= this->maxSubSpace_)
      CErr("Number of initial guess vectors exceeds maximum dimension of iterative subspace");
  } // checkValid
 
  // Resize the Eigen Maps to fit new vectors
  template<typename T>
  void QuasiNewton<T>::resizeMaps(const int NTrial, const int NOld, const int NNew){
    new (&this->SigmaR)   TCMMap(this->SigmaRMem,  this->N_,NTrial);
    new (&this->XTSigmaR) TCMMap(this->XTSigmaRMem,NTrial,  NTrial);
    new (&this->UR)       TCMMap(this->URMem,      this->N_,NTrial);
    new (&this->ResR)     TCMMap(this->ResRMem,    this->N_,NTrial);
    new (&this->NewSR)    TCMMap(this->SigmaRMem+NOld*this->N_,this->N_,NNew);
    new (&this->NewVecR)  TCMMap(this->TVecRMem+ NOld*this->N_,this->N_,NNew);
    if(!this->isHermetian_ || this->symmetrizedTrial_){
      new (&this->RhoR)     TCMMap(this->RhoRMem,    this->N_,NTrial);
      new (&this->XTRhoR)   TCMMap(this->XTRhoRMem,  NTrial,  NTrial);
      new (&this->SigmaL)   TCMMap(this->SigmaLMem,  this->N_,NTrial);
      new (&this->XTSigmaL) TCMMap(this->XTSigmaLMem,NTrial,  NTrial);
      new (&this->RhoL)     TCMMap(this->RhoLMem,    this->N_,NTrial);
      new (&this->XTRhoL)   TCMMap(this->XTRhoLMem,  NTrial,  NTrial);
      new (&this->UL)       TCMMap(this->ULMem,      this->N_,NTrial);
      new (&this->ResL)     TCMMap(this->ResLMem,    this->N_,NTrial);
      new (&this->ASuper)   TCMMap(this->ASuperMem, 2*NTrial,2*NTrial);
      new (&this->SSuper)   TCMMap(this->SSuperMem, 2*NTrial,2*NTrial);

      new (&this->NewRhoR)  TCMMap(this->RhoRMem + NOld*this->N_,this->N_,NNew);
      new (&this->NewRhoL)  TCMMap(this->RhoLMem + NOld*this->N_,this->N_,NNew);
      new (&this->NewSL)    TCMMap(this->SigmaLMem+NOld*this->N_,this->N_,NNew);
      new (&this->NewVecL)  TCMMap(this->TVecLMem+ NOld*this->N_,this->N_,NNew);
    }
  } // resizeMaps

  /*
   *  Compute the linear transformation of the matrix (σ) [and possibly the 
   *  metric (ρ)] onto the basis vectors (b)
   *
   *  For Hermetian matricies in general (viz. Davidson J. Comput. Phys. 17 (1975))
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
  void QuasiNewton<T>::linearTrans(){
    if(this->sdr_ != NULL){
      if(this->sdr_->iMeth() == SDResponse::CIS || this->sdr_->iMeth() == SDResponse::RPA){
        // Linear transformation onto right / gerade
        this->sdr_->formRM3(this->NewVecR,this->NewSR,this->NewRhoL); 
        if(this->sdr_->iMeth() == SDResponse::RPA)   
          // Linear trasnformation onto left / ungerade
          this->sdr_->formRM3(this->NewVecL,this->NewSL,this->NewRhoR);
      }
    } else this->NewSR = (*this->A_) * this->NewVecR;
  } // linearTrans


  /*
   *  Full projection of the matrix (and the metric) onto the reduced subspace
   *
   *  For Hermetian matricies in general (viz. Davidson J. Comput. Phys. 17 (1975))
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
  void QuasiNewton<T>::fullProjection(){
    this->XTSigmaR = this->TrialVecR.transpose()*this->SigmaR; // E(R) or E(R)_gg
    if(!this->isHermetian_ || this->symmetrizedTrial_){
      this->XTRhoR   = this->TrialVecR.transpose()*this->RhoR;   // S(R)_gu
      this->XTSigmaL = this->TrialVecL.transpose()*this->SigmaL; // E(R)_uu
      this->XTRhoL   = this->TrialVecL.transpose()*this->RhoL;   // S(R)_ug
    }
  } // fullProjection

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
    ASuper.setZero();
    SSuper.setZero();
    ASuper.block(0,     0,     NTrial,NTrial) = XTSigmaR;
    ASuper.block(NTrial,NTrial,NTrial,NTrial) = XTSigmaL;
    SSuper.block(0,     NTrial,NTrial,NTrial) = XTRhoR;
    SSuper.block(NTrial,0,     NTrial,NTrial) = XTRhoL;
  } // buildSuperMat

  /*
   *  Reconstruct the approximate eigenvectors
   *
   *  For Hermetian matricies in general (viz. Davidson J. Comput. Phys. 17 (1975))
   *
   *  | X > = | b_i > X(R)_i
   *
   *  For RPA (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 90,91)):
   *
   *  | X_g > = | {b_g}_i > {X(R)_g}_i
   *  | X_u > = | {b_u}_i > {X(R)_u}_i
   */ 
  template<typename T>
  void QuasiNewton<T>::reconstructSolution(){
    this->UR = this->TrialVecR * this->XTSigmaR;
    if(this->symmetrizedTrial_) this->UL = this->TrialVecL * this->XTSigmaL;
    // Stash away current approximation of eigenvalues and eigenvectors (NSek)
    (*this->solutionValues_) = this->ER.block(0,0,this->nSek_,1);
    (*this->solutionVector_) = this->UR.block(0,0,this->N_,this->nSek_); 
    // | X > = | X_g > + | X_u > (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 80))
    if(this->symmetrizedTrial_)(*this->solutionVector_) += this->UL.block(0,0,this->N_,this->nSek_); 
  } // reconstructSolution

  /*
   * Construct the residual vector
   *
   *  For Hermetian matricies in general (viz. Davidson J. Comput. Phys. 17 (1975))
   *
   *  R = A| X > - | X > * ω = | σ_i > * X(R)_i - | X > ω
   *
   *  For RPA (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 92,93)):
   *
   *  R_g = | {σ_g}_i > * {X(R)_g}_i - | {ρ_g}_i > * {X(R)_u}_i * ω
   *  R_u = | {σ_u}_i > * {X(R)_g}_i - | {ρ_u}_i > * {X(R)_g}_i * ω
   */ 
  template<typename T>
  void QuasiNewton<T>::genRes(){
    if(this->isHermetian_ && !this->symmetrizedTrial_) 
      this->ResR = this->SigmaR*this->XTSigmaR - this->UR*this->ER.asDiagonal();
    if(this->symmetrizedTrial_) {
      this->ResR = this->SigmaR*this->XTSigmaR - 
                   this->RhoR*this->XTSigmaL*this->ER.asDiagonal();
      this->ResL = this->SigmaL*this->XTSigmaL - 
                   this->RhoL*this->XTSigmaR*this->ER.asDiagonal();
    }
  } //genRes

  /**
   * Check for convergence
   */ 
  template<typename T>
  std::vector<bool> QuasiNewton<T>::checkConv(int & NNotConv,ostream &output) {
    // Vector to store convergence info
    std::vector<bool> resConv;
    NNotConv = 0;

    // Loop over NSek residual vectors. Decide from which residuals
    // will be made perturbed guess vectors
    for(auto k = 0; k < this->nSek_; k++) {
      double NORM = this->ResR.col(k).norm();
      if(!this->isHermetian_ || this->symmetrizedTrial_) 
        NORM = std::max(NORM,this->ResL.col(k).norm());
      if(NORM < this->resTol_) resConv.push_back(true);
      else {
        resConv.push_back(false); NNotConv++;
      }
    }
    output << "  Checking Quasi-Newton Convergence:" << endl;
    output << "    " << std::setw(8)  << " " << std::setw(32) << std::left << "    Roots at Current Iteration:";
    output << std::setw(32) << std::left << "    (Max) Norm of Residual(s):" << endl;
    for(auto k = 0 ; k < this->nSek_; k++){
      double NORM = this->ResR.col(k).norm();
      if(!this->isHermetian_ || this->symmetrizedTrial_) 
        NORM = std::max(NORM,this->ResL.col(k).norm());

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
