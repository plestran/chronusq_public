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
#include <davidson.h>

namespace ChronusQ {
template<>
void Davidson<double>::runMicro(ostream &output ) {

  /* Explicit Memory Allocation */

  int LenScr  = 0;
  int LenSigma   = this->n_ * this->maxSubSpace_;
  int LenRho     = this->n_ * this->maxSubSpace_;
  int LenXTSigma = this->maxSubSpace_ * this->maxSubSpace_;
  int LenXTRho   = this->maxSubSpace_ * this->maxSubSpace_;
  int LenSuper   = 4 * this->maxSubSpace_ * this->maxSubSpace_;
  int LenU    = this->n_ * this->maxSubSpace_;
  int LenRes  = this->n_ * this->maxSubSpace_;
  int LenTVec = this->n_ * this->maxSubSpace_;
  int LEN_LAPACK_SCR = 0;
  int LWORK = 0;

  /* *******************************************************
   * *********Determine length of scratch space*************
   * *******************************************************
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
   * (18) Local copy of the real part of the eigenvalues (reused for Tau storage for QR)
   *
   * (19) Space for the paired / imaginary part of the eigenvalues
   *
   * (20) Length of LAPACK workspace (used in all LAPACK Calls)
   *
   * (21) Total double precision words required for LAPACK
   */

  LenScr += LenSigma;   // 1
  LenScr += LenXTSigma; // 2
  LenScr += LenU;       // 3
  LenScr += LenRes;     // 4 
  LenScr += LenTVec;    // 5


  if(!this->hermetian_ || this->symmetrized_){
    LenScr += LenRho;     // 7
    LenScr += LenXTRho;   // 8 
    LenScr += LenSigma;   // 9
    LenScr += LenRho;     // 10
    LenScr += LenXTSigma; // 11
    LenScr += LenXTRho;   // 12
    LenScr += LenU;       // 13
    LenScr += LenRes;     // 14
    LenScr += LenTVec;    // 15
    LenScr += LenSuper;   // 16
    LenScr += LenSuper;   // 17
  }

  LWORK = 6*this->n_;
  LEN_LAPACK_SCR += this->maxSubSpace_; // 18
  if(!this->hermetian_ || this->symmetrized_)
    LEN_LAPACK_SCR += this->maxSubSpace_; // 19
  LEN_LAPACK_SCR += LWORK;              // 20
  LenScr += LEN_LAPACK_SCR;             // 21

  double * SCR         = NULL; 
  double * SigmaRMem   = NULL; 
  double * SigmaLMem   = NULL; 
  double * XTSigmaRMem = NULL; 
  double * XTSigmaLMem = NULL; 
  double * RhoRMem     = NULL; 
  double * RhoLMem     = NULL; 
  double * XTRhoRMem   = NULL; 
  double * XTRhoLMem   = NULL; 
  double * ASuperMem   = NULL;
  double * SSuperMem   = NULL;
  double * URMem       = NULL; 
  double * ULMem       = NULL; 
  double * ResRMem     = NULL; 
  double * ResLMem     = NULL; 
  double * TVecRMem    = NULL;         
  double * TVecLMem    = NULL;         
  double * LAPACK_SCR  = NULL;

  // Allocate Scratch Space
  SCR = new double [LenScr]; 

  // Partition scratch space
  SigmaRMem     = SCR;
  XTSigmaRMem   = SigmaRMem   + LenSigma;
  URMem         = XTSigmaRMem + LenXTSigma;
  ResRMem       = URMem       + LenU; 
  TVecRMem      = ResRMem     + LenRes;
  LAPACK_SCR    = TVecRMem    + LenTVec;
  if(!this->hermetian_ || this->symmetrized_){
    RhoRMem       = LAPACK_SCR  + LEN_LAPACK_SCR;
    XTRhoRMem     = RhoRMem     + LenRho;
    SigmaLMem     = XTRhoRMem   + LenXTRho;
    XTSigmaLMem   = SigmaLMem   + LenSigma;
    RhoLMem       = XTSigmaLMem + LenXTSigma;
    XTRhoLMem     = RhoLMem     + LenRho;
    ULMem         = XTRhoLMem   + LenXTRho;
    ResLMem       = ULMem       + LenU;
    TVecLMem      = ResLMem     + LenRes;
    ASuperMem     = TVecLMem    + LenTVec;
    SSuperMem     = ASuperMem   + LenSuper;
  }


  RealCMMap SigmaR(   SigmaRMem,  0,0);
  RealCMMap NewSR(    SigmaRMem,  0,0);
  RealCMMap XTSigmaR( XTSigmaRMem,0,0);
  RealCMMap RhoR(     RhoRMem,    0,0);
  RealCMMap NewRhoR(  RhoRMem,    0,0);
  RealCMMap XTRhoR(   XTRhoRMem,  0,0);
  RealCMMap UR(       URMem,      0,0);
  RealCMMap ResR(     ResRMem,    0,0);
  RealCMMap TrialVecR(TVecRMem,   0,0);
  RealCMMap NewVecR(  TVecRMem,   0,0);

  RealCMMap SigmaL(   SigmaLMem,  0,0);
  RealCMMap NewSL(    SigmaLMem,  0,0);
  RealCMMap XTSigmaL( XTSigmaLMem,0,0);
  RealCMMap RhoL(     RhoLMem,    0,0);
  RealCMMap NewRhoL(  RhoLMem,    0,0);
  RealCMMap XTRhoL(   XTRhoLMem,  0,0);
  RealCMMap UL(       ULMem,      0,0);
  RealCMMap ResL(     ResLMem,    0,0);
  RealCMMap TrialVecL(TVecLMem,   0,0);
  RealCMMap NewVecL(  TVecLMem,   0,0);

  RealCMMap ASuper(ASuperMem,0,0);
  RealCMMap SSuper(SSuperMem,0,0);

  RealCMMap QR(TVecRMem,0,0); 
  RealCMMap QL(TVecLMem,0,0);
  RealCMMap RR(ResRMem ,0,0);
  RealCMMap RL(ResLMem ,0,0);


  RealVecMap ER(LAPACK_SCR,0);
  RealVecMap EI(LAPACK_SCR,0);
  RealMap    VR(LAPACK_SCR,0,0);
  RealMap    VL(LAPACK_SCR,0,0);

  char JOBVR = 'V';
  char JOBVL = 'N';
  char UPLO = 'L';
  int INFO;

  int NTrial = this->nGuess_;
  int NOld   = 0;
  int NNew   = this->nGuess_;

  new (&TrialVecR) RealCMMap(TVecRMem,this->n_,NTrial);
  if(!this->hermetian_ || this->symmetrized_){
    new (&TrialVecL) RealCMMap(TVecLMem,this->n_,NTrial);
  }

  if(this->symmetrized_){
    TrialVecR = (*this->guess_);
    TrialVecL = (*this->guess_);
    TrialVecR.block(this->n_/2,0,this->n_/2,this->nGuess_)
      = TrialVecR.block(0,0,this->n_/2,this->nGuess_);
    TrialVecL.block(this->n_/2,0,this->n_/2,this->nGuess_)
      = -TrialVecL.block(0,0,this->n_/2,this->nGuess_);
    TrialVecR *= std::sqrt(0.5);
    TrialVecL *= std::sqrt(0.5);
  } else {
    TrialVecR = (*this->guess_);
  }

  for(auto iter = 0; iter < this->maxIter_; iter++){
    std::chrono::high_resolution_clock::time_point start,finish;
    std::chrono::duration<double> elapsed;
    output << "Starting Davidson Micro Iteration " << iter + 1 << endl;
    start = std::chrono::high_resolution_clock::now();

    // Resize the Eigen Maps to fit new vectors
    new (&SigmaR)   RealCMMap(SigmaRMem,  this->n_,NTrial);
    new (&XTSigmaR) RealCMMap(XTSigmaRMem,NTrial,  NTrial);
    new (&UR)       RealCMMap(URMem,      this->n_,NTrial);
    new (&ResR)     RealCMMap(ResRMem,    this->n_,NTrial);
    new (&NewSR)    RealCMMap(SigmaRMem+NOld*this->n_,this->n_,NNew);
    new (&NewVecR)  RealCMMap(TVecRMem+ NOld*this->n_,this->n_,NNew);
    if(!this->hermetian_ || this->symmetrized_){
      new (&RhoR)     RealCMMap(RhoRMem,    this->n_,NTrial);
      new (&XTRhoR)   RealCMMap(XTRhoRMem,  NTrial,  NTrial);
      new (&SigmaL)   RealCMMap(SigmaLMem,  this->n_,NTrial);
      new (&XTSigmaL) RealCMMap(XTSigmaLMem,NTrial,  NTrial);
      new (&RhoL)     RealCMMap(RhoLMem,    this->n_,NTrial);
      new (&XTRhoL)   RealCMMap(XTRhoLMem,  NTrial,  NTrial);
      new (&UL)       RealCMMap(ULMem,      this->n_,NTrial);
      new (&ResL)     RealCMMap(ResLMem,    this->n_,NTrial);
      new (&ASuper)   RealCMMap(ASuperMem, 2*NTrial,2*NTrial);
      new (&SSuper)   RealCMMap(SSuperMem, 2*NTrial,2*NTrial);

      new (&NewRhoR)  RealCMMap(RhoRMem + NOld*this->n_,this->n_,NNew);
      new (&NewRhoL)  RealCMMap(RhoLMem + NOld*this->n_,this->n_,NNew);
      new (&NewSL)    RealCMMap(SigmaLMem+NOld*this->n_,this->n_,NNew);
      new (&NewVecL)  RealCMMap(TVecLMem+ NOld*this->n_,this->n_,NNew);
    }

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
    if(this->method_ == SDResponse::CIS || this->method_ == SDResponse::RPA){ 
      this->sdr_->formRM3(NewVecR,NewSR,NewRhoL); // Linear transforms onto right / gerade
      if(this->method_ == SDResponse::RPA)
        this->sdr_->formRM3(NewVecL,NewSL,NewRhoR); // Linear transforms onto left / ungerade
      if(this->debug_){
        prettyPrint(output,SigmaR,"Sigma (g)   ITER: "+std::to_string(iter));
        prettyPrint(output,SigmaL,"Sigma (u)   ITER: "+std::to_string(iter));
        prettyPrint(output,RhoR,"Rho (g)   ITER: "+std::to_string(iter));
        prettyPrint(output,RhoL,"Rho (u)   ITER: "+std::to_string(iter));
      }
    }
    else if(this->AX_==NULL) NewSR = (*this->mat_) * NewVecR;  
    else NewSR = this->AX_(*this->mat_,NewVecR);
   
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
    XTSigmaR = TrialVecR.transpose()*SigmaR; // E(R) or E(R)_gg
    if(this->debug_)
      prettyPrint(output,XTSigmaR,"E(R) (gg)   ITER: "+std::to_string(iter));
    if(!this->hermetian_ || this->symmetrized_){
      XTRhoR   = TrialVecR.transpose()*RhoR;   // S(R)_gu
      XTSigmaL = TrialVecL.transpose()*SigmaL; // E(R)_uu
      XTRhoL   = TrialVecL.transpose()*RhoL;   // S(R)_ug
      if(this->debug_){
        prettyPrint(output,XTRhoR,  "S(R) (gu)   ITER: "+std::to_string(iter));
        prettyPrint(output,XTSigmaL,"E(R) (uu)   ITER: "+std::to_string(iter));
        prettyPrint(output,XTRhoL  ,"S(R) (ug)   ITER: "+std::to_string(iter));
      }

      /*
       * Set up the reduced dimensional supermatricies
       * viz. Kauczor et al. JCTC p. 1610  (Eq 88)
       *
       * E(R) = [ E(R)_gg   0  ]    S(R) = [   0  S(R)_gu ]
       *        [   0  E(R)_uu ]           [ S(R)_ug   0  ]
       *
       */  

      ASuper.setZero();
      SSuper.setZero();
      ASuper.block(0,     0,     NTrial,NTrial) = XTSigmaR;
      ASuper.block(NTrial,NTrial,NTrial,NTrial) = XTSigmaL;
      SSuper.block(0,     NTrial,NTrial,NTrial) = XTRhoR;
      SSuper.block(NTrial,0,     NTrial,NTrial) = XTRhoL;
    }


    /* Diagonalize the subspace */
    int IOff = NTrial; // Space for eigenvalues
    if(!this->hermetian_ || this->symmetrized_){
      IOff += NTrial; // Space for paired eigenvalues
      new (&ER) RealVecMap(LAPACK_SCR,2*NTrial);
    } else{ new (&ER) RealVecMap(LAPACK_SCR,NTrial);}
/*
    if(!this->hermetian_) {
      new (&EI) Eigen::Map<Eigen::VectorXd>(LAPACK_SCR+IOff,NTrial);
      IOff += NTrial;
      new (&VR) RealMap(LAPACK_SCR+IOff,NTrial,NTrial);
      IOff += NTrial*NTrial;
    }
*/
    if(this->hermetian_) {
      if(!this->symmetrized_){
        // Solve E(R)| X(R) > = | X(R) > ω
        dsyev_(&JOBVR,&UPLO,&NTrial,XTSigmaR.data(),&NTrial,ER.data(),
               LAPACK_SCR+IOff,&LWORK,&INFO); 
        if(INFO!=0) CErr("DSYEV failed to converge in Davison Iterations",output);
      } else {
        int iType = 1; // Flag for the type of GEP being solved
        int TwoNTrial = 2*NTrial; // Dimension of the supermatricies
        RealCMMatrix SCPY(SSuper); // Copy of original matrix to use for re-orthogonalization
        /*
         * Solve S(R)| X(R) > = E(R)| X(R) > (1/ω)
         *
         * | X(R) > = | X(R)_g >
         *            | X(R)_u >
         *
         * The opposite (1/ω vs ω) is solved because the metric is not positive definite
         * and can therefore not be solved using DSYGV because of the involved Cholesky
         * decomposition.
         *
         * This must be revisited for Spinor orbital hessians! FIXME
         */ 
        dsygv_(&iType,&JOBVR,&UPLO,&TwoNTrial,SSuper.data(),&TwoNTrial,
               ASuper.data(),&TwoNTrial,ER.data(),LAPACK_SCR+IOff,&LWORK,
               &INFO);
        if(INFO!=0) CErr("DSYGV failed to converge in Davison Iterations",output);

        // Grab the "positive paired" roots (throw away other element of the pair)
        new (&ER)     RealVecMap(LAPACK_SCR+NTrial,NTrial);
        new (&SSuper) RealCMMap (SSuperMem+2*NTrial*NTrial,2*NTrial,NTrial);

        // Swap the ordering because we solve for (1/ω)
        for(auto i = 0 ; i < NTrial; i++) ER(i) = 1.0/ER(i);
        for(auto i = 0 ; i < NTrial/2; i++){
          SSuper.col(i).swap(SSuper.col(NTrial - i - 1));
          double tmp = ER(i);
          ER(i) = ER(NTrial - i - 1);
          ER(NTrial - i - 1) = tmp;
        }

        // Re-orthogonalize the eigenvectors with respect to the metric S(R)
        // because DSYGV orthogonalzies the vectors with respect to E(R)
        // because we solve the opposite problem.
        //
        // Gramm-Schmidt
        for(auto i = 0; i < NTrial; i++){
          double inner = SSuper.col(i).dot(SCPY*SSuper.col(i));
          int sgn = inner / std::abs(inner);
          inner = sgn*std::sqrt(sgn*inner);
          SSuper.col(i) /= inner;
          for(auto j = i+1; j < NTrial; j++){
            SSuper.col(j) -= SSuper.col(i)*(SSuper.col(i).dot(SCPY*SSuper.col(j)));
          }
        }

        XTSigmaR = SSuper.block(0,     0,NTrial,NTrial);
        XTSigmaL = SSuper.block(NTrial,0,NTrial,NTrial);
        if(this->debug_){
          prettyPrint(output,XTSigmaR,"X(R) (g)   ITER: "+std::to_string(iter));
          prettyPrint(output,XTSigmaL,"X(R) (u)   ITER: "+std::to_string(iter));
        }
      }
    } 
/*
    else {
      dgeev_(&JOBVL,&JOBVR,&NTrial,XTSigmaR.data(),&NTrial,ER.data(),
             EI.data(),VL.data(),&NTrial,VR.data(),&NTrial,
             LAPACK_SCR+IOff,&LWORK,&INFO);
      if(INFO!=0) CErr("DGEEV failed to converge in Davison Iterations",output);
      VR.transposeInPlace(); // Convert ro RowMajor...
//    eigSort(false,ER,VR,VL);
    }
*/
    
   
    
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
    if(this->hermetian_)   UR = TrialVecR * XTSigmaR;
    if(this->symmetrized_) UL = TrialVecL * XTSigmaL;
    if(this->symmetrized_ && this->debug_){
      prettyPrint(output,UR,"X (g)   ITER: "+std::to_string(iter));
      prettyPrint(output,UL,"X (u)   ITER: "+std::to_string(iter));
    }
//  else                 UR = TrialVecR * VR;

    // Stash away current approximation of eigenvalues and eigenvectors (NSek)
    (*this->eigenvalues_) = ER.block(0,0,this->nSek_,1);
    (*this->eigenvector_) = UR.block(0,0,this->n_,this->nSek_); 
    // | X > = | X_g > + | X_u > (viz. Kauczor et al. JCTC 7 (2010) p. 1610  (Eq 80))
    if(this->symmetrized_)(*this->eigenvector_) += UL.block(0,0,this->n_,this->nSek_); 
    
    // Construct the residual vector 
    // R = A*U - S*U*E = (AX)*c - (SX)*c*E
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
    if(this->hermetian_ && !this->symmetrized_) ResR = SigmaR*XTSigmaR - UR*ER.asDiagonal();
    if(this->symmetrized_) {
      ResR = SigmaR*XTSigmaR - RhoR*XTSigmaL*ER.asDiagonal();
      ResL = SigmaL*XTSigmaL - RhoL*XTSigmaR*ER.asDiagonal();
      if(this->debug_){
        prettyPrint(output,ResR,"Res (g)   ITER: " + std::to_string(iter));
        prettyPrint(output,ResL,"Res (u)   ITER: " + std::to_string(iter));
      }
   }

/*
    else {
      int NImag = 0;
      for(auto k = 0; k < NTrial; k++) if(std::abs(EI(k))<1e-10) NImag++;
      ResR.resize(this->n_,NTrial+NImag);
      for(auto k = 0; k < NTrial+NImag; k++) {
        ResR.col(k) = SigmaR*VR.col(k) - ER(k)*UR.col(k);
        if(std::abs(EI(k))<1e-10) {
          ResR.col(k+1) = SigmaR*VR.col(k) - EI(k)*UR.col(k);
          k++;
        }
      }
    }
*/
    

    // Vector to store convergence info
    std::vector<bool> resConv;
    int NNotConv = 0;

    // Loop over NSek residual vectors. Decide from which residuals
    // will be made perturbed guess vectors
    for(auto k = 0; k < this->nSek_; k++) {
      double NORM = ResR.col(k).norm();
      if(!this->hermetian_ || this->symmetrized_) 
        NORM = std::max(NORM,ResL.col(k).norm());
      if(NORM < 5e-6) resConv.push_back(true);
      else {
        resConv.push_back(false); NNotConv++;
      }
    }
    output << "  Checking Davidson Convergence:" << endl;
    output << "    " << std::setw(8)  << " " << std::setw(32) << std::left << "    Roots at Current Iteration:";
    output << std::setw(32) << std::left << "    (Max) Norm of Residual(s):" << endl;
    for(auto k = 0 ; k < this->nSek_; k++){
      double NORM = ResR.col(k).norm();
      if(!this->hermetian_ || this->symmetrized_) 
        NORM = std::max(NORM,ResL.col(k).norm());

      output << "    " << std::setw(12) << "State " + std::to_string(k+1) + ":";
      output << std::setw(32) << std::left << std::fixed << (*this->eigenvalues_)(k,0);
      output << std::setw(32) << std::left << std::scientific << NORM;
      if(resConv[k]) output << "     Root has converged" << endl;
      else output << "     Root has not converged" << endl;
    }

    // If NSek (lowest) residuals are converged, exit, else, create
    // perturbed guess vectors for next iteration
    this->converged_ = (NNotConv == 0);
    if(this->converged_) {
      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish - start;
      output << "  Davidson Micro Iteration took " << std::fixed 
             << elapsed.count() << " secs" << endl << endl;
      break;
    }

    // Resize the trial vector dimension to contain the new perturbed
    // guess vectors
    new (&TrialVecR) RealCMMap(TVecRMem,this->n_,NTrial+NNotConv);
    if(this->symmetrized_ || !this->hermetian_){
      new (&TrialVecL) RealCMMap(TVecLMem,this->n_,NTrial+NNotConv);
    }
    int INDX = 0;
    for(auto k = 0; k < this->nSek_; k++) {
      // If the residual for root "k" is not converged, construct
      // a perturbed guess vector according to Davidson's diagonal
      // preconditioning scheme.
      //
      // **WARNING** This is only valid for strongly diagonal dominant
      //             matricies. Convergence will be slow (maybe infinitely)
      //             if this criteria is not met.
      if(!resConv[k]) {
        if(this->method_ == SDResponse::CIS || this->method_ == SDResponse::RPA){
          new (&QR) RealCMMap(TVecRMem+(NTrial+INDX)*this->n_,this->n_,1);
          new (&RR) RealCMMap(ResRMem + k*this->n_,this->n_,1);
          if(this->method_ == SDResponse::RPA){
            new (&QL) RealCMMap(TVecLMem+(NTrial+INDX)*this->n_,this->n_,1);
            new (&RL) RealCMMap(ResLMem + k*this->n_,this->n_,1);
          }
          this->sdr_->formPerturbedGuess(ER(k),RR,QR,RL,QL);
        } else {
          for(auto i = 0; i < this->n_; i++) {
            TrialVecR(i,NTrial+INDX) = - ResR.col(k)(i) / ((*this->mat_)(i,i) - ER(k));
          }
        }
        INDX++;
      }
    }
    // Normalize and orthogonalize the new guess vectors to the
    // existing set using QR factorization
    int N = TrialVecR.cols();
    int M = TrialVecR.rows();
    int LDA = TrialVecR.rows();

    dgeqrf_(&M,&N,TrialVecR.data(),&LDA,LAPACK_SCR,LAPACK_SCR+N,&LWORK,&INFO);
    dorgqr_(&M,&N,&N,TrialVecR.data(),&LDA,LAPACK_SCR,LAPACK_SCR+N,&LWORK,&INFO);
    if(this->symmetrized_ || !this->hermetian_){
      dgeqrf_(&M,&N,TrialVecL.data(),&LDA,LAPACK_SCR,LAPACK_SCR+N,&LWORK,&INFO);
      dorgqr_(&M,&N,&N,TrialVecL.data(),&LDA,LAPACK_SCR,LAPACK_SCR+N,&LWORK,&INFO);
    }
    NOld = NTrial;
    NNew = NNotConv;
    NTrial += NNew;

    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    output << "Davidson Micro Iteration took " << std::fixed 
           << elapsed.count() << " secs" << endl << endl;
  } 
  delete [] SCR; // Cleanup scratch space
}
} // namespace ChronusQ
