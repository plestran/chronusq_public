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

template<typename T>
void QuasiNewton2<T>::run(){
  time_t currentClockTime;

  time(&currentClockTime);

  (*this->out_) << "Quasi-Newton Calculation Started: " << 
                ctime(&currentClockTime);

  if(this->qnObj_->specialAlgorithm_ == FULL_SOLVE){
    this->runFull();
  } else {
    this->runIter();
  }
/*
  this->iniScratchFiles();
  auto start = std::chrono::high_resolution_clock::now();
  auto finish = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = finish - start;
*/
  time(&currentClockTime);
  (*this->out_) << "Quasi-Newton Calculation Completed: " << 
                ctime(&currentClockTime);


}; // QuasiNewton2<T>::run

template <typename T>
void QuasiNewton2<T>::runIter() {
  this->allocIterScr();
  this->iniScratchFiles();
  for(auto iter = 0; iter < this->maxMacroIter_; iter++){
    this->runMicro();
    this->nMacroIter_++;
    this->nTotalIter_ += this->nMicroIter_;
    if(this->isConverged_) break;
  };
  this->cleanupIterScr();
}

template <typename T>
void QuasiNewton2<T>::runFull() {
  if(this->qnObj_->problemType_ == DIAGONALIZATION) {
    cout << "Is Diagonalization" << endl;
    if(this->qnObj_->matrixType_ == HERMETIAN) {
       std::copy_n(this->qnObj_->fullMatrix_,
         this->qnObj_->nSingleDim_*this->qnObj_->nSingleDim_,
         this->qnObj_->solutionVecR_);
       this->stdHermetianDiag('V','L',this->qnObj_->nSingleDim_,
         this->qnObj_->solutionVecR_,this->qnObj_->omega_);
    } else {
      T* E;
      if(std::is_same<double,T>::value)
        E = this->memManager_->template malloc<T>(2*this->qnObj_->nSingleDim_);
      else
        E = this->memManager_->template malloc<T>(this->qnObj_->nSingleDim_);

      this->stdNonHermetianDiag('V','N',this->qnObj_->nSingleDim_,
        this->qnObj_->fullMatrix_,E,this->qnObj_->solutionVecR_,
        this->qnObj_->solutionVecL_);

      std::copy_n(reinterpret_cast<double*>(E),this->qnObj_->nSingleDim_,
        this->qnObj_->omega_);

      RealVecMap Eig(this->qnObj_->omega_,this->qnObj_->nSingleDim_);
      TMap       VR(this->qnObj_->solutionVecR_,
        this->qnObj_->nSingleDim_,this->qnObj_->nSingleDim_); 

      this->eigSrt(VR,Eig);
    }
  }
};

template<typename T>
void QuasiNewton2<T>::runMicro(){
  int NTrial = this->qnObj_->nGuess();
  int NOld   = 0;
  int NNew   = this->qnObj_->nGuess();


  this->readGuess();
  if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL && 
     this->nMacroIter_ == 0) 
    this->symmetrizeTrial(); 

  this->nMicroIter_ = 0;
  for(auto iter = 0; iter < this->maxMicroIter_; iter++){
    this->checkOrthogonality(NTrial);
    this->formLinearTrans(NOld,NNew);
    this->fullProjection(NTrial);
    this->reducedDimDiag(NTrial);
    this->reconstructSolution(NTrial);
    this->generateResiduals(NTrial);

    this->nMicroIter_++;
    int NNotConv;
    auto resConv = this->checkConvergence(NTrial,NNotConv);
    this->isConverged_ = (NNotConv == 0);
    bool restart       = (NTrial + NNotConv > this->maxSubSpace_);
    if(this->isConverged_) {
      break;
    }
    if(restart) CErr("Restart NYI",(*this->out_));

    this->formNewGuess(resConv,NTrial,NNotConv,NOld,NNew);
  };// for iter in [0, maxMicroIter)
}; // QuasiNewton2<T>::runMicro

template<typename T>
void QuasiNewton2<T>::checkOrthogonality(int &NTrial){
//this->checkLinearDependence(NTrial);
  this->orthogonalize(NTrial);
}; // QuasiNewton2<T>::checkOrthogonality

template<typename T>
void QuasiNewton2<T>::formLinearTrans(const int NOld, const int NNew){

  (*this->out_) << "Performing Linear Transformation in QuasiNewton" << endl;
  auto N = this->qnObj_->nSingleDim();

  TMap NewSR  (this->SigmaRMem_ + (NOld*N),N,NNew);
  TMap NewVecR(this->TRMem_     + (NOld*N),N,NNew);

  TMap NewSL  (this->SigmaLMem_,0,0);
  TMap NewVecL(this->TLMem_    ,0,0);
  TMap NewRhoR(this->RhoRMem_  ,0,0);
  TMap NewRhoL(this->RhoLMem_  ,0,0);

  if(this->qnObj_->matrixType_ == HERMETIAN_GEP)
    new (&NewRhoR) TMap(this->RhoRMem_   + (NOld*N),N,NNew);

  if(this->qnObj_->needsLeft_ or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    new (&NewSL  ) TMap(this->SigmaLMem_ + (NOld*N),N,NNew);
    new (&NewVecL) TMap(this->TLMem_     + (NOld*N),N,NNew);
    if(this->qnObj_->matrixType_ == HERMETIAN_GEP)
      new (&NewRhoL) TMap(this->RhoLMem_   + (NOld*N),N,NNew);
  }

//prettyPrint(cout,NewVecR,"TR");
//prettyPrint(cout,NewSR,"SR");
  this->qnObj_->linearTrans(NewVecR,NewVecL,NewSR,NewSL,NewRhoR,NewRhoL);
}; // QuasiNewston<T>::formLinearTrans

template<typename T>
void QuasiNewton2<T>::fullProjection(const int NTrial){

  (*this->out_) << "Performing Full Projection in QuasiNewton" << endl;

  auto N = this->qnObj_->nSingleDim();

  TMap SigmaR  (this->SigmaRMem_  , N     , NTrial);
  TMap XTSigmaR(this->XTSigmaRMem_, NTrial, NTrial);
  TMap TVecR   (this->TRMem_      , N     , NTrial);

  TMap SigmaL  (this->SigmaLMem_  , 0, 0);
  TMap XTSigmaL(this->XTSigmaLMem_, 0, 0);
  TMap TVecL   (this->TLMem_      , 0, 0);
  TMap RhoR    (this->RhoRMem_    , 0, 0);
  TMap XTRhoR  (this->XTRhoRMem_  , 0, 0);
  TMap RhoL    (this->RhoLMem_    , 0, 0);
  TMap XTRhoL  (this->XTRhoLMem_  , 0, 0);

  if(this->qnObj_->matrixType_ == HERMETIAN_GEP){
    new (&RhoR    ) TMap(this->RhoRMem_    , N     , NTrial);
    new (&XTRhoR  ) TMap(this->XTRhoRMem_  , NTrial, NTrial);
  }

  if(this->qnObj_->needsLeft_ or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    new (&SigmaL  ) TMap(this->SigmaLMem_  , N     , NTrial);
    new (&XTSigmaL) TMap(this->XTSigmaLMem_, NTrial, NTrial);
    new (&TVecL   ) TMap(this->TLMem_      , N     , NTrial);
    if(this->qnObj_->matrixType_ == HERMETIAN_GEP){
      new (&RhoL    ) TMap(this->RhoLMem_    , N     , NTrial);
      new (&XTRhoL  ) TMap(this->XTRhoLMem_  , NTrial, NTrial);
    }
  }

  XTSigmaR = TVecR.adjoint() * SigmaR;
  if(this->qnObj_->matrixType_ == HERMETIAN_GEP)
    XTRhoR   = TVecR.adjoint() * RhoR;

  if(this->qnObj_->needsLeft_ or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    XTSigmaL = TVecL.adjoint() * SigmaL;
    if(this->qnObj_->matrixType_ == HERMETIAN_GEP)
      XTRhoL   = TVecL.adjoint() * RhoL;
  }

  if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL)
    this->buildSuperMatricies(NTrial);
}; // QuasiNewton2<T>::fullProjection

template<typename T>
void QuasiNewton2<T>::reconstructSolution(const int NTrial){

  (*this->out_) << "Reconstructing Solution Vectors in QuasiNewton" << endl;

  auto N = this->qnObj_->nSingleDim();

  TMap XTSigmaR(this->XTSigmaRMem_, NTrial, NTrial);
  TMap TVecR   (this->TRMem_      , N     , NTrial);
  TMap UR      (this->URMem_      , N     , NTrial);
  TMap VR      (this->qnObj_->solutionVecR_, N, this->qnObj_->nSek_);

  RealVecMap W (this->qnObj_->omega_,this->qnObj_->nSek_);

  TMap XTSigmaL(this->XTSigmaLMem_, 0, 0);
  TMap TVecL   (this->TLMem_      , 0, 0);
  TMap UL      (this->ULMem_      , 0, 0);
  TMap VL      (this->qnObj_->solutionVecL_,0, 0);

  RealVecMap ER(this->EPersist_,NTrial);
  if(this->qnObj_->needsLeft() or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    new (&XTSigmaL) TMap(this->XTSigmaLMem_, NTrial, NTrial);
    new (&TVecL   ) TMap(this->TLMem_      , N     , NTrial);
    new (&UL      ) TMap(this->ULMem_      , N     , NTrial);
    new (&VL      ) TMap(this->qnObj_->solutionVecL_, N, this->qnObj_->nSek_);
  }


  UR = TVecR * XTSigmaR;
  if(this->qnObj_->needsLeft() or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL) UL = TVecL * XTSigmaL;

  VR = UR.block(0,0,N,this->qnObj_->nSek());
  W  = ER.head(this->qnObj_->nSek());

  if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL)
    VR += UL.block(0,0,N,this->qnObj_->nSek());
  else if(this->qnObj_->needsLeft())
    VL = UL.block(0,0,N,this->qnObj_->nSek());
  
}; // QuasiNewton2<T>::reconstructSolution

template<typename T>
void QuasiNewton2<T>::generateResiduals(const int NTrial){

  (*this->out_) << "Generating Residuals in QuasiNewton" << endl;

  auto N = this->qnObj_->nSingleDim();

  TMap SigmaR  (this->SigmaRMem_  , N     , NTrial);
  TMap XTSigmaR(this->XTSigmaRMem_, NTrial, NTrial);
  TMap TVecR   (this->TRMem_      , N     , NTrial);
  TMap UR      (this->URMem_      , N     , NTrial);
  TMap ResR    (this->ResRMem_    , N     , NTrial);

  RealVecMap E(this->EPersist_,NTrial);

  TMap SigmaL  (this->SigmaLMem_  , 0, 0);
  TMap XTSigmaL(this->XTSigmaLMem_, 0, 0);
  TMap TVecL   (this->TLMem_      , 0, 0);
  TMap RhoR    (this->RhoRMem_    , 0, 0);
  TMap RhoL    (this->RhoLMem_    , 0, 0);
  TMap UL      (this->ULMem_      , 0, 0);
  TMap ResL    (this->ResLMem_    , 0, 0);


  if(this->qnObj_->matrixType_ == HERMETIAN_GEP)
    new (&RhoR    ) TMap(this->RhoRMem_    , N     , NTrial);

  if(this->qnObj_->needsLeft() or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    new (&SigmaL  ) TMap(this->SigmaLMem_  , N     , NTrial);
    new (&XTSigmaL) TMap(this->XTSigmaLMem_, NTrial, NTrial);
    new (&TVecL   ) TMap(this->TLMem_      , N     , NTrial);
    new (&UL      ) TMap(this->ULMem_      , N     , NTrial);
    new (&ResL    ) TMap(this->ResLMem_    , N     , NTrial);

    if(this->qnObj_->matrixType_ == HERMETIAN_GEP)
      new (&RhoL    ) TMap(this->RhoLMem_    , N     , NTrial);
  }

  if(this->qnObj_->matrixType_ == QNMatrixType::HERMETIAN)
    ResR = SigmaR * XTSigmaR - UR * E.asDiagonal();
  else if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    ResR = SigmaR * XTSigmaR - RhoR * XTSigmaL * E.asDiagonal();
    ResL = SigmaL * XTSigmaL - RhoL * XTSigmaR * E.asDiagonal();
  }
}; // QuasiNewton2<T>::generateResiduals

template<typename T>
std::vector<bool> QuasiNewton2<T>::checkConvergence(const int NTrial, 
  int &NNotConv) {

  (*this->out_) << "Checking Convergence in QuasiNewton" << endl;

  auto N = this->qnObj_->nSingleDim();

  TMap     ResR   (this->ResRMem_    , N     , NTrial);
  RealVecMap ER   (this->EPersist_      , NTrial        );

  TMap ResL    (this->ResLMem_    , 0, 0);
  if(this->qnObj_->needsLeft() or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    new (&ResL    ) TMap(this->ResLMem_    , N     , NTrial);
  }

  // Vector to store convergence info
  std::vector<bool> resConv;
  NNotConv = 0;

  // Print Header for convergence check
  (*this->out_) << std::fixed << std::setprecision(12);
  (*this->out_) << "  Checking Quasi-Newton Convergence:" << endl;
  (*this->out_) << "    " << std::setw(8)  << " " << std::setw(32) 
                << std::left << "    Roots at Current Iteration:";
  (*this->out_) << std::setw(32) << std::left 
                << "    (Max) Norm of Residual(s):" 
                << endl;

  // Loop over NSek residual vectors. Decide from which residuals
  // will be made perturbed guess vectors
  for(auto k = 0; k < this->qnObj_->nSek(); k++) {
    double NORM = ResR.col(k).norm();
    
    // If we need the left vectors, decide converence based
    // on maximum of left and right residual norms
    if(this->qnObj_->needsLeft() or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL)
      NORM = std::max(NORM,ResL.col(k).norm());

    // Print current state of k-th root
    (*this->out_) << "    " << std::setw(12) 
                  << "State " + std::to_string(k+1) + ":";
    (*this->out_) << std::setw(32) << std::left << std::fixed 
                  << this->EPersist_[k];
    (*this->out_) << std::setw(32) << std::left << std::scientific << NORM;

    // If the norm meets the tolerence, say that the root is converged
    // else say that it's not (this is used to build new guess vectors
    if(NORM < this->residualTol_) {
      resConv.push_back(true);
      (*this->out_) << "     Root has converged" << endl;
    } else {
      resConv.push_back(false); 
      NNotConv++;
      (*this->out_) << "     Root has not converged" << endl;
    }
  }; // loop k in nSek

  return resConv;
}; // QuasiNewton2<T>::checkConvergence

template<typename T>
void QuasiNewton2<T>::formNewGuess(std::vector<bool> &resConv, int &NTrial,
  int &NNotConv, int &NOld, int &NNew) {

  (*this->out_) << "Forming New Guess Vectors in QuasiNewton" << endl;
  
  auto N = this->qnObj_->nSingleDim();

  TMap TVecR(this->TRMem_  , N, NTrial+NNotConv);
  TMap ResR (this->ResRMem_, N, NTrial         );
  RealVecMap E (this->EPersist_  , NTrial            );   

  TMap TVecL   (this->TLMem_  , 0, 0);
  TMap ResL    (this->ResLMem_, 0, 0);

  if(this->qnObj_->needsLeft() or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
    new (&TVecL) TMap(this->TLMem_  , N, NTrial+NNotConv);
    new (&ResL ) TMap(this->ResLMem_, N, NTrial         );
  }


  TMap QR      (this->TRMem_  , 0, 0);
  TMap RR      (this->ResRMem_, 0, 0);
  TMap QL      (this->TLMem_  , 0, 0);
  TMap RL      (this->ResLMem_, 0, 0);

  int INDX = 0;
  for(auto k = 0; k < this->qnObj_->nSek(); k++) {
    if(!resConv[k]) {
      auto RRSt = this->ResRMem_ + k*N;
      auto QRSt = this->TRMem_   + (NTrial+INDX)*N;
      new (&RR) TMap(RRSt, N, 1);
      new (&QR) TMap(QRSt, N, 1);
      if(this->qnObj_->needsLeft() or this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL){
        auto RLSt = this->ResLMem_ + k*N;
        auto QLSt = this->TLMem_   + (NTrial+INDX)*N;
        new (&RL) TMap(RLSt, N, 1);
        new (&QL) TMap(QLSt, N, 1);
      }

      if(this->qnObj_->guessType_ == QNGuessType::RESIDUAL_DAVIDSON) 
        this->davResidualGuess(E(k),RR,QR,RL,QL);
      INDX++;
    }
  } // loop k in nSek

  // Update the number of vectors
  NOld   =  NTrial;
  NNew   =  NNotConv;
  NTrial += NNew;

  this->writeTrialVectors(NTrial);

}; // QuasiNewton2<T>::formNewGuess

template<typename T>
void QuasiNewton2<T>::davResidualGuess(T Omega, const TMap &RR, TMap &QR,
  const TMap &RL, TMap &QL){

  (*this->out_) << "Generating Davidson Residual Guess in QuasiNewton" << endl;

  if(this->qnObj_->specialAlgorithm_ == NOT_SPECIAL)
    this->davStdResidualGuess(Omega,RR,QR);
  else if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL)
    this->davSymmResidualGuess(Omega,RR,QR,RL,QL);

}; // QuasiNewton2<T>::davResidualGuess

template<typename T>
void QuasiNewton2<T>::davStdResidualGuess(T Omega, const TMap &RR, TMap &QR){

  (*this->out_) << "Using Standard Form for Davidson Guess in QuasiNewton" 
    << endl;

  auto N = this->qnObj_->nSingleDim();
  for(auto i = 0; i < N; i++)
    QR(i) = -RR(i) / (this->qnObj_->diag_[i] - Omega);
}; // QuasiNewton2<T>::davStdResidualGuess

template<typename T>
void QuasiNewton2<T>::davSymmResidualGuess(T Omega, const TMap &RR, TMap &QR,
  const TMap &RL, TMap &QL) {

  (*this->out_) << "Using Symmetrized Form for Davidson Guess in QuasiNewton" 
    << endl;

  auto N = this->qnObj_->nSingleDim();

  for(auto i = 0; i < N; i++){
    QR(i) = RR(i) * this->qnObj_->diag_[i];
    QL(i) = RL(i) * this->qnObj_->diag_[i];
  }
  for(auto i = 0; i < N/2; i++){
    QR(i)       += Omega * RL(i);
    QR(N/2 + i) -= Omega * RL(N/2 + i);
    QL(i)       += Omega * RR(i);
    QL(N/2 + i) -= Omega * RR(N/2 + i);
  }
  for(auto i = 0; i < N; i++){
    QR(i) /= std::pow(this->qnObj_->diag_[i],2.0) - std::pow(Omega,2.0);
    QL(i) /= std::pow(this->qnObj_->diag_[i],2.0) - std::pow(Omega,2.0);
  }
}; // QuasiNewton2<T>::davSymmResidualGuess
