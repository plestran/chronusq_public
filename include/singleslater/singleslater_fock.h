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

/********************************
 * Form Perturbation Tensor (G) *
 ********************************/
// Only available for Libint Integrals
#ifdef USE_LIBINT
template<typename T>
void SingleSlater<T>::formPT(){

  std::vector<std::reference_wrapper<TMap>> mats;
  std::vector<std::reference_wrapper<TMap>> ax;
  std::vector<AOIntegrals::ERI_CONTRACTION_TYPE> contList;
  std::vector<double> scalingFactors;

  mats.emplace_back(*this->onePDMScalar_);
  mats.emplace_back(*this->onePDMScalar_);
  ax.emplace_back(*this->PTScalar_);
  ax.emplace_back(*this->PTScalar_);
  contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::COULOMB);
  contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
  scalingFactors.push_back(1.0);
  scalingFactors.push_back(this->xHF_);

  this->PTScalar_->setZero();

  if(this->nTCS_ == 2 or !this->isClosedShell) {
    mats.emplace_back(*this->onePDMMz_);
    ax.emplace_back(*this->PTMz_);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    scalingFactors.push_back(this->xHF_);
    this->PTMz_->setZero();
  }

  if(this->nTCS_ == 2){
    mats.emplace_back(*this->onePDMMy_);
    mats.emplace_back(*this->onePDMMx_);
    ax.emplace_back(*this->PTMy_);
    ax.emplace_back(*this->PTMx_);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    scalingFactors.push_back(this->xHF_);
    scalingFactors.push_back(this->xHF_);
    this->PTMy_->setZero();
    this->PTMx_->setZero();
  }

  this->aointegrals_->newTwoEContract(mats,ax,contList,scalingFactors);

/*
  if(this->nTCS_ == 1 && this->isClosedShell) {
    (*this->PTA_) *= 0.5; // This is effectively a gather operation
  } else {
    std::vector<std::reference_wrapper<TMap>> toGather;
    toGather.emplace_back(*this->PTScalar_);
    toGather.emplace_back(*this->PTMz_);
    if(this->nTCS_ == 1)
      Quantum<T>::spinGather((*this->PTA_),(*this->PTB_),toGather);
    else {
      toGather.emplace_back(*this->PTMy_);
      toGather.emplace_back(*this->PTMx_);
      Quantum<T>::spinGather((*this->PTA_),toGather);
    }

  }
*/

  if(this->printLevel_ >= 3 && getRank() == 0) this->printPT();
//prettyPrint(cout,*this->onePDMA_,"P in PT");
//prettyPrint(cout,*this->PTA_,"PT in PT");
}
#endif

/********************
 * Form Fock Matrix *
 ********************/
template<typename T>
void SingleSlater<T>::formFock(bool increment){
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes before fock build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
//if(!this->haveDensity) this->formDensity();
#ifndef USE_LIBINT
  if(getRank() == 0){
    // Not even sure if these guys still work let alone are MPI
    // capable
    if(!this->haveCoulomb) this->formCoulomb();
    if(!this->haveExchange) this->formExchange();
  }
#else
  // All MPI processes go to FormPT
  this->formPT();
#endif

  if(getRank() == 0) {
    if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

    if (this->isDFT) this->formVXC_new();
  
/*
    if(this->nTCS_ == 1 && this->isClosedShell) {
      if(!increment) {
        cout << "Not Incrementing" << endl;
        this->fockA_->setZero();
        this->fockA_->real() += (*this->aointegrals_->coreH_);
        this->aointegrals_->addElecDipole(*this->fockA_,this->elecField_);
      }
      (*this->fockA_) += (*this->PTA_);
      if(this->isDFT){ 
        (*this->fockA_) += (*this->vXA_);
      }
    } else {
      if(!increment) {
        this->fockScalar_->setZero();
        this->fockMz_->setZero();
        this->fockScalar_->real() += (*this->aointegrals_->coreH_);
        
        if(this->nTCS_ == 2) {
            this->fockMx_->setZero();
            this->fockMy_->setZero();
          if(this->aointegrals_->doX2C){
            this->fockMx_->real() = 2*(*this->aointegrals_->oneEmx_);
            this->fockMy_->real() = 2*(*this->aointegrals_->oneEmy_);
            this->fockMz_->real() = 2*(*this->aointegrals_->oneEmz_);
            // -----------------------------------
            // SCALE SO parts by 'i'
            // -----------------------------------
            Quantum<T>::complexMyScale(*this->fockMx_);
            Quantum<T>::complexMyScale(*this->fockMy_);
            Quantum<T>::complexMyScale(*this->fockMz_);
            // -----------------------------------
          }
        }

        this->aointegrals_->addElecDipole(*this->fockScalar_,this->elecField_);
        (*this->fockScalar_) *= 2.0;
      }
      (*this->fockScalar_)      += (*this->PTScalar_);        
      (*this->fockMz_)          += (*this->PTMz_);

      // FIXME: Needs to ge generalized for the 2C case
      if(this->nTCS_ == 1 && this->isDFT){
          (*this->fockScalar_) += (*this->vXA_);
          (*this->fockScalar_) += (*this->vXB_);
          (*this->fockMz_)     += (*this->vXA_);
          (*this->fockMz_)     -= (*this->vXB_);
      }

      // Transform the Fock for CUHF
      if(this->Ref_ == CUHF) {
        this->formNO();
        this->fockCUHF();
      }

      std::vector<std::reference_wrapper<TMap>> toGather;
      toGather.emplace_back(*this->fockScalar_);
      toGather.emplace_back(*this->fockMz_);
      if(this->nTCS_ == 1)
        Quantum<T>::spinGather(*this->fockA_,*this->fockB_,toGather);
      else {
        (*this->fockMx_) += (*this->PTMx_);
        (*this->fockMy_) += (*this->PTMy_);
        toGather.emplace_back(*this->fockMy_);
        toGather.emplace_back(*this->fockMx_);
        Quantum<T>::spinGather(*this->fockA_,toGather);

      //  prettyPrint(this->fileio_->out,(*this->fockA_)," fockA_ ");
      }
    }
*/

    if(!increment) {
      this->fockScalar_->setZero();

      if(this->nTCS_ == 2 or !this->isClosedShell)
        this->fockMz_->setZero();

      if(this->nTCS_ == 2){
        this->fockMy_->setZero();
        this->fockMx_->setZero();
      }

      (*this->fockScalar_) = this->aointegrals_->coreH_->template cast<T>();
      this->aointegrals_->subElecDipole(*this->fockScalar_,this->elecField_);

      if(this->nTCS_ == 2 and this->aointegrals_->doX2C) {
        (*this->fockMx_) = this->aointegrals_->oneEmx_->template cast<T>();
        (*this->fockMy_) = this->aointegrals_->oneEmy_->template cast<T>();
        (*this->fockMz_) = this->aointegrals_->oneEmz_->template cast<T>();
        // -----------------------------------
        // SCALE SO parts by 'i'
        // -----------------------------------
        /*
        Quantum<T>::complexMyScale(*this->fockMx_);
        Quantum<T>::complexMyScale(*this->fockMy_);
        Quantum<T>::complexMyScale(*this->fockMz_);
        */
        (*this->fockMx_) *= ComplexScale<T>();
        (*this->fockMy_) *= ComplexScale<T>();
        (*this->fockMz_) *= ComplexScale<T>();
        // -----------------------------------
      }
      for(auto iF = fock_.begin(); iF != fock_.end(); iF++)
        *(*iF) *= 2;
    }

/*
      this->fileio_->out << "Before PT" << endl;
      this->printFock();
      for(auto i = 0; i < this->fockScalar_->size(); i++){
        cout << this->fockScalar_->data()[i] << " + " << this->PTScalar_->data()[i] << " = " << this->fockScalar_->data()[i] + this->PTScalar_->data()[i] << endl;
        this->fockScalar_->data()[i] += this->PTScalar_->data()[i];
      }
*/
    for(auto iF = 0; iF < fock_.size(); iF++) {
      *(fock_[iF]) += *(PT_[iF]);
      if(this->isDFT)
        (*fock_[iF]) += (*vXC_[iF]);
    }
/*
    (*this->fockScalar_) += (*this->PTScalar_);
*/
/*
      this->fileio_->out << "After PT" << endl;
      this->printFock();
*/


    if(this->printLevel_ >= 2) this->printFock(); 
  }
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes after fock build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
};



template<typename T>
void SingleSlater<T>::formVXC_new(){
//Timing
//  this->screenVxc = false;

  int nthreads = omp_get_max_threads();

  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point finish;
  std::chrono::duration<double> duration_formVxc;
  std::chrono::duration<double> T1(0.0);
  std::chrono::duration<double> T2(0.0);
  std::chrono::duration<double> T3(0.0);
  std::chrono::duration<double> T4(0.0);
  std::chrono::duration<double> T5(0.0);
  std::chrono::duration<double> T6(0.0);
  
  std::vector<std::chrono::duration<double>> 
    TF(this->dftFunctionals_.size(),std::chrono::duration<double>(0.0));
/*
  if(this->printLevel_ >= 3) {
    start = std::chrono::high_resolution_clock::now();
  }
*/
//  bool isGGA   = true;
//  bool isGGA   = false;
/*
  RealMatrix SCRATCH2(this->nBasis_,this->nBasis_);
  RealMatrix SCRATCH2X(this->nBasis_,this->nBasis_);
  RealMatrix SCRATCH2Y(this->nBasis_,this->nBasis_);
  RealMatrix SCRATCH2Z(this->nBasis_,this->nBasis_);
  VectorXd   SCRATCH1(this->nBasis_);
  VectorXd   SCRATCH1X(this->nBasis_);
  VectorXd   SCRATCH1Y(this->nBasis_);
  VectorXd   SCRATCH1Z(this->nBasis_);
*/
  
  int NDer = 0;
  if(this->isGGA) NDer = 1; 

  // Parallel w mem manager
  std::vector<RealVecMap> SCRATCH1;
  std::vector<RealVecMap> SCRATCH1X;
  std::vector<RealVecMap> SCRATCH1Y;
  std::vector<RealVecMap> SCRATCH1Z;


  for(auto ithread = 0; ithread < nthreads; ithread++){
    SCRATCH1.emplace_back(
      this->memManager_->template malloc<double>(this->nBasis_),
      this->nBasis_);
    SCRATCH1.back().setZero();

    if(this->isGGA) {
      SCRATCH1X.emplace_back(
        this->memManager_->template malloc<double>(this->nBasis_),
        this->nBasis_);
      SCRATCH1Y.emplace_back(
        this->memManager_->template malloc<double>(this->nBasis_),
        this->nBasis_);
      SCRATCH1Z.emplace_back(
        this->memManager_->template malloc<double>(this->nBasis_),
        this->nBasis_);

      SCRATCH1X.back().setZero();
      SCRATCH1Y.back().setZero();
      SCRATCH1Z.back().setZero();
    }
  }

/*
  std::vector<VectorXd> SCRATCH1(nthreads,VectorXd(this->nBasis_));
  std::vector<VectorXd> SCRATCH1X(nthreads,VectorXd(this->nBasis_));
  std::vector<VectorXd> SCRATCH1Y(nthreads,VectorXd(this->nBasis_));
  std::vector<VectorXd> SCRATCH1Z(nthreads,VectorXd(this->nBasis_));
*/

/*
  std::array<double,3>  drhoT = {0.0,0.0,0.0}; ///< array TOTAL density gradient components
  std::array<double,3>  drhoS = {0.0,0.0,0.0}; ///< array SPIN  density gradient components
  std::array<double,3>  drhoA = {0.0,0.0,0.0}; ///< array ALPHA  density gradient components
  std::array<double,3>  drhoB = {0.0,0.0,0.0}; ///< array BETA  density gradient components
  RealVecMap GradRhoT(&drhoT[0],3);
  RealVecMap GradRhoS(&drhoS[0],3);
  RealVecMap GradRhoA(&drhoA[0],3);
  RealVecMap GradRhoB(&drhoB[0],3);

  double rhoA;
  double rhoB;
  double gammaAA = 0.0;
  double gammaBB = 0.0;
  double gammaAB = 0.0;
*/

  std::vector<std::array<double,3>>  drhoT(nthreads,{0.0,0.0,0.0}); ///< array TOTAL density gradient components
  std::vector<std::array<double,3>>  drhoS(nthreads,{0.0,0.0,0.0}); ///< array SPIN  density gradient components
  std::vector<std::array<double,3>>  drhoA(nthreads,{0.0,0.0,0.0}); ///< array ALPHA  density gradient components
  std::vector<std::array<double,3>>  drhoB(nthreads,{0.0,0.0,0.0}); ///< array BETA  density gradient components
  std::vector<RealVecMap> GradRhoT;
  std::vector<RealVecMap> GradRhoS;
  std::vector<RealVecMap> GradRhoA;
  std::vector<RealVecMap> GradRhoB;
  for(auto ithread = 0; ithread < nthreads; ithread++){
    GradRhoT.emplace_back(&drhoT[ithread][0],3);
    GradRhoS.emplace_back(&drhoS[ithread][0],3);
    GradRhoA.emplace_back(&drhoA[ithread][0],3);
    GradRhoB.emplace_back(&drhoB[ithread][0],3);
  }

//std::vector<double> rhoA(nthreads,0.);
//std::vector<double> rhoB(nthreads,0.);
//std::vector<double> gammaAA(nthreads,0.);
//std::vector<double> gammaBB(nthreads,0.);
//std::vector<double> gammaAB(nthreads,0.);


//std::vector<bool> shMap(this->basisset_->nShell()+1);
  int NSkip2(0);
  int NSkip3(0);
  int NSkip4(0);
  int NSkip5(0);
  bool doTimings(false);

  auto Newstart = std::chrono::high_resolution_clock::now();
  auto Newend = std::chrono::high_resolution_clock::now();


//VectorXd OmegaA(this->nBasis_), OmegaB(this->nBasis_);
  // Parallel w mem manager
  std::vector<RealVecMap> OmegaA, OmegaB;
  for(auto ithread = 0; ithread < nthreads; ithread++){
    OmegaA.emplace_back(
      this->memManager_->template malloc<double>(this->nBasis_),
      this->nBasis_);
    OmegaB.emplace_back(
      this->memManager_->template malloc<double>(this->nBasis_),
      this->nBasis_);
    OmegaA.back().setZero();
    OmegaB.back().setZero();
  }
/*
  std::vector<VectorXd> OmegaA(nthreads,VectorXd(this->nBasis_));
  std::vector<VectorXd> OmegaB(nthreads,VectorXd(this->nBasis_));
*/

//  VectorXd DENCOL(this->nBasis_);
  double fact = 2.0;

  std::vector<std::size_t> shSizes;
  for(auto iSh = 0; iSh < this->basisset_->nShell(); iSh++)
    shSizes.push_back(this->basisset_->shells(iSh).size());

  std::vector<std::vector<std::size_t>> closeShells(nthreads);
// If screen in not on loop over all shell for each atoms, no need to
// figure which are the close shell to each point
  if (!this->screenVxc){
    for(auto ithread = 0; ithread < nthreads; ithread++){
      closeShells[ithread].resize(this->basisset_->nShell());
      std::iota(closeShells[ithread].begin(),closeShells[ithread].end(),0);
    }
  }

  int tracker(0);
  std::vector<int> currentCol(nthreads,-1);
  auto valVxc = [&](std::size_t iAtm, ChronusQ::IntegrationPoint &pt, 
    std::vector<KernelIntegrand<T>> &result) -> void {

    int thread_id = omp_get_thread_num();

/*
    SCRATCH1[thread_id].setZero();
    if(this->isGGA) {
      SCRATCH1X[thread_id].setZero();
      SCRATCH1Y[thread_id].setZero();
      SCRATCH1Z[thread_id].setZero();
    }
*/
    double *SCRATCH1DATA = SCRATCH1[thread_id].data();
    double *SCRATCH1XDATA,*SCRATCH1YDATA,*SCRATCH1ZDATA;
    if(this->isGGA) {
      SCRATCH1XDATA = SCRATCH1X[thread_id].data();
      SCRATCH1YDATA = SCRATCH1Y[thread_id].data();
      SCRATCH1ZDATA = SCRATCH1Z[thread_id].data();
    }
    double *OmegaADATA = OmegaA[thread_id].data();
    double *OmegaBDATA = OmegaB[thread_id].data();

    cartGP &GP = pt.pt;


    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();
    // For each new sphere, determine list of close basis shells (if screen ON)
//  if(pt.I == 0 && this->screenVxc) {
    if(pt.J != currentCol[thread_id] && this->screenVxc) {
//    cout << pt.I << " " << pt.J <<endl;
      currentCol[thread_id] = pt.J;
      closeShells[thread_id].clear();
      double DX = (*this->molecule_->cart())(0,iAtm) - bg::get<0>(GP);
      double DY = (*this->molecule_->cart())(1,iAtm) - bg::get<1>(GP);
      double DZ = (*this->molecule_->cart())(2,iAtm) - bg::get<2>(GP);

      double R = DX*DX + DY*DY + DZ*DZ;
      R = std::sqrt(R);

      for(auto jAtm = 0; jAtm < this->molecule_->nAtoms(); jAtm++){
        double RAB = (*this->molecule_->rIJ())(iAtm,jAtm);
        double DIST = std::abs(RAB - R);
        for(auto iSh = 0; iSh < this->basisset_->nShell(); iSh++){
          if(this->basisset_->mapSh2Cen(iSh)-1 != jAtm) continue;
          if(this->basisset_->radCutSh()[iSh] > DIST)
            closeShells[thread_id].push_back(iSh);
        }
      }
    }
    if(doTimings){ 
      Newend = std::chrono::high_resolution_clock::now();
      T1 += Newend - Newstart;
    }


    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();

    // Compute value of "close" basis functions at the given point
    for(auto iShell : closeShells[thread_id]) {
      int shSize= shSizes[iShell];
      int iSt = this->basisset_->mapSh2Bf(iShell);
      double * buff = this->basisset_->basisDEval(NDer,this->basisset_->shells(iShell),&pt.pt);

      std::memcpy(SCRATCH1DATA + iSt,buff,shSize*sizeof(double));
      if(NDer>0){
        double * ds1EvalX = buff + shSize;
        double * ds1EvalY = ds1EvalX + shSize;
        double * ds1EvalZ = ds1EvalY + shSize;
        std::memcpy(SCRATCH1XDATA + iSt,ds1EvalX,shSize*sizeof(double));
        std::memcpy(SCRATCH1YDATA + iSt,ds1EvalY,shSize*sizeof(double));
        std::memcpy(SCRATCH1ZDATA + iSt,ds1EvalZ,shSize*sizeof(double));
      }
    }
    if(doTimings){
      Newend = std::chrono::high_resolution_clock::now();
      T2 += Newend - Newstart;
    }

    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();

    // Determine if we computed Zeros
    double S1Norm = SCRATCH1[thread_id].norm();
  //double S1XNorm = SCRATCH1X[thread_id].norm();
  //double S1YNorm = SCRATCH1Y[thread_id].norm();
  //double S1ZNorm = SCRATCH1Z[thread_id].norm();

    if(doTimings){
      Newend = std::chrono::high_resolution_clock::now();
      T3 += Newend - Newstart;
    }

    if(this->screenVxc && S1Norm < this->epsScreen) {NSkip4++; return;}
//A    else if(!this->screenVxc && S1Norm < this->epsScreen) {cout << "screenOFF 0" <<endl;}


    double rhoT(0.0);
    double rhoS(0.0);
    T Tt(0.0), Ts(0.0);
    T Pt(0.0), Ps(0.0);
    GradRhoT[thread_id].setZero();
    GradRhoS[thread_id].setZero();

    // Evaluate density and optionally density gradient at the given point
      
    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();
    T * DENT, *DENS;
    // Loop over close shells "I"
    for(auto iShell : closeShells[thread_id]) {
      int iSz= shSizes[iShell];
      int iSt = this->basisset_->mapSh2Bf(iShell);
      for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
        Tt = 0.0; Ts = 0.0;
/*
        if(this->nTCS_ == 1 && this->isClosedShell){
          DENT = this->onePDMA_->data() + iBf*this->nBasis_;
        } else {
          DENT = this->onePDMScalar_->data() + iBf*this->nBasis_;
          DENS = this->onePDMMz_->data() + iBf*this->nBasis_;
        }
*/
        DENT = this->onePDMScalar_->data() + iBf*this->nBasis_;
        if(this->nTCS_ == 2 or !this->isClosedShell)
          DENS = this->onePDMMz_->data() + iBf*this->nBasis_;
        

        // Loop over close shells "J"
        for(auto jShell : closeShells[thread_id]) {
          int jSz= shSizes[jShell];
          int jSt = this->basisset_->mapSh2Bf(jShell);
          for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
            Pt = fact * DENT[jBf];
            if (this->screenVxc) {
              if(std::abs(Pt) > this->epsScreen){
                Tt += Pt * SCRATCH1DATA[jBf];
              } 
            } else {
//               cout << "SCREEN OFF 2" <<endl;
               Tt += Pt * SCRATCH1DATA[jBf];
            } //Screening
          } // jBf
        } // jShell      

        if(this->nTCS_ == 1 && !this->isClosedShell){
          for(auto jShell : closeShells[thread_id]) {
            int jSz= shSizes[jShell];
            int jSt = this->basisset_->mapSh2Bf(jShell);
            for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
              Ps = fact * DENS[jBf];
              if (this->screenVxc) {
                if( std::abs(Ps) > this->epsScreen){
                  Ts += Ps * SCRATCH1DATA[jBf];
                }
              } else {
//                cout << "SCREEN OFF 3" <<endl;
                Ts += Ps * SCRATCH1DATA[jBf];
              }//Screening
            } // jBf
          } // jShell      
        } //UKS
        rhoT += std::real(Tt) * SCRATCH1DATA[iBf];
        if(this->nTCS_ == 1 && !this->isClosedShell){
          rhoS += std::real(Ts) * SCRATCH1DATA[iBf];
        }  //UKS
        if(NDer > 0) {
          drhoT[thread_id][0] += std::real(Tt) * SCRATCH1XDATA[iBf];
          drhoT[thread_id][1] += std::real(Tt) * SCRATCH1YDATA[iBf];
          drhoT[thread_id][2] += std::real(Tt) * SCRATCH1ZDATA[iBf];
          if(this->nTCS_ == 1 && !this->isClosedShell){
            drhoS[thread_id][0] += std::real(Ts) * SCRATCH1XDATA[iBf];
            drhoS[thread_id][1] += std::real(Ts) * SCRATCH1YDATA[iBf];
            drhoS[thread_id][2] += std::real(Ts) * SCRATCH1ZDATA[iBf];
          } //UKS
        } //GGA
      } // iBf
    } // iShell
    if(doTimings){
      Newend = std::chrono::high_resolution_clock::now();
      T4 += Newend - Newstart;
    }

    rhoT *= 0.5;
    rhoS *= 0.5;
    double rhoA = 0.5 * (rhoT + rhoS);
    double rhoB = 0.5 * (rhoT - rhoS);
    double gammaAA,gammaBB,gammaAB;
    if( NDer > 0 ){
      GradRhoA[thread_id].noalias() = 0.5 * (GradRhoT[thread_id] + GradRhoS[thread_id]);
      GradRhoB[thread_id].noalias() = 0.5 * (GradRhoT[thread_id] - GradRhoS[thread_id]);
      gammaAA = GradRhoA[thread_id].dot(GradRhoA[thread_id]);           
      gammaBB = GradRhoB[thread_id].dot(GradRhoB[thread_id]);           
      gammaAB = GradRhoA[thread_id].dot(GradRhoB[thread_id]);           
    }
    
//these if statements prevent numerical instability with zero guesses
    if ((std::abs(rhoT)) <= std::numeric_limits<double>::epsilon() ){return;}
    if (rhoT    <= 0.0 ) {CErr("Numerical noise in the density");}
    if(this->screenVxc && rhoT < this->epsScreen){return;} 

    // Evaluate density functional
    DFTFunctional::DFTInfo kernelXC;
    for(auto i = 0; i < this->dftFunctionals_.size(); i++){
      if(doTimings) 
        Newstart = std::chrono::high_resolution_clock::now();
//      if (NDer > 0) {
      kernelXC += this->dftFunctionals_[i]->eval(
          rhoA,rhoB,gammaAA,gammaAB,gammaBB);
//      } else {
//        kernelXC += this->dftFunctionals_[i]->eval(rhoA, rhoB);
//      }
      if(doTimings){
        Newend = std::chrono::high_resolution_clock::now();
        TF[i] += Newend - Newstart;
      }
    } // loop over kernels


    // Evaluate Functional derivatives

    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();
    if(NDer > 0 ) {
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        double GA(GradRhoA[thread_id](iXYZ)), GB(GradRhoB[thread_id](iXYZ));
          GradRhoA[thread_id](iXYZ) = pt.weight * 
            ( 2.0 * GA * kernelXC.ddgammaAA + GB * kernelXC.ddgammaAB);
          GradRhoB[thread_id](iXYZ) = pt.weight * 
            ( 2.0 * GB * kernelXC.ddgammaBB + GA * kernelXC.ddgammaAB);
      } // Grad
      OmegaA[thread_id].setZero(); 
      for(auto iShell : closeShells[thread_id]) {
        int iSz= shSizes[iShell];
        int iSt = this->basisset_->mapSh2Bf(iShell);
     
        // FIXME: This can be vectorized
        for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
          OmegaADATA[iBf] += GradRhoA[thread_id](0) * SCRATCH1XDATA[iBf] +
                             GradRhoA[thread_id](1) * SCRATCH1YDATA[iBf] +
                             GradRhoA[thread_id](2) * SCRATCH1ZDATA[iBf];
        } // iBf
      } // iShell
      if(this->nTCS_ == 1 && !this->isClosedShell){
        OmegaB[thread_id].setZero(); 
        for(auto iShell : closeShells[thread_id]) {
           int iSz= shSizes[iShell];
           int iSt = this->basisset_->mapSh2Bf(iShell);

           // FIXME: This can be vectorized
           for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
              OmegaBDATA[iBf] += GradRhoB[thread_id](0) * SCRATCH1XDATA[iBf] + 
                                 GradRhoB[thread_id](1) * SCRATCH1YDATA[iBf] + 
                                 GradRhoB[thread_id](2) * SCRATCH1ZDATA[iBf];  
           } // iBf
        } // iShell
      } //UKS
    } //GGA
    result[thread_id].Energy += pt.weight * kernelXC.eps;
    kernelXC.ddrhoA *= pt.weight;
    kernelXC.ddrhoB *= pt.weight;

    for(auto iShell : closeShells[thread_id]) {
      int iSz= shSizes[iShell];
      int iSt = this->basisset_->mapSh2Bf(iShell);

      for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
        double Ta = kernelXC.ddrhoA*SCRATCH1DATA[iBf] + OmegaADATA[iBf];
        for(auto jShell : closeShells[thread_id]) {
          int jSz= shSizes[jShell];
          int jSt = this->basisset_->mapSh2Bf(jShell);
          T *VXCDATA = result[thread_id].VXCScalar.data() + iBf*this->nBasis_;
          for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
            if(jBf < iBf) continue;
            VXCDATA[jBf] += Ta*SCRATCH1DATA[jBf] + SCRATCH1DATA[iBf]*OmegaADATA[jBf]; 
          } // jBf
        } // jShell
      } // iBf
    } // iShell

    if(this->nTCS_ == 1 && !this->isClosedShell){
      for(auto iShell : closeShells[thread_id]) {
        int iSz= shSizes[iShell];
        int iSt = this->basisset_->mapSh2Bf(iShell);
  
        for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
          double Tb = kernelXC.ddrhoB*SCRATCH1DATA[iBf] + OmegaBDATA[iBf];
          for(auto jShell : closeShells[thread_id]) {
            int jSz= shSizes[jShell];
            int jSt = this->basisset_->mapSh2Bf(jShell);
            T *VXCDATB = result[thread_id].VXCMz.data() + iBf*this->nBasis_;
            for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
              if(jBf < iBf) continue;
              VXCDATB[jBf] += Tb*SCRATCH1DATA[jBf] + SCRATCH1DATA[iBf]*OmegaBDATA[jBf]; 
            } // jBf
          } // jShell
        } // iBf
      } // iShell
    } // UKS

    if(doTimings){
      Newend = std::chrono::high_resolution_clock::now();
      T5 += Newend - Newstart;
    }
    tracker++;
  };

//  ChronusQ::AtomicGrid AGrid(100,302,ChronusQ::GRID_TYPE::EULERMAC,
//      ChronusQ::GRID_TYPE::LEBEDEV,ChronusQ::ATOMIC_PARTITION::BECKE,
//      this->molecule_->cartArray(),this->molecule_->rIJ(),0,1.0,false);


//if(this->isGGA) cout << "GGA ON " << this->isGGA <<endl ; 
//if(!this->isGGA) cout << "GGA OFF " << this->isGGA <<endl ; 

/*
  ChronusQ::AtomicGrid AGrid(this->nRadDFTGridPts_,this->nAngDFTGridPts_,
      ChronusQ::GRID_TYPE::GAUSSCHEBFST,ChronusQ::GRID_TYPE::LEBEDEV,
      ChronusQ::ATOMIC_PARTITION::BECKE,this->molecule_->cartArray(),
      this->molecule_->rIJ(),0,this->epsScreen,1e6,1.0,false);

*/   
  ChronusQ::AtomicGrid AGrid(this->nRadDFTGridPts_,this->nAngDFTGridPts_,
      static_cast<ChronusQ::GRID_TYPE>(this->dftGrid_),ChronusQ::GRID_TYPE::LEBEDEV,
      static_cast<ChronusQ::ATOMIC_PARTITION>(this->weightScheme_),this->molecule_->cartArray(),
      this->molecule_->rIJ(),0,this->epsScreen,1e6,1.0,false);
//NOMAG  KernelIntegrand<double> res(this->vXA_->cols());
  //Screaning based on cutoff
//  the radCutoof vector is populated in singleSlater_real_guess.cpp

  this->energyExc    = 0.0;
  std::vector<KernelIntegrand<T>> res(nthreads,
    KernelIntegrand<T>(this->vXCScalar_->cols()));

  for(auto i = 0; i < nthreads; i++) {
    res[i].VXCScalar.setZero();
    res[i].Energy = 0;
    if(!this->isClosedShell or this->nTCS_ == 2)
      res[i].VXCMz.setZero();
  }

  std::vector<double> atomRadCutoff(this->molecule_->nAtoms(),0.0);

  if (this->screenVxc) {
    for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
      for(auto iSh = 0; iSh < this->basisset_->nShell(); iSh++){
        if(this->basisset_->mapSh2Cen(iSh)-1 != iAtm) continue;
        if(this->basisset_->radCutSh()[iSh] > atomRadCutoff[iAtm])
          atomRadCutoff[iAtm] = this->basisset_->radCutSh()[iSh];
//A        atomRadCutoff[iAtm] = 1.e100;
      }
    } 
  }


  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    AGrid.center() = iAtm;
//    cout << "Cutoff Rad = " << atomRadCutoff[iAtm] << "for iAt= " << iAtm << endl;
    if (this->screenVxc) {
      AGrid.setRadCutOff(atomRadCutoff[iAtm]);
      AGrid.findNearestNeighbor();
    } else {
      AGrid.SwitchScreen();
    }
    AGrid.scalingFactor()=0.5 *
      elements[this->molecule_->index(iAtm)].sradius/phys.bohr;

    std::function<void(IntegrationPoint&,std::vector<KernelIntegrand<T>>&)>
      wrapper = std::bind(valVxc,iAtm,std::placeholders::_1,
                  std::placeholders::_2);

    AGrid.integrate<std::vector<KernelIntegrand<T>>>(wrapper,res);
  };


  this->vXCScalar_->setZero();   // Set to zero every occurence of the SCF
  if(!this->isClosedShell or this->nTCS_ == 2)
    this->vXCMz_->setZero();   // Set to zero every occurence of the SCF

  this->energyExc = 0.;
  
  for(auto ithread = 0; ithread < nthreads; ithread++){
    if(this->isClosedShell && this->nTCS_ != 2){
      (*this->vXCScalar_) += 8.0*math.pi*res[ithread].VXCScalar;
    } else if(!this->isClosedShell or this->nTCS_ == 2){
      (*this->vXCScalar_) += 4.0*math.pi*(res[ithread].VXCScalar+res[ithread].VXCMz);  //ALPHA LIKE
      (*this->vXCMz_) += 4.0*math.pi*(res[ithread].VXCScalar-res[ithread].VXCMz);      //BETA LIKE
    }
    this->energyExc += 4*math.pi*res[ithread].Energy;
  }
  (*this->vXCScalar_) = this->vXCScalar_->template selfadjointView<Lower>();
  if(!this->isClosedShell or this->nTCS_ == 2)
    (*this->vXCMz_) = this->vXCMz_->template selfadjointView<Lower>();

  // Cleanup Scratch
  for(auto ithread = 0; ithread < nthreads; ithread++) {
    this->memManager_->free(SCRATCH1[ithread].data(),this->nBasis_);
    if(this->isGGA) {
      this->memManager_->free(SCRATCH1X[ithread].data(),this->nBasis_);
      this->memManager_->free(SCRATCH1Y[ithread].data(),this->nBasis_);
      this->memManager_->free(SCRATCH1Z[ithread].data(),this->nBasis_);
    }

    this->memManager_->free(OmegaA[ithread].data(),this->nBasis_);
    this->memManager_->free(OmegaB[ithread].data(),this->nBasis_);
  }

  if(doTimings) {
    cout << "T1 = " << T1.count() << endl;
    cout << "T2 = " << T2.count() << endl;
    cout << "T3 = " << T3.count() << endl;
    cout << "T4 = " << T4.count() << endl;
    cout << "T5 = " << T5.count() << endl;
    cout << "T6 = " << T6.count() << endl;
  }
//cout << "NSkip2 = " << NSkip2 << endl;
//cout << "NSkip3 = " << NSkip3 << endl;
//cout << "NSkip4 = " << NSkip4 << endl;
//cout << "NSkip5 = " << NSkip5 << endl;
//if(doTimings)
//  for(auto i : TF) cout << "TF " << i.count() << endl;

  if(this->printLevel_ >= 3) {
    finish = std::chrono::high_resolution_clock::now();
    duration_formVxc = finish - start;
    prettyPrintSmart(this->fileio_->out,(*this->vXCScalar()),"LDA Vxc Scalar");
    if(!this->isClosedShell && this->nTCS_ != 2){
      prettyPrintSmart(this->fileio_->out,(*this->vXCMz()),"LDA Vxc Mz");
    }
    this->fileio_->out << "VXC Energy= " <<  this->energyExc << endl, 
    this->fileio_->out << endl << "CPU time for VXC integral:  "
                       << duration_formVxc.count() << " seconds." 
                       << endl;
  }
}

