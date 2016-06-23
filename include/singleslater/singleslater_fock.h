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
template<typename T>
class KernelIntegrand {
  typedef Eigen::Matrix<T,Dynamic,Dynamic> TMatrix;
  public:
  TMatrix VXCA;
  TMatrix VXCB;
  double Energy;

  KernelIntegrand(size_t N) : VXCA(N,N), VXCB(N,N), Energy(0.0){ 
    VXCA.setZero();
    VXCB.setZero();
  };
};

/********************************
 * Form Perturbation Tensor (G) *
 ********************************/
// Only available for Libint Integrals
#ifdef USE_LIBINT
template<typename T>
void SingleSlater<T>::formPT(){
  bool doTCS = (this->Ref_ == TCS);
  bool doRHF = (this->isClosedShell && !doTCS);
  bool doKS  = this->isDFT;
  if(!this->haveDensity) this->formDensity();
//  this->sepReImOnePDM();
//  this->comReImOnePDM();
  this->scatterDensity();

  std::vector<std::reference_wrapper<TMap>> mats;
  std::vector<std::reference_wrapper<TMap>> ax;
  std::vector<AOIntegrals::ERI_CONTRACTION_TYPE> contList;
  std::vector<double> scalingFactors;
  double exchFactor = -0.5;
  if(this->isDFT) exchFactor = 0.0;

  if(this->nTCS_ == 1 && this->isClosedShell) {

    mats.emplace_back(*this->onePDMA_);
    mats.emplace_back(*this->onePDMA_);
    ax.emplace_back(*this->PTA_);
    ax.emplace_back(*this->PTA_);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::COULOMB);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    scalingFactors.push_back(1.0);
    scalingFactors.push_back(exchFactor);

    this->PTA_->setZero();

  } else {
    mats.emplace_back(*this->onePDMScalar_);
    mats.emplace_back(*this->onePDMScalar_);
    mats.emplace_back(*this->onePDMMz_);
    ax.emplace_back(*this->PTScalar_);
    ax.emplace_back(*this->PTScalar_);
    ax.emplace_back(*this->PTMz_);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::COULOMB);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
    scalingFactors.push_back(1.0);
    scalingFactors.push_back(exchFactor);
    scalingFactors.push_back(exchFactor);
    this->PTScalar_->setZero();
    this->PTMz_->setZero();
    if(this->nTCS_ == 2){
      mats.emplace_back(*this->onePDMMy_);
      mats.emplace_back(*this->onePDMMx_);
      ax.emplace_back(*this->PTMy_);
      ax.emplace_back(*this->PTMx_);
      contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
      contList.push_back(AOIntegrals::ERI_CONTRACTION_TYPE::EXCHANGE);
      scalingFactors.push_back(exchFactor);
      scalingFactors.push_back(exchFactor);
      this->PTMy_->setZero();
      this->PTMz_->setZero();
    }
  }

  this->aointegrals_->newTwoEContract(mats,ax,contList,scalingFactors);

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

  if(this->printLevel_ >= 3 && getRank() == 0) this->printPT();
}
#endif

/********************
 * Form Fock Matrix *
 ********************/
template<typename T>
void SingleSlater<T>::formFock(){
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes before fock build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  if(!this->haveDensity) this->formDensity();
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

    bool testNew = true;
    bool isGGA   = true;
    if (this->isDFT){
//    Timing
      std::chrono::high_resolution_clock::time_point start;
      std::chrono::high_resolution_clock::time_point finish;
      std::chrono::duration<double> duration_formVxc;
     
      if(this->printLevel_ >= 3) {
        start = std::chrono::high_resolution_clock::now();
      }

//    this->formVXC();
      if(!testNew) this->formVXC_store();
      else {
        RealMatrix SCRATCH2(this->nBasis_,this->nBasis_);
        VectorXd   SCRATCH1(this->nBasis_);
        RealMatrix SCRATCH2X(this->nBasis_,this->nBasis_);
        RealMatrix SCRATCH2Y(this->nBasis_,this->nBasis_);
        RealMatrix SCRATCH2Z(this->nBasis_,this->nBasis_);
        VectorXd   SCRATCH1X(this->nBasis_);
        VectorXd   SCRATCH1Y(this->nBasis_);
        VectorXd   SCRATCH1Z(this->nBasis_);
        std::chrono::duration<double> T1(0.0);
        std::chrono::duration<double> T2(0.0);
        std::chrono::duration<double> T3(0.0);
        std::chrono::duration<double> T4(0.0);
        int NDer = 0;
        if(isGGA) NDer = 1; 
        auto valVxc = [&](ChronusQ::IntegrationPoint pt, 
        KernelIntegrand<T> &result) {

          auto Newstart = std::chrono::high_resolution_clock::now();
          SCRATCH1.setZero();
          SCRATCH1X.setZero();
          SCRATCH1Y.setZero();
          SCRATCH1Z.setZero();
          auto Newend = std::chrono::high_resolution_clock::now();
          T1 += Newend - Newstart;

          std::array<double,3>  drhoT = {0.0,0.0,0.0}; ///< array TOTAL density gradient components
          std::array<double,3>  drhoS = {0.0,0.0,0.0}; ///< array SPIN  density gradient components
          std::array<double,3>  drhoA = {0.0,0.0,0.0}; ///< array ALPHA  density gradient components
          std::array<double,3>  drhoB = {0.0,0.0,0.0}; ///< array BETA  density gradient components
          RealVecMap GradRhoT(&drhoT[0],3);
          RealVecMap GradRhoS(&drhoS[0],3);
          RealVecMap GradRhoA(&drhoA[0],3);
          RealVecMap GradRhoB(&drhoB[0],3);
          cartGP GP = pt.pt;
          double rhoA;
          double rhoB;
          double gammaAA;
          double gammaBB;
          double gammaAB;
          auto shMap = this->basisset_->MapGridBasis(GP); 
          if(shMap[0]) {return 0.0;}
          Newstart = std::chrono::high_resolution_clock::now();
          for(auto iShell = 0; iShell < this->basisset_->nShell(); iShell++){
            if(!shMap[iShell+1]) {continue;}

            int b_s = this->basisset_->mapSh2Bf(iShell);
            int shSize= this->basisset_->shells(iShell).size();

            libint2::Shell shTmp = this->basisset_->shells(iShell);
            double * buff = this->basisset_->basisDEval(NDer,shTmp,&pt.pt);
            RealMap bMap(buff,shSize,1);
            SCRATCH1.block(b_s,0,shSize,1) = bMap;
            if(NDer>0){
              double * ds1EvalX = buff + shSize;
              double * ds1EvalY = ds1EvalX + shSize;
              double * ds1EvalZ = ds1EvalY + shSize;
              RealMap bMapX(ds1EvalX,shSize,1);
              SCRATCH1X.block(b_s,0,shSize,1) = bMapX;
              RealMap bMapY(ds1EvalY,shSize,1);
              SCRATCH1Y.block(b_s,0,shSize,1) = bMapY;
              RealMap bMapZ(ds1EvalZ,shSize,1);
              SCRATCH1Z.block(b_s,0,shSize,1) = bMapZ;
            }

            delete [] buff;
          };
          Newend = std::chrono::high_resolution_clock::now();
          T2 += Newend - Newstart;

//          if(SCRATCH1.norm() < 1e-8) return 0.0;
          Newstart = std::chrono::high_resolution_clock::now();
          SCRATCH2 = SCRATCH1 * SCRATCH1.transpose();
          double rhoT = this->template computeProperty<double,TOTAL>(SCRATCH2);
          double rhoS = this->template computeProperty<double,MZ>(SCRATCH2);
          Newend = std::chrono::high_resolution_clock::now();
          T3 += Newend - Newstart;

          rhoA = 0.5 * (rhoT + rhoS);
          rhoB = 0.5 * (rhoT - rhoS);

/*
//  Handle numerical instability if screening on
          if (this->screenVxc ) {
//    check if are noise
      if(rhoT   <= 0.0 ) {
        if((std::abs(rhoT)) <= 1.0e-10) {
          return;
        }else{ 
          CErr("Numerical noise in the density");
        }
//    skyp points based on small density
      }else if(rhoT < 1.0e-12){
        return;
      }
    }
*/
          if(NDer>0){
            //Closed Shell GGA
            SCRATCH2X = SCRATCH1 * SCRATCH1X.transpose();
            SCRATCH2Y = SCRATCH1 * SCRATCH1Y.transpose();
            SCRATCH2Z = SCRATCH1 * SCRATCH1Z.transpose();
            drhoT[0] = 2.0*this->template computeProperty<double,TOTAL>(SCRATCH2X); 
            drhoT[1] = 2.0*this->template computeProperty<double,TOTAL>(SCRATCH2Y); 
            drhoT[2] = 2.0*this->template computeProperty<double,TOTAL>(SCRATCH2Z); 
            drhoS[0] = 2.0*this->template computeProperty<double,MZ>(SCRATCH2X); 
            drhoS[1] = 2.0*this->template computeProperty<double,MZ>(SCRATCH2Y); 
            drhoS[2] = 2.0*this->template computeProperty<double,MZ>(SCRATCH2Z); 
            GradRhoA = 0.5 * (GradRhoT + GradRhoS);
            GradRhoB = 0.5 * (GradRhoT - GradRhoS);
            gammaAA = GradRhoA.dot(GradRhoA);           
            gammaBB = GradRhoB.dot(GradRhoB);           
            gammaAB = GradRhoA.dot(GradRhoB);           
            SCRATCH2X += SCRATCH1X * SCRATCH1.transpose();
            SCRATCH2Y += SCRATCH1Y * SCRATCH1.transpose();
            SCRATCH2Z += SCRATCH1Z * SCRATCH1.transpose();
          }
          for(auto i = 0; i < this->dftFunctionals_.size(); i++){
            DFTFunctional::DFTInfo kernelXC;
            if  (rhoT < 1.0e-10) continue;
            if (NDer > 0) {
              Newstart = std::chrono::high_resolution_clock::now();
              kernelXC = this->dftFunctionals_[i]->eval(
                  rhoA,rhoB,gammaAA,gammaAB,gammaBB);
              Newend = std::chrono::high_resolution_clock::now();
              T4 += Newend - Newstart;

              result.VXCA.real() += pt.weight *SCRATCH2X  
                * (  2.0 * GradRhoA[0]*kernelXC.ddgammaAA 
                   + GradRhoB[0]* kernelXC.ddgammaAB);  
              result.VXCA.real() += pt.weight *SCRATCH2Y  
                * (  2.0 * GradRhoA[1]*kernelXC.ddgammaAA 
                   + GradRhoB[1]* kernelXC.ddgammaAB);  
              result.VXCA.real() += pt.weight *SCRATCH2Z  
                * (  2.0 * GradRhoA[2]*kernelXC.ddgammaAA 
                   + GradRhoB[2]* kernelXC.ddgammaAB);  
              result.Energy += pt.weight * kernelXC.eps;
	      if(!this->isClosedShell && this->nTCS_ != 2) {
                result.VXCB.real() += pt.weight *SCRATCH2X  
                  * (  2.0 * GradRhoB[0]*kernelXC.ddgammaBB 
                     + GradRhoA[0]* kernelXC.ddgammaAB);  
                result.VXCB.real() += pt.weight *SCRATCH2Y  
                  * (  2.0 * GradRhoB[1]*kernelXC.ddgammaBB 
                     + GradRhoA[1]* kernelXC.ddgammaAB);  
                result.VXCB.real() += pt.weight *SCRATCH2Z  
                  * (  2.0 * GradRhoB[2]*kernelXC.ddgammaBB 
                     + GradRhoA[2]* kernelXC.ddgammaAB);  
	      }
            } else {
                kernelXC = 
               this->dftFunctionals_[i]->eval(rhoA, rhoB);
            result.Energy += pt.weight * (rhoA+rhoB) * kernelXC.eps;
              }
            result.VXCA.real()   += pt.weight * SCRATCH2 * kernelXC.ddrhoA; 
	    if(!this->isClosedShell && this->nTCS_ != 2) 
              result.VXCB.real()   += pt.weight * SCRATCH2 * kernelXC.ddrhoB; 
          }
        };

        ChronusQ::AtomicGrid AGrid(100,302,ChronusQ::GRID_TYPE::GAUSSCHEBFST,
            ChronusQ::GRID_TYPE::LEBEDEV,ChronusQ::ATOMIC_PARTITION::BECKE,
            this->molecule_->cartArray(),0,1.0,false);
         
        KernelIntegrand<T> res(this->vXA_->cols());
        this->basisset_->radcut(1.0e-10, 50, 1.0e-7);
        this->totalEx    = 0.0;
        this->vXA()->setZero();   // Set to zero every occurence of the SCF
        for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
          AGrid.center() = iAtm;
          AGrid.scalingFactor()=0.5 *
            elements[this->molecule_->index(iAtm)].sradius/phys.bohr;
          AGrid.integrate<KernelIntegrand<T>>(valVxc,res);
        };
        (*this->vXA_) = 4*math.pi*res.VXCA;
        this->totalEx = 4*math.pi*res.Energy;
        if(!this->isClosedShell && this->nTCS_ != 2){
          (*this->vXB_) = 4*math.pi*res.VXCB;
	}
        cout << "T1 = " << T1.count() << endl;
        cout << "T2 = " << T2.count() << endl;
        cout << "T3 = " << T3.count() << endl;
        cout << "T4 = " << T4.count() << endl;
        if(this->printLevel_ >= 3) {
          finish = std::chrono::high_resolution_clock::now();
          duration_formVxc = finish - start;
          prettyPrint(this->fileio_->out,(*this->vXA()),"LDA Vxc alpha");
          if(!this->isClosedShell && this->nTCS_ != 2){
            prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vxc beta");
          }
          this->fileio_->out << "VXC Energy= " <<  this->totalEx << endl, 
          this->fileio_->out << endl << "CPU time for VXC integral:  "
                             << duration_formVxc.count() << " seconds." 
                             << endl;
        }
      }
    }
  
    if(this->nTCS_ == 1 && this->isClosedShell) {
      this->fockA_->setZero();
      this->fockA_->real() += (*this->aointegrals_->oneE_);
      this->aointegrals_->addElecDipole(*this->fockA_,this->elecField_);
      (*this->fockA_) += (*this->PTA_);
      if(this->isDFT){ 
        (*this->fockA_) += (*this->vXA_);
        if(!testNew) (*this->fockA_) += (*this->vCorA_);
      }
    } else {

    //this->fockA_->setZero();
    //this->fockA_->real() += (*this->aointegrals_->oneE_);
    //this->aointegrals_->addElecDipole(*this->fockA_,this->elecField_);
    //(*this->fockA_) += (*this->PTA_);
    //if(this->isDFT){ 
    //  (*this->fockA_) += (*this->vXA_);
    //  (*this->fockA_) += (*this->vCorA_);
    //}
    //this->fockB_->setZero();
    //this->fockB_->real() += (*this->aointegrals_->oneE_);
    //this->aointegrals_->addElecDipole(*this->fockB_,this->elecField_);
    //(*this->fockB_) += (*this->PTB_);
    //if(this->isDFT){ 
    //  (*this->fockB_) += (*this->vXB_);
    //  (*this->fockB_) += (*this->vCorB_);
    //}

      this->fockScalar_->setZero();
      this->fockMz_->setZero();
      this->fockScalar_->real() += (*this->aointegrals_->oneE_);
      this->aointegrals_->addElecDipole(*this->fockScalar_,this->elecField_);
      (*this->fockScalar_) *= 2.0;
      (*this->fockScalar_)      += (*this->PTScalar_);        
      (*this->fockMz_)          += (*this->PTMz_);

      // FIXME: Needs to ge generalized for the 2C case
      if(this->nTCS_ == 1 && this->isDFT){
        if(!testNew) {
          (*this->fockScalar_) += (*this->vXA_) + (*this->vCorA_);
          (*this->fockScalar_) += (*this->vXB_) + (*this->vCorB_);
          (*this->fockMz_)     += (*this->vXA_) + (*this->vCorA_);
          (*this->fockMz_)     -= (*this->vXB_) + (*this->vCorB_);
	} else {
          (*this->fockScalar_) += (*this->vXA_);
          (*this->fockScalar_) += (*this->vXB_);
          (*this->fockMz_)     += (*this->vXA_);
          (*this->fockMz_)     -= (*this->vXB_);
	}
      }

      std::vector<std::reference_wrapper<TMap>> toGather;
      toGather.emplace_back(*this->fockScalar_);
      toGather.emplace_back(*this->fockMz_);
      if(this->nTCS_ == 1)
        Quantum<T>::spinGather(*this->fockA_,*this->fockB_,toGather);
      else {
        this->fockMx_->setZero();
        this->fockMy_->setZero();
        (*this->fockMx_) += (*this->PTMx_);
        (*this->fockMy_) += (*this->PTMy_);
        toGather.emplace_back(*this->fockMy_);
        toGather.emplace_back(*this->fockMx_);
        Quantum<T>::spinGather(*this->fockA_,toGather);
      };


    }

    /*
    // Add in the electric field component if they are non-zero
    std::array<double,3> null{{0,0,0}};
    if(this->elecField_ != null){
      //this->fileio_->out << "Adding in Electric Field Contribution" << endl;
      int NB = this->nTCS_*this->nBasis_;
      int NBSq = NB*NB;
      int iBuf = 0;
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        ConstRealMap mu(&this->aointegrals_->elecDipole_->storage()[iBuf],
          NB,NB);
        fockA_->real() += this->elecField_[iXYZ] * mu;
        if(!this->isClosedShell && this->Ref_ != TCS) 
          fockB_->real() += this->elecField_[iXYZ] * mu;
        iBuf += NBSq;
      }
    }
    */
    if(this->printLevel_ >= 2) this->printFock(); 
  }
#ifdef CQ_ENABLE_MPI
  // Syncronize MPI processes after fock build
  MPI_Barrier(MPI_COMM_WORLD);
#endif
};

