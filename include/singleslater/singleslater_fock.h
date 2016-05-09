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

  KernelIntegrand(size_t N) : VXCA(N,N), Energy(0.0){ 
    VXCA.setZero();
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

    bool testNew = false;
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

        auto valVxc = [&](ChronusQ::IntegrationPoint pt, 
            KernelIntegrand<T> &result) {

          SCRATCH1.setZero();
          cartGP GP = pt.pt;
          double rhoA;
          double rhoB;
          auto shMap = this->basisset_->MapGridBasis(GP); 
          if(shMap[0]) {return 0.0;}
          for(auto iShell = 0; iShell < this->basisset_->nShell(); iShell++){
            if(!shMap[iShell+1]) {continue;}

            int b_s = this->basisset_->mapSh2Bf(iShell);
            int size= this->basisset_->shells(iShell).size();

            libint2::Shell shTmp = this->basisset_->shells(iShell);
            double * buff = this->basisset_->basisDEval(0,shTmp,&pt.pt);
            RealMap bMap(buff,size,1);
            SCRATCH1.block(b_s,0,size,1) = bMap;

            delete [] buff;
          };

          if(SCRATCH1.norm() < 1e-8) return 0.0;
          SCRATCH2 = SCRATCH1 * SCRATCH1.transpose();
          double rhoT = this->template computeProperty<double,TOTAL>(SCRATCH2);
          double rhoS = this->template computeProperty<double,MZ>(SCRATCH2);
          rhoA = 0.5 * (rhoT + rhoS);
          rhoB = 0.5 * (rhoT - rhoS);

          for(auto i = 0; i < this->dftFunctionals_.size(); i++){
            DFTFunctional::DFTInfo kernelXC = 
              this->dftFunctionals_[i]->eval(rhoA, rhoB);
            result.VXCA.real()   += pt.weight * SCRATCH2 * kernelXC.ddrhoA; 
            result.Energy += pt.weight * (rhoA+rhoB) * kernelXC.eps;
          }
        };

        ChronusQ::AtomicGrid AGrid(100,302,ChronusQ::GRID_TYPE::GAUSSCHEBFST,
            ChronusQ::GRID_TYPE::LEBEDEV,ChronusQ::ATOMIC_PARTITION::BECKE,
            this->molecule_->cartArray(),0,1.0,false);
       
        KernelIntegrand<T> res(this->vXA_->cols());
        this->basisset_->radcut(1.0e-10, 50, 1.0e-7);
        for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
          AGrid.center() = iAtm;
          AGrid.scalingFactor()=0.5 *
            elements[this->molecule_->index(iAtm)].sradius/phys.bohr;
          AGrid.integrate<KernelIntegrand<T>>(valVxc,res);
        };
        (*this->vXA_) = 4*math.pi*res.VXCA;
        this->totalEx = 4*math.pi*res.Energy;
        if(this->printLevel_ >= 3) {
          finish = std::chrono::high_resolution_clock::now();
          duration_formVxc = finish - start;
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
        (*this->fockScalar_) += (*this->vXA_) + (*this->vCorA_);
        (*this->fockScalar_) += (*this->vXB_) + (*this->vCorB_);
        (*this->fockMz_)     += (*this->vXA_) + (*this->vCorA_);
        (*this->fockMz_)     -= (*this->vXB_) + (*this->vCorB_);
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

