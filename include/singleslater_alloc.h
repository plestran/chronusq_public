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
/*********************
 * Allocate Matricies *
 **********************/
template<typename T>
void SingleSlater<T>::iniSingleSlater(Molecule * molecule, BasisSet * basisset, 
                                   AOIntegrals * aointegrals, FileIO * fileio, 
                                   Controls * controls) {

  this->communicate(*molecule,*basisset,*aointegrals,*fileio,*controls);
  this->initMeta();

  this->elecField_  = controls->field_;
  this->printLevel_ = controls->printLevel;
  this->guess_      = controls->guess;

  this->isClosedShell = (this->multip_ == 1);
  this->isDFT         = controls->DFT;
  this->ExchKernel_   = SLATER;
  this->CorrKernel_   = VWN5;
  if(controls->HF){
    if(this->isClosedShell && !controls->doCUHF
       && !controls->doTCS)                        this->Ref_ = RHF ; // RHF
    else if(!controls->doCUHF && !controls->doTCS) this->Ref_ = UHF ; // UHF
    else if(controls->doCUHF)                      this->Ref_ = CUHF; // CUHF
    else if(controls->doTCS)                       this->Ref_ = TCS ; // TCS
  } else if(this->isDFT) {
    if(this->isClosedShell && !controls->doCUHF
       && !controls->doTCS)                        this->Ref_ = RKS ; // RKS
    else if(!controls->doCUHF && !controls->doTCS) this->Ref_ = UKS ; // UKs
    else if(controls->doCUHF)                      this->Ref_ = CUKS; // CUKS
    else if(controls->doTCS)                       this->Ref_ = GKS ; // GKS
  }

  this->genMethString();


  if(this->Ref_ == TCS) this->nTCS_ = 2;

  // This is the only way via the C++ interface to set this flag (needed
  // for allocDFT)
  this->isDFT = controls->DFT;
  this->alloc();


};

template<typename T>
void SingleSlater<T>::alloc(){
  this->checkMeta();
  this->allocOp();
//if(this->maxMultipole_ > 0) this->allocMultipole(); 
 
  if(getRank() == 0) {
    if (this->isDFT){
//     Timing
      std::chrono::high_resolution_clock::time_point start;
      std::chrono::high_resolution_clock::time_point finish;
      std::chrono::duration<double> duration_formMap;
      if(this->printLevel_ >= 3) {
        start = std::chrono::high_resolution_clock::now();
      }
      this->genSparseBasisMap();
//      this->genSparseRcrosP();
//      CErr();
//   Timing
      if(this->printLevel_ >= 3) {
        finish = std::chrono::high_resolution_clock::now();
        duration_formMap = finish - start;
        this->fileio_->out << endl << "CPU time for MapVXc:  "
          << duration_formMap.count() << " seconds." << endl;
      }
    }
 
    if(this->isPrimary) {
      if(typeid(T).hash_code() == typeid(double).hash_code())
        this->fileio_->iniStdSCFFilesDouble(
          !this->isClosedShell && this->Ref_ != TCS,this->nTCS_*this->nBasis_
        );
      else if(typeid(T).hash_code() == typeid(dcomplex).hash_code())
        this->fileio_->iniStdSCFFilesComplex(
          !this->isClosedShell && this->Ref_ != TCS,this->nTCS_*this->nBasis_
        );
    }
  }
#ifdef CQ_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
/* Leaks memory
  int i,j,ij;
  this->R2Index_ = new int*[nBasis];
  for(i=0;i<nBasis;i++) this->R2Index_[i] = new int[nBasis];
  for(i=0;i<nBasis;i++) for(j=0;j<nBasis;j++) {
    if(i>=j) ij=j*(nBasis)-j*(j-1)/2+i-j;
    else ij=i*(nBasis)-i*(i-1)/2+j-i;
    this->R2Index_[i][j] = ij;
  };
*/
}

template<typename T>
void SingleSlater<T>::allocOp(){
  this->allocDensity(this->nTCS_*this->nBasis_); // Allocate Den -> Quantum
  this->allocAlphaOp();
  if(!this->isClosedShell && this->Ref_ != TCS) 
    this->allocBetaOp();
}

template<typename T>
void SingleSlater<T>::allocAlphaOp(){
  // Alpha / TCS Density Matrix
//try { 
//  this->densityA_  = std::unique_ptr<TMatrix>( 
//    new TMatrix(this->nTCS_*this->nBasis_, this->nTCS_*this->nBasis_));
//} catch (...) { 
//  if(this->Ref_ == TCS) 
//    CErr(std::current_exception(),"TCS Density Matrix Allocation"  ); 
//  else CErr(std::current_exception(),"Alpha Density Matrix Allocation"); 
//}
  if(getRank() != 0) return;
  // Alpha / TCS Fock Matrix
  try { 
    this->fockA_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) 
      CErr(std::current_exception(),"TCS Fock Matrix Allocation"); 
    else CErr(std::current_exception(),"Alpha Fock Matrix Allocation"); 
  }

  // Alpha / TCS Molecular Orbital Coefficients
  try { 
    this->moA_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  } catch (...) { 
    if(this->Ref_ == TCS) 
      CErr(std::current_exception(),"TCS MO Coefficients Allocation");
    else CErr(std::current_exception(),"Alpha MO Coefficients Allocation"); 
  }

  // Alpha / TCS Eigenorbital Energies
  try { 
    this->epsA_ = std::unique_ptr<RealMatrix>(
      new RealMatrix(this->nTCS_*this->nBasis_,1)); 
  } catch (...) { 
    if(this->Ref_ == TCS) 
      CErr(std::current_exception(),"TCS Eigenorbital Energies"); 
    else CErr(std::current_exception(),"Alpha Eigenorbital Energies"); 
  }

#ifndef USE_LIBINT
  // Alpha / TCS Coulomb Matrix
  try { 
    this->coulombA_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  } catch (...) { 
    if(this->Ref_ == TCS) 
      CErr(std::current_exception(),"TCS Coulomb Tensor Allocation"); 
    else CErr(std::current_exception(),"Alpha Coulomb Tensor Allocation"); 
  }

  // Alpha / TCS Exchange Matrix
  try { 
    this->exchangeA_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) 
      CErr(std::current_exception(),"TCS Exchange Tensor Allocation"); 
    else CErr(std::current_exception(),"Alpha Exchange Tensor Allocation"); 
  }

#else
  // Alpha / TCS Perturbation Tensor
  try { 
    this->PTA_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(),"TCS G[P] Allocation"); 
    else CErr(std::current_exception(),"Alpha G[P] Allocation"); 
  }

#endif

  if(this->isDFT) this->allocAlphaDFT();
}

template<typename T>
void SingleSlater<T>::allocBetaOp(){
  // Beta Density Matrix
//try { 
//  this->densityB_ = std::unique_ptr<TMatrix>(
//    new TMatrix(this->nBasis_,this->nBasis_)); 
//} catch (...) { 
//  CErr(std::current_exception(),"Beta Density Matrix Allocation"); 
//}

  if(getRank() != 0) return;
  // Beta Fock Matrix
  try { 
    this->fockB_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nBasis_,this->nBasis_)); 
  } catch (...) { 
    CErr(std::current_exception(),"Beta Fock Matrix Allocation");
  }

  // Beta Molecular Orbital Coefficients
  try { 
    this->moB_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nBasis_,this->nBasis_));
  } catch (...) { 
    CErr(std::current_exception(),"Beta MO Coefficients Allocation"); 
  }

  // Beta Eigenorbital Energies
  try { 
    this->epsB_ = std::unique_ptr<RealMatrix>(
      new RealMatrix(this->nBasis_,1)); 
  } catch (...) { 
    CErr(std::current_exception(),"Beta Eigenorbital Energies");
  }
#ifndef USE_LIBINT
  // Beta Coulomb Matrix
  try { 
    this->coulombB_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nBasis_,this->nBasis_)); 
  } catch (...) { 
    CErr(std::current_exception(),"Beta Coulomb Tensor Allocation"); 
  }
 
  // Beta Exchange Matrix
  try { 
    this->exchangeB_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nBasis_,this->nBasis_));
  } catch (...) { 
    CErr(std::current_exception(),"Beta Exchange Tensor Allocation"); 
  }
#else
  // Beta Perturbation Tensor
  try { 
    this->PTB_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nBasis_,this->nBasis_));
  } catch (...) { 
    CErr(std::current_exception(),"Beta G[P] Allocation"); 
  }
#endif

  if(this->isDFT) this->allocBetaDFT();

}

template<typename T>
void SingleSlater<T>::allocDFT(){
  this->allocAlphaDFT();
  if(!this->isClosedShell && this->Ref_ != TCS) 
    this->allocBetaDFT();
}

template<typename T>
void SingleSlater<T>::allocAlphaDFT(){
  // Alpha / TCS VXC
  try { 
    this->vXA_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
    this->vCorA_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(), "TCS VXC Allocation"); 
    else CErr(std::current_exception(),"Alpha VXC  Allocation"); 
  }
}

template<typename T>
void SingleSlater<T>::allocBetaDFT(){
  // Alpha / TCS VXC
  try { 
    this->vXB_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
    this->vCorB_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    CErr(std::current_exception(),"Beta VXC  Allocation"); 
  }
}

/*
template<typename T>
void SingleSlater<T>::allocMultipole(){
  if(this->maxMultipole_ >= 1)
    this->dipole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,1));
  if(this->maxMultipole_ >= 2){
    this->quadpole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,3));
    this->tracelessQuadpole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,3));
  }
  if(this->maxMultipole_ >= 3)
    this->octpole_  = std::unique_ptr<RealTensor3d>(new RealTensor3d(3,3,3));
}
*/
