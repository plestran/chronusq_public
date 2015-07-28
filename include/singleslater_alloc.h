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
/*********************
 * Allocate Matricies *
 **********************/
template<typename T>
void SingleSlater<T>::iniSingleSlater(Molecule * molecule, BasisSet * basisset, 
                                   AOIntegrals * aointegrals, FileIO * fileio, 
                                   Controls * controls) {
  this->molecule_ = molecule;
  this->basisset_ = basisset;
  this->fileio_   = fileio;
  this->controls_ = controls;
  this->aointegrals_= aointegrals;
  int nTotalE = molecule->nTotalE();
  this->nBasis_  = basisset->nBasis();
  this->nTT_   = this->nBasis_*(this->nBasis_+1)/2;
  this->multip_  = molecule->multip();
  this->nShell_ = basisset->nShell();
  int nSingleE = this->multip_ - 1;
  this->nOccB_ = (nTotalE - nSingleE)/2;
  this->nVirB_ = this->nBasis_ - this->nOccB_;
  this->nOccA_ = this->nOccB_ + nSingleE;
  this->nVirA_ = this->nBasis_ - this->nOccA_;
  this->energyNuclei = molecule->energyNuclei();
  this->isConverged = false;
  this->doCUHF = controls->doCUHF;

  this->isClosedShell = (this->multip_ == 1);
  if(this->isClosedShell && !controls->doCUHF)   this->Ref_ = RHF ; // RHF
  else if(!controls->doCUHF && !controls->doTCS) this->Ref_ = UHF ; // UHF
  else if(controls->doCUHF)                      this->Ref_ = CUHF; // CUHF / ROHF
  else if(controls->doTCS)                       this->Ref_ = TCS ; // TCS


  this->nTCS_ = 1;
  if(this->Ref_ == TCS) this->nTCS_ = 2;
  

  // Alpha / TCS Density
  try { 
    this->densityA_  = 
      std::unique_ptr<TMatrix>( new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS)
      CErr(std::current_exception(),"TCS Density Matrix Allocation"); 
    else
      CErr(std::current_exception(),"Alpha Density Matrix Allocation"); 
  }

  // Alpha / TCS Fock
  try { 
    this->fockA_ = 
      std::unique_ptr<TMatrix>(new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS)
      CErr(std::current_exception(),"TCS Fock Matrix Allocation"); 
    else
      CErr(std::current_exception(),"Alpha Fock Matrix Allocation"); 
  }

#ifndef USE_LIBINT
  // Alpha / TCS Coulomb Matrix
  try { 
    this->coulombA_  = 
      std::unique_ptr<TMatrix>(new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  } catch (...) { 
    if(this->Ref_ == TCS)
      CErr(std::current_exception(),"TCS Coulomb Tensor (R2) Allocation"); 
    else
      CErr(std::current_exception(),"Alpha Coulomb Tensor (R2) Allocation"); 
  }

  // Alpha / TCS Exchange Matrix
  try { 
    this->exchangeA_ = 
      std::unique_ptr<TMatrix>(new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS)
      CErr(std::current_exception(),"TCS Exchange Tensor (R2) Allocation"); 
    else
      CErr(std::current_exception(),"Alpha Exchange Tensor (R2) Allocation"); 
  }
#else // USE_LIBINT
  // Alpha / TCS Perturbation Tensor
  try { 
    this->PTA_  = 
      std::unique_ptr<TMatrix>(new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS)
      CErr(std::current_exception(),"TCS Perturbation Tensor (G[P]) Allocation"); 
    else
      CErr(std::current_exception(),"Alpha Perturbation Tensor (G[P]) Allocation"); 
  }
#endif
  // Alpha / TCS Molecular Orbital Coefficients
  try { 
    this->moA_ = 
      std::unique_ptr<TMatrix>(new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  } catch (...) { 
    if(this->Ref_ == TCS)
      CErr(std::current_exception(),"TCS MO Coefficients Allocation");
    else
      CErr(std::current_exception(),"Alpha MO Coefficients Allocation"); 
  }

  // Alpha / TCS Eigenorbital Energies
  try { 
    this->epsA_ = 
      std::unique_ptr<TMatrix>(new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  } catch (...) { 
    if(this->Ref_ == TCS)
      CErr(std::current_exception(),"TCS Eigenorbital Energies"); 
    else
      CErr(std::current_exception(),"Alpha Eigenorbital Energies"); 
  }
  

  if(!this->isClosedShell && this->Ref_ != TCS) {
    try { this->densityB_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Density
    catch (...) { CErr(std::current_exception(),"Beta Density Matrix Allocation"); }
    try { this->fockB_     = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Fock
    catch (...) { CErr(std::current_exception(),"Beta Fock Matrix Allocation"); }
#ifndef USE_LIBINT
    try { this->coulombB_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Coulomb Integral
    catch (...) { CErr(std::current_exception(),"Beta Coulomb Tensor (R2) Allocation"); }
    try { this->exchangeB_ = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Exchange Integral
    catch (...) { CErr(std::current_exception(),"Beta Exchange Tensor (R2) Allocation"); }
#else // USE_LIBINT
    try { this->PTB_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Perturbation Tensor
    catch (...) { CErr(std::current_exception(),"Beta Perturbation Tensor (G[P]) Allocation"); }
#endif
    try { this->moB_       = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Molecular Orbital Coefficients
    catch (...) { CErr(std::current_exception(),"Beta MO Coefficients Allocation"); }
    try { this->epsB_       = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Eigenorbital Energies
    catch (...) { CErr(std::current_exception(),"Beta Eigenorbital Energies"); }
  };

  this->dipole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,1));
  this->quadpole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,3));
  this->tracelessQuadpole_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,3));
  this->octpole_  = std::unique_ptr<RealTensor3d>(new RealTensor3d(3,3,3));
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

  this->haveCoulomb = false;
  this->haveExchange= false;
  this->haveDensity = false;
  this->haveMO	    = false;
  this->havePT      = false;
};

template<typename T>
void SingleSlater<T>::initMemLen(){
  this->lenX_  = this->nBasis_ * this->nBasis_;
  this->lenXp_ = this->nBasis_ * this->nBasis_;
  this->lenF_  = this->nBasis_ * this->nBasis_;
  this->lenP_  = this->nBasis_ * this->nBasis_;
  this->lenCoeff_ = 7;
  this->lenB_     = this->lenCoeff_ * this->lenCoeff_;
  this->LWORK_ = 4*this->nBasis_;
  this->lenLambda_ = this->nBasis_ * this->nBasis_;
  this->lenDelF_   = this->nBasis_ * this->nBasis_;
  this->lenOccNum_ = this->nBasis_;

  this->lenScr_ = 0;

  this->lenScr_ += this->lenX_;     // Storage for S^(-0.5)
  this->lenScr_ += this->lenF_;     // Storage for Alpha (Total) Fock
  this->lenScr_ += this->lenP_;     // Storage for Alpha (Total) Density
  this->lenScr_ += this->lenCoeff_; // Storage for CDIIS Coefficients
  this->lenScr_ += this->lenB_;     // Storage for CDIIS Metric
  this->lenScr_ += 2*(this->lenCoeff_ - 1) * this->lenF_; // CDIIS Commutator (A) array
  if(!this->isClosedShell) {
    this->lenScr_ += this->lenF_;     // Storage for Beta Fock
    this->lenScr_ += this->lenP_;     // Storage for Beta Density
    this->lenScr_ += 2*(this->lenCoeff_ - 1) * this->lenF_; // CDIIS Commutator (B) array
  }
  if(this->Ref_ == CUHF) {
    this->lenScr_ += this->lenXp_; // Storage for X^(0.5)
    this->lenScr_ += this->lenOccNum_; // Storage for Occupation Numbers (NOs)
    this->lenScr_ += this->lenLambda_; // Storage for Lambda
    this->lenScr_ += this->lenDelF_;   // Stroage for DelF
    this->lenScr_ += this->lenP_;      // Storage for NOs
  }

  this->lenScr_ += this->LWORK_; // LAPACK Scratch space

};

template <typename T>
void SingleSlater<T>::initSCFPtr(){
  this->SCF_SCR        = NULL;
  this->XMem_          = NULL;
  this->FpAlphaMem_    = NULL;
  this->FpBetaMem_     = NULL;
  this->POldAlphaMem_  = NULL;
  this->POldBetaMem_   = NULL;
  this->ErrorAlphaMem_ = NULL;
  this->ErrorBetaMem_  = NULL;
  this->FADIIS_        = NULL;
  this->FBDIIS_        = NULL;
  this->WORK_          = NULL;
  this->XpMem_         = NULL;
  this->lambdaMem_     = NULL;
  this->delFMem_       = NULL;
  this->PNOMem_        = NULL;
  this->occNumMem_     = NULL;
};

template <typename T>
void SingleSlater<T>::initSCFMem(){
  this->initSCFPtr();
  this->initMemLen();

  T* LAST_FOR_SECTION;
  int LEN_LAST_FOR_SECTION;

  this->SCF_SCR = new double[this->lenScr_];
  std::memset(this->SCF_SCR,0.0,this->lenScr_*sizeof(double));

  this->XMem_          = this->SCF_SCR;
  this->FpAlphaMem_    = this->XMem_          + this->lenX_;
  this->POldAlphaMem_  = this->FpAlphaMem_    + this->lenF_;
  this->ErrorAlphaMem_ = this->POldAlphaMem_  + this->lenP_;
  this->FADIIS_        = this->ErrorAlphaMem_ + this->lenF_*(this->lenCoeff_ -1);
  LAST_FOR_SECTION     = this->FADIIS_;
  LEN_LAST_FOR_SECTION = this->lenF_*(this->lenCoeff_ -1);
  if(!this->isClosedShell){
    this->FpBetaMem_     = LAST_FOR_SECTION + LEN_LAST_FOR_SECTION;
    this->POldBetaMem_   = this->FpBetaMem_    + this->lenF_;
    this->ErrorBetaMem_  = this->POldBetaMem_  + this->lenP_;
    this->FBDIIS_        = this->ErrorBetaMem_ + this->lenF_*(this->lenCoeff_ -1);
    LAST_FOR_SECTION     = this->FBDIIS_;
    LEN_LAST_FOR_SECTION = this->lenF_*(this->lenCoeff_ -1);
  }
  if(this->Ref_ == CUHF) {
    this->XpMem_     = LAST_FOR_SECTION + LEN_LAST_FOR_SECTION;
    this->delFMem_   = this->XpMem_     + this->lenX_;
    this->lambdaMem_ = this->delFMem_   + this->lenDelF_;
    this->PNOMem_    = this->lambdaMem_ + this->lenLambda_;
    this->occNumMem_ = this->PNOMem_    + this->lenP_;
    LAST_FOR_SECTION = this->occNumMem_;
    LEN_LAST_FOR_SECTION = this->lenOccNum_;
  }
  
  this->WORK_ = LAST_FOR_SECTION + LEN_LAST_FOR_SECTION;
};
