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
  this->nAE_   = this->nOccA_;
  this->nBE_   = this->nOccB_;
  this->energyNuclei = molecule->energyNuclei();
  this->isConverged = false;
  this->denTol_ = controls->SCFdenTol_;
  this->eneTol_ = controls->SCFeneTol_;
  this->maxSCFIter_ = controls->SCFmaxIter_;

  this->elecField_ = controls->field_;

  this->isClosedShell = (this->multip_ == 1);
  if(controls->HF){
    if(this->isClosedShell && !controls->doCUHF
       && !controls->doTCS)                        this->Ref_ = RHF ; // RHF
    else if(!controls->doCUHF && !controls->doTCS) this->Ref_ = UHF ; // UHF
    else if(controls->doCUHF)                      this->Ref_ = CUHF; // CUHF
    else if(controls->doTCS)                       this->Ref_ = TCS ; // TCS
  } else if(controls->DFT) {
    if(this->isClosedShell && !controls->doCUHF
       && !controls->doTCS)                        this->Ref_ = RKS ; // RKS
    else if(!controls->doCUHF && !controls->doTCS) this->Ref_ = UKS ; // UKs
    else if(controls->doCUHF)                      this->Ref_ = CUKS; // CUKS
    else if(controls->doTCS)                       this->Ref_ = GKS ; // GKS
  }

  this->getAlgebraicField(); 
  this->SCFType_      = this->algebraicField_      + " ";
  this->SCFTypeShort_ = this->algebraicFieldShort_ + "-";
  if(this->Ref_ == RHF) {
    this->SCFType_      += "Resricted Hartree-Fock"; 
    this->SCFTypeShort_ += "RHF";
  } else if(this->Ref_ == UHF) {
    this->SCFType_      += "Unresricted Hartree-Fock"; 
    this->SCFTypeShort_ += "UHF";
  } else if(this->Ref_ == CUHF) {
    this->SCFType_      += "Constrained Unresricted Hartree-Fock"; 
    this->SCFTypeShort_ += "CUHF";
  } else if(this->Ref_ == TCS) {
    this->SCFType_      += "Generalized Hartree-Fock"; 
    this->SCFTypeShort_ += "GHF";
  }


  this->nTCS_ = 1;
  if(this->Ref_ == TCS) this->nTCS_ = 2;
// Comment out to get rid of DFT tests
//this->controls_->DFT = true;
  

  // Alpha / TCS Density
  try { 
    this->densityA_  = std::unique_ptr<TMatrix>( 
      new TMatrix(this->nTCS_*this->nBasis_, this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(),"TCS Density Matrix Allocation"  ); 
    else                  CErr(std::current_exception(),"Alpha Density Matrix Allocation"); 
  }

  // Alpha / TCS Fock
  try { 
    this->fockA_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(),"TCS Fock Matrix Allocation"); 
    else                  CErr(std::current_exception(),"Alpha Fock Matrix Allocation"); 
  }

#ifndef USE_LIBINT
  // Alpha / TCS Coulomb Matrix
  try { 
    this->coulombA_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(),"TCS Coulomb Tensor Allocation"); 
    else                  CErr(std::current_exception(),"Alpha Coulomb Tensor Allocation"); 
  }

  // Alpha / TCS Exchange Matrix
  try { 
    this->exchangeA_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(),"TCS Exchange Tensor Allocation"); 
    else                  CErr(std::current_exception(),"Alpha Exchange Tensor Allocation"); 
  }
#else // USE_LIBINT
  // Alpha / TCS Perturbation Tensor
  try { 
    this->PTA_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(),
      "TCS Perturbation Tensor (G[P]) Allocation"); 
    else CErr(std::current_exception(),"Alpha Perturbation Tensor (G[P]) Allocation"); 
  }
#endif
  if(this->controls_->DFT) {
    // Alpha / TCS VXC
    try { 
      this->vXCA_  = std::unique_ptr<TMatrix>(
        new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
    } catch (...) { 
      if(this->Ref_ == TCS) CErr(std::current_exception(),
        "TCS VXC Allocation"); 
      else CErr(std::current_exception(),"Alpha VXC  Allocation"); 
    }
  }
  // Alpha / TCS Molecular Orbital Coefficients
  try { 
    this->moA_ = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_)); 
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(),"TCS MO Coefficients Allocation");
    else                  CErr(std::current_exception(),"Alpha MO Coefficients Allocation"); 
  }

  // Alpha / TCS Eigenorbital Energies
  try { 
    this->epsA_ = std::unique_ptr<RealMatrix>(new RealMatrix(this->nTCS_*this->nBasis_,1)); 
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(),"TCS Eigenorbital Energies"); 
    else                  CErr(std::current_exception(),"Alpha Eigenorbital Energies"); 
  }
  

  if(!this->isClosedShell && this->Ref_ != TCS) {
    try { this->densityB_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Density
    catch (...) { CErr(std::current_exception(),"Beta Density Matrix Allocation"); }
    try { this->fockB_     = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Fock
    catch (...) { CErr(std::current_exception(),"Beta Fock Matrix Allocation"); }
#ifndef USE_LIBINT
    try { this->coulombB_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Coulomb Integral
    catch (...) { CErr(std::current_exception(),"Beta Coulomb Tensor Allocation"); }
    try { this->exchangeB_ = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Exchange Integral
    catch (...) { CErr(std::current_exception(),"Beta Exchange Tensor Allocation"); }
#else // USE_LIBINT
    try { this->PTB_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Perturbation Tensor
    catch (...) { CErr(std::current_exception(),"Beta Perturbation Tensor (G[P]) Allocation"); }
#endif
    if(this->controls_->DFT){
      try { this->vXCB_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta VXC
      catch (...) { CErr(std::current_exception(),"Beta VXC Allocation"); }
    }
    try { this->moB_       = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Beta Molecular Orbital Coefficients
    catch (...) { CErr(std::current_exception(),"Beta MO Coefficients Allocation"); }
    try { this->epsB_       = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,1)); } // Beta Eigenorbital Energies
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
void SingleSlater<T>::alloc(){
  this->checkMeta();
  this->allocOp();
  if(this->doMultipole) this->allocMultipole(); 
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
  // Alpha / TCS Density Matrix
  try { 
    this->densityA_  = std::unique_ptr<TMatrix>( 
      new TMatrix(this->nTCS_*this->nBasis_, this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) 
      CErr(std::current_exception(),"TCS Density Matrix Allocation"  ); 
    else CErr(std::current_exception(),"Alpha Density Matrix Allocation"); 
  }

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

  if(!this->isClosedShell && this->Ref_ != TCS) {
    // Beta Density Matrix
    try { 
      this->densityB_ = std::unique_ptr<TMatrix>(
        new TMatrix(this->nBasis_,this->nBasis_)); 
    } catch (...) { 
      CErr(std::current_exception(),"Beta Density Matrix Allocation"); 
    }

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

  if(!this->isClosedShell && this->Ref_ != TCS) {
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

  if(!this->isClosedShell && this->Ref_ != TCS) {
    // Beta Perturbation Tensor
    try { 
      this->PTB_  = std::unique_ptr<TMatrix>(
        new TMatrix(this->nBasis_,this->nBasis_));
    } catch (...) { 
      CErr(std::current_exception(),"Beta G[P] Allocation"); 
    }
  }
#endif

  if(this->isDFT) this->allocDFT();
}

template<typename T>
void SingleSlater<T>::allocDFT(){
  // Alpha / TCS VXC
  try { 
    this->vXCA_  = std::unique_ptr<TMatrix>(
      new TMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
  } catch (...) { 
    if(this->Ref_ == TCS) CErr(std::current_exception(), "TCS VXC Allocation"); 
    else CErr(std::current_exception(),"Alpha VXC  Allocation"); 
  }

  if(!this->isClosedShell && this->Ref_ != TCS) {
    // Beta VXC
    try { 
      this->vXCB_  = std::unique_ptr<TMatrix>(
        new TMatrix(this->nBasis_,this->nBasis_));
    } catch (...) { 
      CErr(std::current_exception(),"Beta VXC  Allocation"); 
    }
  }
}

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
