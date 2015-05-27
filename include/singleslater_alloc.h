template<typename T>
void SingleSlater<T>::iniSingleSlater(Molecule * molecule, BasisSet * basisset, 
                                   AOIntegrals * aointegrals, FileIO * fileio, 
                                   Controls * controls) {
  int nTotalE = molecule->nTotalE();
  this->nBasis_  = basisset->nBasis();
  this->nTT_   = this->nBasis_*(this->nBasis_+1)/2;
  this->spin_  = molecule->spin();
  int nSingleE = this->spin_ - 1;
  this->nOccB_ = (nTotalE - nSingleE)/2;
  this->nVirB_ = this->nBasis_ - this->nOccB_;
  this->nOccA_ = this->nOccB_ + nSingleE;
  this->nVirA_ = this->nBasis_ - this->nOccA_;
  this->energyNuclei = molecule->energyNuclei();
  if(this->spin_!=1) this->RHF_ = 0;
  else this->RHF_ = 1;

  // FIXME Nedd try statements for allocation
  try { this->densityA_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Alpha Density
  catch (...) { CErr(std::current_exception(),"Alpha Density Matrix Allocation"); }
  try { this->fockA_     = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Alpha Fock
  catch (...) { CErr(std::current_exception(),"Alpha Fock Matrix Allocation"); }
#ifndef USE_LIBINT
  try { this->coulombA_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Alpha Coulomb Integral
  catch (...) { CErr(std::current_exception(),"Alpha Coulomb Tensor (R2) Allocation"); }
  try { this->exchangeA_ = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); }// Alpha Exchange Integral
  catch (...) { CErr(std::current_exception(),"Alpha Exchange Tensor (R2) Allocation"); }
#else // USE_LIBINT
  try { this->PTA_  = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Alpha Perturbation Tensor
  catch (...) { CErr(std::current_exception(),"Alpha Perturbation Tensor (G[P]) Allocation"); }
#endif
  try { this->moA_       = std::unique_ptr<TMatrix>(new TMatrix(this->nBasis_,this->nBasis_)); } // Alpha Molecular Orbital Coefficients
  catch (...) { CErr(std::current_exception(),"Alpha MO Coefficients Allocation"); }
  

  if(!this->RHF_) {
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
  };

  this->dipole_ = std::unique_ptr<TMatrix>(new TMatrix(3,1));
  this->quadpole_ = std::unique_ptr<TMatrix>(new TMatrix(3,3));
  this->tracelessQuadpole_ = std::unique_ptr<TMatrix>(new TMatrix(3,3));
  this->octpole_  = std::unique_ptr<RealTensor3d>(new RealTensor3d(3,3,3));
  this->molecule_ = molecule;
  this->basisset_ = basisset;
  this->fileio_   = fileio;
  this->controls_ = controls;
  this->aointegrals_= aointegrals;
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
  this->havePT = false;
};

