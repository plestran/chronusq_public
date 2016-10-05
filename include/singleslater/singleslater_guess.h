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
void SingleSlater<T>::formGuess(){
  if(not this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

  if(this->guess_ == ONLY) {
    this->unOrthoDen3();
    return;
  } else if(this->guess_ == SAD) this->SADGuess();
  else if(this->guess_ == CORE) this->COREGuess();
  else if(this->guess_ == READ) this->READGuess();
  else if(this->guess_ == RANDOM) this->RandomGuess();
  else CErr("Guess NYI",this->fileio_->out);

  if(this->printLevel_ > 3) {
    this->fileio_->out << "Initial Fock Matrix" << endl;
    this->printFock();
  }
  this->orthoFock3();

  this->populateMO4Diag();
  this->diagFock2();
//this->mixOrbitalsSCF();
  if(this->printLevel_ > 3) {
    this->fileio_->out << "Initial Density Matrix Before" << endl;
    this->printDensity();
  }
  this->formDensity();
  if(this->printLevel_ > 3) {
    this->fileio_->out << "Initial Density Matrix Before" << endl;
    this->printDensity();
  }
  this->cpyAOtoOrthoDen();
  this->unOrthoDen3();
//this->scaleDen();

  if(this->printLevel_ > 3) {
    this->fileio_->out << "Initial Density Matrix" << endl;
    this->printDensity();
  }

  this->backTransformMOs();
}

template<typename T>
void SingleSlater<T>::COREGuess(){
  if(not this->aointegrals_->haveAOOneE)
    this->aointegrals_->computeAOOneE(); 


  for(auto fock : this->fock_) fock->setZero();

  (*this->fockScalar_) = this->aointegrals_->coreH_->template cast<T>();
  if(this->nTCS_ == 2 and this->aointegrals_->doX2C) {
    (*this->fockMx_) = this->aointegrals_->oneEmx_->template cast<T>();
    (*this->fockMy_) = this->aointegrals_->oneEmy_->template cast<T>();
    (*this->fockMz_) = this->aointegrals_->oneEmz_->template cast<T>();
    (*this->fockMx_) *= ComplexScale<T>();
    (*this->fockMy_) *= ComplexScale<T>();
    (*this->fockMz_) *= ComplexScale<T>();
  }

  for(auto iF = fock_.begin(); iF != fock_.end(); iF++)
    *(*iF) *= 2;
/*
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    this->fockMz_->setZero();
    if(this->nTCS_ ==  2) {
      this->fockMy_->setZero();
      this->fockMx_->setZero();
    } 
  }
*/
};


template<typename T>
void SingleSlater<T>::SADGuess() {

  this->moA_->setZero();
  if(!this->isClosedShell && this->nTCS_ == 1) 
    this->moB_->setZero();
  
  // Set flags to use in the rest of code
  if(this->molecule_->nAtoms() <= 1) return;

  // Determining unique atoms
  std::vector<Atoms> uniqueElement;
  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    if(iAtm == 0){ 
      uniqueElement.push_back(elements[this->molecule_->index(iAtm)]);
    }
    bool uniq = true;
    for(auto iUn = 0; iUn < uniqueElement.size(); iUn++){
      if(uniqueElement[iUn].atomicNumber == 
        elements[this->molecule_->index(iAtm)].atomicNumber){
        uniq = false;
        break;
      }
    }
    if(uniq) {
      uniqueElement.push_back(elements[this->molecule_->index(iAtm)]);
    }
  }
 
  // Generate a map of unique atoms to centers
  std::vector<std::vector<int>> atomIndex(uniqueElement.size());
  for(auto iUn = 0; iUn < uniqueElement.size(); iUn++){
    for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
      if(uniqueElement[iUn].atomicNumber == 
        elements[this->molecule_->index(iAtm)].atomicNumber){
        
        atomIndex[iUn].push_back(iAtm);
      }
    }
  }
 
  this->fileio_->out << "Running " << uniqueElement.size() << 
                      " atomic SCF calculations to form the initial guess" 
                     << endl;
 
  // Loop and perform CUHF on each atomic center
  for(auto iUn = 0; iUn < uniqueElement.size(); iUn++){
    // Local objects to be constructed and destructed at every loop
    AOIntegrals aointegralsAtom;
    SingleSlater<double> hartreeFockAtom;
    BasisSet basisSetAtom;
    BasisSet dfBasisSetAtom;
    Molecule uniqueAtom(uniqueElement[iUn],this->fileio_->out);
 
    // FIXME: This only makes sense for neutral molecules
    uniqueAtom.setCharge(0);
    uniqueAtom.setMultip(uniqueElement[iUn].defaultMult);
 
    // Construct atomic basis set from the reference
    this->basisset_->constructExtrn(&uniqueAtom,&basisSetAtom);

    // Generate basis maps
    basisSetAtom.makeMaps(&uniqueAtom);
    basisSetAtom.renormShells(); // Libint throws a hissy fit without this
 
 
    // Initialize the local integral and SS classes
    aointegralsAtom.isPrimary = false;
    hartreeFockAtom.isNotPrimary();
    
    // Replaces iniAOIntegrals
    aointegralsAtom.communicate(uniqueAtom,basisSetAtom,*this->fileio_,
      *this->aointegrals_->memManager());
    aointegralsAtom.initMeta();
    aointegralsAtom.integralAlgorithm = AOIntegrals::INTEGRAL_ALGORITHM::DIRECT;
    aointegralsAtom.alloc();

    hartreeFockAtom.communicate(uniqueAtom,basisSetAtom,aointegralsAtom,
      *this->fileio_,*this->memManager_);
  
    hartreeFockAtom.setGuess(GUESS::CORE);
    hartreeFockAtom.initMeta();
    hartreeFockAtom.setField(this->elecField_);
    hartreeFockAtom.isClosedShell = (hartreeFockAtom.multip() == 1); 
    hartreeFockAtom.doDIIS = false;
    hartreeFockAtom.isDFT = false;
    hartreeFockAtom.isHF = true;
    hartreeFockAtom.setRef(UHF);
    hartreeFockAtom.genMethString();
    hartreeFockAtom.alloc();

    if(this->printLevel_ < 4) hartreeFockAtom.setPrintLevel(0);
    else hartreeFockAtom.setPrintLevel(this->printLevel_);
 
    // Zero out the MO coeff for local SS object
    if(getRank() == 0){
      hartreeFockAtom.moA()->setZero();
      if(!hartreeFockAtom.isClosedShell) hartreeFockAtom.moB()->setZero();
    }
 
    // Prime and perform the atomic SCF
    hartreeFockAtom.formGuess();
    hartreeFockAtom.formFock();
    hartreeFockAtom.computeEnergy();
    hartreeFockAtom.SCF3();
    hartreeFockAtom.computeProperties();
    hartreeFockAtom.printProperties();
    
    
    // Place Atomic Densities into Total Densities
    this->placeAtmDen(atomIndex[iUn],hartreeFockAtom);
 
  } // Loop iUn

//prettyPrintSmart(cout,*this->onePDMScalar_,"PS");
  this->scaleDen();
//prettyPrintSmart(cout,*this->onePDMScalar_,"PS");
  this->formFock();
//prettyPrintSmart(cout,*this->fockScalar_,"FS");

//if(not this->isClosedShell)
//  prettyPrintSmart(cout,*this->fockMz_,"FZ");

};

template <typename T>
void SingleSlater<T>::RandomGuess() {
/*
  (*this->onePDMScalar_) = TMatrix::Random(this->nBasis_,this->nBasis_);
  (*this->onePDMScalar_) = 
    this->onePDMScalar_->template selfadjointView<Lower>();
  if(this->nTCS_ == 2 and !this->isClosedShell){
    (*this->onePDMMz_) = TMatrix::Random(this->nBasis_,this->nBasis_);
    (*this->onePDMMz_) = this->onePDMMz_->template selfadjointView<Lower>();
  }
*/
  (*this->moA_) = 
    TMatrix::Random(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_);
  if(this->nTCS_ == 1 and !this->isClosedShell)
    (*this->moB_) = TMatrix::Random(this->nBasis_,this->nBasis_);

  this->diagFock2();
  this->formDensity();
  this->cpyAOtoOrthoDen();
  this->unOrthoDen3();
  this->scaleDen();
  this->formFock();
}

template<typename T>
void SingleSlater<T>::placeAtmDen(std::vector<int> atomIndex, 
  SingleSlater<double> &hfA){
  // Place atomic SCF densities in the right place of the total density
  // ** Note: ALWAYS spin average, even for UHF **
  for(auto iAtm : atomIndex){
    auto iBfSt = this->basisset_->mapCen2Bf(iAtm)[0];
    auto iSize = this->basisset_->mapCen2Bf(iAtm)[1]; 

    this->onePDMScalar_->block(iBfSt,iBfSt,iSize,iSize) =
      hfA.onePDMScalar()->template cast<T>();
    if(this->nTCS_ == 2 or !this->isClosedShell)
      this->onePDMMz_->block(iBfSt,iBfSt,iSize,iSize) =
        hfA.onePDMMz()->template cast<T>();
  } // loop iAtm
}

template<typename T>
void SingleSlater<T>::scaleDen(){
  if(this->nTCS_ == 1 and this->isClosedShell) {
//  prettyPrint(cout,*this->aointegrals_->overlap_,"S");
//  cout <<this->template computeProperty<double,DENSITY_TYPE::TOTAL>(*this->aointegrals_->overlap_) << endl;
    (*this->onePDMScalar_) *= 
      T(this->molecule_->nTotalE()) / 
      T(this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
        *this->aointegrals_->overlap_));
  } else {
    // SCR  = PA
    // SCR2 = PB
    /*
    (*this->NBSqScratch_) = (*this->onePDMScalar_) + (*this->onePDMMz_);
    (*this->NBSqScratch2_) = (*this->onePDMScalar_) - (*this->onePDMMz_);
    (*this->NBSqScratch_) *= 0.5; 
    (*this->NBSqScratch2_) *= 0.5; 

    double TS = this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
      *this->aointegrals_->overlap_); 
    double TZ = this->template computeProperty<double,DENSITY_TYPE::MZ>(
      *this->aointegrals_->overlap_); 

    double TA = 0.5 * (TS + TZ);
    double TB = 0.5 * (TS - TZ);
    cout << "TA " << TA << endl;
    cout << "TB " << TB << endl;


    (*this->NBSqScratch_) *= T(this->nOA_) / T(TA);
    (*this->NBSqScratch2_) *= T(this->nOB_) / T(TB);

    (*this->onePDMScalar_) = (*this->NBSqScratch_) + (*this->NBSqScratch2_);
    (*this->onePDMMz_)     = (*this->NBSqScratch_) - (*this->NBSqScratch2_);
    */
    
    double TS = this->template computeProperty<double,DENSITY_TYPE::TOTAL>(
      *this->aointegrals_->overlap_); 
    double TZ = this->template computeProperty<double,DENSITY_TYPE::MZ>(
      *this->aointegrals_->overlap_); 

   (*this->onePDMScalar_) *= T(this->nOA_ + this->nOB_) / TS;
   (*this->onePDMMz_)     *= T(this->nOA_ - this->nOB_) / TZ;
  }
}

template <typename T>
void SingleSlater<T>::READGuess(){
  SCFDensityScalar_->read(this->onePDMScalar_->data(),H5PredType<T>());
  if(this->nTCS_ == 2 or !this->isClosedShell)
    SCFDensityMz_->read(this->onePDMMz_->data(),H5PredType<T>());
  if(this->nTCS_ == 2){
    SCFDensityMy_->read(this->onePDMMy_->data(),H5PredType<T>());
    SCFDensityMx_->read(this->onePDMMx_->data(),H5PredType<T>());
  }
  this->formFock();
};
