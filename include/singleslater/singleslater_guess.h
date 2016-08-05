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
  if(this->guess_ == SAD) this->SADGuess();
  else if(this->guess_ == CORE) this->COREGuess();
  else if(this->guess_ == READ) this->READGuess();
  else if(this->guess_ == RANDOM) this->RandomGuess();
  else CErr("Guess NYI",this->fileio_->out);
}

template<typename T>
void SingleSlater<T>::COREGuess(){
  this->aointegrals_->computeAOOneE(); 
  if(this->nTCS_ == 1 && this->isClosedShell){
    this->fockA_->real() = (*this->aointegrals_->coreH_);
  } else {
    this->fockScalar_->real() = 2*(*this->aointegrals_->coreH_);
    this->fockMz_->setZero();

    std::vector<std::reference_wrapper<TMap>> toGather;
    toGather.emplace_back(*this->fockScalar_);
    toGather.emplace_back(*this->fockMz_);

    if(this->nTCS_ ==  2) {
      this->fockMy_->setZero();
      this->fockMx_->setZero();
      toGather.emplace_back(*this->fockMy_);
      toGather.emplace_back(*this->fockMx_);

      Quantum<T>::spinGather(*this->fockA_,toGather);
    } else 
      Quantum<T>::spinGather(*this->fockA_,*this->fockB_,toGather);
  }
  this->haveMO = true;
  this->haveDensity = true;

  //this->orthoFock();
  this->orthoFock3();
  if(this->printLevel_ > 3){
    prettyPrint(this->fileio_->out,*this->fockA_,"Initial FA");
    prettyPrint(this->fileio_->out,*this->fockOrthoA_,"Initial FOrthoA");
  }
  this->diagFock2();
  this->mixOrbitalsSCF();
  this->formDensity();
  this->cpyAOtoOrthoDen();
  this->unOrthoDen3();
  if(this->printLevel_ > 3){
    prettyPrint(this->fileio_->out,*this->onePDMA_,"Initial PA");
    prettyPrint(this->fileio_->out,*this->onePDMOrthoA_,"Initial POrthoA");
  }
//this->computeEnergy();
  this->backTransformMOs();
  if(this->printLevel_ > 3) {
    prettyPrint(this->fileio_->out,*this->moA_,"Initial MOA");
  }

};


template<typename T>
void SingleSlater<T>::SADGuess() {

  this->moA_->setZero();
  if(!this->isClosedShell && this->nTCS_ == 1) 
    this->moB_->setZero();
  
  // Set flags to use in the rest of code
  this->haveMO = true;
  if(this->molecule_->nAtoms() > 1) this->haveDensity = true;

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
    aointegralsAtom.integralAlgorithm = 
      AOIntegrals::INTEGRAL_ALGORITHM::DIRECT;
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
 
    // Zero out the MO coeff for local SS object
    if(getRank() == 0){
      hartreeFockAtom.moA()->setZero();
      if(!hartreeFockAtom.isClosedShell) hartreeFockAtom.moB()->setZero();
    }
    hartreeFockAtom.haveMO = true;
 
    // Prime and perform the atomic SCF
    hartreeFockAtom.formFock();
    hartreeFockAtom.computeEnergy();
    hartreeFockAtom.SCF3();
    
    // Place Atomic Densities into Total Densities
    this->placeAtmDen(atomIndex[iUn],hartreeFockAtom);
 
  } // Loop iUn

  this->scaleDen();
  this->scatterDensity();
  this->formFock();
  this->orthoFock3();
  this->diagFock2();
  this->mixOrbitalsSCF();
  this->formDensity();
  this->cpyAOtoOrthoDen();
  this->unOrthoDen3();
  this->backTransformMOs();
};
