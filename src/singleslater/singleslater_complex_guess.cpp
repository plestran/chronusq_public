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
#include <singleslater.h>
#include <aointegrals.h>
#include <basisset.h>
#include <workers.h>
#ifdef USE_LIBINT
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::Molecule;
using ChronusQ::HashNAOs;
namespace ChronusQ {
template<>
void SingleSlater<dcomplex>::placeAtmDen(std::vector<int> atomIndex, SingleSlater<double> &hfA){
  // Place atomic SCF densities in the right place of the total density
  // ** Note: ALWAYS spin average, even for UHF **
  for(auto iAtm : atomIndex){
    auto iBfSt = this->basisset_->mapCen2Bf(iAtm)[0];
    auto iSize = this->basisset_->mapCen2Bf(iAtm)[1]; 
    if(this->Ref_ != TCS){
/*
      this->densityA_->block(iBfSt,iBfSt,iSize,iSize)= (*hfA.densityA_);
      if(!this->isClosedShell){
        if(hfA.isClosedShell)
          this->densityB_->block(iBfSt,iBfSt,iSize,iSize)= 2*(*hfA.densityA_);
        else
          this->densityB_->block(iBfSt,iBfSt,iSize,iSize)= 
            (*hfA.densityB_) + (*hfA.densityA_);
      } else {
        if(!hfA.isClosedShell){
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.densityB_);
        }
      }
*/
      this->densityA_->block(iBfSt,iBfSt,iSize,iSize).real()      = (*hfA.densityA());
      if(this->isClosedShell){
        if(hfA.isClosedShell) 
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityA());
        else
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityB());
      } else {
        this->densityB_->block(iBfSt,iBfSt,iSize,iSize).real()    = (*hfA.densityA());
        if(hfA.isClosedShell){
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityA());
          this->densityB_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityA());
        } else {
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityB());
          this->densityB_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityB());
        }
      }
    } else {
      for(auto I = iBfSt, i = 0; I < (iBfSt +iSize); I += 2, i++)
      for(auto J = iBfSt, j = 0; J < (iBfSt +iSize); J += 2, j++){
        (*this->densityA_)(I,J)     = dcomplex((*hfA.densityA())(i,j) + (*hfA.densityB())(i,j),0.0);
        (*this->densityA_)(I+1,J+1) = dcomplex((*hfA.densityA())(i,j) + (*hfA.densityB())(i,j),0.0);
      }
    }
  } // loop iAtm
}
template<>
void SingleSlater<dcomplex>::scaleDen(){
  // Scale UHF densities according to desired multiplicity
  if(!this->isClosedShell && this->Ref_ != TCS){
    int nE = this->molecule_->nTotalE();
    (*this->densityA_) *= dcomplex((double)this->nAE_/(double)nE,0.0);
    (*this->densityB_) *= dcomplex((double)this->nBE_/(double)nE,0.0);
  } else if(this->Ref_ == TCS) {
    int nE = this->molecule_->nTotalE();
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i += 2)
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j += 2){
      (*this->densityA_)(i,j)      *= dcomplex((double)this->nAE_/(double)nE,0.0);
      (*this->densityA_)(i+1,j+1)  *= dcomplex((double)this->nBE_/(double)nE,0.0);
    }
/*
    double theta = math.pi / 8.0;
    double c = std::cos(theta);
    double s = std::sin(theta);
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i += 2)
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j += 2){
      double Paa = (*this->densityA_)(i,j);
      double Pbb = (*this->densityA_)(i+1,j+1);
      (*this->densityA_)(i,j)     = c*c*Paa + s*s*Pbb;
      (*this->densityA_)(i+1,j+1) = c*c*Pbb + s*s*Paa;
      (*this->densityA_)(i+1,j)   = c*s*(Paa - Pbb);
      (*this->densityA_)(i,j+1)   = c*s*(Paa - Pbb);
     
    }
*/
    
//  (*this->densityA_) *= (double)(this->nAE_+this->nBE_)/(double)nE ;
  }
//CErr();
}; // SingleSlater::scaleDen [T=dcomplex]
//--------------------------------//
// form the initial guess of MO's //
//--------------------------------//
template<>
void SingleSlater<dcomplex>::formGuess() {
  
  int readNPGTO,L, nsize;
  this->moA_->setZero();
  if(!this->isClosedShell && this->Ref_ != TCS) this->moB_->setZero();

  if(this->molecule_->nAtoms() > 1) {
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
                          " atomic SCF calculations to form the initial guess" << endl;
 
    // Loop and perform CUHF on each atomic center
    for(auto iUn = 0; iUn < uniqueElement.size(); iUn++){
      // Local objects to be constructed and destructed at every loop
      AOIntegrals aointegralsAtom;
      SingleSlater<double> hartreeFockAtom;
      Controls controlAtom;
      BasisSet basisSetAtom;
      BasisSet dfBasisSetAtom;
      Molecule uniqueAtom(uniqueElement[iUn],this->fileio_->out);
 
      // FIXME: This only makes sense for neutral molecules
      uniqueAtom.setCharge(0);
      uniqueAtom.setMultip(uniqueElement[iUn].defaultMult);
 
      // Construct atomic basis set from the reference
      this->basisset_->constructExtrn(&uniqueAtom,&basisSetAtom);
      // Generate basis maps
      basisSetAtom.makeMapSh2Bf(1);
      basisSetAtom.makeMapSh2Cen(&uniqueAtom);
      basisSetAtom.renormShells(); // Libint throws a hissy fit without this
 
      controlAtom.iniControls();
      controlAtom.doCUHF = true; // Can set to false too if UHF guess is desired
 
      // Initialize the local integral and SS classes
      aointegralsAtom.iniAOIntegrals(&uniqueAtom,&basisSetAtom,this->fileio_,&controlAtom,
        &dfBasisSetAtom);
      hartreeFockAtom.iniSingleSlater(&uniqueAtom,&basisSetAtom,&aointegralsAtom,
        this->fileio_,&controlAtom);
 
      // Zero out the MO coeff for local SS object
      hartreeFockAtom.moA()->setZero();
      if(!hartreeFockAtom.isClosedShell) hartreeFockAtom.moB()->setZero();
      hartreeFockAtom.haveMO = true;
 
      // Prime and perform the atomic SCF
      hartreeFockAtom.formFock();
      hartreeFockAtom.computeEnergy();
      hartreeFockAtom.SCF();
      
      // Place Atomic Densities into Total Densities
      this->placeAtmDen(atomIndex[iUn],hartreeFockAtom);
 
    } // Loop iUn
 
    this->scaleDen();
  }

  // Set flags to use in the rest of code
  this->haveMO = true;
  if(this->molecule_->nAtoms() > 1) this->haveDensity = true;
};
}; //namespace ChronusQ
#endif
