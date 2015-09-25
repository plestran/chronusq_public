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
void SingleSlater<double>::placeAtmDen(std::vector<int> atomIndex, SingleSlater<double> &hfA){
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
      this->densityA_->block(iBfSt,iBfSt,iSize,iSize)= (*hfA.densityA_);
      if(this->isClosedShell){
        if(hfA.isClosedShell) 
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.densityA_);
        else
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.densityB_);
      } else {
        this->densityB_->block(iBfSt,iBfSt,iSize,iSize)= (*hfA.densityA_);
        if(hfA.isClosedShell){
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.densityA_);
          this->densityB_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.densityA_);
        } else {
          this->densityA_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.densityB_);
          this->densityB_->block(iBfSt,iBfSt,iSize,iSize) += (*hfA.densityB_);
        }
      }
    } else {
      for(auto I = iBfSt, i = 0; I < (iBfSt +iSize); I += 2, i++)
      for(auto J = iBfSt, j = 0; J < (iBfSt +iSize); J += 2, j++){
        (*this->densityA_)(I,J)     = (*hfA.densityA_)(i,j) + (*hfA.densityB_)(i,j);
        (*this->densityA_)(I+1,J+1) = (*hfA.densityA_)(i,j) + (*hfA.densityB_)(i,j);
      }
    }
  } // loop iAtm
}
template<>
void SingleSlater<double>::scaleDen(){
  // Scale UHF densities according to desired multiplicity
  if(!this->isClosedShell && this->Ref_ != TCS){
    int nE = this->molecule_->nTotalE();
    (*this->densityA_) *= (double)this->nAE_/(double)nE ;
    (*this->densityB_) *= (double)this->nBE_/(double)nE ;
  } else if(this->Ref_ == TCS) {
    int nE = this->molecule_->nTotalE();
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i += 2)
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j += 2){
      (*this->densityA_)(i,j)      *= (double)this->nAE_/(double)nE ;
      (*this->densityA_)(i+1,j+1)  *= (double)this->nBE_/(double)nE ;
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
}; // SingleSlater::scaleDen [T=double]
//--------------------------------//
// form the initial guess of MO's //
//--------------------------------//
template<>
void SingleSlater<double>::formGuess() {
  
  std::vector<RealMatrix> atomMO;
  std::vector<RealMatrix> atomMOB;
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
 
    this->fileio_->out << "Found " << uniqueElement.size() << 
                          " unique atoms in molecule" << endl;
    
    this->fileio_->out << endl << "Atomic SCF Starting............." << endl << endl;
 
    // Loop and perform CUHF on each atomic center
    for(auto iUn = 0; iUn < uniqueElement.size(); iUn++){
      // Local objects to be constructed and destructed at every loop
      AOIntegrals aointegralsAtom;
      SingleSlater<double> hartreeFockAtom;
      Controls controlAtom;
      BasisSet basisSetAtom;
      BasisSet dfBasisSetAtom;
      Molecule uniqueAtom(uniqueElement[iUn],this->fileio_);
 
      // FIXME: This only makes sense for neutral molecules
      uniqueAtom.readCharge(0);
      uniqueAtom.readMultip(uniqueElement[iUn].defaultMult);
 
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
      hartreeFockAtom.moA_->setZero();
      if(!hartreeFockAtom.isClosedShell) hartreeFockAtom.moB_->setZero();
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
//------------------------------------------//
// form the initial guess of MOs from input //
//------------------------------------------//
template<>
void SingleSlater<double>::readGuessIO() {
  int i,j;
  this->fileio_->in.clear();
  std::string readString;
  this->fileio_->in.seekg(0,ios::beg);
  this->fileio_->in>>readString;
  readString=stringupper(readString);
  while(!this->fileio_->in.eof()&&readString.compare("$GUESS")) {
    this->fileio_->in>>readString;
    readString=stringupper(readString);
  };
  this->fileio_->in >> readString;
  for(i=0;i<this->nBasis_;i++) for(j=0;j<this->nBasis_;j++){
    this->fileio_->in>>(*(this->moA_))(j,i);
  };
  this->fileio_->in.close();
  if(this->controls_->printLevel>=3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(this->Ref_ != RHF) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
};
//-----------------------------------------------------------------------//
// form the initial guess of MOs from Gaussian raw matrix element file   //
//-----------------------------------------------------------------------//
template<>
void SingleSlater<double>::readGuessGauMatEl(GauMatEl& matEl){
  this->fileio_->out << "Reading MO coefficients from " <<matEl.fname()<< endl;
  if(matEl.nBasis()!=this->nBasis_) CErr("Basis Set mismatch",this->fileio_->out);
  double *scr = matEl.readRec(GauMatEl::moa); 
  for(auto i = 0; i < this->nBasis_*this->nBasis_; i++)
    this->moA_->data()[i] = scr[i];
  this->moA_->transposeInPlace(); // Row Major
  delete [] scr;
  if(this->Ref_ != RHF) {
    scr = matEl.readRec(GauMatEl::mob); 
    for(auto i = 0; i < this->nBasis_*this->nBasis_; i++)
      this->moB_->data()[i] = scr[i];
    this->moB_->transposeInPlace(); // Row Major
    delete [] scr;
  }
  this->matchord();
  if(this->controls_->printLevel>=3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(this->Ref_ != RHF) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
}
//-----------------------------------------------------------------------//
// form the initial guess of MOs from Gaussian formatted checkpoint file //
//-----------------------------------------------------------------------//
template<>
void SingleSlater<double>::readGuessGauFChk(std::string &filename) {
  this->fileio_->out<<"Reading formatted checkpoint file "<<filename<<endl;
  std::string readString;
  int i,j,nBasis;
  double data;
  std::unique_ptr<ifstream> fchk(new ifstream(filename));
  if(fchk->fail()) 
    CErr("Could not find "+filename,this->fileio_->out);

  *fchk>>readString;
  while((!(fchk->eof()))&&(readString.compare("basis"))) *fchk>>readString;
  if(!readString.compare("basis")) {
    *fchk>>readString;
    *fchk>>readString;
    *fchk>>nBasis;
    if(nBasis!=this->nBasis_) {
      this->fileio_->out<<"Basis set mismatch! "<<nBasis<<"!="<<this->nBasis_<<endl;
      exit(1);
    };
  } else {
    this->fileio_->out<<"No basis set found in the formatted checkpoint file! "<<endl;
  };

  while(!(fchk->eof())&&readString.compare("MO")) *fchk>>readString;
  if(!readString.compare("MO")) {
    *fchk>>readString;
    *fchk>>readString;
    *fchk>>readString;
    *fchk>>readString;
    for(i=0;i<nBasis;i++) 
      for(j=0;j<nBasis;j++){
        *fchk>>data;
        (*(this->moA_))(j,i) = data;
    };
  } else {
    this->fileio_->out<<"No basis set found in the formatted checkpoint file! "<<endl;
  };

  fchk->close();
  if(this->controls_->printLevel>=3) {
    prettyPrint(this->fileio_->out,(*this->moA_),"Alpha MO Coeff");
    if(this->Ref_ != RHF) prettyPrint(this->fileio_->out,(*this->moB_),"Beta MO Coeff");
  };
  this->haveMO = true;
};
}; //namespace ChronusQ
#endif