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
#include <singleslater.h>
#include <aointegrals.h>
#include <basisset.h>
#include <workers.h>
#ifdef USE_LIBINT
using ChronusQ::BasisSet;
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
    if(this->nTCS_ == 1){
      this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize).real()      = (*hfA.densityA());
      if(this->isClosedShell){
        if(hfA.isClosedShell) 
          this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityA());
        else
          this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityB());
      } else {
        this->onePDMB_->block(iBfSt,iBfSt,iSize,iSize).real()    = (*hfA.densityA());
        if(hfA.isClosedShell){
          this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityA());
          this->onePDMB_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityA());
        } else {
          this->onePDMA_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityB());
          this->onePDMB_->block(iBfSt,iBfSt,iSize,iSize).real() += (*hfA.densityB());
        }
      }
    } else {
      for(auto I = iBfSt, i = 0; I < (iBfSt +iSize); I += 2, i++)
      for(auto J = iBfSt, j = 0; J < (iBfSt +iSize); J += 2, j++){
        (*this->onePDMA_)(I,J)     = dcomplex((*hfA.densityA())(i,j) + (*hfA.densityB())(i,j),0.0);
        (*this->onePDMA_)(I+1,J+1) = dcomplex((*hfA.densityA())(i,j) + (*hfA.densityB())(i,j),0.0);
      }
    }
  } // loop iAtm
}
template<>
void SingleSlater<dcomplex>::scaleDen(){
  // Scale UHF densities according to desired multiplicity
  if(!this->isClosedShell && this->nTCS_ == 1){
    int nE = this->molecule_->nTotalE();
    (*this->onePDMA_) *= dcomplex((double)this->nAE_/(double)nE,0.0);
    (*this->onePDMB_) *= dcomplex((double)this->nBE_/(double)nE,0.0);
  } else if(this->nTCS_ == 2) {
    int nE = this->molecule_->nTotalE();
    for(auto i = 0; i < this->nTCS_*this->nBasis_; i += 2)
    for(auto j = 0; j < this->nTCS_*this->nBasis_; j += 2){
      (*this->onePDMA_)(i,j)      *= dcomplex((double)this->nAE_/(double)nE,0.0);
      (*this->onePDMA_)(i+1,j+1)  *= dcomplex((double)this->nBE_/(double)nE,0.0);
    }
  }
}; // SingleSlater::scaleDen [T=dcomplex]

//--------------------------------//
// form the initial guess of MO's //
//--------------------------------//
template<>
void SingleSlater<dcomplex>::RandomGuess() {
//JJG make random guess
  auto NTCSxNBASIS = this->nTCS_ * this->nBasis_;
  this->haveMO = true;
  if(this->molecule_->nAtoms() > 1) this->haveDensity = true;
  // Need to init random number otherwise Eigen is not truly random
  srand((unsigned int) time(0));
  *this->onePDMA_ = ComplexMatrix::Random(NTCSxNBASIS,NTCSxNBASIS);
  *this->onePDMA_ = this->onePDMA_->selfadjointView<Lower>();
  if(!this->isClosedShell && this->nTCS_ == 1){
    *this->onePDMB_ = ComplexMatrix::Random(NTCSxNBASIS,NTCSxNBASIS);
    *this->onePDMB_ = this->onePDMB_->selfadjointView<Lower>();
  }  
};

template<>
void SingleSlater<dcomplex>::READGuess(){
  if(getRank() == 0) {
    this->fileio_->out << "Reading SCF Density from disk" << endl;
    H5::DataSpace dataspace = this->fileio_->alphaSCFDen->getSpace();
    this->fileio_->alphaSCFDen->read(this->onePDMA_->data(),*(this->fileio_->complexType),dataspace,dataspace);
    this->fileio_->alphaMO->read(this->moA_->data(),*(this->fileio_->complexType),dataspace,dataspace);
    if(!this->isClosedShell && this->nTCS_ == 1){
      this->fileio_->betaSCFDen->read(this->onePDMB_->data(),*(this->fileio_->complexType),dataspace,dataspace);
      this->fileio_->betaMO->read(this->moB_->data(),*(this->fileio_->complexType),dataspace,dataspace);
    }
  }
  this->haveMO = true;
  if(this->molecule_->nAtoms() > 1) this->haveDensity = true;
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(this->onePDMA_->data(),
    this->nTCS_*this->nTCS_*this->nBasis_*this->nBasis_,MPI_C_DOUBLE_COMPLEX,0,
    MPI_COMM_WORLD);
  if(!this->isClosedShell && this->nTCS_ == 1)
    MPI_Bcast(this->onePDMB_->data(),
      this->nTCS_*this->nTCS_*this->nBasis_*this->nBasis_,MPI_C_DOUBLE_COMPLEX,0,
      MPI_COMM_WORLD);
#endif
  
}
}; //namespace ChronusQ
#endif
