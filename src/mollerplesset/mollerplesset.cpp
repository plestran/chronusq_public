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
#include <mollerplesset.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
using ChronusQ::MollerPlesset;
using ChronusQ::SingleSlater;


void MollerPlesset::iniMollerPlesset( Molecule * molecule, BasisSet * basisSet, MOIntegrals<double> * mointegrals, 
                                FileIO * fileio, Controls * controls, SingleSlater<double> * singleSlater) {
  this->nBasis_         = basisSet->nBasis();
  this->nTCS_           = singleSlater->nTCS();
  this->molecule_       = molecule;
  this->basisSet_       = basisSet;
  this->fileio_         = fileio;
  this->controls_       = controls;
  this->mointegrals_    = mointegrals;
  this->singleSlater_   = singleSlater;
  this->Ref_            = singleSlater->Ref();
  this->nOA_            = this->singleSlater_->nOccA();
  this->nOB_            = this->singleSlater_->nOccB();
  this->nVA_            = this->singleSlater_->nVirA();
  this->nVB_            = this->singleSlater_->nVirB();
  this->nO_             = this->nOA_ + this->nOB_;
  this->nV_             = this->nVA_ + this->nVB_;
}

double MollerPlesset::MP2(){
  this->mointegrals_->formIAJB(false);
  double EMP2 = 0.0;
  double fact = 1.0;
  if(!this->singleSlater_->isClosedShell) fact = 0.5;

  if(this->Ref_ == SingleSlater<double>::TCS) {
    for(auto i = 0; i < this->nO_; i++)
    for(auto j = 0; j < this->nO_; j++)
    for(auto a = 0; a < this->nV_; a++)
    for(auto b = 0; b < this->nV_; b++){
      EMP2 += (this->mointegrals_->IAJB(i,a,j,b)-this->mointegrals_->IAJB(i,b,j,a))*
              (this->mointegrals_->IAJB(i,a,j,b)-this->mointegrals_->IAJB(i,b,j,a))/
              ( (*this->singleSlater_->epsA())(i) + (*this->singleSlater_->epsA())(j)
              - (*this->singleSlater_->epsA())(a + this->nO_)
              - (*this->singleSlater_->epsA())(b + this->nO_));
    }
    EMP2 *= 0.25;
  } else {
    for(auto i = 0; i < this->nOA_; i++)
    for(auto j = 0; j < this->nOA_; j++)
    for(auto a = 0; a < this->nVA_; a++)
    for(auto b = 0; b < this->nVA_; b++){
      EMP2 += fact*this->mointegrals_->IAJB(i,a,j,b,"AAAA")*
                  (this->mointegrals_->IAJB(i,a,j,b,"AAAA") - 
                   this->mointegrals_->IAJB(i,b,j,a,"AAAA")) /
              ( (*this->singleSlater_->epsA())(i) + (*this->singleSlater_->epsA())(j)
              - (*this->singleSlater_->epsA())(a + this->nOA_)
              - (*this->singleSlater_->epsA())(b + this->nOA_));
    }
 
    if(!this->singleSlater_->isClosedShell) {
      for(auto i = 0; i < this->nOB_; i++)
      for(auto j = 0; j < this->nOB_; j++)
      for(auto a = 0; a < this->nVB_; a++)
      for(auto b = 0; b < this->nVB_; b++){
        EMP2 += 0.5*this->mointegrals_->IAJB(i,a,j,b,"BBBB")*
                (this->mointegrals_->IAJB(i,a,j,b,"BBBB") - 
                 this->mointegrals_->IAJB(i,b,j,a,"BBBB")) /
                ( (*this->singleSlater_->epsB())(i) + (*this->singleSlater_->epsB())(j)
                - (*this->singleSlater_->epsB())(a + this->nOB_)
                - (*this->singleSlater_->epsB())(b + this->nOB_));
      }
    }
 
    for(auto i = 0; i < this->nOA_; i++)
    for(auto j = 0; j < this->nOB_; j++)
    for(auto a = 0; a < this->nVA_; a++)
    for(auto b = 0; b < this->nVB_; b++){
      if(this->singleSlater_->isClosedShell)
        EMP2 += this->mointegrals_->IAJB(i,a,j,b,"AABB")*
                this->mointegrals_->IAJB(i,a,j,b,"AABB")/
                ( (*this->singleSlater_->epsA())(i) + (*this->singleSlater_->epsA())(j)
                - (*this->singleSlater_->epsA())(a + this->nOA_)
                - (*this->singleSlater_->epsA())(b + this->nOA_));
      else
        EMP2 += this->mointegrals_->IAJB(i,a,j,b,"AABB")*
                this->mointegrals_->IAJB(i,a,j,b,"AABB")/
                ( (*this->singleSlater_->epsA())(i) + (*this->singleSlater_->epsB())(j)
                - (*this->singleSlater_->epsA())(a + this->nOA_)
                - (*this->singleSlater_->epsB())(b + this->nOB_));
    }
  }
  cout << "EMP2 " << EMP2 +this->singleSlater_->totalEnergy<< endl;
  cout << "EMP2 " << EMP2 << endl;
  return EMP2;
}
