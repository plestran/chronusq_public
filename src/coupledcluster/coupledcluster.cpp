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
#include <coupledcluster.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
using ChronusQ::SingleSlater;
using ChronusQ::CoupledCluster;

void CoupledCluster::iniCoupledCluster( Molecule * molecule, BasisSet * basisSet, MOIntegrals<double> * mointegrals, 
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

double CoupledCluster::CCSD(){
  this->mointegrals_->formIJAB(true);
  this->mointegrals_->formIJKL(true);
  this->mointegrals_->formABCD(true);
  this->mointegrals_->formIAJB(true);
  this->mointegrals_->formIABC(true);
  this->mointegrals_->formIJKA(true);
  
  double ECorr = 0.0;

  // Intermediates
  RealTensor2d Hba, Hbj, Hij, Gca, Gik; 
  RealTensor4d Aijkl, Babcd, Hicak;
  // Amplitudes
  RealTensor2d Tia, Wia;
  RealTensor4d Tiajb, Tau, Wiajb;

  auto NTCSxNBASIS = this->nTCS_*this->nBasis_;

  if(this->Ref_ == SingleSlater<double>::TCS){
    // Intermediates
    Hba     = RealTensor2d(this->nV_,this->nV_);
    Hbj     = RealTensor2d(this->nV_,this->nO_);
    Hij     = RealTensor2d(this->nO_,this->nO_);
    Gik     = RealTensor2d(this->nO_,this->nO_);
    Gca     = RealTensor2d(this->nV_,this->nV_);
    Aijkl   = RealTensor4d(this->nO_,this->nO_,this->nO_,this->nO_);
    Babcd   = RealTensor4d(this->nV_,this->nV_,this->nV_,this->nV_);
    Hicak   = RealTensor4d(this->nO_,this->nV_,this->nV_,this->nO_);
    // Amplitudes 
    Tia     = RealTensor2d(this->nO_,this->nV_);
    Wia     = RealTensor2d(this->nO_,this->nV_);
    Tiajb   = RealTensor4d(this->nO_,this->nV_,this->nO_,this->nV_);
    Tau     = RealTensor4d(this->nO_,this->nV_,this->nO_,this->nV_);
    Wiajb   = RealTensor4d(this->nO_,this->nV_,this->nO_,this->nV_);
  } else {
  }


};
