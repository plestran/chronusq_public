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
#include <sdresponse.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
using ChronusQ::SDResponse;
using ChronusQ::SingleSlater;
using std::cout;
using std::setw;
//------------------------------//
// allocate memory for matrices //
//------------------------------//
void SDResponse::iniSDResponse( Molecule * molecule, BasisSet * basisSet, MOIntegrals<double> * mointegrals, 
                                FileIO * fileio, Controls * controls, SingleSlater<double> * singleSlater) {
  this->nBasis_         = basisSet->nBasis();
  this->nTCS_           = singleSlater->nTCS();
  this->molecule_       = molecule;
  this->basisSet_       = basisSet;
  this->fileio_         = fileio;
  this->controls_       = controls;
  this->mointegrals_    = mointegrals;
  this->singleSlater_   = singleSlater;
  this->aoERI_          = singleSlater->aointegrals()->aoERI_.get();
  this->elecDipole_     = singleSlater->aointegrals()->elecDipole_.get();
  this->Ref_            = singleSlater->Ref();
  this->haveDag_        = false;
  this->nOA_            = this->singleSlater_->nOccA();
  this->nOB_            = this->singleSlater_->nOccB();
  this->nVA_            = this->singleSlater_->nVirA();
  this->nVB_            = this->singleSlater_->nVirB();
  this->nOAVA_          = this->nOA_*this->nVA_;
  this->nOBVB_          = this->nOB_*this->nVB_;
  this->nOAVB_          = this->nOA_*this->nVB_;
  this->nOBVA_          = this->nOB_*this->nVA_;
  this->nVAVA_SLT_      = this->nVA_*(this->nVA_-1)/2;
  this->nVBVB_SLT_      = this->nVB_*(this->nVB_-1)/2;
  this->nVAVA_LT_       = this->nVA_*(this->nVA_+1)/2;
  this->nVAVA_          = this->nVA_*this->nVA_;
  this->nVBVB_          = this->nVB_*this->nVB_;
  this->nOAOA_SLT_      = this->nOA_*(this->nOA_-1)/2;
  this->nOBOB_SLT_      = this->nOB_*(this->nOB_-1)/2;
  this->nOAOA_LT_       = this->nOA_*(this->nOA_+1)/2;
  this->nOAOA_          = this->nOA_*this->nOA_;
  this->nOBOB_          = this->nOB_*this->nOB_;
  this->nVAVB_          = this->nVA_*this->nVB_;
  this->nOAOB_          = this->nOA_*this->nOB_;
  this->nO_             = this->nOA_ + this->nOB_;
  this->nV_             = this->nVA_ + this->nVB_;
  this->nOV_            = this->nO_  * this->nV_;
  this->nVV_SLT_        = this->nV_*(this->nV_-1)/2;
  this->nVV_LT_         = this->nV_*(this->nV_+1)/2;
  this->nVV_            = this->nV_*this->nV_;
  this->nOO_SLT_        = this->nO_*(this->nO_-1)/2;
  this->nOO_LT_         = this->nO_*(this->nO_+1)/2;
  this->nOO_            = this->nO_*this->nO_;
  this->setNSek(this->controls_->SDNSek);
  this->setMeth(this->controls_->SDMethod);
  this->omega_ = std::unique_ptr<VectorXd>(new VectorXd(this->nSek_));
  this->transDen_ = std::unique_ptr<RealCMMatrix>(new RealCMMatrix(this->nSingleDim_,this->nSek_));
  this->oscStrength_ = std::unique_ptr<RealMatrix>(new RealMatrix(this->nSek_+1,this->nSek_+1));
  this->transDipole_ = std::unique_ptr<RealTensor3d>(new RealTensor3d(this->nSek_+1,this->nSek_+1,3));
};
//dbwye
