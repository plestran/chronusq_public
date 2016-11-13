/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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
void SDResponse<T>::alloc(){
  this->checkValid();
  this->omega_       = std::unique_ptr<VectorXd>(new VectorXd(this->nSek_));
  this->transDen_    = std::unique_ptr<TMatrix>(new TMatrix(
                       this->nSingleDim_, this->nSek_));
  this->oscStrength_ = std::unique_ptr<RealMatrix>(new RealMatrix(
                       this->nSek_+1, this->nSek_+1));
  this->transDipole_ = std::unique_ptr<RealTensor3d>(new RealTensor3d(
                       this->nSek_+1, this->nSek_+1, 3));
}

template<typename T>
void SDResponse<T>::iniSDResponse( Molecule * molecule, BasisSet * basisSet, 
  MOIntegrals<T> * mointegrals, FileIO * fileio, Controls * controls, 
  SingleSlater<T> * singleSlater) {

  this->communicate(*molecule,*basisSet,*singleSlater,*mointegrals,*fileio,
                    *controls);
  this->initMeta();
  this->setNSek(this->controls_->SDNSek);
  this->setMeth(this->controls_->SDMethod);
  this->initMeth();
  this->alloc();

}
