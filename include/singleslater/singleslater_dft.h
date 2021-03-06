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
/*
template<typename T>
void SingleSlater<T>::addSlater(double x = 1.0){
  this->dftFunctionals_.emplace_back(std::make_shared<SlaterExchange>(x));
};
template<typename T>
void SingleSlater<T>::addB88(double x = 1.0){
  this->dftFunctionals_.emplace_back(std::make_shared<BEightEight>(x));
  this->isGGA = true;
};
template<typename T>
void SingleSlater<T>::addLYP(double x = 1.0){
  this->dftFunctionals_.emplace_back(std::make_shared<lyp>(x)); 
  this->isGGA = true;
};
template<typename T>
void SingleSlater<T>::addVWN5(double x = 1.0){
  this->dftFunctionals_.emplace_back(std::make_shared<VWNV>(x)); 
};
template<typename T>
void SingleSlater<T>::addVWN3(double x = 1.0){
  this->dftFunctionals_.emplace_back(std::make_shared<VWNIII>(x)); 
};
*/

/*** Instantiation of DFT setup functions ***/

template <typename T>
void SingleSlater<T>::addDFTFunctional(const std::string &func, double x){
  if(!func.compare("Slater"))
    this->dftFunctionals_.emplace_back(std::make_shared<SlaterExchange>(x));
  else if(!func.compare("B88"))
    this->dftFunctionals_.emplace_back(std::make_shared<BEightEight>(x));
  else if(!func.compare("PBE"))
    this->dftFunctionals_.emplace_back(std::make_shared<PBE>(x));
  else if(!func.compare("LYP"))
    this->dftFunctionals_.emplace_back(std::make_shared<lyp>(x));
  else if(!func.compare("VWN5"))
    this->dftFunctionals_.emplace_back(std::make_shared<VWNV>(x));
  else if(!func.compare("VWN3"))
    this->dftFunctionals_.emplace_back(std::make_shared<VWNIII>(x));

  if(!func.compare("B88") or !func.compare("LYP") or !func.compare("PBE")) this->isGGA = true;
};


// EXCHANGE ONLY
template<typename T>
void SingleSlater<T>::createSlater(){
  addDFTFunctional("Slater");
  this->xHF_ = 0.0;
  this->DFTKernel_ = SLATER;
};
template<typename T>
void SingleSlater<T>::createB88(){
  addDFTFunctional("Slater");
  addDFTFunctional("B88");
  this->xHF_ = 0.0;
  this->DFTKernel_ = B88;
};

template<typename T>
void SingleSlater<T>::createPBE(){
//  addDFTFunctional("Slater");
  addDFTFunctional("PBE");
  addDFTFunctional("LYP");
  this->xHF_ = 0.0;
  this->DFTKernel_ = pbe;
};
  
  
// Pure Functionals
template<typename T>
void SingleSlater<T>::createLSDA(){
  addDFTFunctional("Slater");
  addDFTFunctional("VWN3");
  this->xHF_ = 0.0;
  this->DFTKernel_ = LSDA;
};
template<typename T>
void SingleSlater<T>::createSVWN5(){
  addDFTFunctional("Slater");
  addDFTFunctional("VWN5");
  this->xHF_ = 0.0;
  this->DFTKernel_ = SVWN5;
};

template<typename T>
void SingleSlater<T>::createBLYP(){
  addDFTFunctional("Slater");
  addDFTFunctional("B88");
  addDFTFunctional("LYP");
  this->xHF_ = 0.0;
  this->DFTKernel_ = BLYP;
};


// Hybrid Functionals
template<typename T>
void SingleSlater<T>::createB3LYP(){
  addDFTFunctional("Slater",0.8);
  addDFTFunctional("B88",0.72);
  addDFTFunctional("LYP",0.81);
  addDFTFunctional("VWN3",0.19);
  this->xHF_ = -0.5 * 0.2;
  this->DFTKernel_ = B3LYP;
};

template<typename T>
void SingleSlater<T>::createBHandH(){
  addDFTFunctional("Slater",0.5);
  addDFTFunctional("LYP");
  this->xHF_ = -0.5 * 0.5;
  this->DFTKernel_ = BHandH;
};
