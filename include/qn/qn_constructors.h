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

QuasiNewton2(){
  this->maxSubSpace_ = 0;
  this->qnObj_       = NULL;
  this->nMicroIter_  = 0;
  this->nMacroIter_  = 0;
  this->nTotalIter_  = 0;
  this->isConverged_ = 0;

  // Default Values
  this->maxMicroIter_     = 128;
  this->maxMacroIter_     = 20;
  this->residualTol_      = 5.0e-6;
  this->problemType_      = DIAGONALIZATION;
  this->matrixType_       = HERMETIAN;
  this->guessType_        = RESIDUAL_DAVIDSON;
  this->specialAlgorithm_ = NOT_SPECIAL;

  this->out_              = &std::cout;

};

QuasiNewton2(QNCallable<T> * obj) : QuasiNewton2(){
  this->qnObj_ = obj;
  this->maxSubSpace_ = std::min(250,obj->nSingleDim()/2);
};

QuasiNewton2(QNCallable<T> * obj, std::function<H5::DataSet*(const H5::CompType&,
  std::string&, std::vector<hsize_t>&)> fileFactory) : QuasiNewton2(obj){
  this->genScrFile_ = fileFactory;
/*
  std::vector<hsize_t> dims;
  dims.push_back(1);
  dims.push_back(1);

  std::string name = "Test";
  auto ptr = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,name,dims);
*/
};
