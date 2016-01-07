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
#include <qn.h>

namespace ChronusQ {
  template<>
  void QuasiNewton2<double>::iniScratchFiles(){
    (*this->out_) << "Initializing Files for QuasiNewton Calculation" << endl;

    auto N      = this->qnObj_->nSingleDim();
    std::vector<hsize_t> dims;
    dims.push_back(this->maxSubSpace_);
    dims.push_back(N);

    time_t t; time(&t); // Time flag for files to be written

    std::string TRName = "Trial Vectors (Right) "    + std::to_string(t);
    std::string SRName = "Sigma Vectors (Right) "    + std::to_string(t);
    std::string PRName = "Rho Vectors (Right) "      + std::to_string(t);
    std::string RRName = "Residual Vectors (Right) " + std::to_string(t);
    std::string URName = "Solution Vectors (Right) " + std::to_string(t);

    this->TRFile_     = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,TRName,dims);
  //this->SigmaRFile_ = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,SRName,dims);
  //this->ResRFile_   = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,RRName,dims);
  //this->URFile_     = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,URName,dims);

  //if(this->matrixType_ == HERMETIAN_GEP)
  //  this->RhoRFile_  = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,PRName,dims);


    if(this->qnObj_->needsLeft()) {
      std::string TLName = "Trial Vectors (Left) "    + std::to_string(t);
      std::string SLName = "Sigma Vectors (Left) "    + std::to_string(t);
      std::string PLName = "Rho Vectors (Left) "      + std::to_string(t);
      std::string RLName = "Residual Vectors (Left) " + std::to_string(t);
      std::string ULName = "Solution Vectors (Left) " + std::to_string(t);
     
      this->TLFile_     = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,TLName,dims);
    //this->SigmaLFile_ = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,SLName,dims);
    //this->ResLFile_   = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,RLName,dims);
    //this->ULFile_     = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,ULName,dims);
     
    //if(this->matrixType_ == HERMETIAN_GEP)
    //  this->RhoLFile_  = this->genScrFile_(H5::PredType::NATIVE_DOUBLE,PLName,dims);
    }
  };

}; // namespace ChronusQ
