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
#include <qn.h>

namespace ChronusQ {
/*
  template<>
  void QuasiNewton2<double>::iniScratchFiles(){
    (*this->out_) << "Initializing Files for QuasiNewton Calculation" << endl;

    auto N      = this->qnObj_->nSingleDim();
    std::vector<hsize_t> dims;
    dims.push_back(this->maxSubSpace_);
    dims.push_back(N);

    std::vector<hsize_t> dimsSmall;
    dims.push_back(this->maxSubSpace_);
    dims.push_back(this->maxSubSpace_);

    time_t t; time(&t); // Time flag for files to be written

    std::string TRName = "Trial Vectors (Right) "    + std::to_string(t);
    std::string SRName = "Sigma Vectors (Right) "    + std::to_string(t);
    std::string PRName = "Rho Vectors (Right) "      + std::to_string(t);
    std::string RRName = "Residual Vectors (Right) " + std::to_string(t);
    std::string URName = "Solution Vectors (Right) " + std::to_string(t);
    std::string ASName = "A Super Matrix "           + std::to_string(t);
    std::string SSName = "S Super Matrix "           + std::to_string(t);

    this->TRFile_     = this->genScrFile_(H5PredType<double>(),TRName,dims);
  //this->SigmaRFile_ = this->genScrFile_(H5PredType<double>(),SRName,dims);
  //this->ResRFile_   = this->genScrFile_(H5PredType<double>(),RRName,dims);
  //this->URFile_     = this->genScrFile_(H5PredType<double>(),URName,dims);

  //if(this->matrixType_ == HERMETIAN_GEP)
  //  this->RhoRFile_  = this->genScrFile_(H5PredType<double>(),PRName,dims);

 

    if(this->qnObj_->needsLeft()) {
      std::string TLName = "Trial Vectors (Left) "    + std::to_string(t);
      std::string SLName = "Sigma Vectors (Left) "    + std::to_string(t);
      std::string PLName = "Rho Vectors (Left) "      + std::to_string(t);
      std::string RLName = "Residual Vectors (Left) " + std::to_string(t);
      std::string ULName = "Solution Vectors (Left) " + std::to_string(t);
     
      this->TLFile_     = this->genScrFile_(H5PredType<double>(),TLName,dims);
    //this->SigmaLFile_ = this->genScrFile_(H5PredType<double>(),SLName,dims);
    //this->ResLFile_   = this->genScrFile_(H5PredType<double>(),RLName,dims);
    //this->ULFile_     = this->genScrFile_(H5PredType<double>(),ULName,dims);
     
    //if(this->matrixType_ == HERMETIAN_GEP)
    //  this->RhoLFile_  = this->genScrFile_(H5PredType<double>(),PLName,dims);
    }

    if(this->qnObj_->specialAlgorithm_ == SYMMETRIZED_TRIAL) {
      this->ASuperFile_ = this->genScrFile_(H5PredType<double>(),ASName,
        dimsSmall);
      this->SSuperFile_ = this->genScrFile_(H5PredType<double>(),SSName,
        dimsSmall);
    }
  }; // QuasiNewton2<double>::iniScratchFiles

  template<>
  void QuasiNewton2<double>::writeTrialVectors(const int NTrial){
    auto N      = this->qnObj_->nSingleDim();
    hsize_t offset[] = {0,0};
    hsize_t stride[] = {1,1};
    hsize_t block[]  = {1,1};
 
    hsize_t subDim[] = {NTrial,N};
    hsize_t count[]  = {NTrial,N};

    H5::DataSpace memSpace(2,subDim,NULL);
    H5::DataSpace subDataSpace = this->TRFile_->getSpace();
    subDataSpace.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
    this->TRFile_->write(this->TRMem_,H5PredType<double>(),memSpace,
      subDataSpace);
    if(this->qnObj_->needsLeft())
      this->TLFile_->write(this->TLMem_,H5PredType<double>(),memSpace,
        subDataSpace);
  
  }; // QuasiNewton2<double>::writeTrialVectors

  template<>
  void QuasiNewton2<double>::readTrialVectors(const int NTrial){
    auto N      = this->qnObj_->nSingleDim();
    hsize_t offset[] = {0,0};
    hsize_t stride[] = {1,1};
    hsize_t block[]  = {1,1};
 
    hsize_t subDim[] = {NTrial,N};
    hsize_t count[]  = {NTrial,N};

    H5::DataSpace memSpace(2,subDim,NULL);
    H5::DataSpace subDataSpace = this->TRFile_->getSpace();
    subDataSpace.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
    this->TRFile_->read(this->TRMem_,H5PredType<double>(),memSpace,
      subDataSpace);
    if(this->qnObj_->needsLeft())
      this->TLFile_->read(this->TLMem_,H5PredType<double>(),memSpace,
        subDataSpace);
  
  }; // QuasiNewton2<double>::readTrialVectors
*/

}; // namespace ChronusQ
