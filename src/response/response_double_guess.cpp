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
#include <response.h>

namespace ChronusQ{

template<>
void Response<double>::formDiagFOPPA(){

  for(auto iMat = 0; iMat != iMatIter_.size(); iMat++){
    auto diagDim = this->nMatDim_[iMat];
    if(!this->doTDA_) diagDim /= 2;
    this->rmDiag_.push_back(VectorXd(diagDim));

    if(this->iMatIter_[iMat] == FULL){
      for(auto i = 0, ia = 0; i < this->nOA_; i++)
      for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++){
        this->rmDiag_[iMat](ia) = 
          (*this->singleSlater_->epsA())(A) - (*this->singleSlater_->epsA())(i);
      }

      // FIXME: NEED TO DO UHF
      for(auto i = 0, ia = this->nOAVA_; i < this->nOB_; i++)
      for(auto a = 0, A = this->nOB_; a < this->nVB_; a++, A++, ia++){
        this->rmDiag_[iMat](ia) = 
          (*this->singleSlater_->epsA())(A) - (*this->singleSlater_->epsA())(i);
      }
    } //FULL
    else if(this->iMatIter_[iMat] == SINGLETS){
      for(auto i = 0, ia = 0; i < this->nOA_; i++)
      for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++){
        this->rmDiag_[iMat](ia) = 
          (*this->singleSlater_->epsA())(A) - (*this->singleSlater_->epsA())(i);
      }
    } // SINGLETS
    else if(this->iMatIter_[iMat] == TRIPLETS){
      for(auto i = 0, ia = 0; i < this->nOA_; i++)
      for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++){
        this->rmDiag_[iMat](ia) = 
          (*this->singleSlater_->epsA())(A) - (*this->singleSlater_->epsA())(i);
      }
    } // TRIPLETS
    

  } // loop over mats
}; // formDiagFOPPA

template<>
void Response<double>::formGuessFOPPA(){

  this->scratchFile_ = this->fileio_->scr.get();
  for(auto iMat = 0; iMat != iMatIter_.size(); iMat++){
    auto diagDim = this->nMatDim_[iMat];
    if(!this->doTDA_) diagDim /= 2;

    std::string name = "ResponseGuess"+std::to_string(iMat);
    hsize_t dims[] = {this->nMatDim_[iMat],this->nSek_};
    H5::DataSpace dataspace(2,dims);
    this->fileio_->scratchPartitions.push_back(
      FileIO::ScratchPartition(name,dataspace,*this->scratchFile_)
    );
    this->guessFiles_.push_back(&this->fileio_->scratchPartitions.back().data);

    std::vector<int> indx(diagDim,0);
    for(auto i = 0; i < indx.size(); i++){
      indx[i] = i;
    }

    std::sort(indx.begin(),indx.end(),
      [&](const int& a, const int& b){
        return this->rmDiag_[iMat](a) < this->rmDiag_[iMat](b);
      }
    );

    prettyPrint(cout,this->rmDiag_[iMat],"diag");

    for(auto ia = 0; ia < this->nSek_; ia++){
      cout << indx[ia] << " " << this->rmDiag_[iMat](indx[ia]) << endl;
    }
    
  }
     
}; // formGuessFOPPA
} // namespace ChronusQ
