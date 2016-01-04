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
    else if(this->iMatIter_[iMat] == SINGLETS || 
            this->iMatIter_[iMat] == TRIPLETS){

      for(auto i = 0, ia = 0; i < this->nOA_; i++)
      for(auto a = 0, A = this->nOA_; a < this->nVA_; a++, A++, ia++){
        this->rmDiag_[iMat](ia) = 
          (*this->singleSlater_->epsA())(A) - (*this->singleSlater_->epsA())(i);
      }

    } // SINGLETS or TRIPLETS
    
    if(!this->doTDA_)
      this->rmDiag_[iMat].block(diagDim/2,0,diagDim/2,1) = 
      this->rmDiag_[iMat].block(0        ,0,diagDim/2,1);

  } // loop over mats
}; // formDiagFOPPA

template<>
void Response<double>::formDiagPPRPA(){
  cout << "HERE 2" << endl;
  for(auto iMat = 0; iMat != iMatIter_.size(); iMat++){
    auto diagDim = this->nMatDim_[iMat];
    this->rmDiag_.push_back(VectorXd(diagDim));
    cout << "HERE 1" << endl;

    if(this->iMatIter_[iMat] == AA_PPRPA){ 
      auto iOff = this->nVAVA_SLT_;
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < a        ;  b++, ab++){
        this->rmDiag_[iMat](ab) = 
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      } // loop AB

      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j < i         ; j++, ij++){
        this->rmDiag_[iMat](ij) = - 
          ((*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_); 
      } // loop IJ
    } // this->iMatIter_[iMat] == AA_PPRPA
    else if(this->iMatIter_[iMat] == AB_PPRPA){ 
      auto iOff = this->nVAVB_;
      // FIXME: Need to generalize to UHF
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < this->nVB_; b++, ab++){
        this->rmDiag_[iMat](ab) = 
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOB_) - 2*this->rMu_; 
      } // loop AB

      // FIXME: Need to generalize to UHF
      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j < this->nOB_; j++, ij++){
        this->rmDiag_[iMat](ij) = - 
          ((*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_); 
      } // loop IJ
    } // this->iMatIter_[iMat] == AB_PPRPA
    else if(this->iMatIter_[iMat] == BB_PPRPA){ 
    // ************************** //
    // ** BETA-BETA BLOCKS NYI ** //
    // ************************** //
    } // this->iMatIter_[iMat] == BB_PPRPA
    else if(this->iMatIter_[iMat] == AAA_PPTDA){
      cout << "HERE AAA" << endl;
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < a        ;  b++, ab++){
        this->rmDiag_[iMat](ab) = 
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      } // loop AB
      cout << "HERE  AAA" << endl;
    } // this->iMatIter_[iMat] == AAA_PPTDA
    else if(this->iMatIter_[iMat] == AAB_PPTDA){
      cout << "HERE AAB" << endl;
      // FIXME: Need to generalize to UHF
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < this->nVB_; b++, ab++){
        this->rmDiag_[iMat](ab) = 
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOB_) - 2*this->rMu_; 
      } // loop AB
      cout << "HERE AAB" << endl;
    } // this->iMatIter_[iMat] == AAB_PPTDA
    else if(this->iMatIter_[iMat] == ABB_PPTDA){
    // ************************** //
    // ** BETA-BETA BLOCKS NYI ** //
    // ************************** //
    } // this->iMatIter_[iMat] == ABB_PPTDA
    else if(this->iMatIter_[iMat] == CAA_PPTDA){
      for(auto i = 0, ij = 0; i < this->nOA_; i++      )
      for(auto j = 0        ; j < i         ; j++, ij++){
        this->rmDiag_[iMat](ij) = - 
          ((*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_); 
      } // loop IJ
    } // this->iMatIter_[iMat] == CAA_PPTDA
    else if(this->iMatIter_[iMat] == CAB_PPTDA){
      // FIXME: Need to generalize to UHF
      for(auto i = 0, ij = 0; i < this->nOA_; i++      )
      for(auto j = 0        ; j < this->nOB_; j++, ij++){
        this->rmDiag_[iMat](ij) = - 
          ((*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_); 
      } // loop IJ
    } // this->iMatIter_[iMat] == CAB_PPTDA
    else if(this->iMatIter_[iMat] == CBB_PPTDA){
    // ************************** //
    // ** BETA-BETA BLOCKS NYI ** //
    // ************************** //
    } // this->iMatIter_[iMat] == CBB_PPTDA
    else if(this->iMatIter_[iMat] == PPRPA_SINGLETS){
      auto iOff = this->nVAVA_LT_;
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b <= a       ;  b++, ab++){
        this->rmDiag_[iMat](ab) = 
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      } // loop AB

      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j <= i        ; j++, ij++){
        this->rmDiag_[iMat](ij) = - 
          ((*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_); 
      } // loop IJ
    } // this->iMatIter_[iMat] == PPRPA_SINGLETS
    else if(this->iMatIter_[iMat] == A_PPTDA_SINGLETS){
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b <= a       ;  b++, ab++){
        this->rmDiag_[iMat](ab) = 
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      } // loop AB
    } // this->iMatIter_[iMat] == A_PPTDA_SINGLETS
    else if(this->iMatIter_[iMat] == C_PPTDA_SINGLETS){
      for(auto i = 0, ij = 0; i < this->nOA_; i++      )
      for(auto j = 0        ; j <= i        ; j++, ij++){
        this->rmDiag_[iMat](ij) = - 
          ((*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_); 
      } // loop IJ
    } // this->iMatIter_[iMat] == C_PPTDA_SINGLETS
    else if(this->iMatIter_[iMat] == PPRPA_TRIPLETS){
      auto iOff = this->nVAVA_SLT_;
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < a        ;  b++, ab++){
        this->rmDiag_[iMat](ab) = 
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      } // loop AB

      for(auto i = 0, ij = iOff; i < this->nOA_; i++      )
      for(auto j = 0           ; j < i         ; j++, ij++){
        this->rmDiag_[iMat](ij) = - 
          ((*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_); 
      } // loop IJ
    } // this->iMatIter_[iMat] == PPRPA_TRIPLETS
    else if(this->iMatIter_[iMat] == A_PPTDA_TRIPLETS){
      for(auto a = 0, ab = 0; a < this->nVA_; a++      )
      for(auto b = 0        ; b < a        ;  b++, ab++){
        this->rmDiag_[iMat](ab) = 
          (*this->singleSlater_->epsA())(a+this->nOA_) + 
          (*this->singleSlater_->epsA())(b+this->nOA_) - 2*this->rMu_; 
      } // loop AB
    } // this->iMatIter_[iMat] == A_PPTDA_TRIPLETS
    else if(this->iMatIter_[iMat] == C_PPTDA_TRIPLETS){
      for(auto i = 0, ij = 0; i < this->nOA_; i++      )
      for(auto j = 0        ; j < i         ; j++, ij++){
        this->rmDiag_[iMat](ij) = - 
          ((*this->singleSlater_->epsA())(i) + 
          (*this->singleSlater_->epsA())(j) - 2*this->rMu_); 
      } // loop IJ
    } // this->iMatIter_[iMat] == C_PPTDA_TRIPLETS

  } // loop over mats
}; // formDiagPPRPA

template<>
void Response<double>::formGuess(){

  this->scratchFile_ = this->fileio_->scr.get();
  for(auto iMat = 0; iMat != iMatIter_.size(); iMat++){
    auto diagDim = this->nMatDim_[iMat];

    // Determine the dimention of diagonal to sort
    if(!this->doTDA_) {
      if(this->iClass_ == FOPPA)
        diagDim /= 2;
      else if(this->iClass_ == PPPA){
      // This will only do the pp additions part of the
      // pp-RPA guess. FIXME: Need an option to look for
      // the hh additions
        if(this->iMatIter_[iMat] == AA_PPRPA ||
           this->iMatIter_[iMat] == PPRPA_TRIPLETS)
          diagDim = this->nVAVA_SLT_; 
        else if(this->iMatIter_[iMat] == AB_PPRPA)
          diagDim = this->nVAVB_;
        else if(this->iMatIter_[iMat] == PPRPA_TRIPLETS)
          diagDim = this->nVAVA_LT_;
        else if(this->iMatIter_[iMat] == BB_PPRPA)
          diagDim = this->nVBVB_SLT_;
      }
    }

    // Create a new scratch file to store guess vectors
    std::string name = "ResponseGuess"+std::to_string(iMat);
//  hsize_t dims[] = {this->nMatDim_[iMat],this->nSek_};
     
//  HDF5 stores things in RowMajor...
    hsize_t dims[] = {this->nSek_,this->nMatDim_[iMat]};
    H5::DataSpace dataspace(2,dims);
    this->fileio_->scratchPartitions.push_back(
      FileIO::ScratchPartition(name,dataspace,*this->scratchFile_)
    );

    // Keep a pointer to the created files
    this->guessFiles_.push_back(&this->fileio_->scratchPartitions.back().data);

    // Initialize an index vector with increasing ints
    std::vector<int> indx(diagDim,0);
    for(auto i = 0; i < indx.size(); i++){
      indx[i] = i;
    }

    // Sort the index vector based on the diagonal of the RM
    // (uses Lambda expression)
    std::sort(indx.begin(),indx.end(),
      [&](const int& a, const int& b){
        return this->rmDiag_[iMat](a) < this->rmDiag_[iMat](b);
      }
    );

/*
    prettyPrint(cout,this->rmDiag_[iMat],"diag");

    for(auto ia = 0; ia < this->nSek_; ia++){
      cout << indx[ia] << " " << this->rmDiag_[iMat](indx[ia]) << endl;
    }
*/

    this->guessFiles_[iMat]->getSpace();
    // Write the guess vectors to disk
    for(auto iSek = 0; iSek < this->nSek_; iSek++) {
      double one = 1.0;

//    hsize_t offset[] = {indx[iSek],iSek};
  
//    HDF5 Stores things RowMajor...
      hsize_t offset[] = {iSek,indx[iSek]};
      hsize_t count[]  = {1,1};
      hsize_t stride[] = {1,1};
      hsize_t block[]  = {1,1};
      hsize_t subDim[] = {1,1};

      H5::DataSpace memSpace(2,subDim,NULL);
      H5::DataSpace subDataSpace = this->guessFiles_[iMat]->getSpace();
      subDataSpace.selectHyperslab(H5S_SELECT_SET,count,offset,stride,block);
      this->guessFiles_[iMat]->write(
        &one,H5::PredType::NATIVE_DOUBLE,memSpace,subDataSpace
      );
      
    }
    
  }
     
}; // formGuessFOPPA
} // namespace ChronusQ
