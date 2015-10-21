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
#ifndef INCLUDED_SDRESPONSE
#define INCLUDED_SDRESPONSE
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <controls.h>
#include <mointegrals.h>
#include <singleslater.h>
#include <basisset.h>

/****************************/
/* Error Messages 5000-5999 */
/****************************/

namespace ChronusQ {
template<typename T>
class SDResponse {
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TCMMatrix;
  typedef Eigen::Matrix<T,Dynamic,Dynamic,RowMajor> TMatrix;
  typedef Tensor<T,Range3d> TTensor3d;
  typedef Tensor<T,Range4d> TTensor4d;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TVec> TVecMap;
  typedef Eigen::Map<TCMMatrix> TCMMap;

  int       nBasis_;
  int       nTCS_;
  int       Ref_;
  int       nSek_;
  int       nGuess_;
  int       iMeth_;
  int       iPPRPA_;
  bool      haveDag_;
  int       nOA_;
  int       nVA_;
  int       nOB_;
  int       nVB_;
  int       nOAVA_;      // NOA * NVA
  int       nOBVB_;      // NOB * NVB
  int       nOAVB_;      // NOA * NVB
  int       nOBVA_;      // NOB * NVA
  int       nVAVA_SLT_;  // NVA * (NVA - 1) / 2
  int       nVBVB_SLT_;  // NVB * (NVB - 1) / 2
  int       nVAVA_LT_;   // NVA * (NVA + 1) / 2
  int       nVAVA_;      // NVA * NVA
  int       nVBVB_;      // NVB * NVB
  int       nOAOA_SLT_;  // NOA * (NOA - 1) / 2
  int       nOBOB_SLT_;  // NOB * (NOB - 1) / 2
  int       nOAOA_LT_;   // NOA * (NOA + 1) / 2
  int       nOAOA_;      // NOA * NOA
  int       nOBOB_;      // NOB * NOB
  int       nVAVB_;      // NVA * NVB
  int       nOAOB_;      // NOA * NOB
  int       nO_;         // NOA + NOB
  int       nV_;         // NVA + NVB
  int       nOV_;        // NO * NV
  int       nVV_SLT_;    // NV * (NV - 1) / 2
  int       nVV_LT_;     // NV * (NV + 1) / 2
  int       nVV_;        // NV * NV
  int       nOO_SLT_;    // NO * (NO - 1) / 2
  int       nOO_LT_;     // NO * (NO + 1) / 2
  int       nOO_;        // NO * NO

  int       nSingleDim_; // Single dimension of response matrix
  double    rMu_;

  std::unique_ptr<TCMMatrix>  transDen_;
  std::unique_ptr<RealMatrix> oscStrength_;
  std::unique_ptr<VectorXd>   omega_;
  std::unique_ptr<RealTensor3d>  transDipole_;
  BasisSet        * basisSet_;
  Molecule        * molecule_;
  FileIO          * fileio_;
  Controls        * controls_;
  MOIntegrals<T>  * mointegrals_;
  SingleSlater<T> * singleSlater_;
  RealTensor4d       * aoERI_;
  RealTensor3d       * elecDipole_;

  std::unique_ptr<RealCMMatrix> rmDiag_;
  std::unique_ptr<TMatrix>  davGuess_;

  inline void checkWorkers(){
    if(this->fileio_  == NULL) 
      CErr("Fatal: Must initialize SDResponse with FileIO Object");
    if(this->basisSet_ == NULL) 
      CErr("Fatal: Must initialize SDResponse with BasisSet Object",
           this->fileio_->out);
    if(this->molecule_ == NULL) 
      CErr("Fatal: Must initialize SDResponse with Molecule Object",
           this->fileio_->out);
    if(this->singleSlater_ == NULL) 
      CErr("Fatal: Must initialize SDResponse with SingleSlater Object",
           this->fileio_->out);
    if(this->mointegrals_ == NULL) 
      CErr("Fatal: Must initialize SDResponse with MOIntegrals Object",
           this->fileio_->out);
    if(this->controls_ == NULL) 
      CErr("Fatal: Must initialize SDResponse with Controls Object",
           this->fileio_->out);
    
  }

public:
  enum{
    __invalid,
    CIS,
    RPA,
    PPRPA,
    PPATDA,
    PPCTDA,
    STAB
  };
 
  // constructor & destructor
  SDResponse(){
    this->nBasis_     = 0;
    this->nTCS_       = 0;
    this->Ref_        = 0;
    this->nSek_       = 0;
    this->nGuess_     = 0;
    this->iMeth_      = 0;
    this->iPPRPA_     = 0;
    this->nOA_        = 0;
    this->nVA_        = 0;
    this->nOB_        = 0;
    this->nVB_        = 0;
    this->nOAVA_      = 0;      
    this->nOBVB_      = 0;      
    this->nOAVB_      = 0;      
    this->nOBVA_      = 0;      
    this->nVAVA_SLT_  = 0;  
    this->nVBVB_SLT_  = 0;  
    this->nVAVA_LT_   = 0;   
    this->nVAVA_      = 0;      
    this->nVBVB_      = 0;      
    this->nOAOA_SLT_  = 0;  
    this->nOBOB_SLT_  = 0;  
    this->nOAOA_LT_   = 0;   
    this->nOAOA_      = 0;      
    this->nOBOB_      = 0;      
    this->nVAVB_      = 0;      
    this->nOAOB_      = 0;      
    this->nO_         = 0;         
    this->nV_         = 0;         
    this->nOV_        = 0;        
    this->nVV_SLT_    = 0;    
    this->nVV_LT_     = 0;     
    this->nVV_        = 0;        
    this->nOO_SLT_    = 0;    
    this->nOO_LT_     = 0;     
    this->nOO_        = 0;        
    this->nSingleDim_ = 0; 

    this->haveDag_    = false;

    this->rMu_        = 0;

    this->basisSet_     = NULL;
    this->molecule_     = NULL;
    this->fileio_       = NULL;
    this->controls_     = NULL;
    this->mointegrals_  = NULL;
    this->singleSlater_ = NULL;
    this->aoERI_        = NULL;
    this->elecDipole_   = NULL;

    this->transDen_    = nullptr;
    this->oscStrength_ = nullptr;
    this->omega_       = nullptr;
    this->transDipole_ = nullptr;
    this->rmDiag_      = nullptr;
    this->davGuess_    = nullptr;
  };
  ~SDResponse() {;};
  // pseudo-constructor
  void iniSDResponse(Molecule *,BasisSet *,MOIntegrals<T> *,FileIO *,
                     Controls *, SingleSlater<T> *);

  inline void communicate(Molecule &mol, BasisSet &basis, SingleSlater<T> &ss,
    MOIntegrals<T> &moints, FileIO &fileio, Controls &controls) {

    this->molecule_     = &mol;
    this->basisSet_     = &basis;
    this->singleSlater_ = &ss;
    this->mointegrals_  = &moints;
    this->fileio_       = &fileio;
    this->controls_     = &controls;

  }

  inline void initMeta() {
    this->checkWorkers();

    this->nBasis_         = this->basisSet_->nBasis();
    this->nTCS_           = this->singleSlater_->nTCS();
    this->Ref_            = this->singleSlater_->Ref();
    this->nOA_            = this->singleSlater_->nOccA();
    this->nOB_            = this->singleSlater_->nOccB();
    this->nVA_            = this->singleSlater_->nVirA();
    this->nVB_            = this->singleSlater_->nVirB();
    this->aoERI_          = this->singleSlater_->aointegrals()->aoERI_.get();
    this->elecDipole_     = 
      this->singleSlater_->aointegrals()->elecDipole_.get();

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

  }
  void alloc();

  #include <sdresponse_getset.h>
  #include <sdresponse_qnrelated.h>
  #include <sdresponse_io.h>
  #include <sdresponse_misc.h>
  #include <sdresponse_prop.h>

};

#include <sdresponse_alloc.h>
} // namespace ChronusQ
#endif
