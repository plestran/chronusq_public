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
/***********************
 * Form Density Matrix *
 ***********************/
template<typename T>
void SingleSlater<T>::formDensity(){
  if(this->nTCS_ == 1) {
    // Store Pa in Ps
    this->onePDMScalar_->noalias() = 
      this->moA_->block(0,0,this->nBasis_,this->nOA_)*
      this->moA_->block(0,0,this->nBasis_,this->nOA_).adjoint();
    if(!this->isClosedShell) {
      // Store Pb in Scratch
      this->NBSqScratch_->noalias() = 
        this->moB_->block(0,0,this->nBasis_,this->nOB_)*
        this->moB_->block(0,0,this->nBasis_,this->nOB_).adjoint();
      // Overwrite Pz with Pa - Pb
      (*this->onePDMMz_) =     (*this->onePDMScalar_) - (*this->NBSqScratch_);
      // Overwrite Ps with Pa + Pb
      (*this->onePDMScalar_) = (*this->onePDMScalar_) + (*this->NBSqScratch_);
    } else {
      // Factor of 2 for scalar
      (*this->onePDMScalar_) *= 2;
    }
  } else {
    this->NBTSqScratch_->noalias() = 
      this->moA_->block(0,0,this->nTCS_*this->nBasis_,this->nO_)*
      this->moA_->block(0,0,this->nTCS_*this->nBasis_,this->nO_).adjoint();
    std::vector<std::reference_wrapper<TMap>> scatter;
    for(auto iD = 0; iD < this->onePDM_.size(); iD++)
      scatter.emplace_back(*this->onePDM_[iD]);
    this->spinScatter((*this->NBTSqScratch_),scatter); 
  }
}

template<typename T>
void SingleSlater<T>::formFP(){
  // FP(S) = F(S)P(S)
  (*this->NBSqScratch_) = 
    (*this->fockOrthoScalar_) * (*this->onePDMOrthoScalar_);

  // FP(S) += F(z)P(z)
  if(this->nTCS_ == 2 or !this->isClosedShell)
    (*this->NBSqScratch_) += (*this->fockOrthoMz_) * (*this->onePDMOrthoMz_);

  // FP(S) += F(y)P(y) + F(x)P(x)
  if(this->nTCS_ == 2) {
    (*this->NBSqScratch_) += (*this->fockOrthoMy_) * (*this->onePDMOrthoMy_);
    (*this->NBSqScratch_) += (*this->fockOrthoMx_) * (*this->onePDMOrthoMx_);
  }

  this->FPScalar_->write(this->NBSqScratch_->data(),H5PredType<T>());
  
  if(this->nTCS_ == 2 or !this->isClosedShell) {
    // FP(z) = F(S)P(z) + F(z)P(S)
    (*this->NBSqScratch_) = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoMz_);
    (*this->NBSqScratch_) += 
      (*this->fockOrthoMz_) * (*this->onePDMOrthoScalar_);

    // FP(z) += i * (F(x)P(y) - F(y)P(x))
    if(this->nTCS_ == 2) {
      (*this->NBSqScratch_) += DIISComplexScale<T>() * 
        (*this->fockOrthoMx_) * (*this->onePDMOrthoMy_);
      (*this->NBSqScratch_) -= DIISComplexScale<T>() * 
        (*this->fockOrthoMy_) * (*this->onePDMOrthoMx_);
    }
    this->FPMz_->write(this->NBSqScratch_->data(),H5PredType<T>());
  }


  if(this->nTCS_ == 2) {
    // FP(y) = F(S)P(y) + F(y)P(S)
    (*this->NBSqScratch_) = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoMy_);
    (*this->NBSqScratch_) += 
      (*this->fockOrthoMy_) * (*this->onePDMOrthoScalar_);

    // FP(y) += i * (F(z)P(x) - F(x)P(z))
    (*this->NBSqScratch_) += DIISComplexScale<T>() * 
      (*this->fockOrthoMz_) * (*this->onePDMOrthoMx_);
    (*this->NBSqScratch_) -= DIISComplexScale<T>() * 
      (*this->fockOrthoMx_) * (*this->onePDMOrthoMz_);

    this->FPMy_->write(this->NBSqScratch_->data(),H5PredType<T>());

    // FP(x) = F(S)P(x) + F(x)P(S)
    (*this->NBSqScratch_) = 
      (*this->fockOrthoScalar_) * (*this->onePDMOrthoMx_);
    (*this->NBSqScratch_) += 
      (*this->fockOrthoMx_) * (*this->onePDMOrthoScalar_);

    // FP(x) += i * (F(y)P(z) - F(z)P(y))
    (*this->NBSqScratch_) += DIISComplexScale<T>() * 
      (*this->fockOrthoMy_) * (*this->onePDMOrthoMz_);
    (*this->NBSqScratch_) -= DIISComplexScale<T>() * 
      (*this->fockOrthoMz_) * (*this->onePDMOrthoMy_);

    this->FPMx_->write(this->NBSqScratch_->data(),H5PredType<T>());
  }

};

template<typename T>
void SingleSlater<T>::fockCUHF() {
  TMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
  TMap DelF(this->delFMem_,this->nBasis_,this->nBasis_);
  TMap Lambda(this->lambdaMem_,this->nBasis_,this->nBasis_);

  int activeSpace  = this->molecule_->multip() - 1;
  int coreSpace    = (this->molecule_->nTotalE() - activeSpace) / 2;
  int virtualSpace = this->nBasis_ - coreSpace - activeSpace;

  // DelF = X * (F(A) - F(B)) * X

  this->aointegrals_->Ortho1Trans(*this->fockMz_,DelF);
  (*this->fockMz_) *= 0.5;

  // DelF = C(NO)^\dagger * DelF * C(NO) (Natural Orbitals)
  (*this->NBSqScratch_) = P.transpose() * DelF;
  DelF = (*this->NBSqScratch_) * P;

  Lambda.setZero();
  for(auto i = activeSpace + coreSpace; i < this->nBasis_; i++)
  for(auto j = 0                      ; j < coreSpace    ; j++){
    Lambda(i,j) = -DelF(i,j);
    Lambda(j,i) = -DelF(j,i);
  }

  (*this->NBSqScratch_) = P * Lambda;
  Lambda = (*this->NBSqScratch_) * P.transpose();

  this->aointegrals_->Ortho2Trans(Lambda,Lambda);

  (*this->fockMz_) += 2*Lambda;
};

template <typename T>
void SingleSlater<T>::setRef(const std::string &methStr) {
//std::vector<std::string> methComp;
//boost::split(methComp,methStr,boost::is_any_of("- "));

  // Kohn-Sham lists
  std::vector<std::string> KS {
    "SLATER","B88","PBE","LSDA","SVWN5","BLYP","B3LYP","BHANDH"
  };

  std::vector<std::string> RKSL,UKSL,GKSL,X2CKSL;
  for(auto f : KS) {
    RKSL.emplace_back("R" + f);
    UKSL.emplace_back("U" + f);
    GKSL.emplace_back("G" + f);
    X2CKSL.emplace_back("X2C-" + f);
  }

  std::map<std::string,std::function<void()>> createMaps {
    {"LSDA"  ,[&]() -> void {this->createLSDA();}},
    {"SLATER",[&]() -> void {this->createSlater();}},
    {"SVWN5" ,[&]() -> void {this->createSVWN5();}},
    {"BLYP"  ,[&]() -> void {this->createBLYP();}},
    {"B88"   ,[&]() -> void {this->createB88();}},
    {"PBE"   ,[&]() -> void {this->createPBE();}},
    {"B3LYP" ,[&]() -> void {this->createB3LYP();}},
    {"BHANDH",[&]() -> void {this->createBHandH();}}
  };

  // Populate Ref_
  if(!methStr.compare("RHF"))        this->Ref_ = REFERENCE::RHF;
  else if(!methStr.compare("UHF"))   this->Ref_ = REFERENCE::UHF;
  else if(!methStr.compare("GHF"))   this->Ref_ = REFERENCE::GHF;
  else if(!methStr.compare("X2C"))   this->Ref_ = REFERENCE::X2C;

  else if(std::find(RKSL.begin(),RKSL.end(),methStr) != RKSL.end())   
    this->Ref_ = REFERENCE::RKS;
  else if(std::find(UKSL.begin(),UKSL.end(),methStr) != UKSL.end())   
    this->Ref_ = REFERENCE::UKS;
  else if(std::find(GKSL.begin(),GKSL.end(),methStr) != GKSL.end())   
    this->Ref_ = REFERENCE::GKS;
  else if(std::find(X2CKSL.begin(),X2CKSL.end(),methStr) != X2CKSL.end())   
    this->Ref_ = REFERENCE::X2C;

  else CErr("Fatal Error: Reference "+methStr+" is not recognized",
         this->fileio_->out);

  // Decide if restricted
  if(this->Ref_ == RHF or this->Ref_ == RKS)
     this->isClosedShell = true;
  else
     this->isClosedShell = false;

  // Decide if 2C
  if(this->Ref_ == GHF or this->Ref_ == GKS or this->Ref_ == X2C)
    this->nTCS_ = 2;

  // Decide if DFT
  this->isHF = this->Ref_ == RHF or this->Ref_ == UHF or this->Ref_ == GHF or
     !methStr.compare("X2C");
  this->isDFT = !this->isHF;

  // Construct DFT Functional
  if(this->isDFT) {
    for( auto f : KS ) 
      if(methStr.find(f) != std::string::npos) createMaps[f]();
  }


  // If 2C, turn on level shifting
//if(this->nTCS_ == 2) this->doLevelShift = true;
  
  // Generate the Method String
  this->genMethString();
}

template <typename T>
void SingleSlater<T>::setupRef() {
  if( this->Ref_ == X2C ) {
    this->aointegrals_->doX2C                = true;
    this->aointegrals_->useFiniteWidthNuclei = true;
  }
}

template <typename T>
void SingleSlater<T>::genMethString(){
  if(this->Ref_ == _INVALID) 
    CErr("Fatal: SingleSlater reference not set!",this->fileio_->out);


  std::map<DFT,std::string> funcMap {
    {LSDA  ,"LSDA"},
    {SLATER  ,"SLATER"},
    {SVWN5  ,"SVWN5"},
    {BLYP  ,"BLYP"},
    {B88   ,"B88"},
    {pbe   ,"PBE"},
    {B3LYP ,"B3LYP"},
    {BHandH,"BHandH"},
  };


  this->getAlgebraicField(); 
  this->SCFType_      = this->algebraicField_      + " ";
  this->SCFTypeShort_ = this->algebraicFieldShort_ + "-";
  
  std::string generalReference;
  std::string generalRefShort;
  if(this->isHF){
    generalReference = "Hartree-Fock";
    generalRefShort  = "HF";
  } else if(this->isDFT) {
    generalReference = "Kohn-Sham (" + funcMap[this->DFTKernel_] + ")";
    generalRefShort  = funcMap[this->DFTKernel_];
/*
    if(this->DFTKernel_ == USERDEFINED)
      generalRefShort  = "KS";
    else if(this->DFTKernel_ == LSDA) {
      generalReference += " (LSDA)";
      generalRefShort  = "LSDA";
*/
  }

  if(this->isClosedShell) {
    this->SCFType_      += "Restricted " + generalReference; 
    this->SCFTypeShort_ += "R" + generalRefShort;
  } else if(this->nTCS_ == 1) {
    this->SCFType_      += "Unrestricted " + generalReference; 
    this->SCFTypeShort_ += "U" + generalRefShort;
/*
  } else if(this->Ref_ == CUHF) {
    this->SCFType_      += "Constrained Unrestricted " + generalReference; 
    this->SCFTypeShort_ += "CU" + generalRefShort;
*/
  } else if(this->Ref_ == X2C) {
    this->SCFType_      += "Exact Two-Component " + generalReference; 
    this->SCFTypeShort_ += "X2C-"+generalRefShort;
  } else if(this->nTCS_ == 2) {
    this->SCFType_      += "Generalized " + generalReference; 
    this->SCFTypeShort_ += "G" + generalRefShort;
  }
}

template <typename T>
void SingleSlater<T>::McWeeny(std::vector<TMap*> OD, int MAX){
//prettyPrintSmart(cout,*OD[0],"Before");
  (*OD[0]) /= 2.0;
  for(auto I = 0; I < MAX; I++){
    (*NBSqScratch_) = (*OD[0]) * (*OD[0]);

    (*NBSqScratch2_) = (*OD[0]) - (*NBSqScratch_);
    if(NBSqScratch2_->norm() < 1e-10) break;

    (*NBSqScratch2_) = (*NBSqScratch_) *(*OD[0]);
    (*OD[0]) = 3* (*NBSqScratch_) - 2 * (*NBSqScratch2_);

  }
  (*OD[0]) *= 2.0;
//prettyPrintSmart(cout,*OD[0],"After");
  
} 
