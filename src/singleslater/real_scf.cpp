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
#include <singleslater.h>
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;

//----------------------------------------//
// do the SCF                             //
// Sajan                                  //
//----------------------------------------//
namespace ChronusQ {
template<>
void SingleSlater<double>::printDensityInfo(double PAlphaRMS,double EDelta){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<std::scientific<<EDelta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Density = "<<std::setw(15)<<std::scientific<<PAlphaRMS<<endl;
};
template<>
void SingleSlater<double>::printDensityInfo(double PAlphaRMS, double PBetaRMS, double EDelta){
  this->fileio_->out<<"\nSCF Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"    Delta-E = "<<std::setw(15)<<std::scientific<<EDelta<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Alpha Density = "<<std::setw(15)<<std::scientific<<PAlphaRMS<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"RMS Beta Density = "<<std::setw(15)<<std::scientific<<PBetaRMS<<endl;
};

template<>
void SingleSlater<double>::SCF(){
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();
  double EOld;
  double EDelta;
  double PAlphaRMS;
  double PBetaRMS;

  int maxIter    = 128; 
  int n = this->nBasis_; 
  double Dtol = 1e-10;
  double Etol = 1e-8;
  int iter; 
  
  //RealMatrix memory allocation
  int lenX = n*n; // X
  int lenXp= n*n; //Xp
  int lenF = n*n; // Fock
  int lenP = n*n; // Density
  int lenB = 49;  // B
  int lenCoeff = 7;   // Coeff
  int lwork = 4*n; // LAPACK Workspace
  int lenLambda = n*n; //CUHF lambda
  int lenPsort  = n*n; //CUHF Psort
  int lenDelF   = n*n; //CUHF DelF 
  int lenOccNum = n;   //CUHF occupation number
  int LenScr = lenX + lenF + lenP + lenB + lenCoeff + lwork + 2*lenF*(lenCoeff-1);
  if(this->Ref_ != RHF) LenScr += lenF + lenP + 2*lenF*(lenCoeff-1);
  if(doCUHF) LenScr += lenXp + lenOccNum + lenPsort + lenLambda + lenDelF + lenP;

  double *SCR            = NULL;
  double *XMem           = NULL;
  double *FpAlphaMem     = NULL;
  double *POldAlphaMem   = NULL;
  double *FpBetaMem      = NULL;
  double *POldBetaMem    = NULL;
//double *BMem           = NULL;
//double *coef           = NULL;
  double *ErrorAlphaMem  = NULL;  
  double *ErrorBetaMem   = NULL;
  double *FADIIS         = NULL;
  double *FBDIIS         = NULL;
  double *work           = NULL;
  double *XpMem    	 = NULL;
  double *PsortMem       = NULL;
  double *lambdaMem      = NULL;
  double *delFMem        = NULL;
  double *Pmem           = NULL;
  double *occNum         = NULL;
  int info;

  SCR = new double [LenScr];
  XMem    = SCR;
  if(this->Ref_ == RHF) {
    FpAlphaMem   = XMem + lenX;
    POldAlphaMem = FpAlphaMem + lenF;
    ErrorAlphaMem = POldAlphaMem + lenP;
    FADIIS = ErrorAlphaMem + lenF*(lenCoeff-1);
    work = FADIIS + lenF*(lenCoeff-1);
  } else {
    FpAlphaMem   = XMem + lenX;
    FpBetaMem    = FpAlphaMem + lenF;
    POldAlphaMem = FpBetaMem  + lenF;
    POldBetaMem  = POldAlphaMem + lenP;
    ErrorAlphaMem = POldBetaMem + lenP;
    ErrorBetaMem = ErrorAlphaMem + lenF*(lenCoeff-1);
    FADIIS = ErrorBetaMem + lenF*(lenCoeff-1);
    FBDIIS = FADIIS + lenF*(lenCoeff-1);
    work = FBDIIS + lenF*(lenCoeff-1);
  }
  if (doCUHF){
    XpMem = work + lwork;
    delFMem  = XpMem + lenXp;
    lambdaMem= delFMem  + lenDelF;
    PsortMem = lambdaMem+ lenLambda;
    Pmem     = PsortMem + lenPsort;
    occNum   = Pmem     + lenP;
  }
  RealMap X(XMem,n,n);
  RealMap Xp(XpMem,n,n);
  RealMap FpAlpha(FpAlphaMem,n,n);
  RealMap POldAlpha(POldAlphaMem,n,n);
  RealMap FpBeta(FpBetaMem,n,n);
  RealMap POldBeta(POldBetaMem,n,n);
  RealMap DelF(delFMem,n,n);
  RealMap lambda(lambdaMem,n,n);
  RealMap Psort(PsortMem,n,n);
  RealMap P(Pmem,n,n);

  char j='V';
  char u='U';
  int numE;
  int actSpace;
  int coreSpace;
  int virSpace;

  X=(*this->aointegrals_->overlap_).pow(-0.5);
  if (doCUHF){
    Xp= (*this->aointegrals_->overlap_).pow(0.5);
    numE     = this->molecule_->nTotalE();
    actSpace = this->molecule_->multip() - 1;
    coreSpace= (numE - actSpace)/2;
    virSpace = n - coreSpace - actSpace;
  }
  for (iter=0; iter<maxIter; iter++){
    this->fileio_->out << endl << endl << bannerTop <<endl;  
    this->fileio_->out << "SCF iteration:"<< iter+1 <<endl;  
    this->fileio_->out << bannerEnd <<endl;  
    if (doCUHF){
      if (this->Ref_ != RHF){
        P = (*this->densityA_ + *this->densityB_)/2 ;
      }
      else{
        P = *this->densityA_/2;
      }
      P = Xp * P * Xp;
      dsyev_(&j, &u, &n, P.data(), &n, occNum, work, &lwork, &info); //diagonalizing P
      P.transposeInPlace();
    
      for (auto i=0;i<n;i++){
        Psort.col(i)=P.col(n-i-1);
      } 
      P = Psort;
    
      if (this->Ref_ != RHF){
        DelF = (*this->fockA_ - *this->fockB_) / 2; 
      }
      else{
        DelF = (*this->fockA_)/2;
      }
      DelF = X * DelF * X;
      DelF = P.transpose() * DelF * P;
      lambda.setZero();
      for (auto i=coreSpace+actSpace ; i < n ; i++){
        for (auto j = 0; j<coreSpace; j++){
	  lambda(i,j)= -1.0 * DelF(i,j);
          lambda(j,i)= -1.0 * DelF(j,i);
        };
      };
      lambda = P * lambda * P.transpose();
      lambda = Xp *  lambda * Xp;
   
      *this->fockA_ = *this->fockA_ + lambda;
      if (this->Ref_ != RHF){
        *this->fockB_ = *this->fockB_ - lambda;
        POldBeta   = (*this->densityB_);
	FpBeta = X.transpose() * (*this->fockB_) * X;
        dsyev_(&j,&u,&n, FpBeta.data(), &n, occNum, work, &lwork, &info);
        FpBeta.transposeInPlace();
        (*this->moB_)= X * FpBeta;
      }
      POldAlpha = (*this->densityA_);
      EOld = this->totalEnergy;
      FpAlpha = X.transpose() * (*this->fockA_) * X;
      dsyev_(&j,&u,&n, FpAlpha.data(), &n, occNum, work, &lwork, &info);
      FpAlpha.transposeInPlace();
      (*this->moA_)= X * FpAlpha;
    }
    else{
      POldAlpha = (*this->densityA_);
      FpAlpha   = X.transpose()*(*this->fockA_)*X;
      if(this->Ref_ != RHF) {
        POldBeta   = (*this->densityB_);
        FpBeta     = X.transpose()*(*this->fockB_)*X;
      }
      EOld = this->totalEnergy;
      dsyev_(&j,&u,&n, FpAlpha.data(), &n, this->epsA_->data(), work, &lwork, &info);
      FpAlpha.transposeInPlace(); // bc Row Major
      (*this->moA_) = X*FpAlpha;
      if(this->Ref_ != RHF) {
        dsyev_(&j,&u,&n, FpBeta.data(), &n, this->epsB_->data(), work, &lwork, &info);
        FpBeta.transposeInPlace(); // bc Row Major
        (*this->moB_) = X*FpBeta;
      }
    }

    this->formDensity();
    this->formFock();
    this->computeEnergy();
    
    
    // DIIS 
    if (!doCUHF){
      RealMap ErrA(ErrorAlphaMem + (iter%(lenCoeff-1))*lenF,n,n);
      ErrA = (*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_);
      ErrA -= (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_);
      std::memcpy(FADIIS+(iter%(lenCoeff-1))*lenF,this->fockA_->data(),lenF*sizeof(double));
      if(this->Ref_ != RHF){
        RealMap ErrB(ErrorBetaMem + (iter%(lenCoeff-1))*lenF,n,n);
        ErrB = (*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_);
        ErrB -= (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_);
        std::memcpy(FBDIIS+(iter%(lenCoeff-1))*lenF,this->fockB_->data(),lenF*sizeof(double));
      }
    
      if(iter % (lenCoeff-1) == (lenCoeff-2) && iter != 0) 
        this->CDIIS(lenCoeff,ErrorAlphaMem,FADIIS,ErrorBetaMem,FBDIIS);
    }
    PAlphaRMS=((*this->densityA_)-POldAlpha).norm();
    if(this->Ref_ != RHF) PBetaRMS = ((*this->densityB_) - POldBeta).norm();
    EDelta= this->totalEnergy-EOld;
    if(this->Ref_ == RHF) this->printDensityInfo(PAlphaRMS,EDelta);     
    else this->printDensityInfo(PAlphaRMS,PBetaRMS,EDelta);     
    
    if (doCUHF){
      if(this->Ref_ == RHF){
        if(PAlphaRMS < Dtol && pow((this->totalEnergy-EOld),2)<Etol) break;
      } else {
        if((PAlphaRMS < Dtol || PBetaRMS < Dtol) && pow((this->totalEnergy-EOld),2)<Etol) break;
      } 
    }
    else{
      if(this->Ref_ == RHF){
        if(PAlphaRMS < Dtol && pow((this->totalEnergy-EOld),2)<Etol) break;
      } else {
        if(PAlphaRMS < Dtol && PBetaRMS < Dtol && pow((this->totalEnergy-EOld),2)<Etol) break;
      }
    }
  };
  
  //freeing the memory
  delete [] SCR;
  
  if(iter >= maxIter)
    this->fileio_->out << "SCF Failed to converge within maximum number of iterations" << endl << endl;
  this->fileio_->out <<"\n"<<endl; 
  this->fileio_->out << bannerEnd <<endl<<std::fixed;
  this->fileio_->out << "\nRequested convergence on RMS density matrix = " <<std::setw(5)<<Dtol <<"  within  " <<maxIter <<"  cycles."<<endl;
  this->fileio_->out << "Requested convergence on             energy = " <<Etol << endl;
  if(maxIter > iter){
    if(doCUHF) this->fileio_->out << "\nSCF Done: E(CUHF) = "<< this->totalEnergy << "  Eh  after  "<< iter+1 << "  cycles" <<endl;
    else{
      if(this->Ref_ == RHF) this->fileio_->out << "\nSCF Done: E(RHF) = "<< this->totalEnergy << "  Eh  after  "<< iter+1 << "  cycles" <<endl;
      else this->fileio_->out << "\nSCF Done: E(UHF) = "<< this->totalEnergy << "  Eh  after  "<< iter+1 << "  cycles" <<endl;
  }
  }
  this->fileio_->out << bannerEnd <<endl;
}; 

template<>
void SingleSlater<double>::formX(){
  RealMap X(this->XMem_,this->nBasis_,this->nBasis_);
  X = (*this->aointegrals_->overlap_).pow(-0.5); // Make this more efficient... FIXME

  if(this->Ref_ == CUHF){
    RealMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);
    Xp = (*this->aointegrals_->overlap_).pow(0.5); // Make this more efficient... FIXME
  }
}

template<>
void SingleSlater<double>::formNO(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'L';

  RealMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
  RealMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);

  P = 0.5 * Xp * (*this->densityA_) * Xp;
  if(!this->isClosedShell)
    P += 0.5 * Xp * (*this->densityB_) * Xp;

  dsyev_(&JOBZ,&UPLO,&this->nBasis_,this->PNOMem_,&this->nBasis_,this->occNumMem_,
         this->WORK_,&this->LWORK_,&INFO);
  if(INFO != 0) CErr("DSYEV Failed in FormNO",this->fileio_->out);
  P.transposeInPlace();

  // Swap Ordering
  for(auto i = 0; i < this->nBasis_/2; i++) P.col(i).swap(P.col(this->nBasis_ - i- 1));

}

template<>
void SingleSlater<double>::diagFock(){
  int INFO;
  char JOBZ = 'V';
  char UPLO = 'U';

  RealMap X(this->XMem_,this->nBasis_,this->nBasis_);
  RealMap POldAlpha(this->POldAlphaMem_,this->nBasis_,this->nBasis_);
  RealMap FpAlpha(this->FpAlphaMem_,this->nBasis_,this->nBasis_);
  RealMap POldBeta(this->POldBetaMem_,0,0);
  RealMap FpBeta(this->FpBetaMem_,0,0);
  if(!this->isClosedShell){
    new (&POldBeta)  RealMap(this->POldBetaMem_, this->nBasis_,this->nBasis_);
    new (&FpBeta) RealMap(this->FpBetaMem_,this->nBasis_,this->nBasis_);
  }


  if(this->Ref_ == CUHF){
    RealMap P(this->PNOMem_,this->nBasis_,this->nBasis_);
    RealMap Xp(this->XpMem_,this->nBasis_,this->nBasis_);
    RealMap DelF(this->delFMem_,this->nBasis_,this->nBasis_);
    RealMap Lambda(this->lambdaMem_,this->nBasis_,this->nBasis_);

    int activeSpace  = this->molecule_->multip() - 1;
    int coreSpace    = (this->molecule_->nTotalE() - activeSpace) / 2;
    int virtualSpace = this->nBasis_ - coreSpace - activeSpace;

    DelF = 0.5 * X * (*this->fockA_) * X;
    if(!this->isClosedShell)
      DelF -= 0.5 * X * (*this->fockB_) * X;
 
    DelF = P.transpose() * DelF * P;
 
    Lambda.setZero();
    for(auto i = activeSpace + coreSpace; i < this->nBasis_; i++)
    for(auto j = 0                      ; j < coreSpace    ; j++){
      Lambda(i,j) = -DelF(i,j);
      Lambda(j,i) = -DelF(j,i);
    }
    Lambda = P  * Lambda * P.transpose();
    Lambda = Xp * Lambda * Xp;  
 
    (*this->fockA_) += Lambda;
    if(!this->isClosedShell) (*this->fockB_) -= Lambda;
  }

  POldAlpha = (*this->densityA_);
  if(!this->isClosedShell) POldBeta = (*this->densityB_);
  FpAlpha = X.transpose() * (*this->fockA_) * X;
  dsyev_(&JOBZ,&UPLO,&this->nBasis_,this->FpAlphaMem_,&this->nBasis_,this->epsA_->data(),
         this->WORK_,&this->LWORK_,&INFO);
  if(INFO != 0) CErr("DSYEV Failed Fock Alpha",this->fileio_->out);
  FpAlpha.transposeInPlace(); // bc row major
  (*this->moA_) = X * FpAlpha;

  if(!this->isClosedShell){
    FpBeta = X.transpose() * (*this->fockB_) * X;
    dsyev_(&JOBZ,&UPLO,&this->nBasis_,this->FpBetaMem_,&this->nBasis_,this->epsB_->data(),
           this->WORK_,&this->LWORK_,&INFO);
    if(INFO != 0) CErr("DSYEV Failed Fock Beta",this->fileio_->out);
    FpBeta.transposeInPlace(); // bc row major
    (*this->moB_) = X * FpBeta;
  }

}

template<>
void SingleSlater<double>::evalConver(){
  double EOld;
  double EDelta;
  double PAlphaRMS;
  double PBetaRMS;
  double Dtol = 1e-10;
  double Etol = 1e-8;

  RealMap POldAlpha(this->POldAlphaMem_,this->nBasis_,this->nBasis_);
  RealMap POldBeta(this->POldBetaMem_,0,0);
  if(!this->isClosedShell){
    new (&POldBeta)  RealMap(this->POldBetaMem_, this->nBasis_,this->nBasis_);
  }

  EOld = this->totalEnergy;
  this->computeEnergy();
  EDelta = this->totalEnergy - EOld;

  PAlphaRMS = ((*this->densityA_) - POldAlpha).norm();
  if(!this->isClosedShell) PBetaRMS = ((*this->densityB_) - POldBeta).norm();

  if(this->isClosedShell) this->printDensityInfo(PAlphaRMS,EDelta);
  else                    this->printDensityInfo(PAlphaRMS,PBetaRMS,EDelta);

  this->isConverged = (PAlphaRMS < Dtol) && (std::pow(EDelta,2) < Etol);
  if(!this->isClosedShell)
    this->isConverged = this->isConverged && (PBetaRMS < Dtol);
}

template<>
void SingleSlater<double>::SCF2(){
  if(!this->aointegrals_->haveAOOneE) this->aointegrals_->computeAOOneE();

  int maxIter    = 128; 
  double Dtol = 1e-10;
  double Etol = 1e-8;
  int n = this->nBasis_; 
  int iter; 

  this->initSCFMem();
  this->formX();
  for (iter=0; iter<maxIter; iter++){
    this->fileio_->out << endl << endl << bannerTop <<endl;  
    this->fileio_->out << "SCF iteration:"<< iter+1 <<endl;  
    this->fileio_->out << bannerEnd <<endl;  

    if(this->Ref_ == CUHF) this->formNO();
    this->diagFock();
    this->formDensity();
    this->formFock();
  //this->computeEnergy();

    if (this->Ref_ != CUHF){ // DIIS NYI for CUHF
      RealMap ErrA(this->ErrorAlphaMem_ + (iter % (this->lenCoeff_-1)) * this->lenF_,
                   this->nBasis_,this->nBasis_);

      ErrA = (*this->fockA_) * (*this->densityA_) * (*this->aointegrals_->overlap_);
      ErrA -= (*this->aointegrals_->overlap_) * (*this->densityA_) * (*this->fockA_);

      std::memcpy(this->FADIIS_+(iter % (this->lenCoeff_-1)) * this->lenF_,
                  this->fockA_->data(),this->lenF_ * sizeof(double));

      if(!this->isClosedShell){
        RealMap ErrB(this->ErrorBetaMem_ + (iter % (this->lenCoeff_-1)) * this->lenF_,
                     this->nBasis_,this->nBasis_);

        ErrB = (*this->fockB_) * (*this->densityB_) * (*this->aointegrals_->overlap_);
        ErrB -= (*this->aointegrals_->overlap_) * (*this->densityB_) * (*this->fockB_);

        std::memcpy(this->FBDIIS_ + (iter % (this->lenCoeff_-1)) * this->lenF_,
                    this->fockB_->data(),this->lenF_ * sizeof(double));
      }
    
      if(iter % (this->lenCoeff_-1) == (this->lenCoeff_-2) && iter != 0) 
        this->CDIIS(this->lenCoeff_,this->ErrorAlphaMem_,this->FADIIS_,this->ErrorBetaMem_,
                    this->FBDIIS_);
    }
    this->evalConver();
    if(this->isConverged) break;

  }; // SCF Loop
  delete [] this->SCF_SCR;

  if(!this->isConverged)
    CErr("SCF Failed to converge within maximum number of iterations",this->fileio_->out);
  this->fileio_->out <<"\n"<<endl; 
  this->fileio_->out << bannerEnd <<endl<<std::fixed;
  this->fileio_->out << "\nRequested convergence on RMS density matrix = " <<std::setw(5)<<Dtol <<"  within  " <<maxIter <<"  cycles."<<endl;
  this->fileio_->out << "Requested convergence on             energy = " <<Etol << endl;
  if(this->isConverged){
    this->fileio_->out << endl << "SCF Completed: E(";
    if(this->Ref_ == RHF)  this->fileio_->out << "RHF";
    if(this->Ref_ == UHF)  this->fileio_->out << "UHF";
    if(this->Ref_ == CUHF) this->fileio_->out << "CUHF";
    this->fileio_->out << ") = ";
    this->fileio_->out << this->totalEnergy << "  Eh after  " << iter + 1 << "  SCF Iterations" << endl;
  }
  this->fileio_->out << bannerEnd <<endl;
}
} // namespace ChronusQ
