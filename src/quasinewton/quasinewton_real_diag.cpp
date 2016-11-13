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
#include <quasinewton.h>

namespace ChronusQ {
  template<>
  void QuasiNewton<double>::stdHerDiag(int NTrial, ostream &output){
    // Solve E(R)| X(R) > = | X(R) > ω
    char JOBV = 'V';
    char UPLO = 'L';
    int INFO;
    RealMap A(this->XTSigmaRMem,NTrial,NTrial);
  //cout << "HERE" << endl;
  //cout << endl << A << endl;
    dsyev_(&JOBV,&UPLO,&NTrial,this->XTSigmaRMem,&NTrial,
           this->ERMem,this->WORK,&this->LWORK,&INFO); 
    if(INFO!=0) CErr("DSYEV failed to converge in Davison Iterations",output);
  } // stdHerDiag

  template<>
  void QuasiNewton<double>::symmHerDiag(int NTrial, ostream &output){
    /*
     * Solve S(R)| X(R) > = E(R)| X(R) > (1/ω)
     *
     * | X(R) > = | X(R)_g >
     *            | X(R)_u >
     *
     * The opposite (1/ω vs ω) is solved because the metric is not positive definite
     * and can therefore not be solved using DSYGV because of the involved Cholesky
     * decomposition.
     *
     */ 
    char JOBV = 'V';
    char UPLO = 'L';
    int iType = 1;
    int TwoNTrial = 2*NTrial;
    int INFO;

    RealMap SSuper(this->SSuperMem, 2*NTrial,2*NTrial);
    RealMap    SCPY(this->SCPYMem,   TwoNTrial,TwoNTrial);

    SCPY = SSuper; // Copy of original matrix to use for re-orthogonalization

    // Perform diagonalization of reduced subspace using DSYGV
    dsygv_(&iType,&JOBV,&UPLO,&TwoNTrial,this->SSuperMem,&TwoNTrial,
           this->ASuperMem,&TwoNTrial,this->ERMem,this->WORK,&this->LWORK,
           &INFO);
    if(INFO!=0) CErr("DSYGV failed to converge in Davison Iterations",output);
    
    // Grab the "positive paired" roots (throw away other element of the pair)
    this->ERMem += NTrial;
    RealVecMap ER    (this->ERMem,NTrial);
    new (&SSuper) RealMap(this->SSuperMem+2*NTrial*NTrial,2*NTrial,NTrial);

    // Swap the ordering because we solve for (1/ω)
    for(auto i = 0 ; i < NTrial; i++) ER(i) = 1.0/ER(i);
    for(auto i = 0 ; i < NTrial/2; i++){
      SSuper.col(i).swap(SSuper.col(NTrial - i - 1));
      double tmp = ER(i);
      ER(i) = ER(NTrial - i - 1);
      ER(NTrial - i - 1) = tmp;
    }

    /*
     * Re-orthogonalize the eigenvectors with respect to the metric S(R)
     * because DSYGV orthogonalzies the vectors with respect to E(R)
     * because we solve the opposite problem.
     *
     * Gramm-Schmidt
     */
    this->metBiOrth(SSuper,SCPY);

    // Separate the eigenvectors into gerade and ungerade parts
    RealMap XTSigmaR(this->XTSigmaRMem,NTrial,NTrial);
    RealMap XTSigmaL(this->XTSigmaLMem,NTrial,NTrial);
    XTSigmaR = SSuper.block(0,     0,NTrial,NTrial);
    XTSigmaL = SSuper.block(NTrial,0,NTrial,NTrial);
  } // symmHerDiag

  template<>
  void QuasiNewton<double>::symmNonHerDiag(int NTrial, ostream &output){
    char JOBVL = 'N';
    char JOBVR = 'V';
    int TwoNTrial = 2*NTrial;
    int *IPIV = new int[TwoNTrial];
    int INFO;

    RealMap  SSuper(this->SSuperMem, TwoNTrial,TwoNTrial);
    RealMap  ASuper(this->ASuperMem, TwoNTrial,TwoNTrial);
    RealMap    SCPY(this->SCPYMem,   TwoNTrial,TwoNTrial);
    RealMap NHrProd(this->NHrProdMem,TwoNTrial,TwoNTrial);

    SCPY = SSuper; // Copy of original matrix to use for re-orthogonalization

    // Invert the metric (maybe not needed?)
    dgetrf_(&TwoNTrial,&TwoNTrial,this->SSuperMem,&TwoNTrial,IPIV,&INFO);
    dgetri_(&TwoNTrial,this->SSuperMem,&TwoNTrial,IPIV,this->WORK,&this->LWORK,&INFO);
    delete [] IPIV;

    NHrProd = SSuper * ASuper;
  //cout << endl << "PROD" << endl << NHrProd << endl;

    dgeev_(&JOBVL,&JOBVR,&TwoNTrial,NHrProd.data(),&TwoNTrial,this->ERMem,this->EIMem,
           this->SSuperMem,&TwoNTrial,this->SSuperMem,&TwoNTrial,this->WORK,&this->LWORK,
           &INFO);
    // Sort eigensystem using Bubble Sort
    RealVecMap ER(this->ERMem,TwoNTrial);
    RealVecMap EI(this->EIMem,TwoNTrial);
    RealMap  VR(this->SSuperMem,TwoNTrial,TwoNTrial);
//  cout << endl << ER << endl;
    this->eigSrt(VR,ER);
//  cout << endl << ER << endl;
  
    // Grab the "positive paired" roots (throw away other element of the pair)
    this->ERMem += NTrial;
    new (&ER    ) RealVecMap(this->ERMem,NTrial);
    new (&SSuper) RealMap(this->SSuperMem+2*NTrial*NTrial,2*NTrial,NTrial);

    /*
     * Re-orthogonalize the eigenvectors with respect to the metric S(R)
     * because DSYGV orthogonalzies the vectors with respect to E(R)
     * because we solve the opposite problem.
     *
     * Gramm-Schmidt
     */
    this->metBiOrth(SSuper,SCPY);

    // Separate the eigenvectors into gerade and ungerade parts
    RealMap XTSigmaR(this->XTSigmaRMem,NTrial,NTrial);
    RealMap XTSigmaL(this->XTSigmaLMem,NTrial,NTrial);
    XTSigmaR = SSuper.block(0,     0,NTrial,NTrial);
    XTSigmaL = SSuper.block(NTrial,0,NTrial,NTrial);
  //cout << endl << "ER" << endl << ER << endl << endl;
  //cout << endl << "CR" << endl << XTSigmaR << endl << endl;
  //cout << endl << "CR" << endl << XTSigmaL << endl << endl;
//  CErr();
  }

  template<>
  void QuasiNewton<double>::diagMem(int NTrial){
    int IOff = NTrial;
    if(!this->isHermitian_ || this->symmetrizedTrial_){
      IOff += NTrial; // Space for paired eigenvalues or imaginary part
    }
    if(!this->isHermitian_ && this->symmetrizedTrial_){
      IOff += 2*NTrial; // Space for paired eigenvalues or imaginary part
    }
    this->ERMem = LAPACK_SCR;
    if(!this->isHermitian_){
      if(this->symmetrizedTrial_) this->EIMem = this->ERMem + 2*NTrial;
      else                        this->EIMem = this->ERMem +   NTrial;
    }
    this->WORK  = this->ERMem + IOff;
  } // diagMem

  /** Diagoanlize Reduced Problem **/
  template<>
  void QuasiNewton<double>::redDiag(int NTrial,ostream &output){
    this->diagMem(NTrial); 
  //cout << "HERE" << endl;
    if(this->isHermitian_) {
      if(!this->symmetrizedTrial_) this->stdHerDiag(NTrial,output);
      else                         this->symmHerDiag(NTrial,output);
    } else {
      if(this->symmetrizedTrial_) this->symmNonHerDiag(NTrial,output);
    } 
  } // redDiag

}; // namespace ChronusQ
