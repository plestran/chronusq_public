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


  void formAOTDen(const TVecMap &, TMatrix &, TMatrix &);
  void formMOTDen(TVecMap &, const TMatrix &, const TMatrix &);
  void placeVOOV(const TVecMap &,TMatrix &, TMatrix&);
  void placeVVOO(const TVecMap &,TMatrix &);
  void retrvVOOV(TVecMap &,const TMatrix &,const TMatrix&);
  void retrvVVOO(TVecMap &,const TMatrix &);
  void initRMu();
  void scaleDagPPRPA(bool,TVecMap &,TVecMap &,TVecMap *AX=NULL); 
  void initMeth();
  void checkValid();

  // Get the running dimension of the response matrix

  /**
   * Get the Single Dimension for the A Matrix for TDA
   */ 
  inline void nSingleDimCIS(){
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->nSingleDim_ = this->nOV_;
    else 
      this->nSingleDim_ = this->nOAVA_ + this->nOBVB_;
  };

  /**
   * Get the dimension for the First Order Polarization Propagator (FOPP)
   */
  inline void nSingleDimFOPP(){
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->nSingleDim_ = 2*this->nOV_;
    else 
      this->nSingleDim_ = 2*this->nOAVA_ + 2*this->nOBVB_;
  };

  /**
   * Get the dimension for the Particle-Particle Random Phase Approximation (PPRPA)
   *
   *   IPPRPA toggles the spin separation for the response matrix:
   *     IPPRPA = 0 ... All Alpha Block
   *     IPPRPA = 1 ... Mixed Alpha-Beta Block
   *     IPPRPA = 2 ... All Beta Block
   */ 
  inline void nSingleDimPPRPA(){
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->nSingleDim_ = this->nVV_SLT_ + this->nOO_SLT_;
    else {
      if(this->iPPRPA_ == 0)      this->nSingleDim_ = this->nVAVA_SLT_ + this->nOAOA_SLT_;
      else if(this->iPPRPA_ == 1) this->nSingleDim_ = this->nVAVB_     + this->nOAOB_;
      else if(this->iPPRPA_ == 2) this->nSingleDim_ = this->nVBVB_SLT_ + this->nOBOB_SLT_;
    }
  };

  /**
   * Get the dimension for the Particle-Particle Tamm-Danckoff Approximation (PPTDA)
   *
   *   (N+2) -> N
   *
   *   IPPRPA toggles the spin separation for the response matrix:
   *     IPPRPA = 0 ... All Alpha Block
   *     IPPRPA = 1 ... Mixed Alpha-Beta Block
   *     IPPRPA = 2 ... All Beta Block
   */
  inline void nSingleDimPPATDA(){
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->nSingleDim_ = this->nVV_SLT_;
    else {
      if(this->iPPRPA_ == 0)      this->nSingleDim_ = this->nVAVA_SLT_;
      else if(this->iPPRPA_ == 1) this->nSingleDim_ = this->nVAVB_    ; 
      else if(this->iPPRPA_ == 2) this->nSingleDim_ = this->nVBVB_SLT_;
    }
  };

  /**
   * Get the dimension for the Particle-Particle Tamm-Danckoff Approximation (PPTDA)
   *
   *   (N-2) -> N
   *
   *   IPPRPA toggles the spin separation for the response matrix:
   *     IPPRPA = 0 ... All Alpha Block
   *     IPPRPA = 1 ... Mixed Alpha-Beta Block
   *     IPPRPA = 2 ... All Beta Block
   */
  inline void nSingleDimPPCTDA(){
    if(this->Ref_ == SingleSlater<double>::TCS)
      this->nSingleDim_ = this->nVV_SLT_;
    else {
      if(this->iPPRPA_ == 0)      this->nSingleDim_ = this->nOAOA_SLT_;
      else if(this->iPPRPA_ == 1) this->nSingleDim_ = this->nOAOB_    ; 
      else if(this->iPPRPA_ == 2) this->nSingleDim_ = this->nOBOB_SLT_;
    }
  };


  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagSDResponse);
  void mpiRecv(int,int tag=tagSDResponse);

  // Proof of concept routines
  void incorePPRPA();
  void incorePPRPAnew();
  void incoreCIS();
  void incoreRPA();
  void formRM();
  TMatrix formRM2(TMatrix &XMO);
