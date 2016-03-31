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
template <typename T>
Eigen::VectorXd NumericalDifferentiation<T>::ES2GSNACME_CIS(
  SingleSlater<T> &ss_p1, SingleSlater<T> &ss_m1, Response<T> &resp_p1,
  Response<T> &resp_m1, TMatrix &SAO_0_p1, TMatrix &SAO_0_m1,
  TMatrix &SMO_0_p1, TMatrix &SMO_0_m1){

  Eigen::VectorXd NACME(this->responseNRoots_);

  NACME.setZero();


  int nThreads = omp_get_max_threads();
  int nOcc     = this->singleSlater_undisplaced_->nOccA();
  int nVir     = this->singleSlater_undisplaced_->nVirA();

  std::vector<RealMatrix> Prod_0_p1(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA())));
  std::vector<RealMatrix> Prod_0_m1(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA())));

  RealMatrix T_0 = this->response_undisplaced_->transDen()[0].block(0,0,
    this->response_undisplaced_->nMatDim()[0],this->responseNRoots_);

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < this->responseNRoots_; iSt++){
    cout << "Working on iSt = " << iSt << " ES->GS" << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto i = 0, ia = 0; i < nOcc; i++)
      for(auto a = 0; a < nVir; a++, ia++){
        if(ia % nThreads != thread_id) continue;
        if(std::abs(T_0(ia,iSt)) < 1e-8) continue;
  
        Prod_0_p1[thread_id] = SMO_0_p1;
        Prod_0_p1[thread_id].row(i).swap(Prod_0_p1[thread_id].row(nOcc+a)); 
        Prod_0_m1[thread_id] = SMO_0_m1;
        Prod_0_m1[thread_id].row(i).swap(Prod_0_m1[thread_id].row(nOcc+a)); 
        
        double OvLp_IASwap_0_p1 =
          Prod_0_p1[thread_id].block(0,0,nOcc,nOcc).determinant();
        double OvLp_IASwap_0_m1 =
          Prod_0_m1[thread_id].block(0,0,nOcc,nOcc).determinant();
     
        #pragma omp critical
        {
          NACME(iSt) += T_0(ia,iSt) *
            (OvLp_IASwap_0_p1 - OvLp_IASwap_0_m1) /
            (2 * this->step);
        }
      } // loop ove ia
    }
  } // loop over states

  
//for(auto iTh = 0; iTh < nThreads; iTh++) NACME += TMPNACME[iTh];
  NACME = 2 * std::sqrt(0.5) * NACME;
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "GS -> ES NACME took " << elapsed.count() << " s" << endl;
//cout << NACME;
//cout << endl;
//cout << endl;
  return NACME;

};


// This is still the slow version
template <typename T>
RealMatrix NumericalDifferentiation<T>::ES2ESNACME_CIS(
  SingleSlater<T> &ss_p1, SingleSlater<T> &ss_m1, Response<T> &resp_p1,
  Response<T> &resp_m1, TMatrix &SAO_0_p1, TMatrix &SAO_0_m1,
  TMatrix &SMO_0_p1, TMatrix &SMO_0_m1){


  RealMatrix NACME(this->responseNRoots_,this->responseNRoots_);

  NACME.setZero();

  int nThreads = omp_get_max_threads();
  int nOcc     = this->singleSlater_undisplaced_->nOccA();
  int nVir     = this->singleSlater_undisplaced_->nVirA();

  this->checkPhase(resp_p1,resp_m1);
  this->checkPhase((*this->response_undisplaced_),resp_p1);
  this->checkPhase((*this->response_undisplaced_),resp_m1);

  std::vector<RealMatrix> SWAPPED_IA_0(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));
  std::vector<RealMatrix> SWAPPED_JB_p1(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));
  std::vector<RealMatrix> SWAPPED_JB_m1(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));

  std::vector<RealMatrix> Prod_0_p1(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));
  std::vector<RealMatrix> Prod_0_m1(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));

  RealMatrix T_0 = this->response_undisplaced_->transDen()[0].block(0,0,
    this->response_undisplaced_->nMatDim()[0],this->responseNRoots_);
  RealMatrix T_p1 = resp_p1.transDen()[0].block(0,0,
    resp_p1.nMatDim()[0],this->responseNRoots_);
  RealMatrix T_m1 = resp_m1.transDen()[0].block(0,0,
    resp_m1.nMatDim()[0],this->responseNRoots_);


  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < this->responseNRoots_; iSt++)
  for(auto jSt = 0; jSt < this->responseNRoots_; jSt++){
    cout << "Working on iSt = " << iSt <<" jSt = " << jSt << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto i = 0, iajb = 0, ia = 0; i < nOcc; i++)
      for(auto a = 0; a < nVir; a++, ia++){

        SWAPPED_IA_0[thread_id].setZero();
        SWAPPED_IA_0[thread_id] = (*this->singleSlater_undisplaced_->moA());
     
        SWAPPED_IA_0[thread_id].col(i).swap(SWAPPED_IA_0[thread_id].col(nOcc+a)); 
     
        for(auto j = 0, jb = 0; j < nOcc; j++)
        for(auto b = 0; b < nVir; b++, jb++, iajb++){
          if(iajb % nThreads != thread_id) continue;
          if(
            (std::abs(T_p1(jb,jSt)) < 1e-8) &&
            (std::abs(T_m1(jb,jSt)) < 1e-8) 
          ) continue;
       
          SWAPPED_JB_p1[thread_id].setZero();
          SWAPPED_JB_m1[thread_id].setZero();
       
          Prod_0_p1[thread_id].setZero();
          Prod_0_m1[thread_id].setZero();
          
          SWAPPED_JB_p1[thread_id] = (*ss_p1.moA());
          SWAPPED_JB_m1[thread_id] = (*ss_m1.moA());
       
          SWAPPED_JB_p1[thread_id].col(j).swap(SWAPPED_JB_p1[thread_id].col(nOcc+b)); 
          SWAPPED_JB_m1[thread_id].col(j).swap(SWAPPED_JB_m1[thread_id].col(nOcc+b)); 
          
          Prod_0_p1[thread_id] = 
            SWAPPED_IA_0[thread_id].transpose() * SAO_0_p1 * SWAPPED_JB_p1[thread_id];
          Prod_0_m1[thread_id] = 
            SWAPPED_IA_0[thread_id].transpose() * SAO_0_m1 * SWAPPED_JB_m1[thread_id];
       
          
          double OvLp_IASwap_0_p1 =
            Prod_0_p1[thread_id].block(0,0,nOcc,nOcc).determinant();
          double OvLp_IASwap_0_m1 =
            Prod_0_m1[thread_id].block(0,0,nOcc,nOcc).determinant();
     
       
          #pragma omp critical
          {
            NACME(iSt,jSt) += T_0(ia,iSt)*(
              T_p1(jb,jSt)*OvLp_IASwap_0_p1 - 
              T_m1(jb,jSt)*OvLp_IASwap_0_m1
            ) / (2 * this->step);
          }
        } // loop ove jb
      } // loop ove ia
    }
  } // loop over states
//NACME *= 2 * std::sqrt(0.5);
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "ES -> ES NACME took " << elapsed.count() << " s" << endl;

  cout << NACME << endl << endl;
  return NACME;

};

template <typename T>
Eigen::VectorXd NumericalDifferentiation<T>::ES2GSNACME_CIS(
  SingleSlater<T> &ss_p1, SingleSlater<T> &ss_m1, TMatrix &T_0, TMatrix &T_p1,
  TMatrix &T_m1, TMatrix &SAO_0_p1, TMatrix &SAO_0_m1,
  TMatrix &SMO_0_p1, TMatrix &SMO_0_m1){

  Eigen::VectorXd NACME(this->responseNRoots_);

  NACME.setZero();


  int nThreads = omp_get_max_threads();
  int nOcc     = this->singleSlater_undisplaced_->nOccA();
  int nVir     = this->singleSlater_undisplaced_->nVirA();

  std::vector<RealMatrix> Prod_0_p1(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA())));
  std::vector<RealMatrix> Prod_0_m1(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA())));

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < this->responseNRoots_; iSt++){
    cout << "Working on iSt = " << iSt << " ES->GS" << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto i = 0, ia = 0; i < nOcc; i++)
      for(auto a = 0; a < nVir; a++, ia++){
        if(ia % nThreads != thread_id) continue;
        if(std::abs(T_0(ia,iSt)) < 1e-8) continue;
  
        Prod_0_p1[thread_id] = SMO_0_p1;
        Prod_0_p1[thread_id].row(i).swap(Prod_0_p1[thread_id].row(nOcc+a)); 
        Prod_0_m1[thread_id] = SMO_0_m1;
        Prod_0_m1[thread_id].row(i).swap(Prod_0_m1[thread_id].row(nOcc+a)); 
        
        double OvLp_IASwap_0_p1 =
          Prod_0_p1[thread_id].block(0,0,nOcc,nOcc).determinant();
        double OvLp_IASwap_0_m1 =
          Prod_0_m1[thread_id].block(0,0,nOcc,nOcc).determinant();
     
        #pragma omp critical
        {
          NACME(iSt) += T_0(ia,iSt) *
            (OvLp_IASwap_0_p1 - OvLp_IASwap_0_m1) /
            (2 * this->step);
        }
      } // loop ove ia
    }
  } // loop over states

  
  NACME = 2 * std::sqrt(0.5) * NACME;
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "GS -> ES NACME took " << elapsed.count() << " s" << endl;
  prettyPrint(cout,NACME,"GS->ES");

  return NACME;

};

template <typename T>
RealMatrix NumericalDifferentiation<T>::ES2ESNACME_CIS(
  SingleSlater<T> &ss_p1, SingleSlater<T> &ss_m1, TMatrix &T_0, TMatrix &T_p1,
  TMatrix &T_m1, TMatrix &SAO_0_p1, TMatrix &SAO_0_m1,
  TMatrix &SMO_0_p1, TMatrix &SMO_0_m1){


  RealMatrix NACME(this->responseNRoots_,this->responseNRoots_);

  NACME.setZero();

  int nThreads = omp_get_max_threads();
  int nOcc     = this->singleSlater_undisplaced_->nOccA();
  int nVir     = this->singleSlater_undisplaced_->nVirA();


  std::vector<TMatrix> Prod_0_p1(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA()))
  );
  std::vector<TMatrix> Prod_0_m1(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA()))
  );

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < this->responseNRoots_; iSt++)
  for(auto jSt = 0; jSt < this->responseNRoots_; jSt++){
    cout << "Working on iSt = " << iSt <<" jSt = " << jSt << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto i = 0, iajb = 0, ia = 0; i < nOcc; i++)
      for(auto a = 0; a < nVir; a++, ia++)
      for(auto j = 0, jb = 0; j < nOcc; j++)
      for(auto b = 0; b < nVir; b++, jb++, iajb++){
        if(iajb % nThreads != thread_id) continue;
        if(
          (std::abs(T_p1(jb,jSt)) < 1e-8) &&
          (std::abs(T_m1(jb,jSt)) < 1e-8) 
        ) continue;
      
        Prod_0_p1[thread_id] = SMO_0_p1;
        Prod_0_p1[thread_id].row(i).swap(Prod_0_p1[thread_id].row(nOcc+a)); 
        Prod_0_p1[thread_id].col(j).swap(Prod_0_p1[thread_id].col(nOcc+b)); 
        Prod_0_m1[thread_id] = SMO_0_m1;
        Prod_0_m1[thread_id].row(i).swap(Prod_0_m1[thread_id].row(nOcc+a)); 
        Prod_0_m1[thread_id].col(j).swap(Prod_0_m1[thread_id].col(nOcc+b)); 
      
        
        double OvLp_IASwap_0_p1 =
          Prod_0_p1[thread_id].block(0,0,nOcc,nOcc).determinant();
        double OvLp_IASwap_0_m1 =
          Prod_0_m1[thread_id].block(0,0,nOcc,nOcc).determinant();
     
      
        #pragma omp critical
        {
          NACME(iSt,jSt) += T_0(ia,iSt)*(
            T_p1(jb,jSt)*OvLp_IASwap_0_p1 - 
            T_m1(jb,jSt)*OvLp_IASwap_0_m1
          ) / (2 * this->step);
        } // OMP
      } // loop ove jb
    } // OMP
  } // loop over states

  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "ES -> ES NACME took " << elapsed.count() << " s" << endl;

  prettyPrint(cout,NACME,"ES->ES");
  return NACME;

};

template <typename T>
RealMatrix NumericalDifferentiation<T>::ES2ESNACME_PPTDA(
  SingleSlater<T> &ss_p1, SingleSlater<T> &ss_m1, TMatrix &T_0, TMatrix &T_p1,
  TMatrix &T_m1, TMatrix &SAO_0_p1, TMatrix &SAO_0_m1,
  TMatrix &SMO_0_p1, TMatrix &SMO_0_m1){

  RealMatrix NACME(this->responseNRoots_,this->responseNRoots_);

  NACME.setZero();

  int nThreads = omp_get_max_threads();
  int nOcc     = this->singleSlater_undisplaced_->nOccA();
  int nVir     = this->singleSlater_undisplaced_->nVirA();

  std::vector<TMatrix> Prod_0_p1(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA()))
  );
  std::vector<TMatrix> Prod_0_m1(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA()))
  );

  
  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < this->responseNRoots_; iSt++)
  for(auto jSt = 0; jSt < this->responseNRoots_; jSt++){
    cout << "Working on iSt = " << iSt <<" jSt = " << jSt << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto a = 0, ab = 0, abcd = 0; a < nVir; a++      )
      for(auto b = 0        ; b <= a ;  b++, ab++)
      for(auto c = 0, cd = 0; c < nVir; c++      )
      for(auto d = 0        ; d <= c ;  d++, cd++, abcd++){
        if(abcd % nThreads != thread_id) continue; 
        if(
          (std::abs(T_p1(cd,jSt)) < 1e-8) &&
          (std::abs(T_m1(cd,jSt)) < 1e-8)
        ) continue;
     
	// Copy matricies
        Prod_0_p1[thread_id] = SMO_0_p1;
        Prod_0_m1[thread_id] = SMO_0_m1;
     
        // AC OVerlap

        Prod_0_p1[thread_id].row(nOcc).swap(Prod_0_p1[thread_id].row(nOcc+a)); 
        Prod_0_p1[thread_id].col(nOcc).swap(Prod_0_p1[thread_id].col(nOcc+c)); 
        Prod_0_m1[thread_id].row(nOcc).swap(Prod_0_m1[thread_id].row(nOcc+a)); 
        Prod_0_m1[thread_id].col(nOcc).swap(Prod_0_m1[thread_id].col(nOcc+c)); 
        
        double OvLp_AC_0_p1 =
          Prod_0_p1[thread_id].block(0,0,nOcc+1,nOcc+1).determinant();
        double OvLp_AC_0_m1 =
          Prod_0_m1[thread_id].block(0,0,nOcc+1,nOcc+1).determinant();

	// Swap Back
        Prod_0_p1[thread_id].row(nOcc).swap(Prod_0_p1[thread_id].row(nOcc+a)); 
        Prod_0_p1[thread_id].col(nOcc).swap(Prod_0_p1[thread_id].col(nOcc+c)); 
        Prod_0_m1[thread_id].row(nOcc).swap(Prod_0_m1[thread_id].row(nOcc+a)); 
        Prod_0_m1[thread_id].col(nOcc).swap(Prod_0_m1[thread_id].col(nOcc+c)); 
     
     
        // BD Overlap

        Prod_0_p1[thread_id].row(nOcc).swap(Prod_0_p1[thread_id].row(nOcc+b)); 
        Prod_0_p1[thread_id].col(nOcc).swap(Prod_0_p1[thread_id].col(nOcc+d)); 
        Prod_0_m1[thread_id].row(nOcc).swap(Prod_0_m1[thread_id].row(nOcc+b)); 
        Prod_0_m1[thread_id].col(nOcc).swap(Prod_0_m1[thread_id].col(nOcc+d)); 
        
        double OvLp_BD_0_p1 =
          Prod_0_p1[thread_id].block(0,0,nOcc+1,nOcc+1).determinant();
        double OvLp_BD_0_m1 =
          Prod_0_m1[thread_id].block(0,0,nOcc+1,nOcc+1).determinant();

	// Swap Back
        Prod_0_p1[thread_id].row(nOcc).swap(Prod_0_p1[thread_id].row(nOcc+b)); 
        Prod_0_p1[thread_id].col(nOcc).swap(Prod_0_p1[thread_id].col(nOcc+d)); 
        Prod_0_m1[thread_id].row(nOcc).swap(Prod_0_m1[thread_id].row(nOcc+b)); 
        Prod_0_m1[thread_id].col(nOcc).swap(Prod_0_m1[thread_id].col(nOcc+d)); 
     
        // AD Overlap

        Prod_0_p1[thread_id].row(nOcc).swap(Prod_0_p1[thread_id].row(nOcc+a)); 
        Prod_0_p1[thread_id].col(nOcc).swap(Prod_0_p1[thread_id].col(nOcc+d)); 
        Prod_0_m1[thread_id].row(nOcc).swap(Prod_0_m1[thread_id].row(nOcc+a)); 
        Prod_0_m1[thread_id].col(nOcc).swap(Prod_0_m1[thread_id].col(nOcc+d)); 
        
        double OvLp_AD_0_p1 =
          Prod_0_p1[thread_id].block(0,0,nOcc+1,nOcc+1).determinant();
        double OvLp_AD_0_m1 =
          Prod_0_m1[thread_id].block(0,0,nOcc+1,nOcc+1).determinant();

	// Swap Back
        Prod_0_p1[thread_id].row(nOcc).swap(Prod_0_p1[thread_id].row(nOcc+a)); 
        Prod_0_p1[thread_id].col(nOcc).swap(Prod_0_p1[thread_id].col(nOcc+d)); 
        Prod_0_m1[thread_id].row(nOcc).swap(Prod_0_m1[thread_id].row(nOcc+a)); 
        Prod_0_m1[thread_id].col(nOcc).swap(Prod_0_m1[thread_id].col(nOcc+d)); 
     
     
        // BC Overlap

        Prod_0_p1[thread_id].row(nOcc).swap(Prod_0_p1[thread_id].row(nOcc+b)); 
        Prod_0_p1[thread_id].col(nOcc).swap(Prod_0_p1[thread_id].col(nOcc+c)); 
        Prod_0_m1[thread_id].row(nOcc).swap(Prod_0_m1[thread_id].row(nOcc+b)); 
        Prod_0_m1[thread_id].col(nOcc).swap(Prod_0_m1[thread_id].col(nOcc+c)); 
       
        
        double OvLp_BC_0_p1 =
          Prod_0_p1[thread_id].block(0,0,nOcc+1,nOcc+1).determinant();
        double OvLp_BC_0_m1 =
          Prod_0_m1[thread_id].block(0,0,nOcc+1,nOcc+1).determinant();
     
	// Swap Back ** Dont actually need this one for now as it get overwritten by copy at
	// the next iteration **
      //Prod_0_p1[thread_id].row(nOcc).swap(Prod_0_p1[thread_id].row(nOcc+b)); 
      //Prod_0_p1[thread_id].col(nOcc).swap(Prod_0_p1[thread_id].col(nOcc+c)); 
      //Prod_0_m1[thread_id].row(nOcc).swap(Prod_0_m1[thread_id].row(nOcc+b)); 
      //Prod_0_m1[thread_id].col(nOcc).swap(Prod_0_m1[thread_id].col(nOcc+c)); 
     
        double OvLp_Swap_0_p1 = OvLp_AC_0_p1*OvLp_BD_0_p1 + OvLp_AD_0_p1*OvLp_BC_0_p1;
        double OvLp_Swap_0_m1 = OvLp_AC_0_m1*OvLp_BD_0_m1 + OvLp_AD_0_m1*OvLp_BC_0_m1;
     
        double fact = 1;
        if( a == b ) fact *= std::sqrt(0.5);
        if( c == d ) fact *= std::sqrt(0.5);
     
        #pragma omp critical
        {
          NACME(iSt,jSt) += fact*T_0(ab,iSt)*(
            T_p1(cd,jSt)*OvLp_Swap_0_p1 - 
            T_m1(cd,jSt)*OvLp_Swap_0_m1
          ) / (2 * this->step);
        }
      }
    }
  }
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "ES -> ES NACME took " << elapsed.count() << " s" << endl;

  prettyPrint(cout,NACME,"ES->ES");
  prettyPrint(cout,NACME+NACME.transpose(),"SUM");
  return NACME;
};
