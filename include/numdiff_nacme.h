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

  std::vector<RealMatrix> SWAPPED_0(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA())));
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

  std::vector<RealMatrix> SWAPPED_0(nThreads,
    TMatrix((*this->singleSlater_undisplaced_->moA())));
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

/*
  this->checkPhase(resp_p1,resp_m1);
  this->checkPhase((*this->response_undisplaced_),resp_p1);
  this->checkPhase((*this->response_undisplaced_),resp_m1);
*/

  std::vector<RealMatrix> SWAPPED_IA_0(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));
  std::vector<RealMatrix> SWAPPED_JB_p1(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));
  std::vector<RealMatrix> SWAPPED_JB_m1(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));

  std::vector<RealMatrix> Prod_0_p1(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));
  std::vector<RealMatrix> Prod_0_m1(nThreads,RealMatrix((*this->singleSlater_undisplaced_->moA())));

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

//cout << NACME << endl << endl;
  prettyPrint(cout,NACME,"ES->ES");
  return NACME;

};
