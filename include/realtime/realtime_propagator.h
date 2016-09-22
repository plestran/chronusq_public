template <typename T>
void RealTime<T>::formUTrans() {

  auto NB  = ssPropagator_->basisset()->nBasis();
  auto NBT = NB * ssPropagator_->nTCS();

  ComplexMap UTransScalar(UTransScalar_,NB,NB);
  ComplexMap UTransMz(UTransMz_,0,0);
  ComplexMap UTransMy(UTransMy_,0,0);
  ComplexMap UTransMx(UTransMx_,0,0);

  if(ssPropagator_->nTCS() == 2 or !ssPropagator_->isClosedShell){
    new (&UTransMz) ComplexMap(UTransMz_,NB,NB);
  }
  if(ssPropagator_->nTCS() == 2){
    new (&UTransMy) ComplexMap(UTransMy_,NB,NB);
    new (&UTransMx) ComplexMap(UTransMx_,NB,NB);
  }

  if( iMethFormU_ == EigenDecomp ) {
//  prettyPrintSmart(cout,*this->ssPropagator_->fockOrtho()[0],"OFS");
//  this->ssPropagator_->fockOrtho()[0]->printMATLAB(cout);
    ssPropagator_->populateMO4Diag();
    ssPropagator_->diagFock2();
    
    ComplexMap S(NBTSqScratch_,NBT,NBT);

    std::copy_n(ssPropagator_->moA()->data(),NBT*NBT,NBTSqScratch_);

//  prettyPrintSmart(cout,(*ssPropagator_->epsA()),"EPSA");

    for(auto i = 0; i < NBT; i++) {
      double arg = deltaT_ * (*ssPropagator_->epsA())(i);
      S.col(i) *= dcomplex(std::cos(arg),-std::sin(arg));
    }

    if(this->ssPropagator_->nTCS() == 2) {
      ComplexMap S2(NBTSqScratch_,NBT,NBT);
      S2.noalias() = S * ssPropagator_->moA()->adjoint();

      std::vector<std::reference_wrapper<ComplexMap>> Us;
      Us.emplace_back(UTransScalar);
      Us.emplace_back(UTransMz);
      Us.emplace_back(UTransMy);
      Us.emplace_back(UTransMx);

      Quantum<dcomplex>::spinScatter(S2,Us);
    } else {
      // Temporarily store UA in UScalar
      UTransScalar.noalias() = S * ssPropagator_->moA()->adjoint();

      if(!ssPropagator_->isClosedShell) {
        std::copy_n(ssPropagator_->moB()->data(),NBT*NBT,NBTSqScratch_);
       
        for(auto i = 0; i < NBT; i++) {
          double arg = deltaT_ * (*ssPropagator_->epsB())(i);
          S.col(i) *= dcomplex(std::cos(arg),-std::sin(arg));
        }

        // Temporarily store UB in UMz
        UTransMz.noalias() = S * ssPropagator_->moB()->adjoint();

        // Reconstruct UScalar and UMz from UA and UB
        S = UTransScalar;
        UTransScalar.noalias() = S + UTransMz;
        UTransMz.noalias()     = S - UTransMz; 
      } else {
        // UScalar = 2*UA for restricted
        UTransScalar *= 2;
      }
//    cout << deltaT_ << endl;
//    prettyPrintSmart(cout,UTransScalar,"UT");
//    prettyPrintSmart(cout,UTransScalar.adjoint() * UTransScalar,"UUT");
//    CErr();
    }
  }
};

template <typename T>
void RealTime<T>::propDen() {
  auto NB  = ssPropagator_->basisset()->nBasis();

  ComplexMap UTransScalar(UTransScalar_,NB,NB);
  ComplexMap UTransMz(UTransMz_,0,0);
  ComplexMap UTransMy(UTransMy_,0,0);
  ComplexMap UTransMx(UTransMx_,0,0);

  if(ssPropagator_->nTCS() == 2 or !ssPropagator_->isClosedShell){
    new (&UTransMz) ComplexMap(UTransMz_,NB,NB);
  }
  if(ssPropagator_->nTCS() == 2){
    new (&UTransMy) ComplexMap(UTransMy_,NB,NB);
    new (&UTransMx) ComplexMap(UTransMx_,NB,NB);
  }

  ComplexMap S(NBSqScratch_,NB,NB);

  // FIXME: This only works for RHF
  S.noalias() = UTransScalar * (*ssPropagator_->onePDMOrthoScalar());
  ssPropagator_->onePDMOrthoScalar()->noalias() = S * UTransScalar.adjoint();
  (*ssPropagator_->onePDMOrthoScalar()) /= 4.0;
  
};

