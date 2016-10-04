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
//    prettyPrintSmart(cout,UTransMz,"UT");
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

/*
  // FIXME: This only works for RHF
  S.noalias() = UTransScalar * (*ssPropagator_->onePDMOrthoScalar());
  ssPropagator_->onePDMOrthoScalar()->noalias() = S * UTransScalar.adjoint();
  (*ssPropagator_->onePDMOrthoScalar()) /= 4.0;
*/
  
  // S = U**H(S) * PO(S) + U**H(k) * PO(k)
  S.noalias() = UTransScalar * (*ssPropagator_->onePDMOrthoScalar());
  if(ssPropagator_->nTCS() == 2 or !ssPropagator_->isClosedShell)
    S.noalias() += UTransMz * (*ssPropagator_->onePDMOrthoMz());
  if(ssPropagator_->nTCS() == 2){
    S.noalias() += UTransMy * (*ssPropagator_->onePDMOrthoMy());
    S.noalias() += UTransMx * (*ssPropagator_->onePDMOrthoMx());
  }

  // P(S) = S * U(S)
  // P(k) = S * U(k)
  ssPropagator_->onePDMScalar()->noalias() = S * UTransScalar.adjoint();
  if(ssPropagator_->nTCS() == 2 or !ssPropagator_->isClosedShell)
    ssPropagator_->onePDMMz()->noalias() = S * UTransMz.adjoint();
  if(ssPropagator_->nTCS() == 2){
    ssPropagator_->onePDMMy()->noalias() = S * UTransMy.adjoint();
    ssPropagator_->onePDMMx()->noalias() = S * UTransMx.adjoint();
  }

  // S = U**H(S) * PO(z) + U**H(z) * PO(S)
  if(ssPropagator_->nTCS() == 2 or !ssPropagator_->isClosedShell) {
    S.noalias() =  UTransScalar * (*ssPropagator_->onePDMOrthoMz());
    S.noalias() += UTransMz     * (*ssPropagator_->onePDMOrthoScalar());
    // FIXME: Need to add the i * Eps(z,i,j) U**H(i) * PO(j)

    // P(S) += S * U(z)
    // P(z) += S * U(S)
    ssPropagator_->onePDMScalar()->noalias() += S * UTransMz.adjoint();
    ssPropagator_->onePDMMz()->noalias() += S * UTransScalar.adjoint();
  }

  // Copy P -> PO
  ssPropagator_->cpyAOtoOrthoDen(); 

  // Correct factors of 2
  for(auto iODen = 0; iODen < POSav_.size(); iODen++){
    *ssPropagator_->onePDMOrtho()[iODen] *= 0.25;
  }
  
};

