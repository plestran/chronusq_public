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

  ComplexMap S(NBTSqScratch_,NBT,NBT);
  ComplexMap S2(NBTSqScratch2_,NBT,NBT);
  ComplexMap S3(NBTSqScratch3_,NBT,NBT);

  // Scaling Parameters for Polynomial Expansion
  double gamma = (*groundState_->epsA())(NBT-1) - (*groundState_->epsA())(0);
  double EMin = (*groundState_->epsA())(0);

  if(ssPropagator_->nTCS() == 1 and !ssPropagator_->isClosedShell){
    gamma = std::max(gamma,
      (*groundState_->epsB())(NBT-1) - (*groundState_->epsB())(0));
    EMin = std::min(EMin,(*groundState_->epsB())(0));
  }

  gamma *= 3. / 2;

  double alpha = gamma * deltaT_ / 2;

  std::vector<double> CK(nPolyExpMax_);
  int nPolyKeep(0);

  // Copy things into full dimension
  ssPropagator_->populateMO4Diag();

  // Functions for polynomial
  using boost::math::factorial;

  auto binomial = [](int n, int m) -> double {
    return factorial<double>(n) / 
           factorial<double>(m) / factorial<double>(n-m);
  };

  if( iMethFormU_ == EigenDecomp ) {
    ssPropagator_->diagFock2();

    std::copy_n(ssPropagator_->moA()->data(),NBT*NBT,NBTSqScratch_);

    for(auto i = 0; i < NBT; i++) {
      double arg = deltaT_ * (*ssPropagator_->epsA())(i);
      S.col(i) *= dcomplex(std::cos(arg),-std::sin(arg));
    }

    if(this->ssPropagator_->nTCS() == 2) {
      S2.noalias() = S * ssPropagator_->moA()->adjoint();

      std::vector<std::reference_wrapper<ComplexMap>> Us;
      Us.emplace_back(UTransScalar);
      Us.emplace_back(UTransMz);
      Us.emplace_back(UTransMy);
      Us.emplace_back(UTransMx);

      Quantum<dcomplex>::spinScatter(S2,Us);
    } else {
      // Temporarily store UA in UTransScalar
      UTransScalar.noalias() = S * ssPropagator_->moA()->adjoint();

      if(!ssPropagator_->isClosedShell) {
        std::copy_n(ssPropagator_->moB()->data(),NBT*NBT,NBTSqScratch_);
       
        for(auto i = 0; i < NBT; i++) {
          double arg = deltaT_ * (*ssPropagator_->epsB())(i);
          S.col(i) *= dcomplex(std::cos(arg),-std::sin(arg));
        }

        // Temporarily store UB in UMz
        UTransMz.noalias() = S * ssPropagator_->moB()->adjoint();

        // Reconstruct UTransScalar and UMz from UA and UB
        S = UTransScalar;
        UTransScalar.noalias() = S + UTransMz;
        UTransMz.noalias()     = S - UTransMz; 
      } else {
        // UTransScalar = 2*UA for restricted
        UTransScalar *= 2;
      }
    }
  } else if( iMethFormU_ == Taylor ) {

    for(auto i = 0; i < nPolyExpMax_; i++)
      CK[i] = (std::pow(gamma * deltaT_ / 2,i) / factorial<double>(i));

  } else if( iMethFormU_ == Chebyshev ) {

    for(auto iCheb = 0; iCheb < nPolyExpMax_; iCheb++){
      double tmp = 0.0;
      for(auto jCheb = iCheb; jCheb < nPolyExpMax_; jCheb += 2){
        if(iCheb == 0) {
          if(jCheb == 0) tmp += boost::math::cyl_bessel_j(jCheb,alpha);
          else           tmp += 2.0 * boost::math::cyl_bessel_j(jCheb,alpha);
        } else {
          tmp += 2 * std::pow(2.0,iCheb - 1) * jCheb / iCheb *
            binomial((iCheb + jCheb)/2 - 1, iCheb -1) *
            boost::math::cyl_bessel_j(jCheb,alpha);
        }
      }

      CK[iCheb] = tmp;
    }

  }

  if( iMethFormU_ == Taylor or iMethFormU_ == Chebyshev ){
    nPolyKeep = 0;
    for(auto x : CK) { nPolyKeep++; if(x < polyEps_) break; }
    if(nPolyKeep % 2 != 0) nPolyKeep++;

/*
    // to Alpha
    (*ssPropagator_->moA()) *= 0.5;
*/
 
    // Scale matrix for eigenvalues
    (*ssPropagator_->moA()).noalias() -= 
      (gamma/2 + EMin)*ComplexMatrix::Identity(NBT,NBT);
 
    (*ssPropagator_->moA()) *= 2 / gamma;

    // Zeroth and setup for first term
    UTransScalar.noalias() = CK[0] * ComplexMatrix::Identity(NBT,NBT);
    S3.noalias()   = -dcomplex(0,1) * (*ssPropagator_->moA());
 
    // First term
    UTransScalar.noalias() += CK[1] * S3;
    S.noalias()    =  CK[nPolyKeep/2 + 1] * S3;
 
    // Loop over half the kept points
    for(auto iT = 2; iT <= (nPolyKeep/2); iT++) {
      if(iT % 2 == 0) {
        S2.noalias() = dcomplex(0,-1) * (*ssPropagator_->moA()) * S3;
        UTransScalar.noalias() += CK[iT] * S2;
        S.noalias()    += CK[iT + nPolyKeep/2] * S2;
      } else {
        S3.noalias() = dcomplex(0,-1) * (*ssPropagator_->moA()) * S2;
        UTransScalar.noalias() += CK[iT] * S3;
        S.noalias()    += CK[iT + nPolyKeep/2] * S3;
      }
    }
 
    if((nPolyKeep / 2) % 2 == 0) UTransScalar.noalias() += S2 * S;
    else                         UTransScalar.noalias() += S3 * S;
 
    UTransScalar *= std::exp(-dcomplex(0,gamma/2 + EMin)*deltaT_ );

    if(ssPropagator_->nTCS() == 1 and !ssPropagator_->isClosedShell) {
      // Scale matrix for eigenvalues
      (*ssPropagator_->moB()).noalias() -= 
        (gamma/2 + EMin)*ComplexMatrix::Identity(NBT,NBT);
 
      (*ssPropagator_->moB()) *= 2 / gamma;

      // Zeroth and setup for first term
      UTransMz.noalias() = CK[0] * ComplexMatrix::Identity(NBT,NBT);
      S3.noalias()   = -dcomplex(0,1) * (*ssPropagator_->moB());
 
      // First term
      UTransMz.noalias() += CK[1] * S3;
      S.noalias()    =  CK[nPolyKeep/2 + 1] * S3;
 
      // Loop over half the kept points
      for(auto iT = 2; iT <= (nPolyKeep/2); iT++) {
        if(iT % 2 == 0) {
          S2.noalias() = dcomplex(0,-1) * (*ssPropagator_->moB()) * S3;
          UTransMz.noalias() += CK[iT] * S2;
          S.noalias()    += CK[iT + nPolyKeep/2] * S2;
        } else {
          S3.noalias() = dcomplex(0,-1) * (*ssPropagator_->moB()) * S2;
          UTransMz.noalias() += CK[iT] * S3;
          S.noalias()    += CK[iT + nPolyKeep/2] * S3;
        }
      }
    }
 
    if(ssPropagator_->nTCS() == 2) {
      // NYI
    } else if(!ssPropagator_->isClosedShell) {
      // Reconstruct UTransScalar and UMz from UA and UB
      S = UTransScalar;
      UTransScalar.noalias() = S + UTransMz;
      UTransMz.noalias()     = S - UTransMz; 
    } else {
      UTransScalar *= 2;
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

  ComplexMap S(NBTSqScratch_,NB,NB);

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

