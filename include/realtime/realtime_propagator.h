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

  ssPropagator_->populateMO4Diag();

  ComplexMap S(NBTSqScratch_,NBT,NBT);
  ComplexMap S2(NBTSqScratch2_,NBT,NBT);
  ComplexMap S3(NBTSqScratch3_,NBT,NBT);

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
    }
  } else if( iMethFormU_ == Taylor ) {
  }

  using boost::math::factorial;
  ComplexMatrix TEMP(NBT,NBT);
  TEMP.setZero();

  int NTaylor = 100;
  
//prettyPrintSmart(cout,*ssPropagator_->fockOrtho()[0],"FO After");

  double gamma = (*groundState_->epsA())(NBT-1) - (*groundState_->epsA())(0);
  gamma *= 2.5 / 2;
  double EMin = (*groundState_->epsA())(0);

  double alpha = gamma * deltaT_ / 2;

  std::vector<double> CK(NTaylor);
/*
  for(auto i = 0; i < NTaylor; i++)
    CK[i] = (std::pow(gamma * deltaT_ / 2,i) / factorial<double>(i));
*/

  auto binomial = [](int n, int m) -> double {
    return factorial<double>(n) / 
           factorial<double>(m) / factorial<double>(n-m);
  };

  for(auto iCheb = 0; iCheb < NTaylor; iCheb++){
    double tmp = 0.0;
    for(auto jCheb = iCheb; jCheb < NTaylor; jCheb += 2){
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


  NTaylor = 0;
  for(auto x : CK) { NTaylor++; if(x < 1e-15) break; }
  if(NTaylor % 2 != 0) NTaylor++;
  cout << "Keeping " << NTaylor << " terms" << endl;

  // to Alpha
  (*ssPropagator_->fockOrtho()[0]) *= 0.5;

  // Scale matrix for eigenvalues
  (*ssPropagator_->fockOrtho()[0]).noalias() -= 
    (gamma/2 + EMin)*ComplexMatrix::Identity(NBT,NBT);

  (*ssPropagator_->fockOrtho()[0]) *= 2 / gamma;

/*
  TEMP = ComplexMatrix::Identity(NBT,NBT);
  S2 = ComplexMatrix::Identity(NBT,NBT);
  S = S2;

  (*ssPropagator_->fockOrtho()[0]) *= -dcomplex(0,1);
  for(auto iT = 1; iT <= NTaylor; iT++) {
//  double alpha = (std::pow(gamma * deltaT_ / 2,iT) / factorial<double>(iT));
    if(CK[iT] < 1e-15) {break;}
    S2.noalias() = S * (*ssPropagator_->fockOrtho()[0]);
    S = S2;
    TEMP.noalias() +=  CK[iT] * S2;
  }
*/

  TEMP.noalias() = CK[0] * ComplexMatrix::Identity(NBT,NBT);
  S3.noalias()   = -dcomplex(0,1) * (*ssPropagator_->fockOrtho()[0]);

  TEMP.noalias() += CK[1] * S3;
  S.noalias()    =  CK[NTaylor/2 + 1] * S3;

  for(auto iT = 2; iT <= (NTaylor/2); iT++) {
    if(iT % 2 == 0) {
      S2.noalias() = dcomplex(0,-1) * (*ssPropagator_->fockOrtho()[0]) * S3;
      TEMP.noalias() += CK[iT] * S2;
      S.noalias()    += CK[iT + NTaylor/2] * S2;
    } else {
      S3.noalias() = dcomplex(0,-1) * (*ssPropagator_->fockOrtho()[0]) * S2;
      TEMP.noalias() += CK[iT] * S3;
      S.noalias()    += CK[iT + NTaylor/2] * S3;
    }
  }

  if((NTaylor / 2) % 2 == 0) TEMP.noalias() += S2 * S;
  else                       TEMP.noalias() += S3 * S;

  TEMP *= std::exp(-dcomplex(0,gamma/2 + EMin)*deltaT_ );

  TEMP *= 2;
//prettyPrintSmart(this->fileio_->out,UTransScalar,"True");
//prettyPrintSmart(this->fileio_->out,TEMP,"Taylor");
  prettyPrintSmart(cout,UTransScalar - TEMP,"DIFF");
//prettyPrintSmart(cout,UTransScalar.cwiseQuotient(TEMP).cwiseAbs(),"Q");
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

