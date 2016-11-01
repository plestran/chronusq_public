template <typename T>
void MOIntegrals<T>::testMOInts(){
  T EMP2 = 0;
  int NO = this->wfn_->nO();
  int NV = this->wfn_->nV();
  int NOA = this->wfn_->nOA();
  int NVA = this->wfn_->nVA();
  int NOB = this->wfn_->nOB();
  int NVB = this->wfn_->nVB();

  double eI, eJ, eA, eB;

  bool doMP2 = true;
  bool doMP3 = true;
  bool doSpacial = false;

  bool isR = this->wfn_->isClosedShell and this->wfn_->nTCS() == 1;
  bool isU = !isR and this->wfn_->nTCS() == 1;
  bool isG = this->wfn_->nTCS() == 2;

  double * EPSA = this->wfn_->epsA()->data();
  double * EPSB = isR or isG ? EPSA : this->wfn_->epsB()->data();

/*
  if(!doMP2) return;

  if(doSpacial)
    this->formVOVO();

  if(not doSpacial or doMP3)
    this->formFullVOVO();
//return;
/
*/
  this->formFullVOVO();

  if(isR and doSpacial) {
    cout << "Doing RHF Spacial MP2" << endl;
    // RHF Spacial orbital MP2 energy
    for(auto j = 0; j < NOA; j++)
    for(auto i = 0; i < NOA; i++)
    for(auto b = 0; b < NVA; b++)
    for(auto a = 0; a < NVA; a++){
      EMP2 += VOVOAAAA_[a + NVA*i + NVA*NOA*b + NVA*NVA*NOA*j] *
        ( 2.0*VOVOAAAA_[a + NVA*i + NVA*NOA*b + NVA*NVA*NOA*j]
            - VOVOAAAA_[b + NVA*i + NVA*NOA*a + NVA*NVA*NOA*j]) /
      ((*this->wfn_->epsA())(i) + (*this->wfn_->epsA())(j)
      -(*this->wfn_->epsA())(a+NOA) - (*this->wfn_->epsA())(b+NOA));
    }
  } else if(isU and doSpacial) {
    cout << "Doing UHF Spacial MP2" << endl;
    // UHF Spacial orbital MP2 Energy
    // AAAA Coulomb
    EMP2 = 0;
    for(auto j = 0; j < NOA; j++)
    for(auto b = 0; b < NVA; b++)
    for(auto i = 0; i < NOA; i++)
    for(auto a = 0; a < NVA; a++){
      eI = EPSA[i]      ; 
      eJ = EPSA[j]      ; 
      eA = EPSA[a + NOA]; 
      eB = EPSA[b + NOA]; 

      // dIJAB = e(I) + e(J) - e(A) - e(B)
      double deltaIJAB = eI + eJ - eA - eB;

      T CoulIJAB = VOVOAAAA_[a + NVA*i + NVA*NOA*b + NVA*NVA*NOA*j];
 
      EMP2 += CoulIJAB*CoulIJAB / deltaIJAB;

    }

    // AABB Coulomb
    for(auto j = 0; j < NOB; j++)
    for(auto b = 0; b < NVB; b++)
    for(auto i = 0; i < NOA; i++)
    for(auto a = 0; a < NVA; a++){
      eI = EPSA[i]      ; 
      eA = EPSA[a + NOA]; 
      eJ = EPSB[j]      ; 
      eB = EPSB[b + NOB]; 

      // dIJAB = e(I) + e(J) - e(A) - e(B)
      double deltaIJAB = eI + eJ - eA - eB;

      T CoulIJAB = VOVOAABB_[a + NVA*i + NVA*NOA*b + NVA*NVB*NOA*j];
 
      EMP2 += CoulIJAB*CoulIJAB / deltaIJAB;

    }

    // BBAA Coulomb
    for(auto j = 0; j < NOB; j++)
    for(auto b = 0; b < NVB; b++)
    for(auto i = 0; i < NOA; i++)
    for(auto a = 0; a < NVA; a++){
      eI = EPSA[i]      ; 
      eA = EPSA[a + NOA]; 
      eJ = EPSB[j]      ; 
      eB = EPSB[b + NOB]; 

      // dIJAB = e(I) + e(J) - e(A) - e(B)
      double deltaIJAB = eI + eJ - eA - eB;

      T CoulIJAB = VOVOAABB_[a + NVA*i + NVA*NOA*b + NVA*NVB*NOA*j];
 
      EMP2 += CoulIJAB*CoulIJAB / deltaIJAB;

    }

    // BBBB Coulomb
    for(auto j = 0; j < NOB; j++)
    for(auto b = 0; b < NVB; b++)
    for(auto i = 0; i < NOB; i++)
    for(auto a = 0; a < NVB; a++){
      eI = EPSB[i]      ; 
      eJ = EPSB[j]      ; 
      eA = EPSB[a + NOB]; 
      eB = EPSB[b + NOB]; 

      // dIJAB = e(I) + e(J) - e(A) - e(B)
      double deltaIJAB = eI + eJ - eA - eB;

      T CoulIJAB = VOVOBBBB_[a + NVB*i + NVB*NOB*b + NVB*NVB*NOB*j];
 
      EMP2 += CoulIJAB*CoulIJAB / deltaIJAB;

    }

    // AAAA Exchange
    for(auto j = 0; j < NOA; j++)
    for(auto b = 0; b < NVA; b++)
    for(auto i = 0; i < NOA; i++)
    for(auto a = 0; a < NVA; a++){
      eI = EPSA[i]      ; 
      eJ = EPSA[j]      ; 
      eA = EPSA[a + NOA]; 
      eB = EPSA[b + NOA]; 

      // dIJAB = e(I) + e(J) - e(A) - e(B)
      double deltaIJAB = eI + eJ - eA - eB;

      T CoulIJAB = VOVOAAAA_[a + NVA*i + NVA*NOA*b + NVA*NVA*NOA*j];
      T ExchIJAB = VOVOAAAA_[a + NVA*j + NVA*NOA*b + NVA*NVA*NOA*i];
 
      EMP2 -= CoulIJAB*ExchIJAB / deltaIJAB;

    }

    // BBBB Exchange
    for(auto j = 0; j < NOB; j++)
    for(auto b = 0; b < NVB; b++)
    for(auto i = 0; i < NOB; i++)
    for(auto a = 0; a < NVB; a++){
      eI = EPSB[i]      ; 
      eJ = EPSB[j]      ; 
      eA = EPSB[a + NOB]; 
      eB = EPSB[b + NOB]; 

      // dIJAB = e(I) + e(J) - e(A) - e(B)
      double deltaIJAB = eI + eJ - eA - eB;

      T CoulIJAB = VOVOBBBB_[a + NVB*i + NVB*NOB*b + NVB*NVB*NOB*j];
      T ExchIJAB = VOVOBBBB_[a + NVB*j + NVB*NOB*b + NVB*NVB*NOB*i];
 
      EMP2 -= CoulIJAB*ExchIJAB / deltaIJAB;

    }

    EMP2 *= 0.5;
  } else {
    cout << "Doing Spin-Orbital MP2" << endl;
    // Spin-orbital MP2
    for(int j = 0; j < NO; j++)
    for(int i = 0; i < NO; i++)
    for(int b = 0; b < NV; b++)
    for(int a = 0; a < NV; a++){
      
      eI = isG ? EPSA[i] : (i % 2 == 0) ? EPSA[i/2]       : EPSB[i/2];
      eJ = isG ? EPSA[j] : (j % 2 == 0) ? EPSA[j/2]       : EPSB[j/2];
      eA = isG ? EPSA[a+NO] : (a % 2 == 0) ? EPSA[a/2 + NOA] : EPSB[a/2 + NOB];
      eB = isG ? EPSA[b+NO] : (b % 2 == 0) ? EPSA[b/2 + NOA] : EPSB[b/2 + NOB];


      // < IJ || AB > = (AI | BJ) - (BI | AJ)
      T DiracIJAB = 
        VOVO_[a + i*NV + b*NO*NV + j*NO*NV*NV] -
        VOVO_[b + i*NV + a*NO*NV + j*NO*NV*NV];

      // dIJAB = e(I) + e(J) - e(A) - e(B)
      double deltaIJAB = eI + eJ - eA - eB;

      // E(2) += |< IJ || AB >|^2 / dIJAB
      EMP2 += DiracIJAB * std::conj(DiracIJAB) / deltaIJAB;
    }
    
    // E(2) = E(2) / 4
    EMP2 *= 0.25;
  }
  cout << "EMP2 = " << EMP2 << endl;

  if(not doMP3) return;
  cout << "Doing Spin-Orbital MP3" << endl;

  this->formFullVVOO();
  this->formFullVVVV();
  this->formFullOOOO();

  T EMP3_1 = 0;
  T EMP3_2 = 0;
  T EMP3_3 = 0;

  for(auto i = 0; i < NO; i++)
  for(auto j = 0; j < NO; j++)
  for(auto k = 0; k < NO; k++)
  for(auto l = 0; l < NO; l++)
  for(auto a = 0; a < NV; a++)
  for(auto b = 0; b < NV; b++){

    double eI = isG ? (*this->wfn_->epsA())(i)      :  (*this->wfn_->epsA())(i/2);
    double eJ = isG ? (*this->wfn_->epsA())(j)      :  (*this->wfn_->epsA())(j/2);
    double eK = isG ? (*this->wfn_->epsA())(k)      :  (*this->wfn_->epsA())(k/2);
    double eL = isG ? (*this->wfn_->epsA())(l)      :  (*this->wfn_->epsA())(l/2);
    double eA = isG ? (*this->wfn_->epsA())(a + NO) : (*this->wfn_->epsA())(a/2 + NOA);
    double eB = isG ? (*this->wfn_->epsA())(b + NO) : (*this->wfn_->epsA())(b/2 + NOA);

    // < IJ || AB > = (AI | BJ) - (BI | AJ)
    T DiracIJAB = 
      VOVO_[a + i*NV + b*NO*NV + j*NV*NV*NO] - 
      VOVO_[b + i*NV + a*NO*NV + j*NV*NV*NO]; 

    // < IJ || KL > = (IK | JL) - (IL | JK)
    T DiracIJKL = 
      OOOO_[i + k*NO + j*NO*NO + l*NO*NO*NO] - 
      OOOO_[i + l*NO + j*NO*NO + k*NO*NO*NO]; 

    // < KL || AB > = (AK | BL) - (BK | AL)
    T DiracKLAB = 
      VOVO_[a + k*NV + b*NO*NV + l*NV*NV*NO] - 
      VOVO_[b + k*NV + a*NO*NV + l*NV*NV*NO]; 


    // dIJAB = e(I) + e(J) - e(A) - e(B)
    double deltaIJAB = eI + eJ - eA - eB;

    // dKLAB = e(K) + e(L) - e(A) - e(B)
    double deltaKLAB= eK + eL - eA - eB;

//  cout << DiracIJAB << " " << DiracIJKL << " " << DiracKLAB << endl;

    // E(3,1) += < IJ || AB > * < IJ || KL > * < KL || AB > / dIJAB / dKLAB 
    EMP3_1 += DiracIJAB * DiracIJKL * DiracKLAB / deltaIJAB / deltaKLAB;
//  cout << i << " " << j << " " << k << " " << l << " " << a << " " << b <<endl;
//  cout << DiracIJAB << " " << DiracIJKL << " " << DiracKLAB << endl;
//  cout << deltaIJAB << " " << deltaKLAB << endl;
//  cout << EMP3_1 << endl;
  }

  for(auto i = 0; i < NO; i++)
  for(auto j = 0; j < NO; j++)
  for(auto a = 0; a < NV; a++)
  for(auto b = 0; b < NV; b++) 
  for(auto c = 0; c < NV; c++)
  for(auto d = 0; d < NV; d++){

    double eI = isG ? (*this->wfn_->epsA())(i)      : (*this->wfn_->epsA())(i/2);
    double eJ = isG ? (*this->wfn_->epsA())(j)      : (*this->wfn_->epsA())(j/2);
    double eA = isG ? (*this->wfn_->epsA())(a + NO) : (*this->wfn_->epsA())(a/2 + NOA);
    double eB = isG ? (*this->wfn_->epsA())(b + NO) : (*this->wfn_->epsA())(b/2 + NOA);
    double eC = isG ? (*this->wfn_->epsA())(c + NO) : (*this->wfn_->epsA())(c/2 + NOA);
    double eD = isG ? (*this->wfn_->epsA())(d + NO) : (*this->wfn_->epsA())(d/2 + NOA);

    // < IJ || AB > = (AI | BJ) - (BI | AJ)
    T DiracIJAB = 
      VOVO_[a + i*NV + b*NO*NV + j*NV*NV*NO] - 
      VOVO_[b + i*NV + a*NO*NV + j*NV*NV*NO]; 

    // < AB || CD > = (AC | BD) - (AD | BC)
    T DiracABCD = 
      VVVV_[a + c*NV + b*NV*NV + d*NV*NV*NV] - 
      VVVV_[a + d*NV + b*NV*NV + c*NV*NV*NV]; 

    // < IJ || CD > = (CI | DJ) - (DI | CJ)
    T DiracIJCD = 
      VOVO_[c + i*NV + d*NO*NV + j*NV*NV*NO] - 
      VOVO_[d + i*NV + c*NO*NV + j*NV*NV*NO]; 


    // dIJAB = e(I) + e(J) - e(A) - e(B)
    double deltaIJAB = eI + eJ - eA - eB;

    // dIJCD = e(I) + e(J) - e(C) - e(D)
    double deltaIJCD = eI + eJ - eC - eD; 

    // E(3,2) += < IJ || AB > * < AB || CD > * < IJ || CD > / dIJAB / dIJCD 
    EMP3_2 += DiracIJAB * DiracABCD * DiracIJCD / deltaIJAB / deltaIJCD;
    
  }

  for(auto i = 0; i < NO; i++)
  for(auto j = 0; j < NO; j++)
  for(auto k = 0; k < NO; k++)
  for(auto a = 0; a < NV; a++)
  for(auto b = 0; b < NV; b++) 
  for(auto c = 0; c < NV; c++){

    double eI = isG ?(*this->wfn_->epsA())(i)     : (*this->wfn_->epsA())(i/2);       
    double eJ = isG ?(*this->wfn_->epsA())(j)     : (*this->wfn_->epsA())(j/2);
    double eK = isG ?(*this->wfn_->epsA())(k)     : (*this->wfn_->epsA())(k/2);
    double eA = isG ?(*this->wfn_->epsA())(a + NO) : (*this->wfn_->epsA())(a/2 + NOA);
    double eB = isG ?(*this->wfn_->epsA())(b + NO) : (*this->wfn_->epsA())(b/2 + NOA);
    double eC = isG ?(*this->wfn_->epsA())(c + NO) : (*this->wfn_->epsA())(c/2 + NOA);

    // < IJ || AB > = (AI | BJ) - (BI | AJ)
    T DiracIJAB = 
      VOVO_[a + i*NV + b*NO*NV + j*NV*NV*NO] - 
      VOVO_[b + i*NV + a*NO*NV + j*NV*NV*NO]; 

    // < BK || CJ > = (BC | KJ) - (BJ | CK)
    T DiracBKCJ = 
      VVOO_[b + c*NV + k*NV*NV + j*NV*NV*NO] -
      VOVO_[b + j*NV + c*NO*NV + k*NO*NV*NV];

    // < IK || AC > = (AI | CK) - (CI | AK)
    T DiracIKAC = 
      VOVO_[a + i*NV + c*NO*NV + k*NV*NV*NO] - 
      VOVO_[c + i*NV + a*NO*NV + k*NV*NV*NO]; 


    // dIJAB = e(I) + e(J) - e(A) - e(B)
    double deltaIJAB = eI + eJ - eA - eB;

    // dIKAC = e(I) + e(K) - e(A) - e(C)
    double deltaIKAC = eI + eK - eA - eC; 

    // E(3,3) -= < IJ || AB > * < BK || CJ > * < IK || AC > / dIJAB / dIKAC
    EMP3_3 -= DiracIJAB * DiracBKCJ * DiracIKAC / deltaIJAB / deltaIKAC;
  }

  // E(3) = (1/8) * ( E(3,1) + E(3,2) ) + E(3,3)
  T EMP3 = (0.125)*(EMP3_1 + EMP3_2) + EMP3_3;

  cout << "EMP3_1 = " << EMP3_1 << endl;
  cout << "EMP3_2 = " << EMP3_2 << endl;
  cout << "EMP3_3 = " << EMP3_3 << endl;
  cout << "EMP3 = " << EMP3 << endl;

//this->formFullVOOO();
//this->formFullVVVO();
};

