#include "aointegrals.h"
static int hRR_p[6][6]={
1,0,0,0,0,0, //x
0,0,0,1,0,0, 
0,1,0,0,0,0, //y
0,0,0,0,1,0,
0,0,1,0,0,0, //z
0,0,0,0,0,1
};

//
// overlap horizontal recursion     
// (a|b) = (A-B)(a|b-1) + (a+1|b-1) 
//   LA >= LB                       
//
double AOIntegrals::hRRSab(ShellPair *ijSP,int LA,int *lA,int LB,int *lB) {
  int i,j;
  cout<<"hRR "<<LA<<" "<<LB<<endl;
  double tmpVal = 0.0;
  if(LA+LB==0) {
    // (s|s)
    for(i=0;i<ijSP->nPGTOPair;i++) tmpVal+= ijSP->ss[i];
    return tmpVal;
  } else if(LB==0) {
    // (|s)
    for(i=0;i<ijSP->nPGTOPair;i++) tmpVal+= ijSP->ss[i]*this->vRRSa0(ijSP,LA,lA,i);
    return tmpVal;
  };

  int iWork,lAp1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAp1[iWork]=lA[iWork];     
    lBm1[iWork]=lB[iWork];     
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lAp1[iWork]++;
  lBm1[iWork]--;

  tmpVal = this->hRRSab(ijSP,LA+1,lAp1,LB-1,lBm1);
  tmpVal+= ijSP->AB[iWork]*this->hRRSab(ijSP,LA,lA,LB-1,lBm1);
  return tmpVal;
};
//----------------------------------------------------------//
// overlap vertical recursion                               //
// (a|0) = (P-A)(a-1|0) + 1/(2*(alpha+beta))*N_(a-1)*(a-2|0)//
//----------------------------------------------------------//
double AOIntegrals::vRRSa0(ShellPair *ijSP,int LA,int *lA,int iPGTOPair) {
  cout<<"vRR "<<LA<<endl;
  double tmpVal=0.0;
  int iWork;
  int lAm1[3];
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0) iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;

  if(LA==1) return ijSP->PA[iWork][iPGTOPair];

  lAm1[iWork]--;
  tmpVal = ijSP->PA[iWork][iPGTOPair]*this->vRRSa0(ijSP,LA-1,lAm1,iPGTOPair);

  if(LA==2&&lA[iWork]==2) tmpVal += ijSP->halfInverseZeta[iPGTOPair];
  else if (lA[iWork]>=2) {
    lAm1[iWork]--;
    tmpVal += (lAm1[iWork]+1)*ijSP->halfInverseZeta[iPGTOPair]*this->vRRSa0(ijSP,LA-2,lAm1,iPGTOPair);
  };

  return tmpVal;
};

// kinetic vertical recursion
// (a|T|b) = (P-B)(a|T|b-1) + 1/2*1/Zeta*N_(b-1)(a|T|b-2)
//          +1/2*1/Zeta*N_(a)(a-1|T|b-1)
double AOIntegrals::vRRTab(ShellPair *ijShellPair,int LA,int *lA,int LB,int *lB,int *i,int *j) {
  double tmpVal = 0.0;
  if(LA+LB==0) return pairConstants_->TabPar1[*i][*j];
  else if(LB==0) return this->oneevRRTa0(ijShellPair,LA,lA,i,j);
  int iWork,lAm1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAm1[iWork]=lA[iWork];
    lBm1[iWork]=lB[iWork];
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lBm1[iWork]--;
  tmpVal = pairConstants_->TabPar2[*i][*j]*this->oneehRRTSab(ijShellPair,LA,lA,LB,lB,i,j);
  if(LB>2) {
    if(abs(ijShellPair->deltaPB[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPB[iWork][*i][*j]*this->oneevRRTab(ijShellPair,LA,lA,LB-1,lBm1,i,j);
    if (lA[iWork]>0) {
      lAm1[iWork]--;
      tmpVal += (lAm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTab(ijShellPair,LA-1,lAm1,LB-1,lBm1,i,j);
    };
    if(LB==2&&lB[iWork]==2) tmpVal += pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA,lA,i,j) - pairConstants_->TabPar3[*i][*j]*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
    else if (lB[iWork]>=2) {
      lBm1[iWork]--;
      tmpVal += (lBm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTab(ijShellPair,LA,lA,LB-2,lBm1,i,j);
      tmpVal -= (lBm1[iWork]+1)*pairConstants_->TabPar3[*i][*j]*this->oneehRRTSab(ijShellPair,LA,lA,LB-2,lBm1,i,j);
    };
  } else if (LB==2) {
    int iWork2;
    if (lBm1[0]>0)      iWork2=0;
    else if (lBm1[1]>0) iWork2=1;
    else if (lBm1[2]>0) iWork2=2;
    double tmpVal2 = 0.0;
    if(abs(ijShellPair->deltaPB[iWork][*i][*j])>this->controls_->thresholdS) {
      tmpVal2 = pairConstants_->TabPar2[*i][*j]*this->oneehRRTSab(ijShellPair,LA,lA,LB-1,lBm1,i,j);
      if(abs(ijShellPair->deltaPB[iWork2][*i][*j])>this->controls_->thresholdS) tmpVal2 += ijShellPair->deltaPB[iWork2][*i][*j]*this->oneevRRTa0(ijShellPair,LA,lA,i,j);
      if(lA[iWork2]>0) {
	lAm1[iWork2]--;
	tmpVal2 += (lAm1[iWork2]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA-1,lAm1,i,j);
	lAm1[iWork2]++;
      };
      tmpVal += tmpVal2*ijShellPair->deltaPB[iWork][*i][*j];
    };
    if(lA[iWork]>0) {
      lAm1[iWork]--;
      tmpVal2 = pairConstants_->TabPar2[*i][*j]*this->oneehRRTSab(ijShellPair,LA-1,lAm1,LB-1,lBm1,i,j);
      if(abs(ijShellPair->deltaPB[iWork2][*i][*j])>this->controls_->thresholdS) tmpVal2 += ijShellPair->deltaPB[iWork2][*i][*j]*this->oneevRRTa0(ijShellPair,LA-1,lAm1,i,j);
      if(lAm1[iWork2]>0) {
	lAm1[iWork2]--;
	if(LA>2) tmpVal2 += (lAm1[iWork2]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA-2,lAm1,i,j);
	else tmpVal2 += pairConstants_->Sa0Par[*i][*j]*pairConstants_->TabPar1[*i][*j];
      };
      tmpVal += tmpVal2*lA[iWork]*pairConstants_->Sa0Par[*i][*j];
    };
    if(iWork==iWork2) tmpVal += pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA,lA,i,j) - pairConstants_->TabPar3[*i][*j]*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
  } else {
    if(abs(ijShellPair->deltaPB[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPB[iWork][*i][*j]*this->oneevRRTa0(ijShellPair,LA,lA,i,j);
    if (lA[iWork]>0) {
      lAm1[iWork]--;
      tmpVal += (lAm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA-1,lAm1,i,j);
    };
  };
  return tmpVal;
};
// kinetic vertical recursion
//  from (a|T|0) to (0|T|0)
double AOIntegrals::vRRTa0(ShellPair *ijShellPair,int LA,int *lA,int *i,int *j) {
  double tmpVal = 0.0;
  if(LA==0) return pairConstants_->TabPar1[*i][*j];
  int iWork,lAm1[3];
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0)      iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;
  lAm1[iWork]--;
  tmpVal = pairConstants_->TabPar2[*i][*j]*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
  if(LA>2) {
    if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*this->oneevRRTa0(ijShellPair,LA-1,lAm1,i,j);
    if(LA==2&&lA[iWork]==2) {
      tmpVal += pairConstants_->Sa0Par[*i][*j]*pairConstants_->TabPar1[*i][*j];
      tmpVal -= pairConstants_->Ta0Par3[*i][*j];
    } else if (lA[iWork]>=2) {
      lAm1[iWork] -=1;
      tmpVal += (lAm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA-2,lAm1,i,j);
      tmpVal -= (lAm1[iWork]+1)*pairConstants_->Ta0Par3[*i][*j]*this->oneevRRSa0(ijShellPair,LA-2,lAm1,i,j);
    };
  } else if (LA==2) {
    int iWork2;
    if (lAm1[0]>0)      iWork2=0;
    else if (lAm1[1]>0) iWork2=1;
    else if (lAm1[2]>0) iWork2=2;
    double tmpVal2 = 0.0;
    if(iWork!=iWork2) {
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) {
	if(abs(ijShellPair->deltaPA[iWork2][*i][*j])>this->controls_->thresholdS)
	  tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*(pairConstants_->TabPar2[*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j] + ijShellPair->deltaPA[iWork2][*i][*j]*pairConstants_->TabPar1[*i][*j]);
      };
    } else {
      tmpVal += pairConstants_->Sa0Par[*i][*j]*pairConstants_->TabPar1[*i][*j] - pairConstants_->Ta0Par3[*i][*j];
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) 
	tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*(ijShellPair->deltaPA[iWork][*i][*j]*pairConstants_->TabPar1[*i][*j] + pairConstants_->TabPar2[*i][*j]*ijShellPair->deltaPA[iWork][*i][*j]);
    };
  } else if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*pairConstants_->TabPar1[*i][*j];
  return tmpVal;
};




// overlap horizontal recursion
// (a|b) = (A-B)(a|b-1) + (a+1|b-1)
//   LA >= LB
double AOIntegrals::oneehRRSab(ShellPair *ijShellPair,int LA,int *lA,int LB,int *lB) {
  int i,j;
  double tmpVal = 0.0;
  if(LA+LB==0) {
    // (s|s)
    for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) 
      tmpVal+= pairConstants_->ssPair[i][j];
    return tmpVal;
  } else if(LB==0) {
    // (|s)
    for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) 
      if(pairConstants_->ssNonzero[i][j]) tmpVal+= pairConstants_->ssPair[i][j]*this->oneevRRSa0(ijShellPair,LA,lA,&i,&j);
    return tmpVal;
  };
  int iWork,lAp1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAp1[iWork]=lA[iWork];     
    lBm1[iWork]=lB[iWork];     
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  if(LB>2) {
    lAp1[iWork]++;
    lBm1[iWork]--;
    tmpVal = this->oneehRRSab(ijShellPair,LA+1,lAp1,LB-1,lBm1);
    if(abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) 
      tmpVal+= ijShellPair->deltaAB[iWork]*this->oneehRRSab(ijShellPair,LA,lA,LB-1,lBm1);
  } else if(LB==2) {
    // (|d)
    lBm1[iWork]--;
    int iWork2;
    if (lBm1[0]>0)      iWork2=0;
    else if (lBm1[1]>0) iWork2=1;
    else if (lBm1[2]>0) iWork2=2;

    double primitiveVal = 0.0;
    if(iWork==iWork2) {
      if(abs(ijShellPair->deltaAB[iWork2])>this->controls_->thresholdS) {
	double sqrDeltaAB = (ijShellPair->deltaAB[iWork])*(ijShellPair->deltaAB[iWork]);
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++)
          if(pairConstants_->ssNonzero[i][j]){
	    primitiveVal = sqrDeltaAB*this->oneevRRSa0(ijShellPair,LA,lA,&i,&j);
	    lAp1[iWork]++;
	    primitiveVal += math.two*ijShellPair->deltaAB[iWork]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,&i,&j);
	    lAp1[iWork]++;
	    primitiveVal += this->oneevRRSa0(ijShellPair,LA+2,lAp1,&i,&j);
	    lAp1[iWork]-=2;
	    tmpVal += primitiveVal*pairConstants_->ssPair[i][j];
	  };
      } else {
	lAp1[iWork]+=2;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) 
          if(pairConstants_->ssNonzero[i][j]) tmpVal += pairConstants_->ssPair[i][j]*this->oneevRRSa0(ijShellPair,LA+2,lAp1,&i,&j);
      };
    } else {
      if(abs(ijShellPair->deltaAB[iWork2])<this->controls_->thresholdS&&abs(ijShellPair->deltaAB[iWork])<this->controls_->thresholdS) {
	lAp1[iWork]++;
	lAp1[iWork2]++;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) 
          if(pairConstants_->ssNonzero[i][j]) tmpVal += pairConstants_->ssPair[i][j]*this->oneevRRSa0(ijShellPair,LA+2,lAp1,&i,&j);
      } else if(abs(ijShellPair->deltaAB[iWork2])>this->controls_->thresholdS&&abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) {
	double sqrDeltaAB = (ijShellPair->deltaAB[iWork])*(ijShellPair->deltaAB[iWork2]);
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) 
          if(pairConstants_->ssNonzero[i][j]) {
	    primitiveVal = sqrDeltaAB*this->oneevRRSa0(ijShellPair,LA,lA,&i,&j);
	    lAp1[iWork]++;
	    primitiveVal += ijShellPair->deltaAB[iWork2]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,&i,&j);
	    lAp1[iWork2]++;
	    primitiveVal += this->oneevRRSa0(ijShellPair,LA+2,lAp1,&i,&j);
	    lAp1[iWork]--;
	    primitiveVal += ijShellPair->deltaAB[iWork]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,&i,&j);
	    lAp1[iWork2]--;
	    tmpVal += primitiveVal*pairConstants_->ssPair[i][j];
	  };
      } else {
	if(abs(ijShellPair->deltaAB[iWork2])<this->controls_->thresholdS) {
	  int tmpWork = iWork;
	  iWork = iWork2;
	  iWork2 = tmpWork;
	};
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) 
          if(pairConstants_->ssNonzero[i][j]) {
	    lAp1[iWork]++;
	    primitiveVal = ijShellPair->deltaAB[iWork2]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,&i,&j);
	    lAp1[iWork2]++;
	    primitiveVal += this->oneevRRSa0(ijShellPair,LA+2,lAp1,&i,&j);
	    tmpVal += primitiveVal*pairConstants_->ssPair[i][j];
	    lAp1[iWork2]--;
	    lAp1[iWork]--;
	  };
      };
    };
  } else {
    // (|p)
    lAp1[iWork]++;
    if(abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) {
      for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) 
        if(pairConstants_->ssNonzero[i][j]) tmpVal += pairConstants_->ssPair[i][j]*(this->oneevRRSa0(ijShellPair,LA+1,lAp1,&i,&j)+ ijShellPair->deltaAB[iWork]*this->oneevRRSa0(ijShellPair,LA,lA,&i,&j));
    } else {
      for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++)
        if(pairConstants_->ssNonzero[i][j]) tmpVal += pairConstants_->ssPair[i][j]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,&i,&j);
    };
  };
  return tmpVal;
};

// overlap vertical recursion
// (a|0) = (P-A)(a-1|0) + 1/(2*(alpha+beta))*N_(a-1)*(a-2|0)
double AOIntegrals::oneevRRSa0(ShellPair *ijShellPair,int LA,int *lA,int *i,int *j) {
  double tmpVal=0.0;
  int iWork;
  int lAm1[3];
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0) iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;
  if(LA>4) {
    //[>g|s]
    lAm1[iWork]--;
    if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*this->oneevRRSa0(ijShellPair,LA-1,lAm1,i,j);
    if(LA==2&&lA[iWork]==2) tmpVal += pairConstants_->Sa0Par[*i][*j];
    else if (lA[iWork]>=2) {
      lAm1[iWork]--;
      tmpVal += (lAm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRSa0(ijShellPair,LA-2,lAm1,i,j);
    };
  } else if (LA==4) {
    // [g|s]
    if(lA[0]==4||lA[1]==4||lA[2]==4) {
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) {
	double sqrDeltaPA = ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork][*i][*j];
	tmpVal = sqrDeltaPA*sqrDeltaPA + pairConstants_->Sa0Par[*i][*j]*(math.six*sqrDeltaPA+math.three*pairConstants_->Sa0Par[*i][*j]);
      } else tmpVal = math.three*pairConstants_->Sa0Par[*i][*j]*pairConstants_->Sa0Par[*i][*j];
    } else if(lA[0]==3||lA[1]==3||lA[2]==3) {
      int iWork3;
      if(lA[0]==3) iWork3 = 0;
      else if(lA[1]==3) iWork3 = 1;
      else if(lA[2]==3) iWork3 = 2;
      if(lA[0]==1) iWork = 0;
      else if(lA[1]==1) iWork = 1;
      else if(lA[2]==1) iWork = 2;
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) {
	if(abs(ijShellPair->deltaPA[iWork3][*i][*j])>this->controls_->thresholdS) {
	  tmpVal = math.three*pairConstants_->Sa0Par[*i][*j] + ijShellPair->deltaPA[iWork3][*i][*j]*ijShellPair->deltaPA[iWork3][*i][*j];
	  tmpVal *= ijShellPair->deltaPA[iWork3][*i][*j];
	};
      };
    } else if((lA[0]==2||lA[1]==2||lA[2]==2)&&(lA[0]==0||lA[1]==0||lA[2]==0)) {
      int iWork2;
      if(lA[0]==2&&iWork!=0) iWork2 = 0;
      else if(lA[1]==2&&iWork!=1) iWork2 = 1;
      else if(lA[2]==2&&iWork!=2) iWork2 = 2;
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS&&abs(ijShellPair->deltaPA[iWork2][*i][*j])>this->controls_->thresholdS) {
	double sqrDeltaPA1 = ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork][*i][*j];
	double sqrDeltaPA2 = ijShellPair->deltaPA[iWork2][*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j];
	tmpVal = sqrDeltaPA1*sqrDeltaPA2 + pairConstants_->Sa0Par[*i][*j]*(sqrDeltaPA1+sqrDeltaPA2+pairConstants_->Sa0Par[*i][*j]);
      } else if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS){
	double sqrDeltaPA1 = ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork][*i][*j];
	tmpVal = pairConstants_->Sa0Par[*i][*j]*(sqrDeltaPA1+pairConstants_->Sa0Par[*i][*j]);	
      } else if(abs(ijShellPair->deltaPA[iWork2][*i][*j])>this->controls_->thresholdS){
	double sqrDeltaPA2 = ijShellPair->deltaPA[iWork2][*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j];
	tmpVal = pairConstants_->Sa0Par[*i][*j]*(sqrDeltaPA2+pairConstants_->Sa0Par[*i][*j]);	
      } else tmpVal = pairConstants_->Sa0Par[*i][*j]*pairConstants_->Sa0Par[*i][*j];
    } else {
      int iWork1,iWork2;
      if(lA[0]==1&&iWork!=0) iWork1 = 0;
      else if(lA[1]==1&&iWork!=1) iWork1 = 1;
      else if(lA[2]==1&&iWork!=2) iWork1 = 2;
      if(lA[0]==2) iWork2 = 0;
      else if(lA[1]==2) iWork2 = 1;
      else if(lA[2]==2) iWork2 = 2;
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS&&abs(ijShellPair->deltaPA[iWork1][*i][*j])>this->controls_->thresholdS) {
	double sqrDeltaPA2 = ijShellPair->deltaPA[iWork2][*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j];
	tmpVal = ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork1][*i][*j]*(pairConstants_->Sa0Par[*i][*j]+sqrDeltaPA2);
      };
    };
  } else if (LA==3) {
    if(lA[0]==3||lA[1]==3||lA[2]==3) {
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) {
	tmpVal = math.three*pairConstants_->Sa0Par[*i][*j] + ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork][*i][*j];
	tmpVal *= ijShellPair->deltaPA[iWork][*i][*j];
      };
    } else if(lA[0]==lA[1]&&lA[0]==lA[2]) {
      if(abs(ijShellPair->deltaPA[0][*i][*j])>this->controls_->thresholdS&&abs(ijShellPair->deltaPA[1][*i][*j])>this->controls_->thresholdS&&abs(ijShellPair->deltaPA[2][*i][*j]))
	tmpVal = ijShellPair->deltaPA[0][*i][*j]*ijShellPair->deltaPA[1][*i][*j]*ijShellPair->deltaPA[2][*i][*j];
    } else {
      int iWork2;
      if(lA[0]==2) iWork2 = 0;
      else if(lA[1]==2) iWork2 = 1;
      else if(lA[2]==2) iWork2 = 2;
      if(lA[0]==1) iWork = 0;
      else if(lA[1]==1) iWork = 1;
      else if(lA[2]==1) iWork = 2;
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) {
	tmpVal = pairConstants_->Sa0Par[*i][*j];
	if(abs(ijShellPair->deltaPA[iWork2][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPA[iWork2][*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j];
	tmpVal *= ijShellPair->deltaPA[iWork][*i][*j];
      };
    };
  } else if (LA==2) {
    lAm1[iWork]--;
    int iWork2;
    if (lAm1[0]>0)      iWork2=0;
    else if (lAm1[1]>0) iWork2=1;
    else if (lAm1[2]>0) iWork2=2;
    if(iWork!=iWork2) {
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS&&abs(ijShellPair->deltaPA[iWork2][*i][*j])>this->controls_->thresholdS) 
	tmpVal = ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j];
    } else {
      tmpVal = pairConstants_->Sa0Par[*i][*j];
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork][*i][*j];
    };
  } else if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) tmpVal = ijShellPair->deltaPA[iWork][*i][*j];
  return tmpVal;
};
// overlap recursion for kinetic
// (a|b) = (A-B)(a|b-1) + (a+1|b-1)
//   LA >= LB
double AOIntegrals::oneehRRTSab(ShellPair *ijShellPair,int LA,int *lA,int LB,int *lB,int *i, int *j) {
  double tmpVal = 0.0;
  if(LA+LB==0) return 1.0;
  else if(LB==0) return this->oneevRRSa0(ijShellPair,LA,lA,i,j);
  int iWork,lAp1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAp1[iWork]=lA[iWork];     
    lBm1[iWork]=lB[iWork];     
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;

  if(LB>2) {
    lAp1[iWork]++;
    lBm1[iWork]--;
    tmpVal = this->oneehRRTSab(ijShellPair,LA+1,lAp1,LB-1,lBm1,i,j);
    if(abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) tmpVal+= ijShellPair->deltaAB[iWork]*this->oneehRRTSab(ijShellPair,LA,lA,LB-1,lBm1,i,j);
  } else if(LB==2) {
    lBm1[iWork]--;
    int iWork2;
    if (lBm1[0]>0)      iWork2=0;
    else if (lBm1[1]>0) iWork2=1;
    else if (lBm1[2]>0) iWork2=2;

    if(iWork==iWork2) {
      if(abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) {
	double sqrDeltaAB = (ijShellPair->deltaAB[iWork])*(ijShellPair->deltaAB[iWork]);
        tmpVal = sqrDeltaAB*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
        lAp1[iWork]++;
        tmpVal += math.two*ijShellPair->deltaAB[iWork]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,i,j);
        lAp1[iWork]++;
        tmpVal += this->oneevRRSa0(ijShellPair,LA+2,lAp1,i,j);
        lAp1[iWork]-=2;
      } else {
	lAp1[iWork]+=2;
	tmpVal += this->oneevRRSa0(ijShellPair,LA+2,lAp1,i,j);
      };
    } else {
      if(abs(ijShellPair->deltaAB[iWork2])<this->controls_->thresholdS&&abs(ijShellPair->deltaAB[iWork])<this->controls_->thresholdS) {
	lAp1[iWork]++;
	lAp1[iWork2]++;
	tmpVal += this->oneevRRSa0(ijShellPair,LA+2,lAp1,i,j);
      } else if(abs(ijShellPair->deltaAB[iWork2])>this->controls_->thresholdS&&abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) {
	double sqrDeltaAB = (ijShellPair->deltaAB[iWork])*(ijShellPair->deltaAB[iWork2]);
	tmpVal = sqrDeltaAB*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
	lAp1[iWork]++;
	tmpVal += ijShellPair->deltaAB[iWork2]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,i,j);
	lAp1[iWork2]++;
	tmpVal += this->oneevRRSa0(ijShellPair,LA+2,lAp1,i,j);
	lAp1[iWork]--;
	tmpVal += ijShellPair->deltaAB[iWork]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,i,j);
      } else {
	if(abs(ijShellPair->deltaAB[iWork2])<this->controls_->thresholdS) {
	  int tmpWork = iWork;
	  iWork = iWork2;
	  iWork2 = tmpWork;
	};
	lAp1[iWork]++;
	tmpVal = ijShellPair->deltaAB[iWork2]*this->oneevRRSa0(ijShellPair,LA+1,lAp1,i,j);
	lAp1[iWork2]++;
	tmpVal += this->oneevRRSa0(ijShellPair,LA+2,lAp1,i,j);
      };
    };
  } else {
    lAp1[iWork]++;
    tmpVal = this->oneevRRSa0(ijShellPair,LA+1,lAp1,i,j);
    if(abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaAB[iWork]*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
  };
  return tmpVal;
};
// kinetic vertical recursion
//  from (a|T|b) to (a|T|0)
double AOIntegrals::oneevRRTab(ShellPair *ijShellPair,int LA,int *lA,int LB,int *lB,int *i,int *j) {
  double tmpVal = 0.0;
  if(LA+LB==0) return pairConstants_->TabPar1[*i][*j];
  else if(LB==0) return this->oneevRRTa0(ijShellPair,LA,lA,i,j);
  int iWork,lAm1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAm1[iWork]=lA[iWork];
    lBm1[iWork]=lB[iWork];
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lBm1[iWork]--;
  tmpVal = pairConstants_->TabPar2[*i][*j]*this->oneehRRTSab(ijShellPair,LA,lA,LB,lB,i,j);
  if(LB>2) {
    if(abs(ijShellPair->deltaPB[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPB[iWork][*i][*j]*this->oneevRRTab(ijShellPair,LA,lA,LB-1,lBm1,i,j);
    if (lA[iWork]>0) {
      lAm1[iWork]--;
      tmpVal += (lAm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTab(ijShellPair,LA-1,lAm1,LB-1,lBm1,i,j);
    };
    if(LB==2&&lB[iWork]==2) tmpVal += pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA,lA,i,j) - pairConstants_->TabPar3[*i][*j]*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
    else if (lB[iWork]>=2) {
      lBm1[iWork]--;
      tmpVal += (lBm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTab(ijShellPair,LA,lA,LB-2,lBm1,i,j);
      tmpVal -= (lBm1[iWork]+1)*pairConstants_->TabPar3[*i][*j]*this->oneehRRTSab(ijShellPair,LA,lA,LB-2,lBm1,i,j);
    };
  } else if (LB==2) {
    int iWork2;
    if (lBm1[0]>0)      iWork2=0;
    else if (lBm1[1]>0) iWork2=1;
    else if (lBm1[2]>0) iWork2=2;
    double tmpVal2 = 0.0;
    if(abs(ijShellPair->deltaPB[iWork][*i][*j])>this->controls_->thresholdS) {
      tmpVal2 = pairConstants_->TabPar2[*i][*j]*this->oneehRRTSab(ijShellPair,LA,lA,LB-1,lBm1,i,j);
      if(abs(ijShellPair->deltaPB[iWork2][*i][*j])>this->controls_->thresholdS) tmpVal2 += ijShellPair->deltaPB[iWork2][*i][*j]*this->oneevRRTa0(ijShellPair,LA,lA,i,j);
      if(lA[iWork2]>0) {
	lAm1[iWork2]--;
	tmpVal2 += (lAm1[iWork2]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA-1,lAm1,i,j);
	lAm1[iWork2]++;
      };
      tmpVal += tmpVal2*ijShellPair->deltaPB[iWork][*i][*j];
    };
    if(lA[iWork]>0) {
      lAm1[iWork]--;
      tmpVal2 = pairConstants_->TabPar2[*i][*j]*this->oneehRRTSab(ijShellPair,LA-1,lAm1,LB-1,lBm1,i,j);
      if(abs(ijShellPair->deltaPB[iWork2][*i][*j])>this->controls_->thresholdS) tmpVal2 += ijShellPair->deltaPB[iWork2][*i][*j]*this->oneevRRTa0(ijShellPair,LA-1,lAm1,i,j);
      if(lAm1[iWork2]>0) {
	lAm1[iWork2]--;
	if(LA>2) tmpVal2 += (lAm1[iWork2]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA-2,lAm1,i,j);
	else tmpVal2 += pairConstants_->Sa0Par[*i][*j]*pairConstants_->TabPar1[*i][*j];
      };
      tmpVal += tmpVal2*lA[iWork]*pairConstants_->Sa0Par[*i][*j];
    };
    if(iWork==iWork2) tmpVal += pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA,lA,i,j) - pairConstants_->TabPar3[*i][*j]*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
  } else {
    if(abs(ijShellPair->deltaPB[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPB[iWork][*i][*j]*this->oneevRRTa0(ijShellPair,LA,lA,i,j);
    if (lA[iWork]>0) {
      lAm1[iWork]--;
      tmpVal += (lAm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA-1,lAm1,i,j);
    };
  };
  return tmpVal;
};
// kinetic vertical recursion
//  from (a|T|0) to (0|T|0)
double AOIntegrals::oneevRRTa0(ShellPair *ijShellPair,int LA,int *lA,int *i,int *j) {
  double tmpVal = 0.0;
  if(LA==0) return pairConstants_->TabPar1[*i][*j];
  int iWork,lAm1[3];
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0)      iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;
  lAm1[iWork]--;
  tmpVal = pairConstants_->TabPar2[*i][*j]*this->oneevRRSa0(ijShellPair,LA,lA,i,j);
  if(LA>2) {
    if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*this->oneevRRTa0(ijShellPair,LA-1,lAm1,i,j);
    if(LA==2&&lA[iWork]==2) {
      tmpVal += pairConstants_->Sa0Par[*i][*j]*pairConstants_->TabPar1[*i][*j];
      tmpVal -= pairConstants_->Ta0Par3[*i][*j];
    } else if (lA[iWork]>=2) {
      lAm1[iWork] -=1;
      tmpVal += (lAm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*this->oneevRRTa0(ijShellPair,LA-2,lAm1,i,j);
      tmpVal -= (lAm1[iWork]+1)*pairConstants_->Ta0Par3[*i][*j]*this->oneevRRSa0(ijShellPair,LA-2,lAm1,i,j);
    };
  } else if (LA==2) {
    int iWork2;
    if (lAm1[0]>0)      iWork2=0;
    else if (lAm1[1]>0) iWork2=1;
    else if (lAm1[2]>0) iWork2=2;
    double tmpVal2 = 0.0;
    if(iWork!=iWork2) {
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) {
	if(abs(ijShellPair->deltaPA[iWork2][*i][*j])>this->controls_->thresholdS)
	  tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*(pairConstants_->TabPar2[*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j] + ijShellPair->deltaPA[iWork2][*i][*j]*pairConstants_->TabPar1[*i][*j]);
      };
    } else {
      tmpVal += pairConstants_->Sa0Par[*i][*j]*pairConstants_->TabPar1[*i][*j] - pairConstants_->Ta0Par3[*i][*j];
      if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) 
	tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*(ijShellPair->deltaPA[iWork][*i][*j]*pairConstants_->TabPar1[*i][*j] + pairConstants_->TabPar2[*i][*j]*ijShellPair->deltaPA[iWork][*i][*j]);
    };
  } else if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) tmpVal += ijShellPair->deltaPA[iWork][*i][*j]*pairConstants_->TabPar1[*i][*j];
  return tmpVal;
};
// potential integral horizontal recursion
double AOIntegrals::oneehRRVab(ShellPair *ijShellPair,int LA,int *lA,int LB,int *lB){
  int i,j,iAtom;
  double tmp,tmpVal=0.0;

  if(LA<=2){
    if(LB==0){
      if(LA==0){
	// (s||s)
	tmp = 0.0;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++)
	  if(pairConstants_->ssNonzero[i][j]) for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++) 
	    tmpVal += molecularConstants_->atomZ[iAtom]*pairConstants_->FmU[0][i][j][iAtom];
	return tmpVal;
      }else if(LA==1){
	// (p||s)
	int iWork;
	double B0;
	if (lA[0]>0) iWork=0;
	else if (lA[1]>0) iWork=1;
	else if (lA[2]>0) iWork=2;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	  B0 = ijShellPair->deltaPA[iWork][i][j];
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++)
	    tmpVal += molecularConstants_->atomZ[iAtom]*(pairConstants_->FmU[0][i][j][iAtom]*B0 - pairConstants_->FmU[1][i][j][iAtom]*pairConstants_->deltaPZ[iWork][i][j][iAtom]);
	};
	return tmpVal;
      }else if(LA==2){
	/**********/
	/* (d||s) */
	/**********/
	int iWork1,iWork2;
	if (lA[0]>0)      iWork1=0;
	else if (lA[1]>0) iWork1=1;
	else if (lA[2]>0) iWork1=2;
	lA[iWork1]--;
	if (lA[0]>0)      iWork2=0;
	else if (lA[1]>0) iWork2=1;
	else if (lA[2]>0) iWork2=2;
	lA[iWork1]++;
	if(iWork1!=iWork2) {
	  double B0,B1,B0B1;
	  for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	    B0 = ijShellPair->deltaPA[iWork1][i][j];
	    B1 = ijShellPair->deltaPA[iWork2][i][j];
	    B0B1 = B0*B1;
	    for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	      tmpVal += molecularConstants_->atomZ[iAtom]*
		(pairConstants_->FmU[0][i][j][iAtom]*B0B1
		 -pairConstants_->FmU[1][i][j][iAtom]*(pairConstants_->deltaPZ[iWork2][i][j][iAtom]*B0+B1*pairConstants_->deltaPZ[iWork1][i][j][iAtom])
		 +pairConstants_->FmU[2][i][j][iAtom]*pairConstants_->deltaPZ[iWork1][i][j][iAtom]*pairConstants_->deltaPZ[iWork2][i][j][iAtom]);
	    };
	  };
	} else {
	  double S,B0,B02,Coeff0;
	  for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	    B0 = ijShellPair->deltaPA[iWork1][i][j];
	    B02 = math.two*B0;
	    S = pairConstants_->Sa0Par[i][j];
	    Coeff0 = B0*B0+S;
	    for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	      tmpVal += molecularConstants_->atomZ[iAtom]*
		(pairConstants_->FmU[0][i][j][iAtom]*Coeff0
		 -pairConstants_->FmU[1][i][j][iAtom]*(B02*pairConstants_->deltaPZ[iWork1][i][j][iAtom]+S)
		 +pairConstants_->FmU[2][i][j][iAtom]*pairConstants_->deltaPZ[iWork1][i][j][iAtom]*pairConstants_->deltaPZ[iWork1][i][j][iAtom]);
	    };
	  };
	};
	return tmpVal;
      };
    }else if(LB==1){
      if(LA==1){
	/**********/
	/* (p||p) */
	/**********/
	int iWork1,iWork2;
	if (lB[0]>0)      iWork1=0;
	else if (lB[1]>0) iWork1=1;
	else if (lB[2]>0) iWork1=2;
	if (lA[0]>0)      iWork2=0;
	else if (lA[1]>0) iWork2=1;
	else if (lA[2]>0) iWork2=2;
	if(iWork1!=iWork2) {
	  /* (x||y) */
	  double A1pB1,B0,Coeff0;
	  for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	    A1pB1 = ijShellPair->deltaAB[iWork1]+ijShellPair->deltaPA[iWork1][i][j];
	    B0 = ijShellPair->deltaPA[iWork2][i][j];
	    Coeff0 = A1pB1*B0;
	    for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	      tmpVal += molecularConstants_->atomZ[iAtom]*
		(pairConstants_->FmU[0][i][j][iAtom]*Coeff0
		 -pairConstants_->FmU[1][i][j][iAtom]*(A1pB1*pairConstants_->deltaPZ[iWork2][i][j][iAtom]+B0*pairConstants_->deltaPZ[iWork1][i][j][iAtom])
		 +pairConstants_->FmU[2][i][j][iAtom]*pairConstants_->deltaPZ[iWork2][i][j][iAtom]*pairConstants_->deltaPZ[iWork1][i][j][iAtom]);
	    };
	    //(A1+B1)*B0
	    //(A1+B1)*C0 + B0*C1
	    //C0*C1
	  };
	} else {
	  /* (x||x) */
	  double A0pB0,S,A0p2B0,Coeff0;
	  for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	    A0pB0 = ijShellPair->deltaAB[iWork1]+ijShellPair->deltaPA[iWork1][i][j];
	    S = pairConstants_->Sa0Par[i][j];
	    Coeff0 = A0pB0*ijShellPair->deltaPA[iWork1][i][j]+S;
	    A0p2B0 = A0pB0+ijShellPair->deltaPA[iWork1][i][j];
	    for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	      tmpVal += molecularConstants_->atomZ[iAtom]*
		(pairConstants_->FmU[0][i][j][iAtom]*Coeff0
		 -pairConstants_->FmU[1][i][j][iAtom]*(A0p2B0*pairConstants_->deltaPZ[iWork1][i][j][iAtom]+S)
		 +pairConstants_->FmU[2][i][j][iAtom]*pairConstants_->deltaPZ[iWork1][i][j][iAtom]*pairConstants_->deltaPZ[iWork1][i][j][iAtom]);
	    };
	  };
	};
	return tmpVal;
      }else if(LA==2){
	/**********/
	/* (d||p) */
	/**********/
	int iWorkB,iWorkA1,iWorkA2;
	if (lB[0]>0) iWorkB=0;
	else if (lB[1]>0) iWorkB=1;
	else if (lB[2]>0) iWorkB=2;
	if (lA[0]>0) iWorkA1=0;
	else if (lA[1]>0) iWorkA1=1;
	else if (lA[2]>0) iWorkA1=2;
	lA[iWorkA1]--;
	if (lA[0]>0) iWorkA2=0;
	else if (lA[1]>0) iWorkA2=1;
	else if (lA[2]>0) iWorkA2=2;
	lA[iWorkA1]++;
	if(iWorkA1==iWorkA2&&iWorkB!=iWorkA1) {
	  /* (xx||y) */
	  double A1pB1,B0B0pS,B0C02pS,C0C0,Coeff0,S,B02;
	  for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	    S = pairConstants_->Sa0Par[i][j];
	    A1pB1  = ijShellPair->deltaAB[iWorkB]+ijShellPair->deltaPA[iWorkB][i][j];
	    B0B0pS = ijShellPair->deltaPA[iWorkA1][i][j]*ijShellPair->deltaPA[iWorkA1][i][j]+S;
	    B02 = math.two*ijShellPair->deltaPA[iWorkA1][i][j];
	    Coeff0 = A1pB1*B0B0pS;
	    for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	      B0C02pS= B02*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]+S;
	      C0C0   = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	      tmpVal += molecularConstants_->atomZ[iAtom]*
		(pairConstants_->FmU[0][i][j][iAtom]*Coeff0
		 -pairConstants_->FmU[1][i][j][iAtom]*(A1pB1*B0C02pS+B0B0pS*pairConstants_->deltaPZ[iWorkB][i][j][iAtom])
		 +pairConstants_->FmU[2][i][j][iAtom]*(A1pB1*C0C0+B0C02pS*pairConstants_->deltaPZ[iWorkB][i][j][iAtom])
		 -pairConstants_->FmU[3][i][j][iAtom]*(C0C0*pairConstants_->deltaPZ[iWorkB][i][j][iAtom]));
	      //(A1+B1)*(B0*B0+S)
	      //(A1+B1)*(2*B0*C0+S) + (B0*B0+S)*C1
	      //(A1+B1)*C0*C0 + (2*B0*C0+S)*C1
	      //C0*C0*C1
	    };
	  };
	  return tmpVal;
	} else if(iWorkA1==iWorkA2&&iWorkB==iWorkA1) {
	  /* (xx||x) */
	  double B0B0,C0C0,A0p3B0,A02p3B0mB0,A0p3B0mS,C0S3,Coeff0,S,S3;
	  for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	    B0B0    = ijShellPair->deltaPA[iWorkB][i][j]*ijShellPair->deltaPA[iWorkB][i][j];
	    A0p3B0  = ijShellPair->deltaAB[iWorkB]+double(3.0)*ijShellPair->deltaPA[iWorkB][i][j];
	    A02p3B0mB0  = (ijShellPair->deltaAB[iWorkB]+A0p3B0)*ijShellPair->deltaPA[iWorkB][i][j];
	    S = pairConstants_->Sa0Par[i][j];
	    S3 = double(3.0)*S;
	    A0p3B0mS= A0p3B0*S;
	    Coeff0  = (ijShellPair->deltaAB[iWorkB]+ijShellPair->deltaPA[iWorkB][i][j])*B0B0+A0p3B0mS;
	    for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	      C0C0    = pairConstants_->deltaPZ[iWorkB][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB][i][j][iAtom];
	      C0S3    = S3*pairConstants_->deltaPZ[iWorkB][i][j][iAtom];
	      tmpVal += molecularConstants_->atomZ[iAtom]*	    
		(pairConstants_->FmU[0][i][j][iAtom]*Coeff0
		 -pairConstants_->FmU[1][i][j][iAtom]*(A02p3B0mB0*pairConstants_->deltaPZ[iWorkB][i][j][iAtom]+A0p3B0mS+C0S3)
		 +pairConstants_->FmU[2][i][j][iAtom]*(A0p3B0*C0C0+C0S3)
		 -pairConstants_->FmU[3][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB][i][j][iAtom]*C0C0);
	      //(A0+B0)*B0*B0 + (A0+3*B0)*S
	      //(2*A0+3*B0)*B0*C0 + (A0+3*B0)*S + 3*C0*S
	      //(A0+3*B0)*C0*C0 + 3*C0*S
	      //C0*C0*C0
	    };
	  };
	  return tmpVal;
	} else if(iWorkA1!=iWorkA2&&(iWorkA1==iWorkB||iWorkA2==iWorkB)){
	  /* (xy||x) */
	  double A0pB0mB0,A0p2B0,A0p2B0mB1,A0p2B0mB1C0,A0p2B0mC0,B1S,C1S,C0C0,S,Coeff0,B1;
	  if(iWorkA1==iWorkB) iWorkA1 = iWorkA2;
	  for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	    S = pairConstants_->Sa0Par[i][j];
	    B1= ijShellPair->deltaPA[iWorkA1][i][j];
	    A0pB0mB0 = (ijShellPair->deltaAB[iWorkB]+ijShellPair->deltaPA[iWorkB][i][j])*ijShellPair->deltaPA[iWorkB][i][j];
	    A0p2B0   = ijShellPair->deltaAB[iWorkB]+math.two*ijShellPair->deltaPA[iWorkB][i][j];
	    A0p2B0mB1= A0p2B0*B1;
	    B1S      = ijShellPair->deltaPA[iWorkA1][i][j]*S;
	    Coeff0   = A0pB0mB0*ijShellPair->deltaPA[iWorkA1][i][j]+B1S;
	    for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	      A0p2B0mB1C0= A0p2B0mB1*pairConstants_->deltaPZ[iWorkB][i][j][iAtom];
	      A0p2B0mC0= A0p2B0*pairConstants_->deltaPZ[iWorkB][i][j][iAtom];
	      C1S      = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*S;
	      C0C0     = pairConstants_->deltaPZ[iWorkB][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB][i][j][iAtom];
	      tmpVal += molecularConstants_->atomZ[iAtom]*
		(pairConstants_->FmU[0][i][j][iAtom]*Coeff0
		 -pairConstants_->FmU[1][i][j][iAtom]*(A0p2B0mB1C0+A0pB0mB0*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]+B1S+C1S)
		 +pairConstants_->FmU[2][i][j][iAtom]*(B1*C0C0+A0p2B0mC0*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]+C1S)
		 -pairConstants_->FmU[3][i][j][iAtom]*C0C0*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]);
	      //(A0+B0)*B0*B1 + B1*S
	      //(A0+2*B0)*C0*B1 + (A0+B0)*B0*C1 + B1*S + C1*S
	      //B1*C0*C0 + (A0+2*B0)*C0*C1 + C1*S
	      //C0*C0*C1
	    };
	  };
	  return tmpVal;
	} else {
	  /* (xy||z) */
	  double B0,B1,B0B1,A2pB2,B1C0pB0C1,C0C1,Coeff0;
	  for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	    B0 = ijShellPair->deltaPA[iWorkA1][i][j];
	    B1 = ijShellPair->deltaPA[iWorkA2][i][j];
	    B0B1      = B0*B1;
	    A2pB2     = ijShellPair->deltaAB[iWorkB]+ijShellPair->deltaPA[iWorkB][i][j];
	    Coeff0    = A2pB2*B0B1;
	    for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	      B1C0pB0C1 = B1*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]+B0*pairConstants_->deltaPZ[iWorkA2][i][j][iAtom];
	      C0C1      = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA2][i][j][iAtom];
	      tmpVal += molecularConstants_->atomZ[iAtom]*
		(pairConstants_->FmU[0][i][j][iAtom]*Coeff0
		 -pairConstants_->FmU[1][i][j][iAtom]*(A2pB2*B1C0pB0C1+pairConstants_->deltaPZ[iWorkB][i][j][iAtom]*B0B1)
		 +pairConstants_->FmU[2][i][j][iAtom]*(A2pB2*C0C1+B1C0pB0C1*pairConstants_->deltaPZ[iWorkB][i][j][iAtom])
		 -pairConstants_->FmU[3][i][j][iAtom]*C0C1*pairConstants_->deltaPZ[iWorkB][i][j][iAtom]);
	      //(A2+B2)*B0*B1
	      //(A2+B2)*(B1*C0+B0*C1) + B0*B1*C2
	      //(A2+B2)*C0*C1 + (B1*C0+B0*C1)*C2
	      //C0*C1*C2
	    };
	  };
	  return tmpVal;
	};
      };
    }else if(LB==2){
      /**********/
      /* (d||d) */
      /**********/
      int iWorkA1,iWorkA2,iWorkB1,iWorkB2;
      if (lB[0]>0)      iWorkB1=0;
      else if (lB[1]>0) iWorkB1=1;
      else if (lB[2]>0) iWorkB1=2;
      lB[iWorkB1]--;
      if (lB[0]>0)      iWorkB2=0;
      else if (lB[1]>0) iWorkB2=1;
      else if (lB[2]>0) iWorkB2=2;
      lB[iWorkB1]++;

      if (lA[0]>0)      iWorkA1=0;
      else if (lA[1]>0) iWorkA1=1;
      else if (lA[2]>0) iWorkA1=2;
      lA[iWorkA1]--;
      if (lA[0]>0)      iWorkA2=0;
      else if (lA[1]>0) iWorkA2=1;
      else if (lA[2]>0) iWorkA2=2;
      lA[iWorkA1]++;
      if(iWorkB1==iWorkB2&&iWorkA1==iWorkA2&&iWorkA1!=iWorkB1){
	/* (xx||yy) */
	double A1pB1,B0B0pS,sqrA1pB1pS,A1pB1m2C1,A1pB1m2C1pS,B0C02,B0C02pS,C0C0,C1C1,Coeff0,S,B0,B02,A1pB1m2,sum1;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  S = pairConstants_->Sa0Par[i][j];
	  B0 = ijShellPair->deltaPA[iWorkA1][i][j];
	  B02 = math.two*B0;
	  A1pB1      = ijShellPair->deltaAB[iWorkB1]+ijShellPair->deltaPA[iWorkB1][i][j];
	  A1pB1m2    = A1pB1*math.two;
	  B0B0pS     = B0*B0 + S;
	  sqrA1pB1pS = A1pB1*A1pB1 + S;
	  Coeff0 = B0B0pS*sqrA1pB1pS;
	  sum1   = (sqrA1pB1pS+B0B0pS)*S;
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++) {
	    A1pB1m2C1  = A1pB1m2*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom];
	    A1pB1m2C1pS= A1pB1m2C1 + S;
	    B0C02      = B02*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    B0C02pS    = B0C02 + S;
	    C0C0       = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    C1C1       = pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom];
	    tmpVal += molecularConstants_->atomZ[iAtom]*
	      (pairConstants_->FmU[0][i][j][iAtom]*Coeff0
	       -pairConstants_->FmU[1][i][j][iAtom]*(sqrA1pB1pS*B0C02+B0B0pS*A1pB1m2C1+sum1)
	       +pairConstants_->FmU[2][i][j][iAtom]*(sqrA1pB1pS*C0C0+B0C02pS*A1pB1m2C1pS+B0B0pS*C1C1)
	       -pairConstants_->FmU[3][i][j][iAtom]*(A1pB1m2C1pS*C0C0+B0C02pS*C1C1)
	       +pairConstants_->FmU[4][i][j][iAtom]*C0C0*C1C1);
	    //(B0^2+S)*((A1+B1)^2+S)
	    //2*((A1+B1)^2+S)*B0*C0 + 2*(B0^2+S)*(A1+B1)*C1 + ((A1+B1)^2+S+B0^2+S)*S
	    //((A1+B1)^2+S)*C0^2 + (2*B0*C0+S)*(2*(A1+B1)*C1+S) + (B0^2+S)*C1^2
	    //(2*(A1+B1)*C1+S)*C0^2 + (2*B0*C0+S)*C1^2
	    //C0^2*C1^2
	  };
	};
	return tmpVal;
      } else if(iWorkB1==iWorkB2&&iWorkA1==iWorkA2&&iWorkA1==iWorkB1){
	/* (xx||xx) */
	double A0A0,B0B0,A0pB0mB0,A0p2B0,A0p2B0m2C0,C0C0,sum1,sum2,S,Coeff0;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  S = pairConstants_->Sa0Par[i][j];
	  A0A0      = ijShellPair->deltaAB[iWorkA1]*ijShellPair->deltaAB[iWorkA1];
	  B0B0      = ijShellPair->deltaPA[iWorkA1][i][j]*ijShellPair->deltaPA[iWorkA1][i][j];
	  A0pB0mB0  = (ijShellPair->deltaAB[iWorkA1]+ijShellPair->deltaPA[iWorkA1][i][j])*ijShellPair->deltaPA[iWorkA1][i][j];
	  A0p2B0    = ijShellPair->deltaAB[iWorkA1]+math.two*ijShellPair->deltaPA[iWorkA1][i][j];
	  sum1      = A0p2B0*A0p2B0+math.two*A0pB0mB0+double(3.0)*pairConstants_->Sa0Par[i][j];
	  Coeff0    = A0pB0mB0*A0pB0mB0+sum1*S;
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	    A0p2B0m2C0= math.two*(ijShellPair->deltaAB[iWorkA1]+math.two*ijShellPair->deltaPA[iWorkA1][i][j])*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    C0C0      = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    sum2      = double(3.0)*(A0p2B0m2C0+pairConstants_->Sa0Par[i][j]);
	    tmpVal += molecularConstants_->atomZ[iAtom]*
	      (pairConstants_->FmU[0][i][j][iAtom]*Coeff0
	       -pairConstants_->FmU[1][i][j][iAtom]*(A0pB0mB0*A0p2B0m2C0+(sum1+sum2)*S)
	       +pairConstants_->FmU[2][i][j][iAtom]*(sum1*C0C0+(sum2+double(3.0)*C0C0)*S)
	       -pairConstants_->FmU[3][i][j][iAtom]*(A0p2B0m2C0+double(6.0)*S)*C0C0
	       +pairConstants_->FmU[4][i][j][iAtom]*(C0C0*C0C0));
	    //(A0+B0)^2*B0^2 + ((A0+2*B0)^2+2*(A0+B0)*B0+3*S)*S
	    //2*(A0+B0)*(A0+2*B0)*B0*C0 + (A0^2+6*A0*B0+6*B0^2+3*S)*S + (6*(A0+2*B0)*C0+3*S)*S
	    //(A0^2+6*A0*B0+6*B0^2+3*S)*C0^2 + (6*(A0+2*B0)*C0+3*S+3*C0^2)*S
	    //(2*(A0+2*B0)*C0+6*S)*C0^2
	    //C0^4
	  };
	};
	return tmpVal;
      } else if(iWorkA1==iWorkA2&&iWorkB1!=iWorkB2&&(iWorkA1==iWorkB1||iWorkA1==iWorkB2)){
	/* (xx||xy) */
	double A0pB0,A0p3B0,A02p3B0,A1pB1,B0B0,sum1,sum2,B0C0,C0C0,C0C1,S,Coeff0;
	if(iWorkA1==iWorkB2) iWorkB2 = iWorkB1;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  S = pairConstants_->Sa0Par[i][j];
	  A0pB0  = ijShellPair->deltaAB[iWorkA1]+ijShellPair->deltaPA[iWorkA1][i][j];
	  A0p3B0 = ijShellPair->deltaAB[iWorkA1]+double(3.0)*ijShellPair->deltaPA[iWorkA1][i][j];
	  A02p3B0= math.two*ijShellPair->deltaAB[iWorkA1]+double(3.0)*ijShellPair->deltaPA[iWorkA1][i][j];
	  A1pB1  = ijShellPair->deltaAB[iWorkB2]+ijShellPair->deltaPA[iWorkB2][i][j];
	  B0B0   = ijShellPair->deltaPA[iWorkA1][i][j]*ijShellPair->deltaPA[iWorkA1][i][j];
	  sum1   = A0pB0*B0B0+A0p3B0*S;
	  Coeff0 = A1pB1*sum1;
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	    sum2   = (A0p3B0+double(3.0)*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom])*S;
	    B0C0   = ijShellPair->deltaPA[iWorkA1][i][j]*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    C0C0   = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    C0C1   = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB2][i][j][iAtom];
	    tmpVal += molecularConstants_->atomZ[iAtom]*
	      (pairConstants_->FmU[0][i][j][iAtom]*Coeff0
	       -pairConstants_->FmU[1][i][j][iAtom]*(A02p3B0*A1pB1*B0C0+pairConstants_->deltaPZ[iWorkB2][i][j][iAtom]*sum1+sum2*A1pB1)
	       +pairConstants_->FmU[2][i][j][iAtom]*((A0p3B0*C0C0+double(3.0)*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom])*A1pB1+(A02p3B0*B0C0+sum2)*pairConstants_->deltaPZ[iWorkB2][i][j][iAtom])
	       -pairConstants_->FmU[3][i][j][iAtom]*((A1pB1*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]+A0p3B0*pairConstants_->deltaPZ[iWorkB2][i][j][iAtom])*C0C0+double(3.0)*C0C1*S)
	       +pairConstants_->FmU[4][i][j][iAtom]*C0C0*C0C1);
	    //(A1+B1)*((A0+B0)*B0^2+(A0+3*B0)*S)
	    //(2*A0+3*B0)*(A1+B1)*B0*C0 + ((A0+B0)*B0^2+(A0+3*B0)*S)*C1 + (A0+3*B0+3*C0)*(A1+B1)*S
	    //((A0+3*B0)*C0^2+3*C0)*S*(A1+B1) + ((2*A0+3*B0)*B0*C0+(A0+3*B0+3*C0)*S)*C1
	    //((A1+B1)*C0+(A0+3*B0)*C1)*C0^2 + 3*C0*C1*S
	    //C0^2*C0*C1
	  };
	};
	return tmpVal;
      } else if(iWorkB1==iWorkB2&&iWorkA1!=iWorkA2&&(iWorkB1==iWorkA1||iWorkB1==iWorkA2)){
	/* (xy||xx) */ 
	double A0pB0,A0p3B0,A02p3B0,A1pB1,B0B0,sum1,sum2,B0C0,C0C0,C0C1,S,Coeff0;
	if(iWorkB1==iWorkA2) iWorkA2 = iWorkA1;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  S = pairConstants_->Sa0Par[i][j];
	  A0pB0  = -ijShellPair->deltaAB[iWorkB1]+ijShellPair->deltaPB[iWorkB1][i][j];
	  A0p3B0 = -ijShellPair->deltaAB[iWorkB1]+double(3.0)*ijShellPair->deltaPB[iWorkB1][i][j];
	  A02p3B0= -math.two*ijShellPair->deltaAB[iWorkB1]+double(3.0)*ijShellPair->deltaPB[iWorkB1][i][j];
	  A1pB1  = -ijShellPair->deltaAB[iWorkA2]+ijShellPair->deltaPB[iWorkA2][i][j];
	  B0B0   = ijShellPair->deltaPB[iWorkB1][i][j]*ijShellPair->deltaPB[iWorkB1][i][j];
	  sum1   = A0pB0*B0B0+A0p3B0*pairConstants_->Sa0Par[i][j];
	  Coeff0 = A1pB1*sum1;
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	    sum2   = (A0p3B0+double(3.0)*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom])*S;
	    B0C0   = ijShellPair->deltaPB[iWorkB1][i][j]*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom];
	    C0C0   = pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom];
	    C0C1   = pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA2][i][j][iAtom];
	    tmpVal += molecularConstants_->atomZ[iAtom]*
	      (pairConstants_->FmU[0][i][j][iAtom]*Coeff0
	       -pairConstants_->FmU[1][i][j][iAtom]*(A02p3B0*A1pB1*B0C0+pairConstants_->deltaPZ[iWorkA2][i][j][iAtom]*sum1+sum2*A1pB1)
	       +pairConstants_->FmU[2][i][j][iAtom]*((A0p3B0*C0C0+double(3.0)*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom])*A1pB1+(A02p3B0*B0C0+sum2)*pairConstants_->deltaPZ[iWorkA2][i][j][iAtom])
	       -pairConstants_->FmU[3][i][j][iAtom]*((A1pB1*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]+A0p3B0*pairConstants_->deltaPZ[iWorkA2][i][j][iAtom])*C0C0+double(3.0)*C0C1*S)
	       +pairConstants_->FmU[4][i][j][iAtom]*C0C0*C0C1);

	    //(A1+B1)*((A0+B0)*B0^2+(A0+3*B0)*S)
	    //(2*A0+3*B0)*(A1+B1)*B0*C0 + ((A0+B0)*B0^2+(A0+3*B0)*S)*C1 + (A0+3*B0+3*C0)*(A1+B1)*S
	    //((A0+3*B0)*C0^2+3*C0)*S*(A1+B1) + ((2*A0+3*B0)*B0*C0+(A0+3*B0+3*C0)*S)*C1
	    //((A1+B1)*C0+(A0+3*B0)*C1)*C0^2 + 3*C0*C1*S
	    //C0^2*C0*C1
	  };
	};
	return tmpVal;
      } else if(iWorkA1==iWorkA2&&iWorkB1!=iWorkB2&&iWorkA1!=iWorkB1&&iWorkA1!=iWorkB2){
	/* (xx||yz) */
	double A1pB1,A2pB2,prod1,B02pS,B0C0pS,sum1,C0C0,C1C2,S,Coeff0;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  S = pairConstants_->Sa0Par[i][j];
	  A1pB1 = ijShellPair->deltaAB[iWorkB1]+ijShellPair->deltaPA[iWorkB1][i][j];
	  A2pB2 = ijShellPair->deltaAB[iWorkB2]+ijShellPair->deltaPA[iWorkB2][i][j];
	  prod1 = A1pB1*A2pB2;
	  B02pS = ijShellPair->deltaPA[iWorkA1][i][j]*ijShellPair->deltaPA[iWorkA1][i][j]+S;
	  Coeff0= prod1*B02pS;
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	    B0C0pS= math.two*ijShellPair->deltaPA[iWorkA1][i][j]*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]+S;
	    sum1  = A2pB2*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]+A1pB1*pairConstants_->deltaPZ[iWorkB2][i][j][iAtom];
	    C0C0  = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    C1C2  = pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB2][i][j][iAtom];
	    tmpVal += molecularConstants_->atomZ[iAtom]*
	      (pairConstants_->FmU[0][i][j][iAtom]*Coeff0
	       -pairConstants_->FmU[1][i][j][iAtom]*(prod1*B0C0pS+sum1*B02pS)
	       +pairConstants_->FmU[2][i][j][iAtom]*(prod1*C0C0+sum1*B0C0pS+B02pS*C1C2)
	       -pairConstants_->FmU[3][i][j][iAtom]*(sum1*C0C0+B0C0pS*C1C2)
	       +pairConstants_->FmU[4][i][j][iAtom]*C0C0*C1C2);

	    //(A1+B1)*(A2+B2)*(B0^2+S)
	    //(A1+B1)*(A2+B2)*(2*B0*C0+S) + ((A2+B2)*C1+(A1+B1)*C2)*(B0^2+S)
	    //(A1+B1)*(A2+B2)*C0^2 + ((A2+B2)*C1+(A1+B1)*C2)*(2*B0*C0+S) + (B0^2+S)*C1*C2
	    //((A2+B2)*C1+(A1+B1)*C2)*C0^2 + (2*B0*C0+S)*C1*C2
	    //C0^2*C1*C2
	  };
	};
	return tmpVal;
      } else if(iWorkB1==iWorkB2&&iWorkA1!=iWorkA2&&iWorkB1!=iWorkA1&&iWorkB1!=iWorkA2){
	/* (yz||xx) */
	double A1pB1,A2pB2,prod1,B02pS,B0C0pS,sum1,C0C0,C1C2,S,Coeff0;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]){
	  S = pairConstants_->Sa0Par[i][j];
	  A1pB1 = -ijShellPair->deltaAB[iWorkA1]+ijShellPair->deltaPB[iWorkA1][i][j];
	  A2pB2 = -ijShellPair->deltaAB[iWorkA2]+ijShellPair->deltaPB[iWorkA2][i][j];
	  prod1 = A1pB1*A2pB2;
	  B02pS = ijShellPair->deltaPB[iWorkB1][i][j]*ijShellPair->deltaPB[iWorkB1][i][j]+S;
	  Coeff0= prod1*B02pS;
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	    B0C0pS= math.two*ijShellPair->deltaPB[iWorkB1][i][j]*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]+S;
	    sum1  = A2pB2*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]+A1pB1*pairConstants_->deltaPZ[iWorkA2][i][j][iAtom];
	    C0C0  = pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom];
	    C1C2  = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA2][i][j][iAtom];
	    tmpVal += molecularConstants_->atomZ[iAtom]*
	      (pairConstants_->FmU[0][i][j][iAtom]*Coeff0
	       -pairConstants_->FmU[1][i][j][iAtom]*(prod1*B0C0pS+sum1*B02pS)
	       +pairConstants_->FmU[2][i][j][iAtom]*(prod1*C0C0+sum1*B0C0pS+B02pS*C1C2)
	       -pairConstants_->FmU[3][i][j][iAtom]*(sum1*C0C0+B0C0pS*C1C2)
	       +pairConstants_->FmU[4][i][j][iAtom]*C0C0*C1C2);

	    //(A1+B1)*(A2+B2)*(B0^2+S)
	    //(A1+B1)*(A2+B2)*(2*B0*C0+S) + ((A2+B2)*C1+(A1+B1)*C2)*(B0^2+S)
	    //(A1+B1)*(A2+B2)*C0^2 + ((A2+B2)*C1+(A1+B1)*C2)*(2*B0*C0+S) + (B0^2+S)*C1*C2
	    //((A2+B2)*C1+(A1+B1)*C2)*C0^2 + (2*B0*C0+S)*C1*C2
	    //C0^2*C1*C2
	  };
	};
	return tmpVal;
      } else if(iWorkB1!=iWorkB2&&iWorkA1!=iWorkA2&&((iWorkB1==iWorkA1&&iWorkB2==iWorkA2)||(iWorkB1==iWorkA2&&iWorkB2==iWorkA1))){
	/* (xy||xy) */
	double A0pB0mB0,A1pB1mB1,sum1,SS,C0mA02B0,C1mA12B1,sum2,C0C0,C1C1,sum3,S,Coeff0;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  S = pairConstants_->Sa0Par[i][j];
	  A0pB0mB0 = (ijShellPair->deltaAB[iWorkA1]+ijShellPair->deltaPA[iWorkA1][i][j])*ijShellPair->deltaPA[iWorkA1][i][j];
	  A1pB1mB1 = (ijShellPair->deltaAB[iWorkA2]+ijShellPair->deltaPA[iWorkA2][i][j])*ijShellPair->deltaPA[iWorkA2][i][j];
	  sum1     = (A0pB0mB0+A1pB1mB1)*S;
	  SS       = S*S;
	  Coeff0   = A0pB0mB0*A1pB1mB1+sum1+SS;
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	    C0mA02B0 = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*(ijShellPair->deltaAB[iWorkA1]+math.two*ijShellPair->deltaPA[iWorkA1][i][j]);
	    C1mA12B1 = pairConstants_->deltaPZ[iWorkA2][i][j][iAtom]*(ijShellPair->deltaAB[iWorkA2]+math.two*ijShellPair->deltaPA[iWorkA2][i][j]);
	    sum2     = (C0mA02B0+C1mA12B1)*S;
	    C0C0     = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    C1C1     = pairConstants_->deltaPZ[iWorkA2][i][j][iAtom]*pairConstants_->deltaPZ[iWorkA2][i][j][iAtom];
	    sum3     = (C0C0+C1C1)*S;
	    tmpVal += molecularConstants_->atomZ[iAtom]*
	      (pairConstants_->FmU[0][i][j][iAtom]*Coeff0
	       -pairConstants_->FmU[1][i][j][iAtom]*(A1pB1mB1*C0mA02B0+A0pB0mB0*C1mA12B1+sum1+sum2+math.two*SS)
	       +pairConstants_->FmU[2][i][j][iAtom]*(A1pB1mB1*C0C0+C0mA02B0*C1mA12B1+A0pB0mB0*C1C1+sum2+sum3+SS)
	       -pairConstants_->FmU[3][i][j][iAtom]*(C1mA12B1*C0C0+C0mA02B0*C1C1+sum3)
	       +pairConstants_->FmU[4][i][j][iAtom]*C0C0*C1C1);
	    //(A0+B0)*(A1+B1)*B0*B1 + ((A0+B0)*B0+(A1+B1)*B1)*S + S^2
	    //(A1+B1)*B1*C0*(A0+2*B0) + (A0+B0)*B0*C1*(A1+2*B1) + ((A0+B0)*B0+(A1+B1)*B1)*S + ((A1+2*B1)*C1+(A0+2*B0)*C0)*S + 2*S^2
	    //(A1+B1)*B1*C0^2 + (A0+2*B0)*(A1+2*B1)*C0*C1 + (A0+B0)*B0*C1^2 + ((A0+2*B0)*C0+(A1+2*B1)*C1)*S + (C0^2+C1^2)*S + S^2
	    //(A1+2*B1)*C0^2*C1 + (A0+2*B0)*C0*C1^2 + (C0^2+C1^2)*S
	    //C0^2*C1^2
	  };
	};
	return tmpVal;
      } else {
	/* (xy||yz) */
	if(iWorkA1==iWorkB1) iWorkA1=iWorkA2;
	else if(iWorkA1==iWorkB2) {
	  iWorkB2 = iWorkB1;
	  iWorkB1 = iWorkA1;
	  iWorkA1 = iWorkA2;
	} else if(iWorkA2==iWorkB2) {
	  iWorkB2 = iWorkB1;
	  iWorkB1 = iWorkA2;
	};
	double A1pB1mB1,A2pB2,A2pB2mB0,A2pB2mC0,B0C2,A1p2B1mC1,C1C1,C0C2,S,Coeff0;
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  S = pairConstants_->Sa0Par[i][j];
	  A1pB1mB1 = (ijShellPair->deltaAB[iWorkB1]+ijShellPair->deltaPA[iWorkB1][i][j])*ijShellPair->deltaPA[iWorkB1][i][j];
	  A2pB2    = ijShellPair->deltaAB[iWorkB2]+ijShellPair->deltaPA[iWorkB2][i][j];
	  A2pB2mB0 = A2pB2*ijShellPair->deltaPA[iWorkA1][i][j];
	  Coeff0   = (A1pB1mB1+pairConstants_->Sa0Par[i][j])*A2pB2mB0;
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	    A2pB2mC0 = A2pB2*pairConstants_->deltaPZ[iWorkA1][i][j][iAtom];
	    B0C2     = ijShellPair->deltaPA[iWorkA1][i][j]*pairConstants_->deltaPZ[iWorkB2][i][j][iAtom];
	    A1p2B1mC1= (ijShellPair->deltaAB[iWorkB1]+math.two*ijShellPair->deltaPA[iWorkB1][i][j])*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom];
	    C1C1     = pairConstants_->deltaPZ[iWorkB1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB1][i][j][iAtom];
	    C0C2     = pairConstants_->deltaPZ[iWorkA1][i][j][iAtom]*pairConstants_->deltaPZ[iWorkB2][i][j][iAtom];
	    tmpVal += molecularConstants_->atomZ[iAtom]*
	      (pairConstants_->FmU[0][i][j][iAtom]*Coeff0
	       -pairConstants_->FmU[1][i][j][iAtom]*(A1pB1mB1*(A2pB2mC0+B0C2)+A1p2B1mC1*A2pB2mB0+(A2pB2mB0+A2pB2mC0+B0C2)*pairConstants_->Sa0Par[i][j])
	       +pairConstants_->FmU[2][i][j][iAtom]*(A2pB2mC0*A1p2B1mC1+A2pB2mB0*C1C1+A1pB1mB1*C0C2+A1p2B1mC1*B0C2+(A2pB2mC0+B0C2+C0C2)*pairConstants_->Sa0Par[i][j])
	       -pairConstants_->FmU[3][i][j][iAtom]*((A2pB2mC0+B0C2)*C1C1+(A1p2B1mC1+pairConstants_->Sa0Par[i][j])*C0C2)
	       +pairConstants_->FmU[4][i][j][iAtom]*C0C2*C1C1);
	    //((A1+B1)*B1+S)*B0*(A2+B2)
	    //(A1+B1)*B1*((A2+B2)*C0+B0*C2) + (A1+2*B1)*(A2+B2)*B0*C1 + ((A2+B2)*B0+(A2+B2)*C0+B0*C2)*S	  
	    //(A2+B2)*C0*(A1+2*B1)*C1 + (A2+B2)*B0*C1^2 + (A1+B1)*B1*(C0*C2) + (A1+2*B1)*C1*(B0*C2) + (A2+B2)*C0*S + (B0*C2+C0*C2)*S
	    //((A2+B2)*C0+B0*C2)*C1^2 + ((A1+2*B1)*C1+S)*(C0*C2)
	    //C0*C1^2*C2
	  };
	};
	return tmpVal;
      };
    };
  } else {
    int iWork;
    int lAp1[3],lBm1[3];
    for(iWork=0;iWork<3;iWork++) {
      lAp1[iWork]=lA[iWork];
      lBm1[iWork]=lB[iWork];
    };
    if (lB[0]>0) iWork=0;
    else if (lB[1]>0) iWork=1;
    else if (lB[2]>0) iWork=2;
    if(LB>2){
      lAp1[iWork]++;
      lBm1[iWork]--;
      tmpVal = this->oneehRRVab(ijShellPair,LA+1,lAp1,LB-1,lBm1);
      if(abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) tmpVal+= ijShellPair->deltaAB[iWork]*this->oneehRRVab(ijShellPair,LA,lA,LB-1,lBm1);
    } else {
      lAp1[iWork]++;
      if(abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) {
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++) if(abs(pairConstants_->FmU[0][i][j][iAtom])>this->controls_->thresholdS)
	    tmpVal += molecularConstants_->atomZ[iAtom]*(this->oneevRRVa0(ijShellPair,&iAtom,0,LA+1,lAp1,&i,&j)+ ijShellPair->deltaAB[iWork]*this->oneevRRVa0(ijShellPair,&iAtom,0,LA,lA,&i,&j));
	};
      } else {
	for(i=0;i<ijShellPair->nPGTOs[0];i++) for(j=0;j<ijShellPair->nPGTOs[1];j++) if(pairConstants_->ssNonzero[i][j]) {
	  for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++) if(abs(pairConstants_->FmU[0][i][j][iAtom])>this->controls_->thresholdS)
	    tmpVal += molecularConstants_->atomZ[iAtom]*this->oneevRRVa0(ijShellPair,&iAtom,0,LA+1,lAp1,&i,&j);
	};
      };
    };
    return tmpVal;
  };
};
/***********************************************/
/* potential vertical horizontal recursion     */
/***********************************************/
double AOIntegrals::oneevRRVa0(ShellPair *ijShellPair,int *iAtom,int m,int LA,int *lA, int *i, int *j){
  if(LA==0) return pairConstants_->FmU[m][*i][*j][*iAtom];
  double tmpVal=0.0;
  int lAm1[3];
  int iWork;
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0) iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;
  if(LA>3){
    lAm1[iWork] -=1;
    if(abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) tmpVal  = ijShellPair->deltaPA[iWork][*i][*j]*this->oneevRRVa0(ijShellPair,iAtom,m,LA-1,lAm1,i,j);
    if(abs(pairConstants_->deltaPZ[iWork][*i][*j][*iAtom])>this->controls_->thresholdS) tmpVal -= pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*this->oneevRRVa0(ijShellPair,iAtom,m+1,LA-1,lAm1,i,j);
    if(LA==2&&lA[iWork]==2)
      tmpVal += pairConstants_->Sa0Par[*i][*j]*(pairConstants_->FmU[m][*i][*j][*iAtom]-pairConstants_->FmU[m+1][*i][*j][*iAtom]);
    else if (lA[iWork]>=2) {
      lAm1[iWork] -=1;
      tmpVal += (lAm1[iWork]+1)*pairConstants_->Sa0Par[*i][*j]*(this->oneevRRVa0(ijShellPair,iAtom,m,LA-2,lAm1,i,j)-this->oneevRRVa0(ijShellPair,iAtom,m+1,LA-2,lAm1,i,j));
    };
  } else if (LA==3) {
    if(lA[0]==lA[1]&&lA[1]==lA[2]) {
      lAm1[iWork]--;
      int iWork2,iWork3;
      if (lAm1[0]>0)      iWork2=0;
      else if (lAm1[1]>0) iWork2=1;
      else if (lAm1[2]>0) iWork2=2;
      lAm1[iWork2]--;
      if (lAm1[0]>0)      iWork3=0;
      else if (lAm1[1]>0) iWork3=1;
      else if (lAm1[2]>0) iWork3=2;

      tmpVal += pairConstants_->FmU[m][*i][*j][*iAtom]*ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j]*ijShellPair->deltaPA[iWork3][*i][*j];
      tmpVal -= pairConstants_->FmU[m+1][*i][*j][*iAtom]*(ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j]*pairConstants_->deltaPZ[iWork3][*i][*j][*iAtom]+ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork3][*i][*j]*pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom]+ijShellPair->deltaPA[iWork2][*i][*j]*ijShellPair->deltaPA[iWork3][*i][*j]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]);
      tmpVal += pairConstants_->FmU[m+2][*i][*j][*iAtom]*(ijShellPair->deltaPA[iWork][*i][*j]*pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork3][*i][*j][*iAtom]+ijShellPair->deltaPA[iWork2][*i][*j]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork3][*i][*j][*iAtom]+ijShellPair->deltaPA[iWork3][*i][*j]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom]);
      tmpVal -= pairConstants_->FmU[m+3][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork3][*i][*j][*iAtom];
    } else if(lA[0]==3||lA[1]==3||lA[2]==3) {
      double sqrPA = ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork][*i][*j];
      double sqrPZ = pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom];
      tmpVal += pairConstants_->FmU[m][*i][*j][*iAtom]*ijShellPair->deltaPA[iWork][*i][*j]*(sqrPA+math.three*pairConstants_->Sa0Par[*i][*j]);
      tmpVal -= pairConstants_->FmU[m+1][*i][*j][*iAtom]*math.three*(pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*sqrPA+pairConstants_->Sa0Par[*i][*j]*(ijShellPair->deltaPA[iWork][*i][*j]+pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]));
      tmpVal += pairConstants_->FmU[m+2][*i][*j][*iAtom]*math.three*(ijShellPair->deltaPA[iWork][*i][*j]*sqrPZ+pairConstants_->Sa0Par[*i][*j]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]);
      tmpVal -= pairConstants_->FmU[m+3][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*sqrPZ;
    } else if(lA[0]==2||lA[1]==2||lA[2]==2) {
      lAm1[iWork]--;
      int iWork2,iWork3;
      if (lAm1[0]>0)      iWork2=0;
      else if (lAm1[1]>0) iWork2=1;
      else if (lAm1[2]>0) iWork2=2;
      lAm1[iWork2]--;
      if (lAm1[0]>0)      iWork3=0;
      else if (lAm1[1]>0) iWork3=1;
      else if (lAm1[2]>0) iWork3=2;

      if(iWork==iWork2) iWork2 = iWork3;
      else if(iWork2==iWork3) {
	iWork2=iWork;
	iWork=iWork3;
      };

      double sqrPA = ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork][*i][*j];
      double sqrPZ = pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom];
      double PAPZ  = ijShellPair->deltaPA[iWork][*i][*j]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom];
      tmpVal += pairConstants_->FmU[m][*i][*j][*iAtom]*ijShellPair->deltaPA[iWork2][*i][*j]*(sqrPA+pairConstants_->Sa0Par[*i][*j]);
      tmpVal -= pairConstants_->FmU[m+1][*i][*j][*iAtom]*(pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom]*(sqrPA+pairConstants_->Sa0Par[*i][*j])+ijShellPair->deltaPA[iWork2][*i][*j]*(pairConstants_->Sa0Par[*i][*j]+math.two*PAPZ));
      tmpVal += pairConstants_->FmU[m+2][*i][*j][*iAtom]*(ijShellPair->deltaPA[iWork2][*i][*j]*sqrPZ+pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom]*(pairConstants_->Sa0Par[*i][*j]+math.two*PAPZ));
      tmpVal -= pairConstants_->FmU[m+3][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom]*sqrPZ;
    }; 
  } else if (LA==2) {
    lAm1[iWork]--;
    int iWork2;
    if (lAm1[0]>0)      iWork2=0;
    else if (lAm1[1]>0) iWork2=1;
    else if (lAm1[2]>0) iWork2=2;
    if(iWork!=iWork2) {
      tmpVal += pairConstants_->FmU[m][*i][*j][*iAtom]*ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork2][*i][*j];
      tmpVal -= pairConstants_->FmU[m+1][*i][*j][*iAtom]*(pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom]*ijShellPair->deltaPA[iWork][*i][*j]+ijShellPair->deltaPA[iWork2][*i][*j]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]);
      tmpVal += pairConstants_->FmU[m+2][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom];
    } else {
      tmpVal += pairConstants_->FmU[m][*i][*j][*iAtom]*(ijShellPair->deltaPA[iWork][*i][*j]*ijShellPair->deltaPA[iWork][*i][*j]+pairConstants_->Sa0Par[*i][*j]);
      tmpVal -= pairConstants_->FmU[m+1][*i][*j][*iAtom]*(math.two*ijShellPair->deltaPA[iWork][*i][*j]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]+pairConstants_->Sa0Par[*i][*j]);
      tmpVal += pairConstants_->FmU[m+2][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork2][*i][*j][*iAtom];
    };
   } else{
    tmpVal  = pairConstants_->FmU[m][*i][*j][*iAtom]*ijShellPair->deltaPA[iWork][*i][*j];
    tmpVal -= pairConstants_->FmU[m+1][*i][*j][*iAtom]*pairConstants_->deltaPZ[iWork][*i][*j][*iAtom];
  };
  return tmpVal; 
};

