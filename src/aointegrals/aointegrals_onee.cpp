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

  if(LA==2&&lA[iWork]==2) tmpVal += ijSP->halfInvZeta[iPGTOPair];
  else if (lA[iWork]>=2) {
    lAm1[iWork]--;
    tmpVal += (lAm1[iWork]+1)*ijSP->halfInvZeta[iPGTOPair]*this->vRRSa0(ijSP,LA-2,lAm1,iPGTOPair);
  };

  return tmpVal;
};

/*
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
*/
