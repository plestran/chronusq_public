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
/*
#include <aointegrals.h>
using ChronusQ::AOIntegrals;
//----------------------------------------------------//
// two-e horizontal recursion from (ab|cd) to (a0|cd) //
// (ab|cd)=(a+1,b-1|cd)+(A-B)*(a,b-1|cd)              //
//----------------------------------------------------//
double AOIntegrals::twoehRRabcd(  int *nPGTOs, ShellPair *ijShellPair, ShellPair *klShellPair,
                                int LA,int *lA,int LB,int *lB,int LC,int *lC,int LD,int *lD) {
  double tmpVal = 0.0, tmpVal1=0.0;
// iWork is used to indicate which Cartesian angular momentum we are reducing (x,y,z)
  int iWork,jWork,kWork,lWork;
  int totalL = LA + LB + LC + LD;
  int i,j,k,l,ij,kl;
  if(totalL==0) {
    return this->twoeSSSS0(nPGTOs, ijShellPair, klShellPair);
  } else if(totalL==1) {
    // (ps|ss) 
    if (lA[0]>0)      iWork=0;  
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
    int nFmT = 0;
    for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
      tmpVal += (this->quartetConstants_->FmT[0][i][j][k][l]*ijShellPair->deltaPA[iWork][i][j]
                +this->quartetConstants_->FmT[1][i][j][k][l]*this->quartetConstants_->deltaWP[iWork][i][j][k][l]);
    };
    return tmpVal;
  } else if (totalL>4||LA==2||LB==2||LC==2||LD==2) {
    if(LB==0) {
      return this->twoehRRa0cd(nPGTOs,ijShellPair,klShellPair,LA,lA,LC,lC,LD,lD);
    };
    int iWork;
    int lAp1[3],lBm1[3];
    for(iWork=0;iWork<3;iWork++){
      lAp1[iWork]=lA[iWork];     
      lBm1[iWork]=lB[iWork];     
    };
    if (lB[0]>0)      iWork=0;   
    else if (lB[1]>0) iWork=1;
    else if (lB[2]>0) iWork=2;
    lAp1[iWork]+=1;
    if(std::abs(ijShellPair->deltaAB[iWork])>this->controls_->thresholdS) {
      if(LB>1) {
        lBm1[iWork]-=1;
        tmpVal =    this->twoehRRabcd(nPGTOs,ijShellPair,klShellPair,LA+1,lAp1,LB-1,lBm1,LC,lC,LD,lD) 
                    +ijShellPair->deltaAB[iWork]*this->twoehRRabcd(nPGTOs,ijShellPair,klShellPair,LA,lA,LB-1,lBm1,LC,lC,LD,lD);
      } else {
        tmpVal =    this->twoehRRa0cd(nPGTOs,ijShellPair,klShellPair,LA+1,lAp1,LC,lC,LD,lD) 
                    +ijShellPair->deltaAB[iWork]*this->twoehRRa0cd(nPGTOs,ijShellPair,klShellPair,LA,lA,LC,lC,LD,lD);
        return tmpVal;
      };
    } else {
      if(LB>1) {
        lBm1[iWork]-=1;
        tmpVal =    this->twoehRRabcd(nPGTOs,ijShellPair,klShellPair,LA+1,lAp1,LB-1,lBm1,LC,lC,LD,lD);
      } else {
        tmpVal =    this->twoehRRa0cd(nPGTOs,ijShellPair,klShellPair,LA+1,lAp1,LC,lC,LD,lD);
      };
    };
    return tmpVal;
  } else if(totalL==2) {
    return this->twoepp00(nPGTOs,ijShellPair,klShellPair,LA,lA,LB,lB,LC,lC,LD,lD);
  } else if(totalL==3) {
    return this->twoeppp0(nPGTOs,ijShellPair,klShellPair,LA,lA,LB,lB,LC,lC,LD,lD);
  } else if(totalL==4) {
    return this->twoepppp(nPGTOs,ijShellPair,klShellPair,LA,lA,LB,lB,LC,lC,LD,lD);
  };
};
//----------------------------------------------------//
// two-e horizontal recursion from (a0|cd) to (a0|c0) //
// (a0|cd)=(a,0|c+1,d-1)+(C-D)*(a,0|c,d-1)            //
//----------------------------------------------------//
double AOIntegrals::twoehRRa0cd(  int *nPGTOs, ShellPair *ijShellPair, ShellPair *klShellPair,
                                int LA,int *lA,int LC,int *lC,int LD,int *lD)  {
  int i=0,j=0,k=0,l=0;
  int iWork;
  int lCp1[3],lDm1[3];  
  double tmpVal=0.0;
  for(iWork=0;iWork<3;iWork++){
    lCp1[iWork]=lC[iWork];
    lDm1[iWork]=lD[iWork];
  };
  if (lD[0]>0) iWork=0;
  else if (lD[1]>0) iWork=1;
  else if (lD[2]>0) iWork=2;
  lCp1[iWork]+=1;
  if(LD==0) {
    for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++)
      tmpVal +=     this->twoevRRa0c0(ijShellPair,klShellPair,0,LA,lA,LC,lC,&i,&j,&k,&l);
  } else {
    if(std::abs(klShellPair->deltaAB[iWork])>this->controls_->thresholdAB) {
      if(LD>1) { // (a0|cd)
        lDm1[iWork]-=1;
        tmpVal =    this->twoehRRa0cd(nPGTOs,ijShellPair,klShellPair,LA,lA,LC+1,lCp1,LD-1,lDm1) 
                    +klShellPair->deltaAB[iWork]*this->twoehRRa0cd(nPGTOs,ijShellPair,klShellPair,LA,lA,LC,lC,LD-1,lDm1);
      } else {
        for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++)
          tmpVal += (this->twoevRRa0c0(ijShellPair,klShellPair,0,LA,lA,LC+1,lCp1,&i,&j,&k,&l) 
                    +klShellPair->deltaAB[iWork]*this->twoevRRa0c0(ijShellPair,klShellPair,0,LA,lA,LC,lC,&i,&j,&k,&l));
      };
    } else {
      if(LD>1) {
        lDm1[iWork]-=1;
        tmpVal =    this->twoehRRa0cd(nPGTOs,ijShellPair,klShellPair,LA,lA,LC+1,lCp1,LD-1,lDm1);
      } else {
        for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++)
          tmpVal +=     this->twoevRRa0c0(ijShellPair,klShellPair,0,LA,lA,LC+1,lCp1,&i,&j,&k,&l);
      };
    };
  };
  return tmpVal;
};

double AOIntegrals::twoepp00( int *nPGTOs, ShellPair *ijShellPair, ShellPair *klShellPair,
                            int LA,int *lA,int LB,int *lB,int LC,int *lC,int LD,int *lD) {
  double tmpVal = math.zero;
  double Par0,Par1,Par2;
  int iWork,jWork,kWork,lWork,i,j,k,l;
  if(LC==0) {
    if (lA[0]>0)      iWork=0;
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
    if (lB[0]>0)      jWork=0;  
    else if (lB[1]>0) jWork=1;
    else if (lB[2]>0) jWork=2;
    if(iWork==jWork){ //(xx|00)
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) {
        Par0 = ijShellPair->inversezeta[i][j]+ijShellPair->deltaPAtPB[iWork][iWork][i][j];
        for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
          Par1    =     -this->quartetConstants_->a000Par2[i][j][k][l]*ijShellPair->inversezeta[i][j]
                        +this->quartetConstants_->deltaWP[iWork][i][j][k][l]*ijShellPair->deltaPApPB[iWork][iWork][i][j];
          Par2    =     this->quartetConstants_->deltaWP[iWork][i][j][k][l]*this->quartetConstants_->deltaWP[iWork][i][j][k][l];
          tmpVal +=     (this->quartetConstants_->FmT[0][i][j][k][l]*Par0
                        +this->quartetConstants_->FmT[1][i][j][k][l]*Par1
                        +this->quartetConstants_->FmT[2][i][j][k][l]*Par2);
        };
      };
      return tmpVal;
    } else { //(xy|00)
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) {
        for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
          Par0    =     ijShellPair->deltaPAtPB[iWork][jWork][i][j];
          Par1    =     this->quartetConstants_->deltaWP[iWork][i][j][k][l]*ijShellPair->deltaPB[jWork][i][j]
                        +this->quartetConstants_->deltaWP[jWork][i][j][k][l]*ijShellPair->deltaPA[iWork][i][j];
          Par2    =     this->quartetConstants_->deltaWP[iWork][i][j][k][l]*this->quartetConstants_->deltaWP[jWork][i][j][k][l];
          tmpVal +=     (this->quartetConstants_->FmT[0][i][j][k][l]*Par0
                        +this->quartetConstants_->FmT[1][i][j][k][l]*Par1
                        +this->quartetConstants_->FmT[2][i][j][k][l]*Par2);
        };
      };
      return tmpVal;
    };
  } else if(LC==1) { //<p0|p0>
    if (lA[0]>0)      iWork=0;  
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
    if (lC[0]>0)      kWork=0;  
    else if (lC[1]>0) kWork=1;
    else if (lC[2]>0) kWork=2;
    if(iWork==kWork){ //(x0|x0)
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) {
        for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {             
          Par0    =     ijShellPair->deltaPA[iWork][i][j]*klShellPair->deltaPA[iWork][k][l];
          Par1    =     klShellPair->deltaPA[iWork][k][l]*this->quartetConstants_->deltaWP[iWork][i][j][k][l]
                        +ijShellPair->deltaPA[iWork][i][j]*this->quartetConstants_->deltaWQ[iWork][i][j][k][l]
                        +this->quartetConstants_->a0c0Par3[i][j][k][l];
          Par2    =     this->quartetConstants_->deltaWP[iWork][i][j][k][l]*this->quartetConstants_->deltaWQ[iWork][i][j][k][l];
          tmpVal +=     (this->quartetConstants_->FmT[0][i][j][k][l]*Par0
                        +this->quartetConstants_->FmT[1][i][j][k][l]*Par1
                        +this->quartetConstants_->FmT[2][i][j][k][l]*Par2);
        };
      };
      return tmpVal;
    } else { //(x0|y0)
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) {
        for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
          Par0 = ijShellPair->deltaPA[iWork][i][j]*klShellPair->deltaPA[kWork][k][l];
          Par1 = klShellPair->deltaPA[kWork][k][l]*this->quartetConstants_->deltaWP[iWork][i][j][k][l]+ijShellPair->deltaPA[iWork][i][j]*this->quartetConstants_->deltaWQ[kWork][i][j][k][l];
          Par2 = this->quartetConstants_->deltaWP[iWork][i][j][k][l]*this->quartetConstants_->deltaWQ[kWork][i][j][k][l];
          tmpVal +=     (this->quartetConstants_->FmT[0][i][j][k][l]*Par0
                        +this->quartetConstants_->FmT[1][i][j][k][l]*Par1
                        +this->quartetConstants_->FmT[2][i][j][k][l]*Par2);
        };
      };
      return tmpVal;
    };
  };
};

double AOIntegrals::twoeppp0( int *nPGTOs,ShellPair *ijShellPair,ShellPair *klShellPair,
                            int LA,int *lA,int LB,int *lB,int LC,int *lC,int LD,int *lD) {
  double tmpVal = math.zero;
  double CO1,CO2,CO3,CO4,CO5;
  double tC1,tC2,tC3,tC4,tC5,tC6,tC7,tC8,tC9,tC0,tCa,tCb,tCc,tCd,tCe,tCf;
  double AB0,AB1,AB2,CD0,CD1,CD2,PA0,PA1,PA2,QC0,QC1,QC2,WP0,WP1,WP2,WQ0,WQ1,WQ2;
  int iWork,jWork,kWork,lWork,i,j,k,l;
  if (lA[0]>0)      iWork=0;  
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;
  if (lB[0]>0)      jWork=0;  
  else if (lB[1]>0) jWork=1;
  else if (lB[2]>0) jWork=2;
  if (lC[0]>0)      kWork=0;  
  else if (lC[1]>0) kWork=1;
  else if (lC[2]>0) kWork=2;
  if(iWork==jWork){
    if(kWork==iWork){ //(xx|x0)
      AB0=ijShellPair->deltaAB[iWork];
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
        CO1=ijShellPair->inversezeta[i][j];
        PA0=ijShellPair->deltaPA[iWork][i][j];
        tC1=ijShellPair->deltaPB[iWork][i][j];
        tC2=PA0*tC1+CO1;
        CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
        CO5=this->quartetConstants_->a000Par2[i][j][k][l];
        QC0=klShellPair->deltaPA[iWork][k][l];
        WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
        WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];
        tC3=QC0*WP0+math.two*CO3;
        tC4=WP0*tC1-CO1*CO5;
        tC5=WQ0*WP0;
        tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(QC0*tC2)
                    +this->quartetConstants_->FmT[1][i][j][k][l]*(PA0*tC3+QC0*tC4+WQ0*tC2+AB0*CO3)
                    +this->quartetConstants_->FmT[2][i][j][k][l]*(WP0*tC3+PA0*tC5+WQ0*tC4)
                    +this->quartetConstants_->FmT[3][i][j][k][l]*(WP0*tC5));
      };
      return tmpVal;
    } else { //(xx|y0)
      AB0=ijShellPair->deltaAB[iWork];
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
        CO1=ijShellPair->inversezeta[i][j];
        PA0=ijShellPair->deltaPA[iWork][i][j];
        tC1=ijShellPair->deltaPB[iWork][i][j];
        tC2=tC1+PA0;
        tC3=PA0*tC1+CO1;
        CO5=this->quartetConstants_->a000Par2[i][j][k][l];
        WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
        QC1=klShellPair->deltaPA[kWork][k][l];
        WQ1=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];
        tC4=WP0*tC2-CO1*CO5;
        tC5=WP0*WP0;
        tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(QC1*tC3)
                    +this->quartetConstants_->FmT[1][i][j][k][l]*(QC1*tC4+WQ1*tC3)
                    +this->quartetConstants_->FmT[2][i][j][k][l]*(QC1*tC5+WQ1*tC4)
                    +this->quartetConstants_->FmT[3][i][j][k][l]*(WQ1*tC5));
      };
      return tmpVal;
    };
  } else {
    if(kWork==iWork){ //(xy|x0)
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
        CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
        PA0=ijShellPair->deltaPA[iWork][i][j];
        QC0=klShellPair->deltaPA[iWork][k][l];
        WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
        WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];
        AB1=ijShellPair->deltaAB[jWork];
        PA1=ijShellPair->deltaPA[jWork][i][j];
        WP1=this->quartetConstants_->deltaWP[jWork][i][j][k][l];                  
        tC1=ijShellPair->deltaPB[jWork][i][j];
        tC2=PA0*QC0;             
        tC3=PA0*WQ0+QC0*WP0+CO3; 
        tC4=WP0*WQ0;                            
        tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC1*tC2)
                    +this->quartetConstants_->FmT[1][i][j][k][l]*(tC1*tC3+WP1*tC2)
                    +this->quartetConstants_->FmT[2][i][j][k][l]*(tC1*tC4+WP1*tC3)
                    +this->quartetConstants_->FmT[3][i][j][k][l]*(WP1*tC4));
      };
      return tmpVal;
    } else if(kWork==jWork){ //(xy|y0)   ===>>  (yx|y0)
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
        CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
        PA0=ijShellPair->deltaPB[jWork][i][j];
        QC0=klShellPair->deltaPA[jWork][k][l];
        WP0=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
        WQ0=this->quartetConstants_->deltaWQ[jWork][i][j][k][l];
        AB1=-ijShellPair->deltaAB[iWork];
        PA1=ijShellPair->deltaPB[iWork][i][j];
        WP1=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
        tC1=ijShellPair->deltaPA[iWork][i][j];
        tC2=PA0*QC0; 
        tC3=PA0*WQ0+QC0*WP0+CO3;
        tC4=WP0*WQ0;
        tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC1*tC2)
                    +this->quartetConstants_->FmT[1][i][j][k][l]*(tC1*tC3+WP1*tC2)
                    +this->quartetConstants_->FmT[2][i][j][k][l]*(tC1*tC4+WP1*tC3)
                    +this->quartetConstants_->FmT[3][i][j][k][l]*(WP1*tC4));
      };
      return tmpVal;
    } else { //(xy|z0)
      for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
        PA0=ijShellPair->deltaPA[iWork][i][j];
        WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
        AB1=ijShellPair->deltaAB[jWork];
        PA1=ijShellPair->deltaPA[jWork][i][j];
        WP1=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
        QC2=klShellPair->deltaPA[kWork][k][l];
        WQ2=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];
        tC1=ijShellPair->deltaPB[jWork][i][j];
        tC2=PA0*QC2;
        tC3=PA0*WQ2+QC2*WP0;
        tC4=WP0*WQ2;
        tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC1*tC2)
                    +this->quartetConstants_->FmT[1][i][j][k][l]*(tC1*tC3+WP1*tC2)
                    +this->quartetConstants_->FmT[2][i][j][k][l]*(WP1*tC3+tC1*tC4)
                    +this->quartetConstants_->FmT[3][i][j][k][l]*(WP1*tC4));
      };
      return tmpVal;
    };
  };
};

double AOIntegrals::twoepppp(int *nPGTOs,ShellPair *ijShellPair,ShellPair *klShellPair,int LA,int *lA,int LB,int *lB,int LC,int *lC,int LD,int *lD) {
  double tmpVal = 0.0;
  double CO1,CO2,CO3,CO4,CO5; //depent on exponent of 'i,j,k,l'
  double tC1,tC2,tC3,tC4,tC5,tC6,tC7,tC8,tC9,tC0,tCa,tCb,tCc,tCd,tCe,tCf;
  double AB0,AB1,AB2,CD0,CD1,CD2,PA0,PA1,PA2,QC0,QC1,QC2,WP0,WP1,WP2,WQ0,WQ1,WQ2;
  int i,j,k,l;
  int iWork,jWork,kWork,lWork;
  int totalL=LA+LB+LC+LD;
  if(LA==1){
    //<pp|pp>
    if (lA[0]>0)      iWork=0;  
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
    if (lB[0]>0)      jWork=0;  
    else if (lB[1]>0) jWork=1;
    else if (lB[2]>0) jWork=2;
    if (lC[0]>0)      kWork=0;  
    else if (lC[1]>0) kWork=1;
    else if (lC[2]>0) kWork=2;
    if (lD[0]>0)      lWork=0;  
    else if (lD[1]>0) lWork=1;
    else if (lD[2]>0) lWork=2;
    //    format 1: <xx|xx>
    //    format 2: <xx|yy>
    //    format 3: <xx|xy>
    //    format 4: <xx|yz>
    //    format 5: <xy|xy>
    //    format 6: <xy|xz>
    if(iWork==jWork){    //<xx|??>
      if(kWork==iWork){  //<xx|x?>
        if(lWork==iWork){//<xx|xx> ==>>format 1         
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
            CO1=ijShellPair->inversezeta[i][j]; 
            CO2=klShellPair->inversezeta[k][l];
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
            CO4=this->quartetConstants_->a0c0Par2[i][j][k][l];
            CO5=this->quartetConstants_->a000Par2[i][j][k][l];
            AB0=ijShellPair->deltaAB[iWork];
            CD0=klShellPair->deltaAB[iWork];
            PA0=ijShellPair->deltaPA[iWork][i][j];
            QC0=klShellPair->deltaPA[iWork][k][l];
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];
//          tC1=AB0+PA0;
            tC1=ijShellPair->deltaPB[iWork][i][j];
            tC2=CD0+QC0;
            tC3=tC1+PA0;
            tC4=tC2+QC0;
            tC5=CO1+PA0*tC1;
            tC6=CO2+QC0*tC2;
            tCa=CO1*CO5;
            tCb=CO2*CO4;
            tC7=WP0*tC3-tCa;
            tC8=WQ0*tC4-tCb;
            tC9=WP0*WP0;
            tC0=WQ0*WQ0;                  
            tCc=math.two*QC0*tC1+math.two*PA0*tC2+AB0*CD0;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC5*tC6)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC5*tC8+tC6*tC7+tCc*CO3)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC5*tC0+tC6*tC9+tCc*WP0*WQ0-tCa*tC8-tCb*tC7-tCa*tCb+math.two*CO3*(WQ0*tC3+WP0*tC4+CO3))
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC8*tC9+tC7*tC0+math.four*CO3*WP0*WQ0)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC9*tC0));
          };
          return tmpVal;
        } else {//<xx|xy> or <xx|xz> ==>format 3          
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
            CO1=ijShellPair->inversezeta[i][j];                  
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
            CO4=this->quartetConstants_->a0c0Par2[i][j][k][l];
            CO5=this->quartetConstants_->a000Par2[i][j][k][l];            
            AB0=ijShellPair->deltaAB[iWork]; 
            PA0=ijShellPair->deltaPA[iWork][i][j];
            QC0=klShellPair->deltaPA[iWork][k][l];
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];             
            CD1=klShellPair->deltaAB[lWork];             
            QC1=klShellPair->deltaPA[lWork][k][l];               
            WQ1=this->quartetConstants_->deltaWQ[lWork][i][j][k][l];
            tC1=ijShellPair->deltaPB[iWork][i][j];
            tC2=CD1+QC1;
            tC3=tC1+PA0;
            tC4=tC1*PA0;
            tCa=tC4+CO1;
            tCb=CO3+QC0*WP0;
            tC5=tCa*QC0;
            tC6=tCa*WQ0 + tCb*tC3 - CO5*QC0*CO1;
            tC7=(CO3+tCb)*WP0+(tC3*WP0-CO1*CO5)*WQ0;
            tCc=WP0*WP0*WQ0;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC5*tC2)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC5*WQ1+tC6*tC2)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC6*WQ1+tC7*tC2)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC7*WQ1+tC2*tCc)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tCc*WQ1));
          };
          return tmpVal;
        };
      } else {              //<xx|y?>
        if(lWork==iWork){//<xx|yx>  ===>format 3    
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
            CO1=ijShellPair->inversezeta[i][j];                  
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
            CO4=this->quartetConstants_->a0c0Par2[i][j][k][l];
            CO5=this->quartetConstants_->a000Par2[i][j][k][l];            
            AB0=ijShellPair->deltaAB[iWork]; 
            PA0=ijShellPair->deltaPA[iWork][i][j];
            QC0=klShellPair->deltaPB[iWork][k][l];              //QD
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];             
            CD1=-klShellPair->deltaAB[kWork];                  //-CD
            QC1=klShellPair->deltaPB[kWork][k][l];              //QD 
            WQ1=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];
            tC1=ijShellPair->deltaPB[iWork][i][j];
            tC2=CD1+QC1;
            tC3=tC1+PA0;
            tC4=tC1*PA0;
            tCa=tC4+CO1;
            tCb=CO3+QC0*WP0;
            tC5=tCa*QC0;
            tC6=tCa*WQ0 + tCb*tC3 - CO5*QC0*CO1;
            tC7=(CO3+tCb)*WP0+(tC3*WP0-CO1*CO5)*WQ0;
            tCc=WP0*WP0*WQ0;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC5*tC2)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC5*WQ1+tC6*tC2)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC6*WQ1+tC7*tC2)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC7*WQ1+tC2*tCc)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tCc*WQ1));
          };
          return tmpVal;
        } else if(lWork==kWork){//<xx|yy> ==>format 2   
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
            CO1=ijShellPair->inversezeta[i][j]; 
            CO2=klShellPair->inversezeta[k][l];           
            CO4=this->quartetConstants_->a0c0Par2[i][j][k][l];
            CO5=this->quartetConstants_->a000Par2[i][j][k][l];
            AB0=ijShellPair->deltaAB[iWork];
            PA0=ijShellPair->deltaPA[iWork][i][j];               
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];             
            CD1=klShellPair->deltaAB[kWork];             
            QC1=klShellPair->deltaPA[kWork][k][l];              
            WQ1=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];
            tC1=CO1+PA0*ijShellPair->deltaPB[iWork][i][j];
            tC2=CO2+QC1*(CD1+QC1);
            tC3=WP0*(AB0+math.two*PA0)-CO1*CO5;
            tC4=WQ1*(CD1+math.two*QC1)-CO2*CO4;
            tC5=WP0*WP0;
            tC6=WQ1*WQ1;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC1*tC2)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC1*tC4+tC2*tC3)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC3*tC4+tC1*tC6+tC2*tC5)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC3*tC6+tC4*tC5)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC5*tC6));
          };
          return tmpVal;
        } else {//<xx|yz> ==>format 4   
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
            CO1=ijShellPair->inversezeta[i][j];                   
            CO5=this->quartetConstants_->a000Par2[i][j][k][l];
            AB0=ijShellPair->deltaAB[iWork];             
            PA0=ijShellPair->deltaPA[iWork][i][j];               
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];              
            QC1=klShellPair->deltaPA[kWork][k][l];               
            WQ1=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];              
            CD2=klShellPair->deltaAB[lWork];              
            QC2=klShellPair->deltaPA[lWork][k][l];               
            WQ2=this->quartetConstants_->deltaWQ[lWork][i][j][k][l];
//          tC1=AB0+PA0;
            tC1=ijShellPair->deltaPB[iWork][i][j];
            tC2=CD2+QC2;
            tC3=tC1+PA0;
            tC4=PA0*tC1+CO1;
            tC5=QC1*tC2;
            tC6=WQ1*WQ2;
            tC7=WP0*WP0;
            tC8=WQ1*tC2+QC1*WQ2;
            tC9=WP0*tC3-CO1*CO5;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC5*tC4)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC8*tC4+tC9*tC5)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC5*tC7+tC9*tC8+tC4*tC6)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC8*tC7+tC9*tC6)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC6*tC7));
          };
          return tmpVal;
        };
      };
    } else {               //<xy|??>
      if(kWork==iWork){ //<xy|x?>
        if(lWork==iWork){//<xy|xx>  ===>format 3        
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){                                                                                                 
            CO1=klShellPair->inversezeta[k][l];      //==>CO2                                                   
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];                                                              
            CO4=this->quartetConstants_->a000Par2[i][j][k][l];     //==>CO5                                                 
            CO5=this->quartetConstants_->a0c0Par2[i][j][k][l];     //==>CO4                                                 
            AB0=klShellPair->deltaAB[kWork];           //==>CD                                              
            PA0=klShellPair->deltaPA[kWork][k][l];     //==>QC                                              
            QC0=ijShellPair->deltaPA[kWork][i][j];     //==>PA                                              
            WP0=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];      //==>WQ                                     
            WQ0=this->quartetConstants_->deltaWP[kWork][i][j][k][l];      //==>WP                                          
            CD1=ijShellPair->deltaAB[jWork];                 //==>AB                                                
            QC1=ijShellPair->deltaPA[jWork][i][j];           //==>PA                                                
            WQ1=this->quartetConstants_->deltaWP[jWork][i][j][k][l];      //==>WP                                   
            tC1=klShellPair->deltaPB[kWork][k][l];
            tC2=CD1+QC1;                                                                        
            tC3=tC1+PA0;                                                                        
            tC4=tC1*PA0;                                                                        
            tCa=tC4+CO1;
            tCb=CO3+QC0*WP0;
            tC5=tCa*QC0;                                                                        
            tC6=tCa*WQ0 + tCb*tC3 - CO5*QC0*CO1;                                                
            tC7=(CO3+tCb)*WP0+(tC3*WP0-CO1*CO5)*WQ0;                                            
            tCc=WP0*WP0*WQ0;                                                                    
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC5*tC2)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC5*WQ1+tC6*tC2)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC6*WQ1+tC7*tC2)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC7*WQ1+tC2*tCc)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tCc*WQ1));
          };
          return tmpVal;
        } else if(lWork==jWork){//<xy|xy> ==>format 5           
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
            AB1=ijShellPair->deltaAB[jWork];
            CD1=klShellPair->deltaAB[jWork];
            PA0=ijShellPair->deltaPA[iWork][i][j];
            QC0=klShellPair->deltaPA[iWork][k][l];                
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];      
            PA1=ijShellPair->deltaPA[jWork][i][j];
            QC1=klShellPair->deltaPA[jWork][k][l];
            WP1=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
            WQ1=this->quartetConstants_->deltaWQ[jWork][i][j][k][l];
//          tC1=AB1+PA1;
            tC1=ijShellPair->deltaPB[jWork][i][j];
            tC2=CD1+QC1;
            tC3=tC1*tC2;
            tC4=CO3+QC0*WP0+PA0*WQ0;
            tC5=CO3+tC1*WQ1+tC2*WP1;
            tC6=PA0*QC0;
            tC7=WP0*WQ0;
            tC8=WP1*WQ1;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC3*tC6)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC3*tC4+tC5*tC6)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC4*tC5+tC3*tC7+tC6*tC8)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC5*tC7+tC4*tC8)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC7*tC8));
          };
          return tmpVal;
        } else {//<xy|xz> ==>format 6   
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){           
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];           
            AB1=ijShellPair->deltaAB[jWork];
            CD2=klShellPair->deltaAB[lWork];
            PA0=ijShellPair->deltaPA[iWork][i][j];
            QC0=klShellPair->deltaPA[iWork][k][l];
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];             
            PA1=ijShellPair->deltaPA[jWork][i][j];
            QC2=klShellPair->deltaPA[lWork][k][l];
            WP1=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
            WQ2=this->quartetConstants_->deltaWQ[lWork][i][j][k][l];             
//          tC1=AB1+PA1;
            tC1=ijShellPair->deltaPB[jWork][i][j];
            tC2=CD2+QC2;
            tC3=tC1*tC2;
            tC4=PA0*QC0;
            tC5=WP1*WQ2;
            tC6=WP0*WQ0;
            tC7=PA0*WQ0+QC0*WP0+CO3;
            tC8=tC1*WQ2+tC2*WP1;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC3*tC4)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC3*tC7+tC4*tC8)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC4*tC5+tC3*tC6+tC7*tC8)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC5*tC7+tC6*tC8)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC5*tC6));
          };
          return tmpVal;
        };
      } else if(kWork==jWork){//<xy|y?>
        if(lWork==iWork){//<xy|yx>  ===>format 5        
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){           
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
            AB1=ijShellPair->deltaAB[jWork];
            CD1=-klShellPair->deltaAB[jWork];      //-CD
            PA0=ijShellPair->deltaPA[iWork][i][j];
            QC0=klShellPair->deltaPB[iWork][k][l]; //QD           
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];      
            PA1=ijShellPair->deltaPA[jWork][i][j]; 
            QC1=klShellPair->deltaPB[jWork][k][l]; //QD
            WP1=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
            WQ1=this->quartetConstants_->deltaWQ[jWork][i][j][k][l];
//          tC1=AB1+PA1;
            tC1=ijShellPair->deltaPB[jWork][i][j];
            tC2=CD1+QC1;
            tC3=tC1*tC2;
            tC4=CO3+QC0*WP0+PA0*WQ0;
            tC5=CO3+tC1*WQ1+tC2*WP1;
            tC6=PA0*QC0;
            tC7=WP0*WQ0;
            tC8=WP1*WQ1;          
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC3*tC6)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC3*tC4+tC5*tC6)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC4*tC5+tC3*tC7+tC6*tC8)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC5*tC7+tC4*tC8)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC7*tC8));
          };
          return tmpVal;
        } else if(lWork==jWork){//<xy|yy> ==>format 3   
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
            CO1=klShellPair->inversezeta[k][l];           
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];
            CO4=this->quartetConstants_->a000Par2[i][j][k][l];       
            CO5=this->quartetConstants_->a0c0Par2[i][j][k][l];     
            AB0=klShellPair->deltaAB[kWork];             
            PA0=klShellPair->deltaPA[kWork][k][l];       
            QC0=ijShellPair->deltaPB[kWork][i][j];       
            WP0=this->quartetConstants_->deltaWQ[kWork][i][j][k][l]; 
            WQ0=this->quartetConstants_->deltaWP[kWork][i][j][k][l];  
            CD1=-ijShellPair->deltaAB[iWork];            
            QC1=ijShellPair->deltaPB[iWork][i][j];       
            WQ1=this->quartetConstants_->deltaWP[iWork][i][j][k][l]; 
//          tC1=AB0+PA0;
            tC1=klShellPair->deltaPB[kWork][k][l];
            tC2=CD1+QC1;
            tC3=tC1+PA0;
            tC4=tC1*PA0;
            tCa=tC4+CO1;
            tCb=CO3+QC0*WP0;
            tC5=tCa*QC0;
            tC6=tCa*WQ0 + tCb*tC3 - CO5*QC0*CO1;
            tC7=(CO3+tCb)*WP0+(tC3*WP0-CO1*CO5)*WQ0;
            tCc=WP0*WP0*WQ0;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC5*tC2)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC5*WQ1+tC6*tC2)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC6*WQ1+tC7*tC2)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC7*WQ1+tC2*tCc)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tCc*WQ1));
          };
          return tmpVal;
        } else {//<xy|yz> ==>format 6    
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){           
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];           
            AB1=-ijShellPair->deltaAB[iWork];
            CD2=klShellPair->deltaAB[lWork];
            PA0=ijShellPair->deltaPB[jWork][i][j];
            QC0=klShellPair->deltaPA[jWork][k][l];
            WP0=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[jWork][i][j][k][l];             
            PA1=ijShellPair->deltaPB[iWork][i][j];
            QC2=klShellPair->deltaPA[lWork][k][l];
            WP1=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ2=this->quartetConstants_->deltaWQ[lWork][i][j][k][l];             
//          tC1=AB1+PA1;
            tC1=ijShellPair->deltaPA[iWork][i][j];
            tC2=CD2+QC2;
            tC3=tC1*tC2;
            tC4=PA0*QC0;
            tC5=WP1*WQ2;
            tC6=WP0*WQ0;
            tC7=PA0*WQ0+QC0*WP0+CO3;
            tC8=tC1*WQ2+tC2*WP1;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC3*tC4)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC3*tC7+tC4*tC8)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC4*tC5+tC3*tC6+tC7*tC8)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC5*tC7+tC6*tC8)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC5*tC6));
          };
          return tmpVal;
        };
      } else {//<xy|z?>
        if(lWork==iWork){//<xy|zx>  ===>format 6
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){           
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];           
            AB1=ijShellPair->deltaAB[jWork];
            CD2=-klShellPair->deltaAB[kWork];
            PA0=ijShellPair->deltaPA[iWork][i][j];
            QC0=klShellPair->deltaPB[iWork][k][l];
            WP0=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[iWork][i][j][k][l];             
            PA1=ijShellPair->deltaPA[jWork][i][j];
            QC2=klShellPair->deltaPB[kWork][k][l];
            WP1=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
            WQ2=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];             
//          tC1=AB1+PA1;
            tC1=ijShellPair->deltaPB[jWork][i][j];
            tC2=CD2+QC2;
            tC3=tC1*tC2;
            tC4=PA0*QC0;
            tC5=WP1*WQ2;
            tC6=WP0*WQ0;
            tC7=PA0*WQ0+QC0*WP0+CO3;
            tC8=tC1*WQ2+tC2*WP1;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC3*tC4)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC3*tC7+tC4*tC8)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC4*tC5+tC3*tC6+tC7*tC8)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC5*tC7+tC6*tC8)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC5*tC6));
          };
          return tmpVal;
        } else if(lWork==jWork){//<xy|zy> ==>format 6
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){           
            CO3=this->quartetConstants_->a0c0Par3[i][j][k][l];           
            AB1=-ijShellPair->deltaAB[iWork];
            CD2=-klShellPair->deltaAB[kWork];
            PA0=ijShellPair->deltaPB[jWork][i][j];
            QC0=klShellPair->deltaPB[jWork][k][l];
            WP0=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
            WQ0=this->quartetConstants_->deltaWQ[jWork][i][j][k][l];             
            PA1=ijShellPair->deltaPB[iWork][i][j];
            QC2=klShellPair->deltaPB[kWork][k][l];
            WP1=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            WQ2=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];             
//          tC1=AB1+PA1;
            tC1=ijShellPair->deltaPA[iWork][i][j];
            tC2=CD2+QC2;
            tC3=tC1*tC2;
            tC4=PA0*QC0;
            tC5=WP1*WQ2;
            tC6=WP0*WQ0;
            tC7=PA0*WQ0+QC0*WP0+CO3;
            tC8=tC1*WQ2+tC2*WP1;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC3*tC4)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC3*tC7+tC4*tC8)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC4*tC5+tC3*tC6+tC7*tC8)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC5*tC7+tC6*tC8)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC5*tC6));
          };return tmpVal;      
        } else {//<xy|zz> ==>Format 4
          for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++){
            CO1=klShellPair->inversezeta[k][l];
            CO5=this->quartetConstants_->a0c0Par2[i][j][k][l];
            AB0=klShellPair->deltaAB[kWork];
            PA0=klShellPair->deltaPA[kWork][k][l];
            WP0=this->quartetConstants_->deltaWQ[kWork][i][j][k][l];
            QC1=ijShellPair->deltaPA[iWork][i][j];
            WQ1=this->quartetConstants_->deltaWP[iWork][i][j][k][l];
            CD2=ijShellPair->deltaAB[jWork];
            QC2=ijShellPair->deltaPA[jWork][i][j];
            WQ2=this->quartetConstants_->deltaWP[jWork][i][j][k][l];
//          tC1=AB0+PA0;
            tC1=klShellPair->deltaPB[kWork][k][l];
            tC2=CD2+QC2;
            tC3=tC1+PA0;
            tC4=PA0*tC1+CO1;
            tC5=QC1*tC2;
            tC6=WQ1*WQ2;
            tC7=WP0*WP0;
            tC8=WQ1*tC2+QC1*WQ2;
            tC9=WP0*tC3-CO1*CO5;
            tmpVal +=   (this->quartetConstants_->FmT[0][i][j][k][l]*(tC5*tC4)
                        +this->quartetConstants_->FmT[1][i][j][k][l]*(tC8*tC4+tC9*tC5)
                        +this->quartetConstants_->FmT[2][i][j][k][l]*(tC5*tC7+tC9*tC8+tC4*tC6)
                        +this->quartetConstants_->FmT[3][i][j][k][l]*(tC8*tC7+tC9*tC6)
                        +this->quartetConstants_->FmT[4][i][j][k][l]*(tC6*tC7));
          };
          return tmpVal;
        };
      };
    };
  };
};

//---------------------------------------------------------------//
// two-e vertical recursion from [a0|c0] to [a0|00]              //
// [a0|c0]^m = (Q-C)*[a0|c-1,0]^m                                //
//           + (W-Q)*[a0|c-1,0]^(m+1)                            //
//           + N(a)/(2*(zeta+eta))*[a-1,0|c-1,0]^(m+1)           //
//           + (N(c)-1)/(2*eta)*[a0|c-2,0]^m                     //
//           - (N(c)-1)/(2*eta)*zeta/(zeta+eta)*[a0|c-2,0]^(m+1) //
//---------------------------------------------------------------//
double AOIntegrals::twoevRRa0c0(  ShellPair *ijShellPair, ShellPair *klShellPair,
                                int m, int LA, int *lA, int LC, int *lC, int *i, int *j, int *k, int *l) {
  if(LC==0) return this->twoevRRa000(ijShellPair,klShellPair,m,LA,lA,i,j,k,l);
  int iWork;
  int lAm1[3],lCm1[3];  
  for(iWork=0;iWork<3;iWork++){
    lAm1[iWork]=lA[iWork];     
    lCm1[iWork]=lC[iWork];
  };
  double tmpVal=0.0;
  if (lC[0]>0) iWork=0;
  else if (lC[1]>0) iWork=1;
  else if (lC[2]>0) iWork=2;
  if(LC>1) {
    lCm1[iWork]-=1;
    if(std::abs(klShellPair->deltaPA[iWork][*k][*l])>this->controls_->thresholdS) 
      tmpVal += klShellPair->deltaPA[iWork][*k][*l]*this->twoevRRa0c0(ijShellPair,klShellPair,m,LA,lA,LC-1,lCm1,i,j,k,l);
    if(std::abs(this->quartetConstants_->deltaWQ[iWork][*i][*j][*k][*l])>this->controls_->thresholdS) 
      tmpVal=this->quartetConstants_->deltaWQ[iWork][*i][*j][*k][*l]*this->twoevRRa0c0(ijShellPair,klShellPair,m+1,LA,lA,LC-1,lCm1,i,j,k,l);
    if (lA[iWork]>0) {
      lAm1[iWork] -= 1;
      tmpVal += (lAm1[iWork]+1)*this->quartetConstants_->a0c0Par3[*i][*j][*k][*l]*this->twoevRRa0c0(ijShellPair,klShellPair,m+1,LA-1,lAm1,LC-1,lCm1,i,j,k,l);
    };
    if(lC[iWork]==2&&LC==2) 
      tmpVal += klShellPair->inversezeta[*k][*l]*(this->twoevRRa000(ijShellPair,klShellPair,m,LA,lA,i,j,k,l)
               -this->quartetConstants_->a0c0Par2[*i][*j][*k][*l]*this->twoevRRa000(ijShellPair,klShellPair,m+1,LA,lA,i,j,k,l));
    else if (lC[iWork]>=2){
      lCm1[iWork] -=1; 
      tmpVal += (lCm1[iWork]+1)*klShellPair->inversezeta[*k][*l]*(this->twoevRRa0c0(ijShellPair,klShellPair,m,LA,lA,LC-2,lCm1,i,j,k,l)
               -this->quartetConstants_->a0c0Par2[*i][*j][*k][*l]*this->twoevRRa0c0(ijShellPair,klShellPair,m+1,LA,lA,LC-2,lCm1,i,j,k,l));
    };
  } else {
    if(std::abs(this->quartetConstants_->deltaWQ[iWork][*i][*j][*k][*l])>this->controls_->thresholdS) 
      tmpVal = this->quartetConstants_->deltaWQ[iWork][*i][*j][*k][*l]*this->twoevRRa000(ijShellPair,klShellPair,m+1,LA,lA,i,j,k,l);
    if(std::abs(klShellPair->deltaPA[iWork][*k][*l])>this->controls_->thresholdS) 
      tmpVal+= klShellPair->deltaPA[iWork][*k][*l]*this->twoevRRa000(ijShellPair,klShellPair,m,LA,lA,i,j,k,l);
    if (lA[iWork]>0) {
      lAm1[iWork]-=1;
      tmpVal += (lAm1[iWork]+1)*this->quartetConstants_->a0c0Par3[*i][*j][*k][*l]*this->twoevRRa000(ijShellPair,klShellPair,m+1,LA-1,lAm1,i,j,k,l);
    };
  };
  return tmpVal;
};
//---------------------------------------------------------------//
// two-e vertical recursion from [a0|c0] to [00|00]              //
// [a0|00]^m = (P-A)*[a-1,0|00]^m                                //
//           + (W-P)*[a-1,0|00]^(m+1)                            //
//           + (N(a)-1)/(2*zeta)*[a-2,0|00]^m                    //
//           - (N(a)-1)/(2*zeta)*eta/(zeta+eta)*[a-2,0|00]^(m+1) //
//---------------------------------------------------------------//
double AOIntegrals::twoevRRa000(ShellPair *ijShellPair,ShellPair *klShellPair,int m,int LA,int *lA,int *i,int *j,int *k,int *l) {
  if(LA==0) return this->quartetConstants_->FmT[m][*i][*j][*k][*l];
  double tmpVal=0.0;
  int iWork;
  int lAm1[3];
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0) iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;
  if(LA>1) {
    lAm1[iWork]-=1;
    if(std::abs(this->quartetConstants_->deltaWP[iWork][*i][*j][*k][*l])>this->controls_->thresholdS) 
      tmpVal = this->quartetConstants_->deltaWP[iWork][*i][*j][*k][*l]*this->twoevRRa000(ijShellPair,klShellPair,m+1,LA-1,lAm1,i,j,k,l);
    if(std::abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) 
      tmpVal+= ijShellPair->deltaPA[iWork][*i][*j]*this->twoevRRa000(ijShellPair,klShellPair,m,LA-1,lAm1,i,j,k,l);
    if(LA==2&&lA[iWork]==2) 
      tmpVal += ijShellPair->inversezeta[*i][*j]*(this->quartetConstants_->FmT[m][*i][*j][*k][*l] 
               -this->quartetConstants_->a000Par2[*i][*j][*k][*l]*this->quartetConstants_->FmT[m+1][*i][*j][*k][*l]);
    else if (lA[iWork]>=2) {
      lAm1[iWork] -=1;
      tmpVal += (lAm1[iWork]+1)*ijShellPair->inversezeta[*i][*j]*(this->twoevRRa000(ijShellPair,klShellPair,m,LA-2,lAm1,i,j,k,l) 
               -this->quartetConstants_->a000Par2[*i][*j][*k][*l]*this->twoevRRa000(ijShellPair,klShellPair,m+1,LA-2,lAm1,i,j,k,l));
    };
  } else {
    if(std::abs(this->quartetConstants_->deltaWP[iWork][*i][*j][*k][*l])>this->controls_->thresholdS) 
      tmpVal = this->quartetConstants_->deltaWP[iWork][*i][*j][*k][*l]*this->quartetConstants_->FmT[m+1][*i][*j][*k][*l];
    if(std::abs(ijShellPair->deltaPA[iWork][*i][*j])>this->controls_->thresholdS) 
      tmpVal+= ijShellPair->deltaPA[iWork][*i][*j]*this->quartetConstants_->FmT[m][*i][*j][*k][*l];
  };
  return tmpVal;
};

double AOIntegrals::twoeSSSS0(int *nPGTOs, ShellPair *ijShellPair, ShellPair *klShellPair){
  // integrate (SS|SS) 
  int i,j,k,l,m;
  double PQ,sqrPQ,SSSS0=0.0,FmT[1],T,Upq,expo1,expo2,expoT;

  for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
    Upq = (ijShellPair->UAB[i][j])*(klShellPair->UAB[k][l]);
    if(std::abs(Upq)>this->controls_->thresholdS) {
      sqrPQ = math.zero;
      for(m=0;m<3;m++) {
        PQ = (ijShellPair->centerP[m][i][j]-klShellPair->centerP[m][k][l]);
        sqrPQ += PQ*PQ;
      };
      expoT = ijShellPair->zeta[i][j]+klShellPair->zeta[k][l];
      if(sqrPQ>this->controls_->thresholdS) {
        T = sqrPQ/(ijShellPair->invzeta[i][j]+klShellPair->invzeta[k][l]);
        this->computeFmTTaylor(FmT,T,0,0);
        SSSS0 += Upq*FmT[0]/sqrt(expoT);
      }else{
        SSSS0 += Upq/sqrt(expoT);
      };
    };
  };
  return SSSS0;
};
*/
