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
#include <aointegrals.h>
using ChronusQ::AOIntegrals;
static double smallT[21]={
1.00000000000000,
0.33333333333333,
0.20000000000000,
0.14285714285714,
0.11111111111111,
0.09090909090909,
0.07692307692308,
0.06666666666667,
0.05882352941176,
0.05263157894737,
0.04761904761905,
0.04347826086957,
0.04000000000000,
0.03703703703704,
0.03448275862069,
0.03225806451613,
0.03030303030303,
0.02857142857143,
0.02702702702703,
0.02564102564103,
0.02439024390244
};

static double factTLarge[21] = {
0.88622692545275800,
0.44311346272637900,
0.66467019408956900,
1.66167548522392000,
5.81586419828372000,
26.1713888922768000,
143.942638907522000,
935.627152898894000,
7017.20364674171000,
59646.2309973045000,
566639.194474393000,
5949711.54198112000,
68421682.7327829000,
855271034.159787000,
11546158961.1571000,
167419304936.778000,
2594999226520.06000,
42817487237581.0000,
749306026657668.000,
13862161493166900.0,
270312149116754000.0
};

//---------------------
// initialize AOIntegrals
//---------------------
void AOIntegrals::iniAOIntegrals(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basisset, 
                                 std::shared_ptr<FileIO> fileio, std::shared_ptr<Controls> controls){
  this->molecule_ = molecule;
  this->basisSet_ = basisset;
  this->fileio_   = fileio;
  this->controls_ = controls;
  this->nBasis_   = basisset->nBasis();
  this->nTT_      = this->nBasis_*(this->nBasis_+1)/2;

  // FIXME need a try statement for alloc
#ifndef USE_LIBINT // We don't need to allocate these if we're using Libint
  this->twoEC_ = std::make_shared<RealMatrix>(this->nTT_,this->nTT_); // Raffenetti Two Electron Coulomb AOIntegrals
  this->twoEX_ = std::make_shared<RealMatrix>(this->nTT_,this->nTT_); // Raffenetti Two Electron Exchange AOIntegrals
#else // Allocate space for all N^4 AO Integrals in BTAS Tensor object (TODO need to set this up to be conditional)
  if(this->controls_->buildn4eri) 
    this->aoERI_ = std::make_shared<Tensor<double>>(this->nBasis_,this->nBasis_,this->nBasis_,this->nBasis_);
#endif
  this->oneE_      = std::make_shared<RealMatrix>(this->nBasis_,this->nBasis_); // One Electron Integral
  this->overlap_   = std::make_shared<RealMatrix>(this->nBasis_,this->nBasis_); // Overlap
  this->kinetic_   = std::make_shared<RealMatrix>(this->nBasis_,this->nBasis_); // Kinetic
  this->potential_ = std::make_shared<RealMatrix>(this->nBasis_,this->nBasis_); // Potential
#ifdef USE_LIBINT
  this->schwartz_ = std::make_shared<RealMatrix>(this->basisSet_->nShell(),this->basisSet_->nShell()); // Schwartz
#endif

  int i,j,ij;
  this->R2Index_ = new int*[this->nBasis_];
  for(i=0;i<this->nBasis_;i++) this->R2Index_[i] = new int[this->nBasis_];
  for(i=0;i<this->nBasis_;i++) for(j=0;j<this->nBasis_;j++) {
    if(i>=j) ij=j*(this->nBasis_)-j*(j-1)/2+i-j;
    else ij=i*(this->nBasis_)-i*(i-1)/2+j-i;
    this->R2Index_[i][j] = ij;
  };

  pairConstants_ = new PairConstants;
  molecularConstants_ = new MolecularConstants;
  quartetConstants_ = new QuartetConstants;

  this->haveAOTwoE = false;
  this->haveAOOneE = false;
#ifdef USE_LIBINT
  this->haveSchwartz = false;
#endif

// initialize the FmT table
// Need to know the max L first
  this->FmTTable_ = new double*[MaxFmTPt];
  for(i=0;i<MaxFmTPt;i++) this->FmTTable_[i] = new double[MaxTotalL];
  this->generateFmTTable();
};

void AOIntegrals::generateFmTTable() {
  double intervalFmT = 0.025;
  double T = 0.0;
  int i,m;
  double critT = 33.0;
  double expT, factor, term, sum, twoT, Tn;
  for(i=0;i<MaxFmTPt;i++){
    if(std::abs(T)<= math.small) {
      for(m=0;m<=MaxTotalL;m++) this->FmTTable_[i][m]=math.one/(math.two*m+1);
    } else if(T>critT) {
      this->FmTTable_[i][0] = math.half*sqrt(math.pi/T);
      twoT = math.two*T;
      Tn=math.one;
      for(m=1;m<=MaxTotalL;m++) {
        Tn*=twoT;
        this->FmTTable_[i][m] = this->FmTTable_[i][0]*dFactorial[m]/Tn;
      };
    } else {
      expT   = exp(-T);
      factor = MaxTotalL + math.half;
      term   = math.half/factor;
      sum    = term;
      while(term>math.small) {
        factor += math.one;
        term   *= T/factor;
        sum    += term;
      };
      this->FmTTable_[i][MaxTotalL] = expT*sum;
      twoT = math.two*T;
      for(m=MaxTotalL-1;m>=0;m--) this->FmTTable_[i][m] = (twoT*this->FmTTable_[i][m+1]+expT)/(2*m+1);
    };
    T += intervalFmT;
  };
};

void AOIntegrals::computeFmTTaylor(double *FmT, double T, int maxM, int minM){
  int m,i;
  double intervalFmT = 0.025;
  double critT = 33.0;
  if(T>critT) {
    double Tn=sqrt(T);
    FmT[0] = factTLarge[0]/Tn;
    for(m=minM+1;m<=maxM;m++) {
      Tn*=T;
      FmT[m] = factTLarge[m]/Tn;
    };
  } else {
    int j=T/intervalFmT;
    double deltaT = j*intervalFmT-T;
    double sum = this->FmTTable_[j][maxM];
    double tmp = math.one;
    double change;
    for(i=1;i<5;i++){
      tmp   *= deltaT/i;
      change = tmp*this->FmTTable_[j][maxM+i];
      sum   += change;
    };
//  down-recursion to obtain m<maxM values
    FmT[maxM]=sum;
    if(minM<maxM) {
      double twoT = math.two*T;
      double expT = exp(-T);
      for(m=maxM-1;m>=minM;m--) FmT[m] = (twoT*FmT[m+1]+expT)*smallT[m];
    };
  };
};

void AOIntegrals::iniQuartetConstants(ShellPair *ijShellPair, ShellPair *klShellPair){
  /* compute four-center info for shell quartet with L >= 1*/
  int i,j,k,l,m,nFmT=0;
  int nPGTOs[4];
  double centerQuartet[3][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];
  nPGTOs[0] = ijShellPair->nPGTOs[0];
  nPGTOs[1] = ijShellPair->nPGTOs[1];
  nPGTOs[2] = klShellPair->nPGTOs[0];
  nPGTOs[3] = klShellPair->nPGTOs[1];
  int totalL = (*ijShellPair).LTotal+(*klShellPair).LTotal;
  double PQ,sqrPQ=0.0;
  double T,Upq,normijkl;
  double *FmT = new double[totalL+1];
  double expo1,expo2,expoT;
  for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
    Upq = ((*ijShellPair).UAB[i][j])*((*klShellPair).UAB[k][l]);
//  Upq = ((*ijShellPair).KAB[i][j])*((*klShellPair).KAB[k][l]);
    if(std::abs(Upq)>this->controls_->thresholdS) {
      expo1=(*ijShellPair).zeta[i][j];
      expo2=(*klShellPair).zeta[k][l];
      expoT=expo1+expo2;
      for(m=0;m<3;m++) {
        centerQuartet[m][i][j][k][l]=((*ijShellPair).centerPZeta[m][i][j]+(*klShellPair).centerPZeta[m][k][l])/expoT;
        this->quartetConstants_->deltaWP[m][i][j][k][l] = centerQuartet[m][i][j][k][l] - (*ijShellPair).centerP[m][i][j];
        this->quartetConstants_->deltaWQ[m][i][j][k][l] = centerQuartet[m][i][j][k][l] - (*klShellPair).centerP[m][k][l];
      };
      if(totalL>1) {
        this->quartetConstants_->a0c0Par2[i][j][k][l] = expo1/expoT;
        this->quartetConstants_->a000Par2[i][j][k][l] = expo2/expoT;
        this->quartetConstants_->a0c0Par3[i][j][k][l] = math.half/expoT;
      };

      sqrPQ = math.zero;
      for(m=0;m<3;m++) {
        PQ = ((*ijShellPair).centerP[m][i][j]-(*klShellPair).centerP[m][k][l]);
        sqrPQ += PQ*PQ;
      };
      Upq = Upq/sqrt(expoT);
      if(sqrPQ>this->controls_->thresholdAB) {
        T = sqrPQ/((*ijShellPair).invzeta[i][j]+(*klShellPair).invzeta[k][l]);
        this->computeFmTTaylor(FmT,T,totalL,0);
        for(m=0;m<=totalL;m++) this->quartetConstants_->FmT[m][i][j][k][l] = Upq*FmT[m];
      }else{
        for(m=0;m<=totalL;m++) this->quartetConstants_->FmT[m][i][j][k][l] = Upq*smallT[m];
      };
    } else {
      for(m=0;m<=totalL;m++) this->quartetConstants_->FmT[m][i][j][k][l] = math.zero;
    };
  };
//  for(i=0;i<nPGTOs[0];i++) for(j=0;j<nPGTOs[1];j++) for(k=0;k<nPGTOs[2];k++) for(l=0;l<nPGTOs[3];l++) {
//    normijkl = ijShellPair->norm[i][j]*klShellPair->norm[k][l];
//    for(m=0;m<=totalL;m++) this->quartetConstants_->FmT[m][i][j][k][l] *= normijkl; 
//  };
};
//-----------------------------//
// initialize pair constants   //
//-----------------------------//
void AOIntegrals::iniPairConstants(ShellPair *ijShellPair){
  int i,j,k;
  double expo[MAXCONTRACTION][2];
  double center[3][2];
  this->pairConstants_->intSmall = controls_->thresholdS;
  /*-------------------------*/
  /* compute one-center info */
  /*-------------------------*/
  int totalL = ijShellPair->LTotal;
  for(i=0;i<2;i++) {
    this->pairConstants_->nPGTOs[i]=ijShellPair->nPGTOs[i];
    this->pairConstants_->aoIndex[i]=ijShellPair->aoIndex[i];
    this->pairConstants_->L[i]=ijShellPair->L[i];
    this->pairConstants_->nBasis[i]=ijShellPair->nBasis[i];
  };
  for(j=0;j<2;j++) {
    for(k=0;k<3;k++) center[k][j]=(*(molecule_->cart()))(k,((*ijShellPair).center[j]));
    for(i=0;i<this->pairConstants_->nPGTOs[j];i++){
      expo[i][j]=(basisSet_->shells[(*ijShellPair).shIndex[j]]).expo[i];
    };
  };
  /*-------------------------*/
  /* compute two-center info */
  /*-------------------------*/
  double expo1,expo2,expoT,constAtom;
  double sqrAB=0.0, sqrPZ;
  double *FmT = new double[totalL+1];
  int m,iAtom;
  for(k=0;k<3;k++) {
    this->pairConstants_->deltaAB[k] = ijShellPair->deltaAB[k];
    sqrAB += (this->pairConstants_->deltaAB[k])*(this->pairConstants_->deltaAB[k]);
  };
  this->pairConstants_->ssPairTotal = 0.0;
  for(i=0;i<this->pairConstants_->nPGTOs[0];i++) for(j=0;j<this->pairConstants_->nPGTOs[1];j++){
    expo1=expo[i][0];
    expo2=expo[j][1];
    expoT=ijShellPair->zeta[i][j];
    for(k=0;k<3;k++) {
      this->pairConstants_->deltaPA[k][i][j] = ijShellPair->deltaPA[k][i][j];
      this->pairConstants_->deltaPB[k][i][j] = ijShellPair->deltaPB[k][i][j];
    };
    this->pairConstants_->Sa0Par[i][j] = ijShellPair->inversezeta[i][j];
    this->pairConstants_->TabPar1[i][j] = math.three*expo1*expo2/expoT - math.two*(expo1*expo2/expoT)*(expo1*expo2/expoT)*sqrAB;
    this->pairConstants_->TabPar2[i][j] = math.two*expo1*expo2/expoT;
    this->pairConstants_->TabPar3[i][j] = this->pairConstants_->TabPar2[i][j]/(math.two*expo2);
    this->pairConstants_->Ta0Par3[i][j] = this->pairConstants_->TabPar2[i][j]/(math.two*expo1);
    this->pairConstants_->ssPair[i][j]  = ijShellPair->norm[i][j]*sqrt((math.pi/expoT)*(math.pi/expoT)*(math.pi/expoT))*exp(-expo1*expo2*sqrAB/expoT);
    this->pairConstants_->ssPairTotal += this->pairConstants_->ssPair[i][j];
    if(std::abs(this->pairConstants_->ssPair[i][j])>this->pairConstants_->intSmall) this->pairConstants_->ssNonzero[i][j] = true;
    else this->pairConstants_->ssNonzero[i][j] = false;
    /*******************************/
    /* compute FmU[m][i][j][iAtom] */
    /*******************************/
    if(this->pairConstants_->ssNonzero[i][j]) {
      constAtom = math.two*sqrt(expoT/math.pi)*this->pairConstants_->ssPair[i][j];
      for(iAtom=0;iAtom<molecularConstants_->nAtom;iAtom++){
	sqrPZ=0.0;
	for(k=0;k<3;k++) {
	  this->pairConstants_->deltaPZ[k][i][j][iAtom] = ijShellPair->centerP[k][i][j]-molecularConstants_->cart[k][iAtom];
	  sqrPZ += (this->pairConstants_->deltaPZ[k][i][j][iAtom])*(this->pairConstants_->deltaPZ[k][i][j][iAtom]);
	};
	this->computeFmTTaylor(FmT,expoT*sqrPZ,totalL,0);
	for(m=0;m<=totalL;m++) this->pairConstants_->FmU[m][i][j][iAtom]=constAtom*FmT[m];
      };
    };
  };
};

void AOIntegrals::iniMolecularConstants(){
  this->molecularConstants_->nAtom=molecule_->nAtoms();
  for(int i=0;i<this->molecularConstants_->nAtom;i++){
    this->molecularConstants_->atomZ[i]=atom[molecule_->index(i)].atomicNumber;
    for(int j=0;j<3;j++) this->molecularConstants_->cart[j][i]=(*(molecule_->cart()))(j,i);
  };
};

void AOIntegrals::printTimings() {
    this->fileio_->out << endl << "Timing Statistics: "<<endl << bannerTop << endl;
    this->fileio_->out << endl << "One Electron Integral Timings" << endl << bannerMid << endl;
    this->fileio_->out << std::left << std::setw(60) << "Wall time for Overlap evaluation:" 
                       << std::left << std::setw(15) << this->SED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << "Wall time for Kinetic evaluation:" 
                       << std::left << std::setw(15) << this->TED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << "Wall time for Nuclear Attraction Potential evaluation:" 
                       << std::left << std::setw(15) << this->VED.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << " "
                       << std::left << std::setw(15) << "---------------" << "----" << endl;
    this->fileio_->out << std::left << std::setw(60) << "Total wall time for one-electron integral evaluation:" 
                       << std::left << std::setw(15) << this->OneED.count() << " sec" << endl;
    this->fileio_->out << endl << endl;
    this->fileio_->out << "Two Electron Integral Timings" << endl << bannerMid << endl;
    this->fileio_->out << std::left << std::setw(60) << "Wall time for Schwartz Bound evaluation:" 
                       << std::left << std::setw(15) << this->SchwartzD.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << "Wall time for Density Shell Block Norm evaluation:" 
                       << std::left << std::setw(15) << this->DenShBlkD.count() << " sec" << endl;
    this->fileio_->out << std::left << std::setw(60) << "Wall time for Perturbation Tensor evaluation:" 
                       << std::left << std::setw(15) << this->PTD.count() << " sec" << endl;
      
    this->fileio_->out << bannerEnd << endl;
}
