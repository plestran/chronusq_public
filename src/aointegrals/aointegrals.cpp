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

void AOIntegrals::generateFmTTable() {
  double intervalFmT = 0.025;
  double T = 0.0;
  int i,m;
  double critT = 33.0;
  double expT, factor, term, sum, twoT, Tn;
  for(i = 0; i < MaxFmTPt; i++){
    if(std::abs(T) <= math.small) {
      for(m = 0; m <= MaxTotalL; m++) 
        this->FmTTable_[i][m] = math.one / (math.two * m + 1);
    } else if(T > critT) {
      this->FmTTable_[i][0] = math.half * sqrt(math.pi / T);
      twoT = math.two * T;
      Tn = math.one;
      for(m = 1; m <= MaxTotalL; m++) {
        Tn *= twoT;
        this->FmTTable_[i][m] = this->FmTTable_[i][0] * dFactorial[m] / Tn;
      };
    } else {
      expT   = exp(-T);
      factor = MaxTotalL + math.half;
      term   = math.half / factor;
      sum    = term;
      while(term > math.small) {
        factor += math.one;
        term   *= T / factor;
        sum    += term;
      };
      this->FmTTable_[i][MaxTotalL] = expT * sum;
      twoT = math.two * T;
      for(m = MaxTotalL - 1; m >= 0; m--) 
        this->FmTTable_[i][m] = (twoT * this->FmTTable_[i][m+1] + expT)/(2*m + 1);
    };
    T += intervalFmT;
  };
};

void AOIntegrals::computeFmTTaylor(double *FmT, double T, int maxM, int minM){
  int m,i;
  double intervalFmT = 0.025;
  double critT = 33.0;
  if(T > critT) {
    double Tn = sqrt(T);
    FmT[0] = factTLarge[0] / Tn;
    for(m = minM + 1; m <= maxM; m++) {
      Tn    *= T;
      FmT[m] = factTLarge[m]/Tn;
    };
  } else {
    int    j      = T / intervalFmT;
    double deltaT = j * intervalFmT - T;
    double sum = this->FmTTable_[j][maxM];
    double tmp = math.one;
    double change;
    for(i = 1; i < 5; i++){
      tmp   *= deltaT / i;
      change = tmp * this->FmTTable_[j][maxM + i];
      sum   += change;
    };
//  down-recursion to obtain m<maxM values
    FmT[maxM] = sum;
    if(minM < maxM) {
      double twoT = math.two * T;
      double expT = exp(-T);
      for(m = maxM - 1; m >= minM; m--) FmT[m] = (twoT * FmT[m+1] + expT) * smallT[m];
    };
  };
};
/*
void AOIntegrals::iniQuartetConstants(ShellPair *ijShellPair, ShellPair *klShellPair){
  // compute four-center info for shell quartet with L >= 1
  int i,j,k,l,m,nFmT=0;
  int nPGTOs[4];
  double centerQuartet[3][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION][MAXCONTRACTION];
  nPGTOs[0] = ijShellPair->nPGTOs[0];
  nPGTOs[1] = ijShellPair->nPGTOs[1];
  nPGTOs[2] = klShellPair->nPGTOs[0];
  nPGTOs[3] = klShellPair->nPGTOs[1];
  int totalL = ijShellPair->LTotal + klShellPair->LTotal;
  double PQ,sqrPQ=0.0;
  double T,Upq,normijkl;
  double *FmT = new double[totalL+1];
  double expo1,expo2,expoT;
  for(i = 0; i < nPGTOs[0]; i++) 
  for(j = 0; j < nPGTOs[1]; j++) 
  for(k = 0; k < nPGTOs[2]; k++) 
  for(l = 0; l < nPGTOs[3]; l++) {
    Upq = ijShellPair->UAB[i][j] * klShellPair->UAB[k][l];
    if(std::abs(Upq) > this->thresholdS_) {
      expo1 = ijShellPair->zeta[i][j];
      expo2 = klShellPair->zeta[k][l];
      expoT = expo1 + expo2;
      for(m = 0; m < 3; m++) {
        centerQuartet[m][i][j][k][l]=
          (ijShellPair->centerPZeta[m][i][j] + klShellPair->centerPZeta[m][k][l]) / expoT;

        this->quartetConstants_->deltaWP[m][i][j][k][l] = 
          centerQuartet[m][i][j][k][l] - ijShellPair->centerP[m][i][j];
        this->quartetConstants_->deltaWQ[m][i][j][k][l] = 
          centerQuartet[m][i][j][k][l] - klShellPair->centerP[m][k][l];
      };
      if(totalL > 1) {
        this->quartetConstants_->a0c0Par2[i][j][k][l] = expo1 / expoT;
        this->quartetConstants_->a000Par2[i][j][k][l] = expo2 / expoT;
        this->quartetConstants_->a0c0Par3[i][j][k][l] = math.half / expoT;
      };

      sqrPQ = math.zero;
      for(m = 0; m < 3; m++) {
        PQ = ijShellPair->centerP[m][i][j] - klShellPair->centerP[m][k][l];
        sqrPQ += PQ*PQ;
      };
      Upq = Upq / sqrt(expoT);
      if(sqrPQ > this->thresholdAB_) {

        T = sqrPQ / (ijShellPair->invzeta[i][j] + klShellPair->invzeta[k][l]);
        this->computeFmTTaylor(FmT,T,totalL,0);

        for(m = 0; m<= totalL; m++) 
          this->quartetConstants_->FmT[m][i][j][k][l] = Upq * FmT[m];

      } else {
        for(m = 0; m <= totalL; m++) 
          this->quartetConstants_->FmT[m][i][j][k][l] = Upq * smallT[m];
      };
    } else {
      for(m = 0; m <= totalL; m++) 
        this->quartetConstants_->FmT[m][i][j][k][l] = math.zero;
    };
  };
};
//-----------------------------//
// initialize pair constants   //
//-----------------------------//
void AOIntegrals::iniPairConstants(ShellPair *ijShellPair){
  int i,j,k;
  double expo[MAXCONTRACTION][2];
  double center[3][2];
  this->pairConstants_->intSmall = this->thresholdS_;
  // compute one-center info 
  int totalL = ijShellPair->LTotal;
  for(i = 0; i < 2; i++) {
    this->pairConstants_->nPGTOs[i]  = ijShellPair->nPGTOs[i];
    this->pairConstants_->aoIndex[i] = ijShellPair->aoIndex[i];
    this->pairConstants_->L[i]       = ijShellPair->L[i];
    this->pairConstants_->nBasis[i]  = ijShellPair->nBasis[i];
  };
  for(j = 0; j < 2; j++) {
    for(k = 0; k < 3; k++)  
      center[k][j] = (*(molecule_->cart()))(k,(ijShellPair->center[j]));

    for(i = 0; i < this->pairConstants_->nPGTOs[j];i++)
      expo[i][j] = (basisSet_->shells_old[ijShellPair->shIndex[j]]).expo[i];
  };
  // compute two-center info
  double expo1,expo2,expoT,constAtom;
  double sqrAB = 0.0, sqrPZ;
  double *FmT = new double[totalL + 1];
  int m,iAtom;
  for(k = 0; k < 3; k++) {
    this->pairConstants_->deltaAB[k] = ijShellPair->deltaAB[k];
    sqrAB += this->pairConstants_->deltaAB[k]*this->pairConstants_->deltaAB[k];
  };
  this->pairConstants_->ssPairTotal = 0.0;
  for(i = 0; i < this->pairConstants_->nPGTOs[0]; i++) 
  for(j = 0; j < this->pairConstants_->nPGTOs[1]; j++){
    expo1 = expo[i][0];
    expo2 = expo[j][1];
    expoT = ijShellPair->zeta[i][j];
    for(k = 0; k < 3; k++) {
      this->pairConstants_->deltaPA[k][i][j] = ijShellPair->deltaPA[k][i][j];
      this->pairConstants_->deltaPB[k][i][j] = ijShellPair->deltaPB[k][i][j];
    };
    this->pairConstants_->Sa0Par[i][j]  = ijShellPair->inversezeta[i][j];
    this->pairConstants_->TabPar1[i][j] = 
      math.three * expo1 * expo2/expoT - 
      math.two   * (expo1 * expo2 / expoT) * (expo1*expo2/expoT) * sqrAB;
    this->pairConstants_->TabPar2[i][j] = math.two * expo1 * expo2 / expoT;
    this->pairConstants_->TabPar3[i][j] = 
      this->pairConstants_->TabPar2[i][j] / (math.two * expo2);
    this->pairConstants_->Ta0Par3[i][j] = 
      this->pairConstants_->TabPar2[i][j] / (math.two * expo1);
    this->pairConstants_->ssPair[i][j]  = 
      ijShellPair->norm[i][j] * 
      sqrt((math.pi/expoT) * (math.pi/expoT) * (math.pi/expoT)) * 
      exp(- expo1 * expo2 * sqrAB / expoT);
    this->pairConstants_->ssPairTotal += this->pairConstants_->ssPair[i][j];
    if(std::abs(this->pairConstants_->ssPair[i][j]) > this->pairConstants_->intSmall) 
      this->pairConstants_->ssNonzero[i][j] = true;
    else this->pairConstants_->ssNonzero[i][j] = false;
    // compute FmU[m][i][j][iAtom]
    if(this->pairConstants_->ssNonzero[i][j]) {
      constAtom = math.two * sqrt(expoT/math.pi) * this->pairConstants_->ssPair[i][j];
      for(iAtom = 0; iAtom < molecularConstants_->nAtom; iAtom++){
	sqrPZ = 0.0;
	for(k = 0; k < 3; k++) {
	  this->pairConstants_->deltaPZ[k][i][j][iAtom] = 
            ijShellPair->centerP[k][i][j]-molecularConstants_->cart[k][iAtom];
	  sqrPZ += 
            (this->pairConstants_->deltaPZ[k][i][j][iAtom]) * 
            (this->pairConstants_->deltaPZ[k][i][j][iAtom]);
	};
	this->computeFmTTaylor(FmT,expoT*sqrPZ,totalL,0);
	for(m = 0; m <= totalL; m++) 
          this->pairConstants_->FmU[m][i][j][iAtom]=constAtom*FmT[m];
      };
    };
  };
};
*/
void AOIntegrals::iniMolecularConstants(){
  this->molecularConstants_->nAtom=molecule_->nAtoms();
  for(int i = 0; i < this->molecularConstants_->nAtom ; i++){
    this->molecularConstants_->atomZ[i] = atom[molecule_->index(i)].atomicNumber;
    for(int j = 0; j < 3; j++) 
      this->molecularConstants_->cart[j][i]=(*(molecule_->cart()))(j,i);
  };
};

void AOIntegrals::printTimings() {
    this->fileio_->out << endl << "Integral Timings: "<<endl << BannerTop 
      << endl;
    this->fileio_->out << endl << "One Electron Integral Timings" << endl 
                       << bannerMid << endl;

/*
    if(this->maxMultipole_ >= 4) {
      this->fileio_->out << std::left << std::setw(60) 
            << "Wall time for Overlap + Dipole + Quadrupole + Octupole + RxP"; 
    } else if(this->maxMultipole_ >= 1) {
      this->fileio_->out << std::left << std::setw(60) 
                         << "Wall time for Overlap + Dipole evaluation:"; 
    } else {
      this->fileio_->out << std::left << std::setw(60) 
                         << "Wall time for Overlap evaluation:"; 
    }
*/

    this->fileio_->out << std::left << std::setw(60) <<
      "Wall time for Overlap + Multipole:";
    this->fileio_->out << std::left << std::setw(15) 
      << this->SED.count() << " sec" << endl;

    this->fileio_->out << std::left << std::setw(60) 
      << "Wall time for Kinetic evaluation:" 
      << std::left << std::setw(15) << this->TED.count() << " sec" << endl;

    this->fileio_->out << std::left << std::setw(60) 
      << "Wall time for Nuclear Attraction Potential evaluation:" 
      << std::left << std::setw(15) << this->VED.count() << " sec" << endl;

    if(this->orthoType == LOWDIN)
      this->fileio_->out << std::left << std::setw(60) 
        << "Wall time for formation of Lowdin transformation:";
    if(this->orthoType == CHOLESKY)
      this->fileio_->out << std::left << std::setw(60) 
        << "Wall time for formation of Cholesky transformation:";

    this->fileio_->out
      << std::left << std::setw(15) << this->OrthoD.count() << " sec" << endl;

    if(this->doX2C)
      this->fileio_->out << std::left << std::setw(60) 
        << "Wall time for formation of X2C transformation:"
        << std::left << std::setw(15) << this->X2CD.count() << " sec" << endl;

    this->fileio_->out << std::left << std::setw(60) 
      << "Wall time for Schwartz Bound evaluation:" 
      << std::left << std::setw(15) << this->SchwartzD.count() << " sec" 
      << endl;
        



    this->fileio_->out << std::left << std::setw(60) << " "
      << std::left << std::setw(15) << "---------------" 
      << "----" << endl;

    this->fileio_->out << std::left << std::setw(60) 
      << "Total wall time for one-electron integral evaluation:" 
      << std::left << std::setw(15) << this->OneED.count() << " sec" 
      << endl;


/*
    this->fileio_->out << endl << endl;
    this->fileio_->out << "Two Electron Integral Timings" << endl << bannerMid 
      << endl;
*/


    this->fileio_->out << BannerEnd << endl;
}

void AOIntegrals::printOneE(){
  std::vector<RealMap> mat;
  int NB = this->nBasis_;
  int NBSq = NB*NB;

  mat.push_back(RealMap(this->overlap_->data(),NB,NB));
  mat.push_back(RealMap(this->kinetic_->data(),NB,NB));
  mat.push_back(RealMap(this->potential_->data(),NB,NB));
  mat.push_back(RealMap(this->coreH_->data(),NB,NB));
  if(this->maxMultipole_ >= 1)
    for(auto i = 0, IOff=0; i < 3; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecDipole_->storage()[IOff],NB,NB));
  if(this->maxMultipole_ >= 2)
    for(auto i = 0, IOff=0; i < 6; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecQuadpole_->storage()[IOff],NB,NB));
  if(this->maxMultipole_ >= 3)
    for(auto i = 0, IOff=0; i < 10; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecOctpole_->storage()[IOff],NB,NB));
  if(this->maxMultipole_ >= 4) 
    for(auto i = 0, IOff=0; i < 3; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->RcrossDel_->storage()[IOff],NB,NB));
  
  /*
  if(this->nTCS_ == 2){
    prettyPrintTCS(this->fileio_->out,(mat[0]),"Overlap");
    if(this->maxMultipole_ >= 1){
      prettyPrintTCS(this->fileio_->out,(mat[4]),"Electric Dipole (x)");
      prettyPrintTCS(this->fileio_->out,(mat[5]),"Electric Dipole (y)");
      prettyPrintTCS(this->fileio_->out,(mat[6]),"Electric Dipole (z)");
    }
    if(this->maxMultipole_ >= 2){
      prettyPrintTCS(this->fileio_->out,(mat[7]), "Electric Quadrupole (xx)");
      prettyPrintTCS(this->fileio_->out,(mat[8]), "Electric Quadrupole (xy)");
      prettyPrintTCS(this->fileio_->out,(mat[9]), "Electric Quadrupole (xz)");
      prettyPrintTCS(this->fileio_->out,(mat[10]),"Electric Quadrupole (yy)");
      prettyPrintTCS(this->fileio_->out,(mat[11]),"Electric Quadrupole (yz)");
      prettyPrintTCS(this->fileio_->out,(mat[12]),"Electric Quadrupole (zz)");
    }
    if(this->maxMultipole_ >= 3){
      prettyPrintTCS(this->fileio_->out,(mat[13]),"Electric Octupole (xxx)");
      prettyPrintTCS(this->fileio_->out,(mat[14]),"Electric Octupole (xxy)");
      prettyPrintTCS(this->fileio_->out,(mat[15]),"Electric Octupole (xxz)");
      prettyPrintTCS(this->fileio_->out,(mat[16]),"Electric Octupole (xyy)");
      prettyPrintTCS(this->fileio_->out,(mat[17]),"Electric Octupole (xyz)");
      prettyPrintTCS(this->fileio_->out,(mat[18]),"Electric Octupole (xzz)");
      prettyPrintTCS(this->fileio_->out,(mat[19]),"Electric Octupole (yyy)");
      prettyPrintTCS(this->fileio_->out,(mat[20]),"Electric Octupole (yyz)");
      prettyPrintTCS(this->fileio_->out,(mat[21]),"Electric Octupole (yzz)");
      prettyPrintTCS(this->fileio_->out,(mat[22]),"Electric Octupole (zzz)");
    }
    if(this->maxMultipole_ >= 4){
      prettyPrintTCS(this->fileio_->out,(mat[23]),"R cross Del (x)");
      prettyPrintTCS(this->fileio_->out,(mat[24]),"R cross Del (y)");
      prettyPrintTCS(this->fileio_->out,(mat[25]),"R cross Del (z)");
    }

    prettyPrintTCS(this->fileio_->out,(mat[1]),"Kinetic");
    prettyPrintTCS(this->fileio_->out,(mat[2]),"Potential");
    prettyPrintTCS(this->fileio_->out,(mat[3]),"Core Hamiltonian");
  } else {
  */

    prettyPrint(this->fileio_->out,(mat[0]),"Overlap");
    if(this->maxMultipole_ >= 1){
      prettyPrint(this->fileio_->out,(mat[4]),"Electric Dipole (x)");
      prettyPrint(this->fileio_->out,(mat[5]),"Electric Dipole (y)");
      prettyPrint(this->fileio_->out,(mat[6]),"Electric Dipole (z)");
    }
    if(this->maxMultipole_ >= 2){
      prettyPrint(this->fileio_->out,(mat[7]), "Electric Quadrupole (xx)");
      prettyPrint(this->fileio_->out,(mat[8]), "Electric Quadrupole (xy)");
      prettyPrint(this->fileio_->out,(mat[9]), "Electric Quadrupole (xz)");
      prettyPrint(this->fileio_->out,(mat[10]),"Electric Quadrupole (yy)");
      prettyPrint(this->fileio_->out,(mat[11]),"Electric Quadrupole (yz)");
      prettyPrint(this->fileio_->out,(mat[12]),"Electric Quadrupole (zz)");
    }
    if(this->maxMultipole_ >= 3){
      prettyPrint(this->fileio_->out,(mat[13]),"Electric Octupole (xxx)");
      prettyPrint(this->fileio_->out,(mat[14]),"Electric Octupole (xxy)");
      prettyPrint(this->fileio_->out,(mat[15]),"Electric Octupole (xxz)");
      prettyPrint(this->fileio_->out,(mat[16]),"Electric Octupole (xyy)");
      prettyPrint(this->fileio_->out,(mat[17]),"Electric Octupole (xyz)");
      prettyPrint(this->fileio_->out,(mat[18]),"Electric Octupole (xzz)");
      prettyPrint(this->fileio_->out,(mat[19]),"Electric Octupole (yyy)");
      prettyPrint(this->fileio_->out,(mat[20]),"Electric Octupole (yyz)");
      prettyPrint(this->fileio_->out,(mat[21]),"Electric Octupole (yzz)");
      prettyPrint(this->fileio_->out,(mat[22]),"Electric Octupole (zzz)");
    }
    if(this->maxMultipole_ >= 4){
      prettyPrint(this->fileio_->out,(mat[23]),"R cross Del (x)");
      prettyPrint(this->fileio_->out,(mat[24]),"R cross Del (y)");
      prettyPrint(this->fileio_->out,(mat[25]),"R cross Del (z)");
    }
    prettyPrint(this->fileio_->out,(mat[1]),"Kinetic");
    prettyPrint(this->fileio_->out,(mat[2]),"Potential");
    prettyPrint(this->fileio_->out,(mat[3]),"Core Hamiltonian");
//  }
}

void AOIntegrals::alloc(){
  this->memManager_->setBlockSize(this->nBasis_*this->nBasis_*sizeof(double)); 
  this->memManager_->allocMem();

  this->checkMeta();
  this->allocOp();
  this->allocOrth();
//  this->memManager_->printSummary(cout);
  if(getRank() == 0) {
    if(this->maxMultipole_ >= 1) this->allocMultipole();
 
    this->pairConstants_ = std::unique_ptr<PairConstants>(new PairConstants);
    this->molecularConstants_ = 
      std::unique_ptr<MolecularConstants>(new MolecularConstants);
    this->quartetConstants_ = 
      std::unique_ptr<QuartetConstants>(new QuartetConstants);
 
/*
    if(this->isPrimary) 
      this->fileio_->iniStdOpFiles(this->basisSet_->nBasis());
*/
  }
#ifdef CQ_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
/* This whole block leaks memory like a siv (~ 8MB leaked for test 4!)
  int i,j,ij;
  this->R2Index_ = new int*[this->nBasis_];
  for(i=0;i<this->nBasis_;i++) this->R2Index_[i] = new int[this->nBasis_];
  for(i=0;i<this->nBasis_;i++) for(j=0;j<this->nBasis_;j++) {
    if(i>=j) ij=j*(this->nBasis_)-j*(j-1)/2+i-j;
    else ij=i*(this->nBasis_)-i*(i-1)/2+j-i;
    this->R2Index_[i][j] = ij;
  };


// initialize the FmT table
// Need to know the max L first
  this->FmTTable_ = new double*[MaxFmTPt];
  for(i=0;i<MaxFmTPt;i++) this->FmTTable_[i] = new double[MaxTotalL];
  this->generateFmTTable();
*/
}

void AOIntegrals::allocOp(){
  if(getRank() == 0) {
#ifndef USE_LIBINT
    try {
      // Raffenetti Two Electron Coulomb AOIntegrals
      this->twoEC_ = std::unique_ptr<RealMatrix>(
          new RealMatrix(this->nTT_,this->nTT_)); 
      // Raffenetti Two Electron Exchange AOIntegrals
      this->twoEX_ = std::unique_ptr<RealMatrix>(
          new RealMatrix(this->nTT_,this->nTT_)); 
    } catch (...) {
      CErr(std::current_exception(),
          "Coulomb and Exchange Tensor(R4) Allocation");
    }
#else
    try {
      if(this->integralAlgorithm == INCORE){ // Allocate R4 ERI Tensor
        this->aoERI_ = std::unique_ptr<RealTensor4d>(
          new RealTensor4d(this->nBasis_,this->nBasis_,this->nBasis_,this->nBasis_)); 
      }
    } catch (...) {
      CErr(std::current_exception(),"N^4 ERI Tensor Allocation");
    }
#endif

//******** BEGIN SUNSC *************
//   number of spherical harmonic gaussian
    int nSphBasis = 0;
    int nCartBasis = 0;
    for (auto mm = 0 ; mm < this->basisSet_->nShell() ; mm++ ) {
       nSphBasis += 2*(this->basisSet_->shellsCQ[mm].l)+1;
       nCartBasis += (this->basisSet_->shellsCQ[mm].l+1)*(this->basisSet_->shellsCQ[mm].l+2)/2;
    }
    this->nSphBasis_ = nSphBasis;
    this->nCartBasis_ = nCartBasis;
//**********************************
    try {
 
      // One Electron Integral
      auto NBSq = this->nBasis_ * this->nBasis_;
      this->coreH_      = 
        std::unique_ptr<RealMap>(
            new RealMap(this->memManager_->malloc<double>(NBSq),this->nBasis_,this->nBasis_)
        ); 
      this->oneEmx_      = 
        std::unique_ptr<RealMap>(
            new RealMap(this->memManager_->malloc<double>(NBSq),this->nBasis_,this->nBasis_)
        ); 
      this->oneEmy_      = 
        std::unique_ptr<RealMap>(
            new RealMap(this->memManager_->malloc<double>(NBSq),this->nBasis_,this->nBasis_)
        ); 
      this->oneEmz_      = 
        std::unique_ptr<RealMap>(
            new RealMap(this->memManager_->malloc<double>(NBSq),this->nBasis_,this->nBasis_)
        ); 
      // Overlap
      this->overlap_   = 
        std::unique_ptr<RealMap>(
            new RealMap(this->memManager_->malloc<double>(NBSq),this->nBasis_,this->nBasis_)
        ); 
      // Kinetic
      this->kinetic_   = 
        std::unique_ptr<RealMap>(
            new RealMap(this->memManager_->malloc<double>(NBSq),this->nBasis_,this->nBasis_)
        ); 
      // Potential
      this->potential_ = 
        std::unique_ptr<RealMap>(
            new RealMap(this->memManager_->malloc<double>(NBSq),this->nBasis_,this->nBasis_)
        ); 

//******** BEGIN SUNSC *************
      // Angular momentum
      this->angular_   = std::unique_ptr<RealTensor3d>(
        new RealTensor3d(this->nBasis_,this->nBasis_,3));
      this->SOCoupling_= std::unique_ptr<RealTensor3d>(
        new RealTensor3d(this->nBasis_,this->nBasis_,3));
      this->pVp_ = std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nBasis_));
      // pV dot p 
//**********************************

 
    } catch(...) {
      CErr(std::current_exception(),"One Electron Integral Tensor Allocation");
    }
    this->coreH_->setZero();
    this->oneEmx_->setZero();
    this->oneEmy_->setZero();
    this->oneEmz_->setZero();
    this->overlap_->setZero();
    this->kinetic_->setZero();
    this->potential_->setZero();
  } // End serial allocation

#ifdef USE_LIBINT
  try { 
    this->schwartz_ = std::unique_ptr<RealMap>(
      new RealMap(this->memManager_->malloc<double>(
          this->basisSet_->nShell()*this->basisSet_->nShell()),
          this->basisSet_->nShell(),this->basisSet_->nShell())
      ); 
  } catch(...) {
    CErr(std::current_exception(),"Schwartx Bound Tensor Allocation");
  } 
  this->schwartz_->setZero();

  if(getRank() == 0) {
    if(this->integralAlgorithm == DENFIT){
 
      try { 
        this->aoRII_ = std::unique_ptr<RealTensor3d>(
          new RealTensor3d(this->nBasis_,this->nBasis_,this->DFbasisSet_->nBasis())); 
 
        this->aoRIS_ = std::unique_ptr<RealTensor2d>(
          new RealTensor2d(
            this->DFbasisSet_->nBasis(),this->DFbasisSet_->nBasis()
          )
        );
 
      } catch (...) { 
        CErr(std::current_exception(),"Density Fitting Tensor Allocation");
      }
    }
  } // End Serial Allocation
#endif
}

void AOIntegrals::allocOrth() {
  auto NBSq = this->nBasis_ * this->nBasis_;
  this->ortho1_ = std::unique_ptr<RealMap>(
      new RealMap(this->memManager_->malloc<double>(NBSq),
        this->nBasis_,this->nBasis_));
  this->ortho2_ = std::unique_ptr<RealMap>(
      new RealMap(this->memManager_->malloc<double>(NBSq),
        this->nBasis_,this->nBasis_));
  this->ortho1_->setZero();
  this->ortho2_->setZero();
};

void AOIntegrals::allocMultipole(){
  try {

    if(this->maxMultipole_ >= 1)
      this->elecDipole_ = std::unique_ptr<RealTensor3d>(
        new RealTensor3d(this->nBasis_,this->nBasis_,3)
      );

    if(this->maxMultipole_ >= 2)
      this->elecQuadpole_ = std::unique_ptr<RealTensor3d>(
        new RealTensor3d(this->nBasis_,this->nBasis_,6)
      );

    if(this->maxMultipole_ >= 3)
      this->elecOctpole_ = std::unique_ptr<RealTensor3d>(
        new RealTensor3d(this->nBasis_,this->nBasis_,10)
      );

    if(this->maxMultipole_ >= 4)
      this->RcrossDel_ = std::unique_ptr<RealTensor3d>(
        new RealTensor3d(this->nBasis_,this->nBasis_,3)
      );

  } catch(...) {
    CErr(std::current_exception(),"Multipole Tensor Allocation");
  }
}

/*
void AOIntegrals::writeOneE(){
  this->fileio_->overlap->write(this->overlap_->data(),H5::PredType::NATIVE_DOUBLE);
  this->fileio_->kinetic->write(this->kinetic_->data(),H5::PredType::NATIVE_DOUBLE);
  this->fileio_->nucRepl->write(this->potential_->data(),H5::PredType::NATIVE_DOUBLE);
  this->fileio_->coreHam->write(this->coreH_->data(),H5::PredType::NATIVE_DOUBLE);
}
*/

//xslis
void AOIntegrals::createShellPairs() {
  int i,j,k,l,m,n,iPP,nPP,iAtom;
  double squareAB,prodexp,sumexp,KAB,pairNorm;
  double sqrt2pi54 = 5.91496717279561; // sqrt(2)*pi^{5/4}
  std::array<double,MaxTotalL> tmpFmT;
  std::array<double,3> tmpP,tmpPA,tmpPB,tmpPZeta,tmpssL,AC,BC,ACxBC;
  ChronusQ::ShellCQ *iS,*jS;
  ChronusQ::ShellPair *ijS;

  int nShell = this->basisSet_->nShell();
  this->nShellPair_ = nShell*(nShell+1)/2;

  if (this->useFiniteWidthNuclei == true) {
    this->molecule_->generateFiniteWidthNuclei();
  }         // generate finite width nuclei 
//cout<<this->molecule_->nucShell(0)<<endl;

//this->shellPairs_ = new ShellPair[this->nShellPair_];
// This does the same thing as a 'new' call when passed a dimension and
// dones't require annoying memory managment
  this->shellPairs_ = std::vector<ChronusQ::ShellPair>(this->nShellPair_);
  n = 0;


  int spherical_lsize[nShell],cart_lsize[nShell];
  int mapSh2spBf[nShell],mapSh2caBf[nShell];
  int tempmapSh2spBf=0,tempmapSh2caBf=0;
  for ( i =0 ; i<nShell ; i++ ) {
    mapSh2spBf[i] = tempmapSh2spBf;
    mapSh2caBf[i] = tempmapSh2caBf;
/*    if (this->basisSet_->shellsCQ[i].cartesian_l.size() == 6){
       spherical_lsize[i]=5;
    }
    else {
      spherical_lsize[i]=this->basisSet_->shellsCQ[i].cartesian_l.size();
    }
*/
    spherical_lsize[i]=2*(this->basisSet_->shellsCQ[i].l)+1;
    cart_lsize[i] = (this->basisSet_->shellsCQ[i].l+1)*(this->basisSet_->shellsCQ[i].l+2)/2;
    tempmapSh2spBf += spherical_lsize[i];
    tempmapSh2caBf += cart_lsize[i];
  }

cout<<"total bumber of spherical function= "<<tempmapSh2spBf<<endl;
cout<<"total number of cartesian function= "<<tempmapSh2caBf<<endl;  

  for(i=0;i<nShell;i++) 
    {      
//  int ibf_s1 = this->basisSet_->mapSh2Bf(i);
//  int n1 = this->basisSet_->shells(i).size();
//  int n1 = this->basisSet_->shellsCQ[i].l;
    
    this->fileio_->out << "ShellIndex= " << i <<endl;
    int n1 = this->basisSet_->shellsCQ[i].cartesian_l.size();
    this->fileio_->out << "Size= "  <<n1 <<endl;
    this->fileio_->out << "Size_of_spherical= " << spherical_lsize[i]<<endl;
    this->fileio_->out << "start of shell function" << this->basisSet_->mapSh2Bf(i)<<endl;
    this->fileio_->out << "start of shell for spherical function" << mapSh2spBf[i]<<endl;
    for (int ii =0; ii < n1; ii++){
    this->fileio_->out << "BasisSet ind " << ii <<endl;
      for(k=0;k<3;k++){
       auto  lxyz = this->basisSet_->shellsCQ[i].cartesian_l[ii][k];
         this->fileio_->out << lxyz<<"\t";
      } 
     this->fileio_->out << endl;
      
    }
  }

  for(i=0;i<nShell;i++) for(j=i;j<nShell;j++) {
    ijS = &(this->shellPairs_[n]);
    if(this->basisSet_->shellsCQ[i].l>this->basisSet_->shellsCQ[j].l) {    
      ijS->iShell = &(this->basisSet_->shellsCQ[i]);
      ijS->jShell = &(this->basisSet_->shellsCQ[j]);
      ijS->ibf_s = this->basisSet_->mapSh2Bf(i);
      ijS->jbf_s = this->basisSet_->mapSh2Bf(j);
//start number for spherical gaussian       
      ijS->isphbf_s = mapSh2spBf[i];
      ijS->jsphbf_s = mapSh2spBf[j];
//start number for cartesian gaussian
      ijS->icarbf_s = mapSh2caBf[i];
      ijS->jcarbf_s = mapSh2caBf[j];
//cout<< "iShell->l= "<<ijS->iShell->l<<endl;
//cout<<ijS->ibf_s<<"\t"<<ijS->icarbf_s<<endl;
//cout<<ijS->jbf_s<<"\t"<<ijS->jcarbf_s<<endl;
     
    } else {
      ijS->iShell = &(this->basisSet_->shellsCQ[j]);
      ijS->jShell = &(this->basisSet_->shellsCQ[i]);
      ijS->ibf_s = this->basisSet_->mapSh2Bf(j);
      ijS->jbf_s = this->basisSet_->mapSh2Bf(i);
      ijS->isphbf_s = mapSh2spBf[j];
      ijS->jsphbf_s = mapSh2spBf[i];
      ijS->icarbf_s = mapSh2caBf[j];
      ijS->jcarbf_s = mapSh2caBf[i];
//cout<<ijS->ibf_s<<"\t"<<ijS->icarbf_s<<endl;
//cout<<ijS->jbf_s<<"\t"<<ijS->jcarbf_s<<endl;
 
   };
//cout<<"isbf_s= "<<ijS->ibf_s<<"\t jbf_s= "<<ijS->jbf_s<<endl;
//cout<<"isphbf_s= "<<ijS->isphbf_s<<"\t jsphbf_s= "<<ijS->jsphbf_s<<endl<<endl;
    ijS->isphsize = 2*(ijS->iShell->l)+1;
    ijS->jsphsize = 2*(ijS->jShell->l)+1;
    ijS->icarsize = (ijS->iShell->l+1)*(ijS->iShell->l+2)/2;
    ijS->jcarsize = (ijS->jShell->l+1)*(ijS->jShell->l+2)/2;

    iS = ijS->iShell;
    jS = ijS->jShell;

    ijS->lTotal = iS->l + jS->l;

    ijS->A = iS->O;
    ijS->B = jS->O;
    squareAB = math.zero;
    for(m=0;m<3;m++) {
      ijS->AB[m] = ijS->A[m]-ijS->B[m];
      AC[m]      = ijS->A[m];
      BC[m]      = ijS->B[m];
      squareAB  += ijS->AB[m]*ijS->AB[m];
    };
    ACxBC[0]     = ijS->A[1]*ijS->B[2]-ijS->A[2]*ijS->B[1];  
    ACxBC[1]     = ijS->A[2]*ijS->B[0]-ijS->A[0]*ijS->B[2];
    ACxBC[2]     = ijS->A[0]*ijS->B[1]-ijS->A[1]*ijS->B[0];   
    
    nPP = 0;
    for(k=0; k<iS->coeff.size();k++) 
    for(l=0; l<jS->coeff.size();l++) {
      sumexp  = iS->alpha[k] + jS->alpha[l];
      prodexp = iS->alpha[k] * jS->alpha[l];
      
      KAB = exp(-squareAB*prodexp/sumexp);
      ijS->KAB.push_back(KAB);
      ijS->beta.push_back(jS->alpha[l]);
      ijS->UAB.push_back(sqrt2pi54*KAB/sumexp);
      ijS->zetaa.push_back(iS->alpha[k]);
      ijS->zetab.push_back(jS->alpha[l]);
//      cout<<iS->alpha[k]<<";"<<ijS->zetaa[nPP]<<endl;
//      cout<<jS->alpha[l]<<";"<<ijS->zetab[nPP]<<endl;
      ijS->Zeta.push_back(sumexp);    // zeta = alpha + beta
      ijS->Xi.push_back(prodexp/sumexp); // xi = alpha*beta/(alpha+beta)
      ijS->invZeta.push_back(math.one/sumexp);
      ijS->halfInvZeta.push_back(math.half/sumexp);
      pairNorm = iS->norm[k]*jS->norm[l];
      ijS->norm.push_back(pairNorm);
      ijS->ss.push_back(pairNorm*sqrt((math.pi/sumexp)*(math.pi/sumexp)*(math.pi/sumexp))*KAB);
      for (m=0;m<3;m++) {
        tmpP[m]     = (iS->alpha[k]*ijS->A[m]+jS->alpha[l]*ijS->B[m])/sumexp;// P=(alpha*A+beta*B)/zeta
        tmpPZeta[m] = tmpP[m] * sumexp;
        tmpPA[m]    = tmpP[m] - ijS->A[m]; // PA = P - A
        tmpPB[m]    = tmpP[m] - ijS->B[m]; // PB = P - B
//        cout<<ijS->B[m]<<"t";
      };
//      cout<<endl;

      ijS->P.push_back(tmpP);
      ijS->PA.push_back(tmpPA);
      ijS->PB.push_back(tmpPB);
      ijS->PZeta.push_back(tmpPZeta);

//      if ( this->useFiniteWidthNuclei == false ) {        
        ijS->ssV.push_back(math.two*sqrt(ijS->Zeta[nPP]/math.pi)*ijS->ss[nPP]);
//      }
      ijS->ssT.push_back((ijS->Xi[nPP])*(math.three-math.two*(ijS->Xi[nPP])*squareAB)*ijS->ss[nPP]);
      for (m=0;m<3;m++) {
        tmpssL[m]    = math.two*(ijS->Xi[nPP])*ACxBC[m]*ijS->ss[nPP]; 
      }; 
      ijS->ssL.push_back(tmpssL);

//print ssL
/*
        std::string xyz[3];
        xyz[0] = "x";
        xyz[1] = "y";
        xyz[2] = "z"; 
      cout<<"ssL"<<"\t"<<i+1<<";"<<j+1<<";"<<"\t"<<nPP+1<<"\t";
         
      for (m=0;m<3;m++) {       
        cout<<xyz[m]<<"\t"<<ijS->ssL[nPP][m]/pairNorm<<","<<"\t";
      };
      cout<<endl;
*/


      nPP++;
    };
    ijS->nPGTOPair = nPP;
    n++;
  }

       

};
//xslie
