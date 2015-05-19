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
#include <basisset.h>
using ChronusQ::atom; 
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
//-------------//
// constructor //
//-------------//
BasisSet::BasisSet(int nBasis, int nShell) {
  if(nBasis>0&&nShell>0) {
    nBasis_ = nBasis;
    nShell_ = nShell;
    nShellPair_ = nShell_*(nShell_+1)/2;
    iniBasisSet();
  };
};
//-------------------//
// initialize memory //
//-------------------//
void BasisSet::iniBasisSet(){
  // FIXME Need to move these over to unique_ptr
  ao = new (nothrow) AOCartesian[this->nBasis_];
  if(ao==NULL) CErr("BasisSet::AOCartesian Allocation");
  shells = new (nothrow) Shell[this->nShell_];
  if(shells==NULL) CErr("BasiSset::Shell Allocation");
  sortedShells = new (nothrow) int[this->nShell_];
  if(sortedShells==NULL) CErr("BasisSet::sortedShells");
  shellPairs = new (nothrow) ShellPair[this->nShellPair_];
  if(shellPairs==NULL) CErr("BasisSet::shellPairs");
};
//---------------------------------//
// print out basis set information //
//---------------------------------//
void BasisSet::printAO(ostream &output){
  int i,j,k,aoIndex;
  output.precision(6);
  output.fill(' ');
  output.setf(ios::right,ios::adjustfield);
  output.setf(ios::fixed,ios::floatfield);
  output<<endl<<"Cartesian Basis Functions:"<<endl;
  output<<bannerTop<<endl;
  output<<std::setw(16)<<"         Center"<<std::setw(4)<<"lx"<<std::setw(4)<<"ly"<<std::setw(4)<<"lz"
	<<std::setw(18)<<"exponent"<<std::setw(18)<<"coefficient"<<endl;
  output<<bannerMid;
  for(i=0;i<this->nShell_;i++) {
    output<<endl;
    output<<std::setw(5)<<i+1<<std::setw(4)<<shells[i].name<<std::setw(5)<<shells[i].center+1<<"  ";
    aoIndex = shells[i].aoIndex;
    for(k=aoIndex;k<(aoIndex+HashNAOs(shells[i].L));k++) {
      if(k!=aoIndex) output<<std::setw(16)<<" ";
      output<<std::setw(4)<<ao[k].lx<<std::setw(4)<<ao[k].ly<<std::setw(4)<<ao[k].lz;
      for(j=0;j<this->shells[i].nPGTOs;j++) {
	if(j!=0) output<<std::setw(28)<<" ";
	output<<std::setw(18)<<shells[i].expo[j]<<std::setw(18)<<shells[i].coef[j]<<endl;
      };
    };
  };
  output<<bannerEnd<<endl;
};
//--------------------------------------------//
// Read basis set information from input file //
//--------------------------------------------//
void BasisSet::readBasisSet(std::shared_ptr<FileIO> fileio, std::shared_ptr<Molecule> mol){
  int i,j,n,k;
  // TODO Have this be read to a string
  // TODO Check for basis in BASIS_PATH (should be set up though cmake)
  std::string readString;
  fileio->in>>readString;
  cout << readString <<endl;
/*
  std::unique_ptr<ifstream> fileBasis;
  if(fexists(BASIS_PATH+"/"+readString)) {
    fileBasis = std::unique_ptr<ifstream>(new ifstream(readString));
  } else {
    cout << "Could not find basis set file" << endl;
    exit(EXIT_FAILURE);
  }
*/
  std::string basis_path = "/" + readString;
  basis_path.insert(0,BASIS_PATH);
  std::unique_ptr<ifstream> fileBasis(new ifstream (basis_path));
  if(fileBasis->fail()){ // Check if file is in BASIS_PATH
    fileBasis.reset();
    fileBasis = std::unique_ptr<ifstream>(new ifstream(readString));
    if(fileBasis->fail()) { // Check if file is in PWD
      CErr("Could not find basis set file",fileio->out);
    } else {
      fileio->out << "Reading Basis Set from: ./" << readString<< endl;
    }
  } else {
    fileio->out << "Reading Basis Set from: " << basis_path<< endl;
  }
  double threePI=math.pi*math.pi*math.pi,readNorm, readExp, readCoeff1, readCoeff2;
  int    nBasis, nShell, readNPGTO, L;
  std::string   atomStr;
  //----------------------------------------------------------//
  // Read cartesian basis functions                           //
  // 1. first round of reading tests the dimension of problem //
  // 2. second round reads information                        //
  //----------------------------------------------------------//
  for(i=0; i<20; i++) this->nLShell_[i] = 0;
  this->nShell_=0;
  for(i=0; i<mol->nAtoms(); i++){
    atomStr = "-";
    atomStr = atomStr+atom[mol->index(i)].symbol;
    while(atomStr.compare(readString)&&(!fileBasis->eof())) *fileBasis>>readString;
    *fileBasis >> readString;
    while(readString.compare("****")) { 
      *fileBasis >> readNPGTO;
      *fileBasis >> readNorm ;
      // TODO Extend this to higher angular momentum (follow though with the rest of
      //      the code)
      if(!(readString.compare("S"))||!(readString.compare("P"))||!(readString.compare("D"))||
         !(readString.compare("F"))||!(readString.compare("G"))){
        //S,P,D,F,G shell
        L = HashL(readString);
        (this->nLShell_[L])++;
        (this->nShell_)++;
        (this->nBasis_)=(this->nBasis_)+HashNAOs(L);
        (this->nPrimitive_)=(this->nPrimitive_)+readNPGTO*HashNAOs(L);
        for(j=0;j<2*readNPGTO;j++) *fileBasis>>readString;
      // q
      // Is this obsolete?
      } else if(!(readString.compare("SP"))){
        //SP shell
        (this->nLShell_[0])++;
        (this->nLShell_[1])++;
        (this->nShell_)++;
        (this->nShell_)++;
        (this->nBasis_)+=4;
        (this->nPrimitive_)=(this->nPrimitive_)+readNPGTO*4;
        for(j=0;j<3*readNPGTO;j++) *fileBasis>>readString;
      } else {
        fileio->out<<"Error: unrecognized shell symbol!"<<endl;
        CErr("Unrecognized Shell symbol");
      };
      *fileBasis>>readString;
    };
    fileBasis->seekg(0);
  };
  this->nShellPair_=this->nShell_*(this->nShell_+1)/2;
  this->iniBasisSet();

  fileBasis->seekg(0);
  nBasis = 0;
  nShell = 0;
  for(i=0; i<mol->nAtoms(); i++){
    atomStr="-";
    atomStr=atomStr+atom[mol->index(i)].symbol;
    while((atomStr.compare(readString))&&(!fileBasis->eof())) *fileBasis>>readString;
    *fileBasis >> readString;
    while(readString.compare("****")) { 
      *fileBasis >> readNPGTO;
      *fileBasis >> readNorm ;
      // TODO Extend this to higher angular momentum (follow though with the rest of
      //      the code)
      if((!readString.compare("S"))||!(readString.compare("P"))||!(readString.compare("D"))||
         !(readString.compare("F"))||!(readString.compare("G"))){
        //S,P,D,F,G shell
        L = HashL(readString);
	strcpy(this->shells[nShell].name,readString.c_str());
        (this->shells[nShell]).center = i;
        (this->shells[nShell]).SP     = false;
        (this->shells[nShell]).L      = L;
        (this->shells[nShell]).nPGTOs = readNPGTO;
        (this->shells[nShell]).aoIndex= nBasis;
        for(n=0;n<HashNAOs(L);n++) {
          (this->ao[nBasis+n]).shIndex = nShell;
          (this->ao[nBasis+n]).divConst = 1;
          if(L==0) {
	    (this->ao[nBasis+n]).lx = 0;
	    (this->ao[nBasis+n]).ly = 0;
	    (this->ao[nBasis+n]).lz = 0;
	    for(k=0;k<3;k++) (this->ao[nBasis+n]).l[k] = 0;
          } else if(L==1) {
	    (this->ao[nBasis+n]).lx = (n==0)?1:0;
	    (this->ao[nBasis+n]).ly = (n==1)?1:0;
	    (this->ao[nBasis+n]).lz = (n==2)?1:0;
	    (this->ao[nBasis+n]).l[0] = (n==0)?1:0;
	    (this->ao[nBasis+n]).l[1] = (n==1)?1:0;
	    (this->ao[nBasis+n]).l[2] = (n==2)?1:0;
          } else if(L==2) {
	    (this->ao[nBasis+n]).lx = (n==3||n==4)?1:0;
	    (this->ao[nBasis+n]).ly = (n==3||n==5)?1:0;
	    (this->ao[nBasis+n]).lz = (n==4||n==5)?1:0;
	    if(n==0) (this->ao[nBasis+n]).lx = 2;
	    if(n==1) (this->ao[nBasis+n]).ly = 2;
	    if(n==2) (this->ao[nBasis+n]).lz = 2;
	    (this->ao[nBasis+n]).l[0] = (n==3||n==4)?1:0;
	    (this->ao[nBasis+n]).l[1] = (n==3||n==5)?1:0;
	    (this->ao[nBasis+n]).l[2] = (n==4||n==5)?1:0;
	    if(n==0) (this->ao[nBasis+n]).l[0] = 2;
	    if(n==1) (this->ao[nBasis+n]).l[1] = 2;
	    if(n==2) (this->ao[nBasis+n]).l[2] = 2;
	    if(n<3 ) (this->ao[nBasis+n]).divConst = sqrt(3.0);
          } else if(L==3) {
          } else if(L==4) {
          };
        };
        for(j=0;j<readNPGTO;++j){
          *fileBasis >> readExp >> readCoeff2;
          if(L==0) readNorm = sqrt(sqrt(8*readExp*readExp*readExp/threePI));
          else if (L==1) readNorm = sqrt(sqrt(128*readExp*readExp*readExp*readExp*readExp/threePI));
          else if (L==2) readNorm = sqrt(sqrt(2048*readExp*readExp*readExp*readExp*readExp*readExp*readExp/threePI));
          (this->shells[nShell]).expo[j] = readExp;
          (this->shells[nShell]).coef[j] = readCoeff2;
          (this->shells[nShell]).norm[j] = readNorm*readCoeff2;
        };
        nShell++;
        nBasis=nBasis+HashNAOs(L);
      } else if(!(readString.compare("SP"))){
        //SP Shell
        strcpy((this->shells[nShell]).name,"S");
        (this->shells[nShell]).center = i;
        (this->shells[nShell]).SP     = false;
        (this->shells[nShell]).L      = 0;
        (this->shells[nShell]).nPGTOs = readNPGTO;
        (this->shells[nShell]).aoIndex= nBasis;
        (this->ao[nBasis]).lx = 0;
        (this->ao[nBasis]).ly = 0;
        (this->ao[nBasis]).lz = 0;
	for(k=0;k<3;k++) (this->ao[nBasis]).l[k] = 0;
        (this->ao[nBasis]).shIndex = nShell;
        (this->ao[nBasis]).divConst = 1;
        strcpy((this->shells[nShell+1]).name,"P");
        (this->shells[nShell+1]).center = i;
        (this->shells[nShell+1]).L      = 1;
        (this->shells[nShell+1]).nPGTOs = readNPGTO;
        (this->shells[nShell+1]).aoIndex= nBasis+1;
        for(n=1;n<4;n++) {
          (this->ao[nBasis+n]).lx = (n==1)?1:0;
          (this->ao[nBasis+n]).ly = (n==2)?1:0;
          (this->ao[nBasis+n]).lz = (n==3)?1:0;
          (this->ao[nBasis+n]).l[0] = (n==1)?1:0;
          (this->ao[nBasis+n]).l[1] = (n==2)?1:0;
          (this->ao[nBasis+n]).l[2] = (n==3)?1:0;
          (this->ao[nBasis+n]).shIndex = nShell;
          (this->ao[nBasis+n]).divConst = 1;
        };
        for(j=0;j<readNPGTO;++j){
          *fileBasis >> readExp >> readCoeff1 >> readCoeff2;
          (this->shells[nShell  ]).expo[j] = readExp;
          (this->shells[nShell  ]).coef[j] = readCoeff1;
          (this->shells[nShell+1]).expo[j] = readExp;
          (this->shells[nShell+1]).coef[j] = readCoeff2;
          readNorm = sqrt(sqrt(8*readExp*readExp*readExp/threePI));
          (this->shells[nShell]).norm[j] = readNorm*readCoeff1;
          readNorm = sqrt(sqrt(128*readExp*readExp*readExp*readExp*readExp/threePI));
          (this->shells[nShell+1]).norm[j] = readNorm*readCoeff2;
        };
        nShell+=2;
        nBasis+=4;
      };
      *fileBasis>>readString;
    };
    fileBasis->seekg(0);
  };
  fileBasis->close();
  this->createShellPair(mol);
};
//------------------------------------//
// print out sorted shell information //
//------------------------------------//
void BasisSet::printShell(ostream &output){
  int i,j,k,shIndex,aoIndex;
  output.precision(6);
  output.fill(' ');
  output.setf(ios::right,ios::adjustfield);
  output.setf(ios::fixed,ios::floatfield);
  output<<endl<<"Atmoic Orbital Shells:"<<endl;
  output<<bannerTop<<endl;
  output<<std::setw(16)<<"         Center"<<std::setw(4)<<"lx"<<std::setw(4)<<"ly"<<std::setw(4)<<"lz"
	<<std::setw(18)<<"exponent"<<std::setw(18)<<"coefficient"<<endl;
  output<<bannerMid;
  for(i=0;i<this->nShell_;i++) {
    output<<endl;
    shIndex = sortedShells[i];
    output<<std::setw(5)<<i+1<<std::setw(4)<<shells[shIndex].name<<std::setw(5)<<shells[shIndex].center+1<<"  ";
    aoIndex = shells[shIndex].aoIndex;
    for(k=aoIndex;k<(aoIndex+HashNAOs(shells[shIndex].L));k++) {
      if(k!=aoIndex) output<<std::setw(16)<<" ";
      output<<std::setw(4)<<ao[k].lx<<std::setw(4)<<ao[k].ly<<std::setw(4)<<ao[k].lz;
      for(j=0;j<this->shells[shIndex].nPGTOs;j++) {
	if(j!=0) output<<std::setw(28)<<" ";
	output<<std::setw(18)<<shells[shIndex].expo[j]<<std::setw(18)<<shells[shIndex].coef[j]<<endl;
      };
    };
  };
  output<<bannerEnd<<endl;
};
//------------------------------------//
// print out shell pair information //
//------------------------------------//
void BasisSet::printShellPair(ostream &output){
  int i,j,k,shIndex,aoIndex;
  output.precision(6);
  output.fill(' ');
  output.setf(ios::right,ios::adjustfield);
  output.setf(ios::fixed,ios::floatfield);
  output<<endl<<"Atmoic Orbital Shell Pairs:"<<endl;
  output<<bannerTop<<endl;
  output<<std::setw(8)<<"Index"
	<<std::setw(12)<<"Center-1"<<std::setw(11)<<"Shell-1"<<std::setw(6)<<"L-1"<<std::setw(6)<<"*"
	<<std::setw(10)<<"Center-2"<<std::setw(11)<<"Shell-2"<<std::setw(6)<<"L-2"<<endl;
  output<<bannerMid<<endl;
  for(i=0;i<this->nShellPair_;i++)
    output<<std::setw(8)<<i+1
	  <<std::setw(12)<<shellPairs[i].center[0]+1<<std::setw(11)<<shellPairs[i].shIndex[0]+1<<std::setw(6)<<shellPairs[i].L[0]<<std::setw(6)<<"*"
	  <<std::setw(10)<<shellPairs[i].center[1]+1<<std::setw(11)<<shellPairs[i].shIndex[1]+1<<std::setw(6)<<shellPairs[i].L[1]<<endl;
  output<<bannerEnd<<endl;
};
//---------------------------------------//
// print a general basis set information //
//---------------------------------------//
void BasisSet::printInfo(std::shared_ptr<FileIO> fileio,std::shared_ptr<Controls> controls) {
  fileio->out<<"\nBasis Function Information:"<<endl;
  fileio->out<< std::setw(15) << "nBasis =" << std::setw(8)<< nBasis_<< std::setw(5) <<" "
	     << std::setw(20) << "nPrimitive =" << std::setw(8)<< nPrimitive_<< endl;
  fileio->out<< std::setw(15) << "nSShell =" << std::setw(8)<< nLShell_[0]<< std::setw(5) <<" "
	     << std::setw(20) << "nPShell =" << std::setw(8)<< nLShell_[1]<< endl;
  fileio->out<< std::setw(15) << "nDShell =" << std::setw(8)<< nLShell_[2]<< std::setw(5) <<" "
	     << std::setw(20) << "nFShell =" << std::setw(8)<< nLShell_[3]<< endl;
  fileio->out<< std::setw(15) << "nGShell =" << std::setw(8)<< nLShell_[4]<< endl;
  if(controls->printLevel>=2) printAtomO(fileio->out);
  //if(controls->printLevel>=3) printShellPair(fileio->out);
};
//--------------------------------------------------------------//
// create and sort shell pairs according to the angular momenta //
//--------------------------------------------------------------//
void BasisSet::createShellPair(std::shared_ptr<Molecule> mol) {
  int i,j,k,l,n=0;
  // sort shells according to:
  //   (1) angular momentum
  for(i=0;i<nShell_;i++) sortedShells[i]=i;
  for(i=1;i<nShell_;i++) for(j=0;j<(nShell_-i);j++) {
    if(shells[sortedShells[j]].L<shells[sortedShells[j+1]].L) {
      n = sortedShells[j];
      sortedShells[j] = sortedShells[j+1];
      sortedShells[j+1] = n;
    };
  };

  n = 0;
  for(i=0;i<nShell_;i++) for(j=i;j<nShell_;j++) {
    shellPairs[n].L[0]=shells[sortedShells[i]].L;
    shellPairs[n].L[1]=shells[sortedShells[j]].L;
    shellPairs[n].shIndex[0]=sortedShells[i];
    shellPairs[n].shIndex[1]=sortedShells[j];
    n++;
  };

  // compute two-e constants related to shellPairs
  double sqrAB;
  double prodexp,expoi,expoj;
  int    sh0, sh1, iMax, jMax, nAOPair;
  for (n=0;n<nShellPair_;n++) {
    sh0 = shellPairs[n].shIndex[0];
    sh1 = shellPairs[n].shIndex[1];
    shellPairs[n].center[0]=shells[sh0].center;
    shellPairs[n].center[1]=shells[sh1].center;
    shellPairs[n].aoIndex[0]=shells[sh0].aoIndex;
    shellPairs[n].aoIndex[1]=shells[sh1].aoIndex;
    shellPairs[n].nPGTOs[0] = shells[sh0].nPGTOs;
    shellPairs[n].nPGTOs[1] = shells[sh1].nPGTOs;
    shellPairs[n].nBasis[0] = HashNAOs(shells[sh0].L);
    shellPairs[n].nBasis[1] = HashNAOs(shells[sh1].L);
    shellPairs[n].LTotal=shellPairs[n].L[0]+shellPairs[n].L[1];

    nAOPair = 0;
    if(shellPairs[n].aoIndex[0]!=shellPairs[n].aoIndex[1]) {
      iMax = shellPairs[n].aoIndex[0]+shellPairs[n].nBasis[0];
      jMax = shellPairs[n].aoIndex[1]+shellPairs[n].nBasis[1];
      for (i=shellPairs[n].aoIndex[0];i<iMax;i++) for (j=shellPairs[n].aoIndex[1];j<jMax;j++) {
        shellPairs[n].aoPairIndex[nAOPair][0]=i;
        shellPairs[n].aoPairIndex[nAOPair][1]=j;
        shellPairs[n].divConst[nAOPair] = ao[i].divConst*ao[j].divConst;
        nAOPair++;
      };
    } else {
      iMax = shellPairs[n].aoIndex[0]+shellPairs[n].nBasis[0];
      for (i=shellPairs[n].aoIndex[0];i<iMax;i++) for (j=i;j<iMax;j++) {
        shellPairs[n].aoPairIndex[nAOPair][0]=i;
        shellPairs[n].aoPairIndex[nAOPair][1]=j;
        shellPairs[n].divConst[nAOPair] = ao[i].divConst*ao[j].divConst;
        nAOPair++;
      };
    };
    shellPairs[n].nAOPair = nAOPair;

    sqrAB=0.0;
    for(k=0;k<3;k++) {
      shellPairs[n].centerA[k] = (*(mol->cart()))(k,(shellPairs[n].center[0]));
      shellPairs[n].centerB[k] = (*(mol->cart()))(k,(shellPairs[n].center[1]));
      shellPairs[n].deltaAB[k] = shellPairs[n].centerA[k]-shellPairs[n].centerB[k];
      sqrAB += shellPairs[n].deltaAB[k]*shellPairs[n].deltaAB[k];
    };

    for(i=0;i<(shells[sh0]).nPGTOs;i++)
      for(j=0;j<(shells[sh1]).nPGTOs;j++) {
        shellPairs[n].norm[i][j] = (shells[sh0]).norm[i]*(shells[sh1]).norm[j];
        expoi = (shells[sh0]).expo[i];
        expoj = (shells[sh1]).expo[j];
        shellPairs[n].zeta[i][j] = expoi + expoj;
        shellPairs[n].inversezeta[i][j] = math.half/shellPairs[n].zeta[i][j];
        shellPairs[n].invzeta[i][j] = math.one/shellPairs[n].zeta[i][j];
        prodexp = expoi * expoj;
        shellPairs[n].KAB[i][j] = math.sqrt2pi54*exp(-sqrAB*prodexp/(shellPairs[n].zeta[i][j]))/(shellPairs[n].zeta[i][j]);
        shellPairs[n].UAB[i][j] = shellPairs[n].norm[i][j]*math.sqrt2pi54*exp(-sqrAB*prodexp/(shellPairs[n].zeta[i][j]))/(shellPairs[n].zeta[i][j]);
        for (k=0;k<3;k++) {
          shellPairs[n].centerP[k][i][j] = (expoi*shellPairs[n].centerA[k]+expoj*shellPairs[n].centerB[k])/shellPairs[n].zeta[i][j];
          shellPairs[n].centerPZeta[k][i][j] = shellPairs[n].centerP[k][i][j]*shellPairs[n].zeta[i][j];
          shellPairs[n].deltaPA[k][i][j] = shellPairs[n].centerP[k][i][j] - shellPairs[n].centerA[k];
          shellPairs[n].deltaPB[k][i][j] = shellPairs[n].centerP[k][i][j] - shellPairs[n].centerB[k];
        }
        for (k=0;k<3;k++) for (l=0;l<3;l++) {
          shellPairs[n].deltaPApPB[k][l][i][j] = shellPairs[n].deltaPA[k][i][j] + shellPairs[n].deltaPB[l][i][j];
          shellPairs[n].deltaPAtPB[k][l][i][j] = shellPairs[n].deltaPA[k][i][j]*shellPairs[n].deltaPB[l][i][j];
        };
    };
  };
};

/*************************/
/* MPI Related Routines  */
/*************************/
/*
void BasisSet::mpiSend(int toID,int tag) {
  int i,j;
  OOMPI_COMM_WORLD[toID].Send(this->nCartBasis_,tag);
  OOMPI_COMM_WORLD[toID].Send(this->nCartPrimitive_,tag);
  for(i=0;i<this->nCartBasis_;i++) {
    OOMPI_COMM_WORLD[toID].Send(this->aoCart[i].center,tag);
    OOMPI_COMM_WORLD[toID].Send(this->aoCart[i].nPGTOs,tag);
    OOMPI_COMM_WORLD[toID].Send(this->aoCart[i].x,tag);
    OOMPI_COMM_WORLD[toID].Send(this->aoCart[i].y,tag);
    OOMPI_COMM_WORLD[toID].Send(this->aoCart[i].z,tag);
    for(j=0;j<3;j++)
      OOMPI_COMM_WORLD[toID].Send(this->aoCart[i].pGTO[j],MAXCONTRACTION,tag);
  }
};
void BasisSet::mpiRecv(int fromID,int tag) {
  int i,j;
  OOMPI_COMM_WORLD[fromID].Recv(this->nCartBasis_,tag);
  OOMPI_COMM_WORLD[fromID].Recv(this->nCartPrimitive_,tag);
  for(i=0;i<this->nCartBasis_;i++) {
    OOMPI_COMM_WORLD[fromID].Recv(this->aoCart[i].center,tag);
    OOMPI_COMM_WORLD[fromID].Recv(this->aoCart[i].nPGTOs,tag);
    OOMPI_COMM_WORLD[fromID].Recv(this->aoCart[i].x,tag);
    OOMPI_COMM_WORLD[fromID].Recv(this->aoCart[i].y,tag);
    OOMPI_COMM_WORLD[fromID].Recv(this->aoCart[i].z,tag);
    for(j=0;j<3;j++)
      OOMPI_COMM_WORLD[fromID].Recv(this->aoCart[i].pGTO[j],MAXCONTRACTION,tag);
  }
};
*/
