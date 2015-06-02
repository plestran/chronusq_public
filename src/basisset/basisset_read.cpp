
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

#ifdef USE_LIBINT
using ChronusQ::BasisSet;
using ChronusQ::Molecule;
using ChronusQ::HashL;
using ChronusQ::HashNAOs;
typedef libint2::Shell CShell;
typedef libint2::Shell LIShell;


void BasisSet::basisSetRead(FileIO * fileio, Molecule * mol){
  std::vector<double> coeff;
  std::vector<double> coeffP;
  std::vector<double> exp;
  std::vector<double> expP;
  std::array<double,3> center;
  std::string readString;
  fileio->in>>readString;
  std::string basis_path = "/" + readString;
  basis_path.insert(0,BASIS_PATH);
  std::unique_ptr<ifstream> fileBasis(new ifstream (basis_path));
  if(fileBasis->fail()){ // Check if file is in BASIS_PATH
    fileBasis.reset();
    fileBasis = std::unique_ptr<ifstream>(new ifstream(readString));
    if(fileBasis->fail()) { // Check if file is in PWD
      CErr("Could not find basis set file",fileio->out);
    } else {
      fileio->out << "Reading Basis Set from:\n ./" << readString<< endl;
    }
  } else {
    fileio->out << "Reading Basis Set from:" << endl << basis_path<< endl;
  }

  std::string atomStr;
  int readNPGTO,L,Ls,Lp;
  auto count=0;
  auto temp=0;
  double readNorm,expval,coeffD,coeffDp;
  std::string coefval, coefvalp;
  std::string nameOfAtom;
  struct basisSet_libint{
    std::string atomName;
    std::vector<libint2::Shell> refShell;
  };
  std::vector <basisSet_libint> allBasis;
  std::vector <libint2::Shell> refShell;
  center={{0,0,0}};
  while (!(fileBasis->eof())){
    while (readString.compare("****")){
      *fileBasis>>readString;
    };
    *fileBasis >> readString;
    nameOfAtom =readString;
    *fileBasis >> readString; 
    *fileBasis >> readString;
    while(readString.compare("****")){
      *fileBasis >> readNPGTO;   
      *fileBasis >> readNorm; 
      for (auto j=0;j<readNPGTO;++j){
        if (readString.size()==1){
          L=HashL(readString);
	  *fileBasis >> expval;
	  *fileBasis >> coefval;
	  exp.push_back(expval);

	  
	  auto e(coefval.find_first_of("Dd"));
	  if (e!=std::string::npos)
	    coefval[e]='E';
          coeffD=std::stod(coefval);
	  coeff.push_back(coeffD);  	
        }
        else if(readString.size()==2){
          Ls=HashL(readString.substr(0,1));
          Lp=HashL(readString.substr(1,1));
	  *fileBasis >> expval;
	  *fileBasis >> coefval;
	  *fileBasis >> coefvalp;
	  exp.push_back (expval);
	  expP.push_back(expval);
	  
	  auto e(coefval.find_first_of("Dd"));
	  if (e!=std::string::npos)
	    coefval[e]='E';
          coeffD=std::stod(coefval);
	  
	  auto f(coefvalp.find_first_of("Dd"));
	  if (f!=std::string::npos)
	    coefvalp[f]='E';
          coeffDp=std::stod(coefvalp);
	  
	  coeff.push_back (coeffD);
	  coeffP.push_back(coeffDp);
      }
     }
      if (readString.size()==1){
        refShell.push_back(
        LIShell{
          exp,
	  {
            {L, false, coeff}
	  },
	    center
          }
        );
        coeff.resize(0);
        exp.resize(0);
      }
      else if(readString.size()==2){
        refShell.push_back(
        LIShell{
          exp,
	  {
            {Ls, false, coeff}
	  },
	  center
          }
        );
        refShell.push_back(
        LIShell{
          expP,
	  {
            {Lp, false, coeffP}
	  },
	  center
          }
        ); 
        coeff.resize(0);
        exp.resize(0);
        expP.resize(0);
        coeffP.resize(0);
      };
     *fileBasis >> readString;
    }
    allBasis.push_back({nameOfAtom,refShell});
    refShell.resize(0);
  }
  count=0;
  for (auto i=0; i<mol->nAtoms();++i){
   atomStr=atom[mol->index(i)].symbol;
   center = {{(*mol->cart())(0,i),
              (*mol->cart())(1,i),
              (*mol->cart())(2,i)}};
   for(auto k=0; k<allBasis.size();++k){
     if (atomStr.compare(allBasis[k].atomName)==0){
       for (auto j=0;j<allBasis[k].refShell.size();++j){
         this->shells_libint.push_back(LIShell{allBasis[k].refShell[j].alpha,{allBasis[k].refShell[j].contr[0]},center});
         shells_libint_unnormal.push_back(LIShell{allBasis[k].refShell[j].alpha,{allBasis[k].refShell[j].contr[0]},center});
	 int L=allBasis[k].refShell[j].contr[0].l;
	 if (this->nLShell_.size()==L) this->nLShell_.push_back(0);
         temp=(this->nLShell_)[L];
         (this->nLShell_)[L]=temp+1;
	 this->shells_libint[count].renorm();
	 count++;
	 readNPGTO=allBasis[k].refShell[j].alpha.size();
	 this->nBasis_=(this->nBasis_)+HashNAOs(L);
	 this->atomNum.push_back(i);
         (this->nPrimitive_)=(this->nPrimitive_)+readNPGTO*HashNAOs(L);
         if(count==1) {
           this->maxPrim = readNPGTO;
           this->maxL = L;
         }   
	 else {
           if(readNPGTO > this->maxPrim)this->maxPrim = readNPGTO;
           if(L > this->maxL) this->maxL = L;
         }
       }
       break;
     }
     if (k==allBasis.size()-1){
       fileio->out<<"Error: unrecognized shell symbol!"<<endl;
       CErr("Unrecognized Shell symbol");
     }
   }
  }
  this->convToLI=true;
  fileBasis->close();
}

void BasisSet::makeMap(Molecule *  mol) {
  if(!this->convToLI) this->convShell(mol);
  int n = 0;
  for( auto shell: this->shells_libint) {
    this->mapSh2Bf.push_back(n);
    n += shell.size();
  }
  this->haveMap = true;
}

void BasisSet::computeShBlkNorm(Molecule * mol, const RealMatrix *D){
  // This will be much easier in Eigen
  //if(!this->convToLI) this->convShell(mol);
  if(!this->haveMap)  this->makeMap(mol);
  int nOfShell=this->shells_libint.size();
  this->shBlkNorm = std::unique_ptr<RealMatrix>(new RealMatrix(nOfShell,nOfShell));
  for(int s1 = 0; s1 < nOfShell; s1++) {
    int bf1 = this->mapSh2Bf[s1];
    int n1  = this->shells_libint[s1].size();
    for(int s2 = 0; s2 < nOfShell; s2++) {
      int bf2 = this->mapSh2Bf[s2];
      int n2  = this->shells_libint[s2].size();
     
      (*this->shBlkNorm)(s1,s2) = D->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
    }
  }
  
}

//---------------------------------------//
// print a general basis set information //
//---------------------------------------//
void BasisSet::printInfo_libint(FileIO * fileio,Controls * controls) {
  fileio->out<<"\nBasis Function Information:"<<endl;
  for(auto i=0;i<this->nLShell_.size();i=i+2){
    fileio->out <<std::setw(15) << "n" <<HashS(i)<<"Shell =" << std::setw(8) << this->nLShell_[i] << std::setw(5) << " "
             << std::setw(20) << "n" <<HashS(i+1)<<"Shell =" << std::setw(8) << this->nLShell_[i+1] << endl;
  }
  if(controls->printLevel>=2) printAtomO(fileio->out);
  //if(controls->printLevel>=3) printShellPair(fileio->out);
  
};
void BasisSet::printAtomO(ostream &output){
    
  output.precision(6);
  output.fill(' ');
  output.setf(ios::right,ios::adjustfield);
  output.setf(ios::fixed,ios::floatfield);
  output<<endl<<"Cartesian Basis Functions:"<<endl;
  output<<bannerTop<<endl;
  output<<std::setw(16)<<"        Shell Type"<<"    Center"<<std::setw(8)<<"L"<<std::setw(26)<<"exponent"<<std::setw(18)<<"coefficient"<<endl;
  output<<bannerMid << endl;
  for (auto i=0;i<this->shells_libint.size();i++){
    
    output<<std::setw(5)<<i+1<<std::setw(8)<<HashS(this->shells_libint[i].contr[0].l)<<std::setw(13)<<this->atomNum[i]+1<<"  ";
    output<<std::setw(8)<<this->shells_libint[i].contr[0].l;
      for(auto j=0;j<this->shells_libint[i].alpha.size();j++) {
	if(j!=0) output<<std::setw(36)<<" ";
	output<<std::setw(26)<<this->shells_libint[i].alpha[j]<<std::setw(18)<<shells_libint_unnormal[i].contr[0].coeff[j]<<endl;
      };
      output << endl;
  } 
  output << bannerEnd <<endl;
}
#endif

