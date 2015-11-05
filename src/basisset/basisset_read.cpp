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
namespace ChronusQ{
/**
 *  Parse a basis file and construct reference and local basis set defintion  
 */
void BasisSet::basisSetRead(FileIO * fileio, Molecule * mol, Controls *controls){

  std::string readString;
  int nTCS = 1;
  if(controls->doTCS) nTCS = 2;
  
  //this->fileio_ = fileio;
  this->communicate(*fileio);
 
  this->fileio_->in >> readString; // read the name of the basis set file
  this->findBasisFile(readString); // Try to find the basis set file
  this->parseGlobal();
  this->constructLocal(mol);
  this->makeMaps(nTCS,mol);
  this->printInfo();
  this->renormShells();

}; // BasisSet::basisSetRead

/**
 * Attempt to find the file containing the basis set definition
 */
void BasisSet::findBasisFile(std::string fName){
  std::string tmpStr;

  tmpStr = "/" + fName;
  tmpStr.insert(0,BASIS_PATH);
  this->setBasisPath(tmpStr);

  this->basisFile_ = std::unique_ptr<ifstream>(new ifstream(this->basisPath_));

  
  if(!this->basisFile_->fail()){ // Check if file is in BASIS_PATH
//  this->fileio_->out << "Reading Basis Set from: " << this->basisPath_ << endl;
  } else {
    this->basisFile_.reset();
    this->setBasisPath(fName);
    this->basisFile_ = std::unique_ptr<ifstream>(new ifstream(this->basisPath_));
    if(!this->basisFile_->fail()){ // Check if file is in PWD
      this->fileio_->out << "Reading Basis Set from: ./" << this->basisPath_ << endl;
    } else CErr("Could not find basis set file \"" + fName + "\"");
  }
}; // BasisSet::findBasisFile

/**
 *  Parse the basis set file and generate a set of reference shells
 *  from which local and external basis set objects are constructed
 */
void BasisSet::parseGlobal(){

  std::string readString;
  std::string nameOfAtom;
  std::string shSymb;
  int         contDepth;
  int atomicNumber;
  int indx;
  std::vector<libint2::Shell> tmpShell;

  bool readRec = false;
  bool newRec  = false;
  bool firstRec = true;
  int nEmpty = 0;
  int nComm  = 0;
  int nRec   = 0;

  while(!this->basisFile_->eof()){
    std::getline(*this->basisFile_,readString);
    if(readString.size() == 0)    nEmpty++;
    else if(readString[0] == '!') nComm++;
    else if(!readString.compare("****")){
      std::getline(*this->basisFile_,readString);
      if(readString.size() == 0) { nEmpty++; readRec = false; continue;}
      nRec++;
      readRec = true;
      newRec  = true;
    }

    if(readRec){
      std::istringstream iss(readString);
      std::vector<std::string> tokens(std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>{});
      if(newRec){
        if(!firstRec) {
          this->refShells_.push_back(ReferenceShell{atomicNumber,indx,tmpShell});
        }
        indx = HashAtom(tokens[0],0);
        atomicNumber = elements[indx].atomicNumber;
        newRec = false;
        firstRec = false;
        tmpShell.clear();
      } else {
        contDepth = std::stoi(tokens[1]);
        shSymb    = tokens[0];
        std::vector<double> exp;
        std::vector<double> contPrimary;
        std::vector<double> contSecondary;

        for(auto i = 0; i < contDepth; i++) {
          std::getline(*this->basisFile_,readString);
          std::istringstream iss2(readString);
          std::vector<std::string> tokens2(std::istream_iterator<std::string>{iss2},
            std::istream_iterator<std::string>{});

          exp.push_back(std::stod(tokens2[0]));
          contPrimary.push_back(std::stod(tokens2[1]));
          if(!shSymb.compare("SP"))
            contSecondary.push_back(std::stod(tokens2[2]));
        }

        if(!shSymb.compare("SP")) {
          tmpShell.push_back(
            libint2::Shell{ exp, {{0,this->doSph_,contPrimary}}, {{0,0,0}} }
          );
          tmpShell.push_back(
            libint2::Shell{ exp, {{1,this->doSph_,contSecondary}}, {{0,0,0}} }
          );
        } else {
          tmpShell.push_back(
            libint2::Shell{ exp, {{HashL(shSymb),this->doSph_,contPrimary}}, {{0,0,0}} }
          );
        }
      }
    }
  }
  // Append the last Rec
  this->refShells_.push_back(ReferenceShell{atomicNumber,indx,tmpShell}); 
 
}; // BasisSet::parseGlobal

template<>
double * BasisSet::basisEval(int iShell, std::array<double,3> center, cartGP *pt){
  auto shSize = this->shells(iShell).size(); 
  auto contDepth = this->shells(iShell).alpha.size(); 
  double * fEVal = new double[shSize];
  
  std::vector<std::array<int,3>> L;
  if(this->shells(iShell).contr[0].l == 0){
    L.push_back({{0,0,0}});
  } else if(this->shells(iShell).contr[0].l == 1){
    L.push_back({{1,0,0}});
    L.push_back({{0,1,0}});
    L.push_back({{0,0,1}});
  } else if(this->shells(iShell).contr[0].l == 2){
    L.push_back({{2,0,0}});
    L.push_back({{1,1,0}});
    L.push_back({{1,0,1}});
    L.push_back({{0,2,0}});
    L.push_back({{0,1,1}});
    L.push_back({{0,0,2}});
  } else CErr("L > 2 NYI");

  std::memset(fEVal,0,shSize*sizeof(double));

  double x = bg::get<0>(*pt) - center[0];
  double y = bg::get<1>(*pt) - center[1];
  double z = bg::get<2>(*pt) - center[2];
  double rSq = x*x + y*y + z*z;
  for(auto i = 0; i < shSize; i++){
    cout << endl << fEVal[i] << endl;
    for(auto k = 0; k < contDepth; k++){
      fEVal[i] += 
        this->shells(iShell).contr[0].coeff[k] *
        std::exp(-this->shells(iShell).alpha[k]*rSq);
    }
    auto l = L[i][0];
    auto m = L[i][1];
    auto n = L[i][2];

    fEVal[i] *= std::pow(x,l);
    fEVal[i] *= std::pow(y,m);
    fEVal[i] *= std::pow(z,n);
  }

  return fEVal;
}
template<>
double * BasisSet::basisEval(int iShell, std::array<double,3> center,sph3GP *ptSph){
  cartGP pt;
  bg::transform(*ptSph,pt);
  auto shSize = this->shells(iShell).size(); 
  auto contDepth = this->shells(iShell).alpha.size(); 
  double * fEVal = new double[shSize];
  
  std::vector<std::array<int,3>> L;
  if(this->shells(iShell).contr[0].l == 0){
    L.push_back({{0,0,0}});
  } else if(this->shells(iShell).contr[0].l == 1){
    L.push_back({{1,0,0}});
    L.push_back({{0,1,0}});
    L.push_back({{0,0,1}});
  } else if(this->shells(iShell).contr[0].l == 2){
    L.push_back({{2,0,0}});
    L.push_back({{1,1,0}});
    L.push_back({{1,0,1}});
    L.push_back({{0,2,0}});
    L.push_back({{0,1,1}});
    L.push_back({{0,0,2}});
  } else CErr("L > 2 NYI");

  std::memset(fEVal,0,shSize*sizeof(double));
  cout << "Point Cart " << bg::get<0>(pt) << " " << bg::get<1>(pt) << " " << bg::get<2>(pt) <<endl;
  cout << "Center     " << center[0] << " " << center[1] << " " << center[2] <<endl;
  double x = bg::get<0>(pt) - center[0];
  double y = bg::get<1>(pt) - center[1];
  double z = bg::get<2>(pt) - center[2];
  cout << "Point Scaled" << x << " " << y << " " << z << endl;
  double rSq = x*x + y*y + z*z;
  cout << " rSq " << rSq << endl;
  cout << "shSize " << shSize << endl;
  cout << "contDepth " << contDepth << endl;
  for(auto i = 0; i < shSize; i++){
//    cout << endl << fEVal[i] << endl;
    for(auto k = 0; k < contDepth; k++){
      fEVal[i] += 
        this->shells(iShell).contr[0].coeff[k] *
        std::exp(-this->shells(iShell).alpha[k]*rSq);
//        cout << "AP " << this->shells(iShell).contr[0].coeff[k] << endl;
//        cout <<  this->shells(iShell).alpha[k] << endl;
    }
    auto l = L[i][0];
    auto m = L[i][1];
    auto n = L[i][2];
    cout << "l= " << l << "m= " << m << "n= " << n << endl;
    fEVal[i] *= std::pow(x,l);
    fEVal[i] *= std::pow(y,m);
    fEVal[i] *= std::pow(z,n);
  }

  return fEVal;
}
template<>
double * BasisSet::basisEval(libint2::Shell &liShell, cartGP *pt){
  auto shSize = liShell.size(); 
  auto contDepth = liShell.alpha.size(); 
  auto center = liShell.O;
  double * fEVal = new double[shSize];
  
  std::vector<std::array<int,3>> L;
  if(liShell.contr[0].l == 0){
    L.push_back({{0,0,0}});
  } else if(liShell.contr[0].l == 1){
    L.push_back({{1,0,0}});
    L.push_back({{0,1,0}});
    L.push_back({{0,0,1}});
  } else if(liShell.contr[0].l == 2){
    L.push_back({{2,0,0}});
    L.push_back({{1,1,0}});
    L.push_back({{1,0,1}});
    L.push_back({{0,2,0}});
    L.push_back({{0,1,1}});
    L.push_back({{0,0,2}});
  } else CErr("L > 2 NYI");

  std::memset(fEVal,0,shSize*sizeof(double));

  double x = bg::get<0>(*pt) - center[0];
  double y = bg::get<1>(*pt) - center[1];
  double z = bg::get<2>(*pt) - center[2];
  double rSq = x*x + y*y + z*z;
  for(auto i = 0; i < shSize; i++){
//    cout << endl << fEVal[i] << endl;
    for(auto k = 0; k < contDepth; k++){
      fEVal[i] += 
        liShell.contr[0].coeff[k] *
        std::exp(-liShell.alpha[k]*rSq);
    }
    auto l = L[i][0];
    auto m = L[i][1];
    auto n = L[i][2];

    fEVal[i] *= std::pow(x,l);
    fEVal[i] *= std::pow(y,m);
    fEVal[i] *= std::pow(z,n);
  }

  return fEVal;
}
template<>
double * BasisSet::basisEval(libint2::Shell &liShell, sph3GP *ptSph){
  cartGP pt;
  bg::transform(*ptSph,pt);
  auto shSize = liShell.size(); 
  auto contDepth = liShell.alpha.size(); 
  auto center = liShell.O;
  double * fEVal = new double[shSize];
  
  std::vector<std::array<int,3>> L;
  if(liShell.contr[0].l == 0){
    L.push_back({{0,0,0}});
  } else if(liShell.contr[0].l == 1){
    L.push_back({{1,0,0}});
    L.push_back({{0,1,0}});
    L.push_back({{0,0,1}});
  } else if(liShell.contr[0].l == 2){
    L.push_back({{2,0,0}});
    L.push_back({{1,1,0}});
    L.push_back({{1,0,1}});
    L.push_back({{0,2,0}});
    L.push_back({{0,1,1}});
    L.push_back({{0,0,2}});
  } else CErr("L > 2 NYI");

  std::memset(fEVal,0,shSize*sizeof(double));

  double x = bg::get<0>(pt) - center[0];
  double y = bg::get<1>(pt) - center[1];
  double z = bg::get<2>(pt) - center[2];
  double rSq = x*x + y*y + z*z;
  for(auto i = 0; i < shSize; i++){
//    cout << endl << fEVal[i] << endl;
    for(auto k = 0; k < contDepth; k++){
      fEVal[i] += 
        liShell.contr[0].coeff[k] *
        std::exp(-liShell.alpha[k]*rSq);
    }
    auto l = L[i][0];
    auto m = L[i][1];
    auto n = L[i][2];

    fEVal[i] *= std::pow(x,l);
    fEVal[i] *= std::pow(y,m);
    fEVal[i] *= std::pow(z,n);
//    cout << "inside " << i << "  " << fEVal[i] << endl;
  }

  return fEVal;
}

template<>
double * BasisSet::basisProdEval(libint2::Shell s1, libint2::Shell s2, cartGP *pt){
          
//  cout << "{" <<bg::get<0>(pt) << ", "<<bg::get<1>(pt)<<", " <<bg::get<2>(pt) <<"}, "<< endl;
  double *fEVal = new double[s1.size()*s2.size()];
  double *s1Eval = basisEval(s1,pt);
  double *s2Eval = basisEval(s2,pt);

  double   temp;
  double   temp2;
  double   zero = 0.0;
  for(auto i = 0, ij = 0; i < s1.size(); i++)
  for(auto j = 0; j < s2.size(); j++, ij++){
    fEVal[ij] = s1Eval[i]*s2Eval[j];
  }
  delete [] s1Eval;
  delete [] s2Eval;
  
  return fEVal;
  
}

template<>
double * BasisSet::basisProdEval(libint2::Shell s1, libint2::Shell s2, sph3GP *pt){
  double *fEVal = new double[s1.size()*s2.size()];
  double *s1Eval = basisEval(s1,pt);
  double *s2Eval = basisEval(s2,pt);

  for(auto i = 0, ij = 0; i < s1.size(); i++)
  for(auto j = 0; j < s2.size(); j++, ij++){
    fEVal[ij] = s1Eval[i]*s2Eval[j];
//    cout << "Print Inside =" << fEVal[ij] <<endl;
}
  delete [] s1Eval;
  delete [] s2Eval;
  
  return fEVal;
  
}

double radcut(int IAtom, double thr){
  double radius = 0.0;
//  double *s1Eval = basisEval(s1,pt);
//  for(auto i = 0; i < s1.size(); i++){
//     s1Eval[i] 
//  }
  return radius;
};


} // namespace ChronusQ

