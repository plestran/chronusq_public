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

  std::string fNameUpper = boost::to_upper_copy<std::string>(fName);
  auto isGBS = fName.find(".gbs");
  auto isKey = this->basisKey.find(fNameUpper);

  if(isGBS != std::string::npos) {
    tmpStr = "/" + fName;
    tmpStr.insert(0,BASIS_PATH);
  } else if(isKey != this->basisKey.end()){
    tmpStr = this->basisMap[this->basisKey[fNameUpper]];
  }
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
            libint2::Shell{ exp, {{0,false,contPrimary}}, {{0,0,0}} }
          );
          tmpShell.push_back(
            libint2::Shell{ exp, {{1,false,contSecondary}}, {{0,0,0}} }
          );
        } else {
          int L = HashL(shSymb);
          bool doSph = (L > 1);
          if(this->forceCart_) doSph = false;

          tmpShell.push_back(
            libint2::Shell{ exp, {{L,doSph,contPrimary}}, {{0,0,0}} }
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
  bool time = false;
  std::chrono::high_resolution_clock::time_point start_1;
  std::chrono::high_resolution_clock::time_point start_2;
  std::chrono::high_resolution_clock::time_point start_3;
  std::chrono::high_resolution_clock::time_point start_4;
  std::chrono::high_resolution_clock::time_point start_5;
  std::chrono::high_resolution_clock::time_point finish_1;
  std::chrono::high_resolution_clock::time_point finish_2;
  std::chrono::high_resolution_clock::time_point finish_3;
  std::chrono::high_resolution_clock::time_point finish_4;
  std::chrono::high_resolution_clock::time_point finish_5;
//T
  if (time) {
     start_1 = std::chrono::high_resolution_clock::now();
     }
//T
  auto shSize = liShell.size(); 
  auto contDepth = liShell.alpha.size(); 
  auto center = liShell.O;
//T
  if (time) {
  finish_1 = std::chrono::high_resolution_clock::now();  
  this->duration_1 += finish_1 - start_1;
  start_4 = std::chrono::high_resolution_clock::now();
  }
//T
  double * fEVal = new double[shSize];
/*
  std::vector<std::array<int,3>> L(shSize);
  if(liShell.contr[0].l == 0){
    L[0] = {0,0,0};
  } else if(liShell.contr[0].l == 1){
    L[0] = {1,0,0};
    L[1] = {0,1,0};
    L[2] = {0,0,1};
  } else if(liShell.contr[0].l == 2){
    L[0] = {2,0,0};
    L[1] = {1,1,0};
    L[2] = {1,0,1};
    L[3] = {0,2,0};
    L[4] = {0,1,1};
    L[5] = {0,0,2};
  } else CErr("L > 2 NYI");
*/

  if (time) {
  finish_4 = std::chrono::high_resolution_clock::now();  
  this->duration_4 += finish_4 - start_4;
  start_5 = std::chrono::high_resolution_clock::now();
  }
  std::memset(fEVal,0,shSize*sizeof(double));
//T
  if (time) {
    finish_5 = std::chrono::high_resolution_clock::now();  
    this->duration_5 += finish_5 - start_5;
    start_2 = std::chrono::high_resolution_clock::now();
    }
//T
  double x = bg::get<0>(*pt) - center[0];
  double y = bg::get<1>(*pt) - center[1];
  double z = bg::get<2>(*pt) - center[2];
  double rSq = x*x + y*y + z*z;
//T
  if (time) {
    finish_2 = std::chrono::high_resolution_clock::now();  
    this->duration_2 += finish_2 - start_2;
    start_3 = std::chrono::high_resolution_clock::now();
    }
//T
/*
  for(auto i = 0; i < shSize; i++){
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
*/

  double expFactor = 0.0;
  for(auto k = 0; k < contDepth; k++){
    expFactor += 
      liShell.contr[0].coeff[k] *
      std::exp(-liShell.alpha[k]*rSq);
  }

  if(liShell.contr[0].l == 0){
    fEVal[0] = expFactor;
  }else if(liShell.contr[0].l == 1){
    fEVal[0] = expFactor*x;
    fEVal[1] = expFactor*y;
    fEVal[2] = expFactor*z;
  } else if(liShell.contr[0].l == 2){
    fEVal[0] = expFactor*x*x;
    fEVal[1] = expFactor*y*x;
    fEVal[2] = expFactor*z*x;
    fEVal[3] = expFactor*y*y;
    fEVal[4] = expFactor*y*z;
    fEVal[5] = expFactor*z*z;
  }
//T
  if (time) {
    finish_3 = std::chrono::high_resolution_clock::now();  
    this->duration_3 += finish_3 - start_3;
  }
//T

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


std::vector<bool> BasisSet::MapGridBasis(cartGP pt){
//  Set map_[ishell] to be avaluated (true) or not (false)
//  note: radCutSh_ has to be already populated by calling before radcut
//bool * map_ = new bool[this->nShell()+1];
  std::vector<bool> map_(this->nShell()+1);
  double x ;
  double y ;
  double z ;
  double r ;
  bool   nodens = true;  //becomes truee if at least one shell hs to evaluated and it is stored in map_[0]
  for(auto s1=0l; s1 < this->nShell(); s1++){  //loop over shells
    auto center = shells(s1).O;
    x = bg::get<0>(pt) - center[0];
    y = bg::get<1>(pt) - center[1];
    z = bg::get<2>(pt) - center[2];
    r = std::pow((x*x + y*y + z*z),(0.5));
    map_[s1+1] = false;        
    if (r < this->radCutSh_[s1]) {
      map_[s1+1] = true;
      nodens = false;
      }
    } //End loop over shells
  map_[0] = nodens;
  return map_;
}


void BasisSet::radcut(double thr, int maxiter, double epsConv){
  this->radCutSh_ = new double[this->nShell()];
  double alphaMin;
//  double *s1Eval = basisEval(s1,pt);
  for(auto s1=0l; s1 < this->nShell(); s1++){
//    find smallest alpha coeff for each shell
      auto contDepth = this->shells(s1).alpha.size(); 
      alphaMin = 1.0e15;
      for(auto k = 0; k < contDepth; k++){
        if (this->shells(s1).alpha[k] <= alphaMin){
        alphaMin = this->shells(s1).alpha[k];
        }
      }
//       cout << "s1 " << s1 << endl;
//       this->fSpAv (2, shells(s1).contr[0].l, alphaMin, 1.0e-5);
//       this->fSpAv (1, shells(s1).contr[0].l, alphaMin, 3);
//     Populate a Vector storing all the cut off radius (Av_xi(r_cut)<thr)
       radCutSh_[s1] = this->fRmax (shells(s1).contr[0].l, alphaMin, thr, epsConv, maxiter);

  }
  return ;
}

double BasisSet::fRmax (int l, double alpha, double thr, double epsConv, int maxiter){
  double root ;
  double root1 ;
       root =  fSpAv (2, l,alpha, thr);
  for (auto i=0; i < maxiter; i++){
       root1  =  - (this->fSpAv(0, l,alpha, root) - thr);
       root1 /=  this->fSpAv (1, l,alpha, root);
       root1 +=  root;
       if(std::abs(root1-root) <= epsConv){
//       cout << "l "<< l << " alpha " << alpha <<endl;
//       cout << "root(n-1)= " << root  << " root(n)= "<<root1 <<" abs_err " << std::abs(root1-root)  << endl;
//       cout << "Root found " << root1 << " It " << i << " froot " << this->fSpAv(0, l,alpha, root) << endl;
         return root1;
       }else{     
           root = root1;
    }
   }
           this->fileio_->out << "Convergence Failure in fRmax, change maxiter or turn off screening " << endl;    
           this->fileio_->out << "root(n-1)= " << root  << " root(n)= "<<root1 <<" abs_err " << std::abs(root1-root)  << endl;
           CErr("Convergence Failure",this->fileio_->out);
}   

double BasisSet::fSpAv (int iop, int l, double alpha, double r){
       double fAv = 0.0;
       double two = 2.0;
       double threeOv2 = 1.5;
       double oneOv2 = 0.5;
       fAv = std::pow((two*alpha),(l+threeOv2)) ;
       fAv /= two*math.pi*boost::math::tgamma(l+threeOv2); 
       fAv  = std::pow(fAv,(oneOv2)) ;

       if (iop == 0){
       fAv *= std::exp(-alpha*r*r) ;
       fAv *= std::pow(r,l) ;
       }else if(iop ==1){
       fAv *= std::exp(-alpha*r*r);
       fAv *= std::pow(r,(l-1)) ;
       fAv *= (l - two*alpha*r*r);
       }else if(iop = 2){
       fAv = -std::log(fAv);
       fAv += std::log(r);
       fAv /= -alpha;
       fAv = std::pow(fAv,0.5);
      }   
//       cout << "l "<< l << " alpha " << alpha << " fAv " << fAv <<endl;
       return fAv;
     }

double * BasisSet::basisonFlyProdEval(libint2::Shell s1, int s1size, libint2::Shell s2, int s2size, double rx, double ry, double rz){
          
  double *fEVal = new double[s1size*s2size];

// Evaluate shell pair new center and i think the distance of r from the new center
// Evaluate the shell pair exponential factor (is going to be the same no matter the angular momentum)
// build the angular momentum part and ad in f[ij]

/*  
  auto contDepth1 = s1.alpha.size(); 
  auto center1 = s1.O;
  auto contDepth2 = s2.alpha.size(); 
  auto center2 = s2.O;
  double x1 = rx - center1[0];
  double y1 = ry - center1[1];
  double z1 = rz - center1[2];
  double rSq1 = x1*x1 + y1*y1 + z1*z1;
  double x2 = rx - center2[0];
  double y2 = ry - center2[1];
  double z2 = rz - center2[2];
  double rSq2 = x2*x2 + y2*y2 + z2*z2;

  double expFactor1 = 0.0;
  double expFactor2 = 0.0;
  for(auto k1 = 0; k1 < contDepth1; k1++){
    expFactor1 += 
      s1.contr[0].coeff[k1] *
      std::exp(-s1.alpha[k1]*rSq1);
  }

  for(auto k2 = 0; k2 < contDepth2; k2++){
    expFactor2 += 
      s2.contr[0].coeff[k2] *
      std::exp(-s2.alpha[k2]*rSq2);
  }
*/
/*
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
*/  
  return fEVal;
  
}


void BasisSet::popExpPairSh(){
  this->expPairSh_ = new double[(this->nShell()*(this->nShell()+1))/2];
     for(auto s1 = 0, ij = 0; s1 < this->nShell(); s1++){
        auto contDepth1 = this->shells(s1).alpha.size(); 
        auto center1 = this->shells(s1).O;
        for(auto s2 = s1; s2 < this->nShell(); s2++, ij++){
          auto contDepth2 = this->shells(s2).alpha.size(); 
          auto center2 = this->shells(s2).O;
          double rpx = center1[0] - center2[0];
          double rpy = center1[1] - center2[1];
          double rpz = center1[2] - center2[2];
          double rpsq = rpx*rpx + rpy*rpy + rpz*rpz;
                 for(auto k1 = 0; k1 < contDepth1; k1++){
                    for(auto k2 = 0; k2 < contDepth2; k1++){
                      this->expPairSh_[ij] = this->shells(s1).contr[0].coeff[k1]*this->shells(s2).contr[0].coeff[k2];
                      this->expPairSh_[ij] *= std::exp(-this->shells(s1).alpha[k1]-this->shells(s2).alpha[k2]*rpsq);
                 }
               }
             }
           }
  return ;
}

template<>
double * BasisSet::basisDEval(int iop, libint2::Shell &liShell, cartGP *pt){
// IOP Derivative

  if (iop >1) CErr("Derivative order NYI in basisDEval");
  auto shSize = liShell.size(); 
  auto contDepth = liShell.alpha.size(); 
  auto center = liShell.O;
  double * fEVal = new double[(3*iop + 1)*shSize];
//double * DfEval = new double[shSize];
  double * f = fEVal;
  double * DfEval = f + shSize;

  double * dx = DfEval;
  double * dy = dx + shSize;
  double * dz = dy + shSize;

  double x = bg::get<0>(*pt) - center[0];
  double y = bg::get<1>(*pt) - center[1];
  double z = bg::get<2>(*pt) - center[2];
  double rSq = x*x + y*y + z*z;
  double alpha = 0.0;
  double expFactor = 0.0;

  for(auto k = 0; k < contDepth; k++){
    expFactor += 
      liShell.contr[0].coeff[k] *
      std::exp(-liShell.alpha[k]*rSq);
      if (iop == 1) alpha += 
      2.0*liShell.alpha[k]*liShell.contr[0].coeff[k] *
      std::exp(-liShell.alpha[k]*rSq);
  }
 
  if(liShell.contr[0].l == 0){
    f[0] = expFactor;
   if (iop == 1) {
      dx[0] =  alpha*x;
      dy[0] =  alpha*y;
      dz[0] =  alpha*z;
      }
  }else if(liShell.contr[0].l == 1){
    f[0] = expFactor*x;
    f[1] = expFactor*y;
    f[2] = expFactor*z;
    if (iop == 1) {
//    dpx/dxi
      dx[0] = alpha*x*x-expFactor;
      dy[0] = alpha*x*y;
      dz[0] = alpha*x*z;
//    dpy/dxi
      dx[1] = alpha*y*x;
      dy[1] = alpha*y*y - expFactor;
      dz[1] = alpha*y*z;
//    dpz/dxi
      dx[2] = alpha*z*x;
      dy[2] = alpha*z*y;
      dz[2] = alpha*z*z - expFactor;
    }
  } else if(liShell.contr[0].l == 2){
    f[0] = expFactor*x*x;
    f[1] = expFactor*y*x;
    f[2] = expFactor*z*x;
    f[3] = expFactor*y*y;
    f[4] = expFactor*y*z;
    f[5] = expFactor*z*z;
    if (iop == 1) {
//    dDx^2/dxi
      dx[0] = alpha*x*x*x-2.0*x*expFactor;
      dy[0] = alpha*x*x*y;
      dz[0] = alpha*x*x*z;
//    dDyx/dxi
      dx[1] = alpha*x*x*y-y*expFactor;
      dy[1] = alpha*x*y*y-x*expFactor;
      dz[1] = alpha*x*y*z;
//    dDzx/dxi
      dx[2] = alpha*x*x*z-z*expFactor;
      dy[2] = alpha*x*y*z;
      dz[2] = alpha*x*z*z-x*expFactor;
//    dDy^2/dxi
      dx[3] = alpha*y*y*x;
      dy[3] = alpha*y*y*y-2.0*y*expFactor;
      dz[3] = alpha*y*y*z;
//    dDyz/dxi
      dx[4] = alpha*x*y*z;
      dy[4] = alpha*y*y*z-z*expFactor;
      dz[4] = alpha*y*z*z-y*expFactor;
//    dDz^2/dxi
      dx[5] = alpha*z*z*x;
      dy[5] = alpha*z*z*y;
      dz[5] = alpha*z*z*z-2.0*z*expFactor;
    }
  }
  return fEVal;
}

} // namespace ChronusQ

