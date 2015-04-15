#include "basisset.h"

#ifdef USE_LIBINT

using ChronusQ::BasisSet;
using ChronusQ::Molecule;
typedef ChronusQ::Shell CShell;
typedef libint2::Shell LIShell;

void BasisSet::convShell(Molecule* mol) {
  std::vector<double> coeff;
  std::vector<double> exp;
  std::array<double,3> center;
  for(auto i = 0; i < this->nShell_; ++i) {
    for(auto iPrim = 0; iPrim < this->shells[i].nPGTOs; ++iPrim){
      coeff.push_back(this->shells[i].coef[iPrim]);
      exp.push_back(this->shells[i].expo[iPrim]);
    }
    center = {(*mol->cart())(0,this->shells[i].center),
              (*mol->cart())(1,this->shells[i].center),
              (*mol->cart())(2,this->shells[i].center)};

    int L = ChronusQ::HashL(&this->shells[i].name[0]);

    this->shells_libint.push_back(
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
    shells_libint[i].renorm();
    if(i==0) {
      this->maxPrim = this->shells[i].nPGTOs;
      this->maxL = L;
    } else {
      if(this->shells[i].nPGTOs > this->maxPrim)this->maxPrim = this->shells[i].nPGTOs;
      if(L > this->maxL) this->maxL = L;
    }
  }
  this->convToLI = true;
}

void BasisSet::makeMap(Molecule * mol) {
  if(!this->convToLI) this->convShell(mol);
  int n = 0;
  for( auto shell: this->shells_libint) {
    this->mapSh2Bf.push_back(n);
    n += shell.size();
  }
  this->haveMap = true;
}

void BasisSet::computeShBlkNorm(Molecule *mol, Matrix<double> *D){
  // This will be much easier in Eigen
  if(!this->convToLI) this->convShell(mol);
  if(!this->haveMap)  this->makeMap(mol);

  this->shBlkNorm = new Matrix<double>(this->nShell(),this->nShell());
  Matrix<double> *tmp;
  for(int s1 = 0; s1 < this->nShell(); s1++) {
    int bf1 = this->mapSh2Bf[s1];
    int n1  = this->shells_libint[s1].size();
    for(int s2 = 0; s2 < this->nShell(); s2++) {
      int bf2 = this->mapSh2Bf[s2];
      int n2  = this->shells_libint[s2].size();
     
      // Grab submat
      tmp = new Matrix<double>(n1,n2);
      tmp->clearAll();
      for(int i = 0; i < n1; i++) {
        for(int j = 0; j < n2; j++) {
           (*tmp)(i,j) = (*D)(bf1+i,bf2+j);
        } 
      } 
      (*this->shBlkNorm)(s1,s2) = tmp->infNorm();
      delete tmp;
    }
  }
}
#endif
