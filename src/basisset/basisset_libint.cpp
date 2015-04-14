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
#endif
