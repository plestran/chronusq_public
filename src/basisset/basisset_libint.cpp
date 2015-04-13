#include "basisset.h"

#ifdef USE_LIBINT

using ChronusQ::BasisSet;
using ChronusQ::Molecule;
typedef ChronusQ::Shell CShell;
typedef libint2::Shell LIShell;

void BasisSet::convShell(Molecule* mol) {
  ChronusQ::Matrix<double> cart = (*mol->cart());
  std::vector<double> coeff;
  std::vector<double> exp;
  std::array<double,3> center;
  for(auto i = 0; i < this->nShell_; ++i) {
    for(auto iPrim = 0; iPrim < this->shells[i].nPGTOs; ++iPrim){
      coeff.push_back(this->shells[i].coef[iPrim]);
      exp.push_back(this->shells[i].expo[iPrim]);
    }
    center = {cart(0,this->shells[i].center),
              cart(1,this->shells[i].center),
              cart(2,this->shells[i].center)};
    int L = 0;
    switch(this->shells[i].name[0]) {
      case 'S':
        L = 0; 
	break;
      case 'P':
        L = 1; 
	break;
      case 'D':
        L = 2; 
	break;
      case 'F':
        L = 3; 
	break;
      case 'G':
        L = 4; 
	break;
      case 'H':
        L = 5; 
	break;
      case 'I':
        L = 6; 
	break;
      case 'J':
        L = 7; 
	break;
      case 'K':
        L = 8; 
	break;
      case 'L':
        L = 9; 
	break;
      case 'M':
        L = 10; 
	break;
      case 'N':
        L = 11; 
	break;
      case 'O':
        L = 12; 
	break;
      case 'Q':
        L = 13; 
	break;
      case 'R':
        L = 14; 
	break;
      case 'T':
        L = 15; 
	break;
      case 'U':
        L = 16; 
	break;
      case 'V':
        L = 17; 
	break;
      case 'W':
        L = 18; 
	break;
      case 'X':
        L = 19; 
	break;
      case 'Y':
        L = 20; 
	break;
      case 'Z':
        L = 21; 
	break;
      default:
        L = 0;
    }

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
//  shells_libint[i].renorm();
    cout << shells_libint[i];
  }

}
#endif
