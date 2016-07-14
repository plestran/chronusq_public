#ifndef INCLUDED_EULERMAC_H
#define INCLUDED_EULERMAC_H

#include <grid2/odgrid.h>

namespace ChronusQ {
class EulerMac : public OneDGrid2 {
  public:
    EulerMac(size_t npts, bool onTheFly = true) : OneDGrid2(npts,onTheFly){
        if(!onTheFly_) this->generateGridPoints();
    };
    ~EulerMac(){ };
    IntegrationPoint operator[](size_t i){ 
      if(!onTheFly_ && haveGPs_) 
        return IntegrationPoint(gPoints_[i],weights_[i]);

/*
 IMplemented as Eq 6. and 7 page 1001.
Christopher W. Murray, Nicholas C. Handy & Gregory J. Laming (1993):
Quadrature schemes for integrals of density functional theory, Molecular Physics: An
International Journal at the Interface Between Chemistry and Physics, 78:4, 997-1014 
Note that the Jacobian r^2 is not included (to added later on in the weights) and
also the equeations are mapped from 0,inf respect to the paper according the following
substitution 1 -> N+1, q -> i and w -> w * (N+1).
The resulting equation including the Jacobian with m=2 are eq 24 and 25 in:
Development, implementation and applications of efficient
methodologies for density functional calculations
Johnson 169 Theor. Comp. Chem. Vol 2 (Jacobian included in the weight) 

 // m=1
      double pt  = (i + 1.0) / 
        ((this->nPts_ - i));
      double wgt = (this->nPts_ + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i) );
*/
// m= 2
      double pt  = (i + 1.0) * (i + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i));
      double wgt = 2.0 * (i + 1.0) * (this->nPts_ + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i) * (this->nPts_ - i));
/*
// m = 3
      double pt  = (i + 1.0) * (i + 1.0) * (i + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i) * (this->nPts_ - i) );
      double wgt = 3.0 * (i + 1.0) * (i + 1.0) * (this->nPts_ + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i) * (this->nPts_ - i) * (this->nPts_ - i) );
*/
      return IntegrationPoint(pt,wgt);
    };
};
}
#endif
