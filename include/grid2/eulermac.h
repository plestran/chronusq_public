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

      double pt  = (i + 1.0) * (i + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i));
      double wgt = 2.0 * (i + 1.0) * (this->nPts_ + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i) * (this->nPts_ - i));

      return IntegrationPoint(pt,wgt);
    };
};
}
#endif
