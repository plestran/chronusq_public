#ifndef INCLUDED_ODGRID_H
#define INCLUDED_ODGRID_H

#include <grid2/grid2_def.h>

namespace ChronusQ {

class OneDGrid2 : public Grid2 {
  public:
    OneDGrid2(size_t npts = 0, double screenTol = 0.0, bool onTheFly = true) : Grid2(npts,screenTol,onTheFly){ };
    virtual ~OneDGrid2(){ };
    virtual IntegrationPoint operator[](size_t) = 0;
};
}
#endif
