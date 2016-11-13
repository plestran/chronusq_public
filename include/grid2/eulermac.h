/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
#ifndef INCLUDED_EULERMAC_H
#define INCLUDED_EULERMAC_H

#include <grid2/odgrid.h>

namespace ChronusQ {
class EulerMac : public OneDGrid2 {
  public:
    EulerMac(size_t npts, double screenTol = 0.0, bool onTheFly = true) : OneDGrid2(npts,screenTol,onTheFly){
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
