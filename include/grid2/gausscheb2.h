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
#ifndef INCLUDED_GAUSSCHEBTWO_H
#define INCLUDED_GAUSSCHEBTWO_H

#include <grid2/odgrid.h>

namespace ChronusQ {
/**
 * Gauss-Chebyshev Quadrature of the second kind
 *
 * \int _{-1}^{1} f(x)\sqrt{1 - x^2} \mathrm{d} x = \sum_n w_i f(x_i)
 *
 * x_i = \cos\left( \frac{i}{n + 1} \pi \right)
 * w_i = \frac{\pi}{n+1} \sin^2\left( \frac{i}{n+1} \pi \right)
 *
 * where "n" is the number of integration points
 *
 */
class GaussChebSnd : public OneDGrid2 {
  public:
    GaussChebSnd(size_t npts = 0, double screenTol = 0., bool onTheFly = true) : 
      OneDGrid2(npts,screenTol,onTheFly){
        if(!onTheFly_) this->generateGridPoints();
      };

    ~GaussChebSnd(){ };
    IntegrationPoint operator[](size_t i){ 
      if(!onTheFly_ && haveGPs_) 
        return IntegrationPoint(gPoints_[i],weights_[i]);

      // Generate points and weights for (-1,1)
      double pt = std::cos(math.pi * i / (this->nPts_ +1));
      double wgt = math.pi / (this->nPts_ + 1) * 
        std::pow(
            std::sin(math.pi * i / (this->nPts_ +1)),
            2.0
        );

      // Map (-1,1) -> (0,\inf)
        
      // As we're integrating f(x) not f(x)\sqrt{1 - x^2}
      // we must factor the \frac{1}{\sqrt{1-x^2}} into the weights
      //
      // Also, we must factor in the Jacobian of the domain shift into
      // the weights 
      //
      // dr = \frac{2}{(1-x)^2} dx
      //
      // Both of these factors lead to
      //
      // w_i \rightarrow w_i(x_i) = \frac{2\sqrt{1-x^2}}{(1-x)^2}
      //
      wgt *= 2.0 / (std::sqrt(1 - pt * pt) * (1 - pt) * (1 - pt) );


      // Perform the coordinate shift
      //
      // r_i = \frac{(1+x_i)}{(1-x_i)}
      pt   = (1 + pt) / (1 - pt);

      return IntegrationPoint(pt,wgt);
    };
};
}
#endif
