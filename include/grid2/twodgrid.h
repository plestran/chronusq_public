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
#ifndef INCLUDED_TWODGRID_H
#define INCLUDED_TWODGRID_H

#include <grid2/grid2_def.h>
#include <grid2/odgrid.h>
#include <grid2/lebedev.h>
#include <grid2/gausscheb1.h>
#include <grid2/gausscheb2.h>
#include <grid2/eulermac.h>

namespace ChronusQ {
class TwoDGrid2 : public Grid2 {
  protected:
    std::unique_ptr<OneDGrid2> GRad;
    std::unique_ptr<OneDGrid2> GAng;
  public:
    TwoDGrid2(size_t nPtsRad, size_t nPtsAng, double screenTol, GRID_TYPE GTypeRad, 
        GRID_TYPE GTypeAng, bool onTheFly = true) : 
      Grid2(nPtsRad*nPtsAng,screenTol,true) {
      if(GTypeRad == GAUSSCHEBFST) {
        GRad = std::unique_ptr<OneDGrid2>(new GaussChebFst(nPtsRad,onTheFly));
      } else if(GTypeRad == GAUSSCHEBSND) {
        GRad = std::unique_ptr<OneDGrid2>(new GaussChebSnd(nPtsRad,onTheFly));
      } else if(GTypeRad == EULERMAC) {
        GRad = std::unique_ptr<OneDGrid2>(new EulerMac(nPtsRad,onTheFly));
      };

      if(GTypeAng == LEBEDEV) {
        GAng = std::unique_ptr<OneDGrid2>(new Lebedev(nPtsAng));
      };

    };

    virtual inline IntegrationPoint operator[](size_t i) {
      // Running Angular Grid as fastest running
      // IJ = i
      // I = IJ % NAng (Angular Point)
      // J = IJ / NAng (Radial Point)

      cartGP totalPoint; 
      std::size_t I = i % GAng->npts();
      std::size_t J = i / GAng->npts();

      IntegrationPoint gpAng = (*GAng)[I];
      IntegrationPoint gpRad = (*GRad)[J];
      
      totalPoint.set<0>(bg::get<0>(gpRad.pt) * bg::get<0>(gpAng.pt));
      totalPoint.set<1>(bg::get<0>(gpRad.pt) * bg::get<1>(gpAng.pt));
      totalPoint.set<2>(bg::get<0>(gpRad.pt) * bg::get<2>(gpAng.pt));

      double totalWeight = gpRad.weight * gpAng.weight;

      IntegrationPoint totalIntPt(totalPoint,totalWeight);
      totalIntPt.I = I;
      totalIntPt.J = J;
      totalIntPt.NI = GAng->npts();
      totalIntPt.NJ = GRad->npts();

      return totalIntPt;
    };

    inline void generateGridPoints() { };

};
}
#endif
