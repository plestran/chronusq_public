/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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
#ifndef INCLUDED_GRID2
#define INCLUDED_GRID2
#include <global.h>
namespace ChronusQ {

enum GRID_TYPE {
  GAUSSCHEBFST,
  GAUSSCHEBSND,
  EULERMAC,
  LEBEDEV
};

struct IntegrationPoint {
  cartGP pt;
  double weight;
};

// Classes
class Grid2 {
protected:
    size_t nPts_;  ///< number of grid points       
public:
    Grid2(size_t npts = 0) : nPts_(npts){ };
    virtual ~Grid2(){ };
    virtual void   printGrid() = 0; ///<function to print the grid points
    virtual IntegrationPoint operator[](size_t) = 0;
    size_t npts(){return this->nPts_;};
}; // class Grid2

class OneDGrid2 : public Grid2 {
  public:
    OneDGrid2(size_t npts = 0) : Grid2(npts){ };
};

class TwoDGrid2 : public Grid2 {
  OneDGrid GRad;
  OneDGrid GAng;

  public:
    TwoDGrid() = delete;
    TwoDGrid(size_t nPtsRad, size_t nPtsAng, GRID_TYPE GTypeRad, 
        GRID_TYPE GTypeAng) : nPts_(nPtsRad*nPtsAng){

      if(GTypeRad == GAUSSCHEBFST) {
        GRad = new GaussChebFst(nPtsRad);
      } else if(GTypeRad == GAUSSCHEBSND) {
        GRad = new GaussChebSnd(nPtsRad);
      } else if(GTypeRad == EULERMAC) {
        GRad = new EulerMac(nPtsRad);
      };

      if(GTypeAng == LEBEDEV) {
        GAng = new Lebedev(nPtsAng);
      };

    };

    inline IntegrationPoint operator[](size_t i) {
      // Running Angular Grid as fastest running
      // IJ = i
      // I = IJ % NAng (Angular Point)
      // J = IJ / NAng (Radial Point)
      IntegrationPoint gpAng = GAng[i % GAng.npts()];
      IntegrationPoint gpRad = GRad[i / GAng.npts()];
      
      carGP totalPoint; 
      totalPoint.set<0>(bg::get<0>(gpRad.pt) * bg::get<0>(gpAng.pt));
      totalPoint.set<1>(bg::get<0>(gpRad.pt) * bg::get<1>(gpAng.pt));
      totalPoint.set<2>(bg::get<0>(gpRad.pt) * bg::get<2>(gpAng.pt));

      double totalWeight = gpRad.weight * gpAng.weight;

      IntegrationPoint totalIntPt;
      totalIntPt.pt = totalPoint;
      totalIntPt.weight = totalWeight;

      return totalIntPt;
    };

};

};
