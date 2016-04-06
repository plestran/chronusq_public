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
  IntegrationPoint(double pt_ = 0, double weight_ = 0) :pt(pt_), weight(weight_){ };
};

// Classes
class Grid2 {
protected:
    size_t nPts_;  ///< number of grid points       
public:
    Grid2(size_t npts = 0) : nPts_(npts){ };
    virtual ~Grid2(){ };
  //virtual void   printGrid() = 0; ///<function to print the grid points
    virtual IntegrationPoint operator[](size_t) = 0;
    size_t npts(){return this->nPts_;};
}; // class Grid2

class OneDGrid2 : public Grid2 {
  public:
    OneDGrid2(size_t npts = 0) : Grid2(npts){ };
    virtual ~OneDGrid2(){ };
    virtual IntegrationPoint operator[](size_t) = 0;
};

class GaussChebFst : public OneDGrid2 {
  public:
    GaussChebFst(size_t npts = 0) : OneDGrid2(npts){ };
    ~GaussChebFst(){ };
    IntegrationPoint operator[](size_t i){ 
      double pt  = std::cos( (2.0*(i+1)-1.0) / (2*this->nPts_) * math.pi );
      double wgt = (math.pi / this->nPts_);

      wgt *= 2.0 * std::sqrt(1 - pt * pt) / ( (1 - pt) * (1 - pt) );
      pt   = (1 + pt) / (1 - pt);

      return IntegrationPoint(pt,wgt);

    };
};

class GaussChebSnd : public OneDGrid2 {
  public:
    GaussChebSnd(size_t npts) : OneDGrid2(npts){ };
    ~GaussChebSnd(){ };
    IntegrationPoint operator[](size_t i){ return IntegrationPoint(0,0);};
};

class EulerMac : public OneDGrid2 {
  public:
    EulerMac(size_t npts) : OneDGrid2(npts){ };
    ~EulerMac(){ };
    IntegrationPoint operator[](size_t i){ return IntegrationPoint(0,0);};
};

class Lebedev : public OneDGrid2 {
  public:
    Lebedev(size_t npts) : OneDGrid2(npts){ };
    ~Lebedev(){ };
    IntegrationPoint operator[](size_t i){ return IntegrationPoint(0,0);};
};

class TwoDGrid2 : public Grid2 {
  std::unique_ptr<OneDGrid2> GRad;
  std::unique_ptr<OneDGrid2> GAng;

  public:
    TwoDGrid2(size_t nPtsRad, size_t nPtsAng, GRID_TYPE GTypeRad, 
        GRID_TYPE GTypeAng) : Grid2(nPtsRad*nPtsAng) {

      if(GTypeRad == GAUSSCHEBFST) {
        GRad = std::unique_ptr<OneDGrid2>(new GaussChebFst(nPtsRad));
      } else if(GTypeRad == GAUSSCHEBSND) {
        GRad = std::unique_ptr<OneDGrid2>(new GaussChebSnd(nPtsRad));
      } else if(GTypeRad == EULERMAC) {
        GRad = std::unique_ptr<OneDGrid2>(new EulerMac(nPtsRad));
      };

      if(GTypeAng == LEBEDEV) {
        GAng = std::unique_ptr<OneDGrid2>(new Lebedev(nPtsAng));
      };

    };

    inline IntegrationPoint operator[](size_t i) {
      // Running Angular Grid as fastest running
      // IJ = i
      // I = IJ % NAng (Angular Point)
      // J = IJ / NAng (Radial Point)
      IntegrationPoint gpAng = (*GAng)[i % GAng->npts()];
      IntegrationPoint gpRad = (*GRad)[i / GAng->npts()];
      
      cartGP totalPoint; 
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
#endif
