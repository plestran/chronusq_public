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
#include <cerr.h>
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
  IntegrationPoint(double pt_ = 0, double weight_ = 0) : 
    pt(pt_), weight(weight_){ };
  IntegrationPoint(cartGP pt_ = cartGP(0.0,0.0,0.0), double weight_ = 0) : 
    pt(pt_), weight(weight_){ };
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
    template<typename T>
    inline T integrate(std::function< T(cartGP) > func) {

      T result = (*this)[0].weight * func((*this)[0].pt);

      for(auto iPt = 1; iPt < this->nPts_; iPt++){
        IntegrationPoint intPoint = (*this)[iPt];
        cartGP gp = intPoint.pt;
        result += intPoint.weight * func(gp);
      }
      return result;
    };

    template <typename T>
    inline void integrate(std::function< T(IntegrationPoint) > func,
        T& result) {

      result = func((*this)[0]);
      for(auto iPt = 1; iPt < this->nPts_; iPt++)
        result += func((*this)[iPt]);

    };

    template <typename T>
    inline void integrate(std::function< void(IntegrationPoint,T&) > func,
        T& result) {
      for(auto iPt = 0; iPt < this->nPts_; iPt++)
        if((*this)[iPt].weight > 1e-6){
          func((*this)[iPt],result);
        }
    };
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
    IntegrationPoint operator[](size_t i){ 
      double pt  = (i + 1.0) * (i + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i));
      double wgt = 2.0 * (i + 1.0) * (this->nPts_ + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i) * (this->nPts_ - i));

      return IntegrationPoint(pt,wgt);
    };
};

enum LEBEDEV_ALGEBRAIC_ORDER {
  LEBEDEV_3, LEBEDEV_5, LEBEDEV_7, LEBEDEV_9, LEBEDEV_11, LEBEDEV_13,
  LEBEDEV_15, LEBEDEV_17, LEBEDEV_19, LEBEDEV_21, LEBEDEV_23,
  LEBEDEV_25, LEBEDEV_27, LEBEDEV_29, LEBEDEV_31, LEBEDEV_35,
  LEBEDEV_41, LEBEDEV_47, LEBEDEV_53, LEBEDEV_59, LEBEDEV_65,
  LEBEDEV_71, LEBEDEV_77, LEBEDEV_83, LEBEDEV_89, LEBEDEV_95,
  LEBEDEV_101, LEBEDEV_107, LEBEDEV_113, LEBEDEV_119, LEBEDEV_125,
  LEBEDEV_131
};

class Lebedev : public OneDGrid2 {
  LEBEDEV_ALGEBRAIC_ORDER algOrder_;
  std::vector<cartGP> sphPoints_;
  std::vector<double> sphWeights_;

  void loadAlgebraicPoints();
  template <LEBEDEV_ALGEBRAIC_ORDER ORDER> void loadLebedev();

  public:
    Lebedev(size_t N) : OneDGrid2(N){
      if(     N ==    6) this->algOrder_ = LEBEDEV_3;  
      else if(N ==   14) this->algOrder_ = LEBEDEV_5;  
      else if(N ==   26) this->algOrder_ = LEBEDEV_7;  
      else if(N ==   38) this->algOrder_ = LEBEDEV_9;  
      else if(N ==   50) this->algOrder_ = LEBEDEV_11; 
      else if(N ==   74) this->algOrder_ = LEBEDEV_13; 
      else if(N ==   86) this->algOrder_ = LEBEDEV_15; 
      else if(N ==  110) this->algOrder_ = LEBEDEV_17; 
      else if(N ==  146) this->algOrder_ = LEBEDEV_19; 
      else if(N ==  170) this->algOrder_ = LEBEDEV_21; 
      else if(N ==  194) this->algOrder_ = LEBEDEV_23; 
      else if(N ==  230) this->algOrder_ = LEBEDEV_25; 
      else if(N ==  266) this->algOrder_ = LEBEDEV_27; 
      else if(N ==  302) this->algOrder_ = LEBEDEV_29; 
      else if(N ==  350) this->algOrder_ = LEBEDEV_31; 
      else if(N ==  434) this->algOrder_ = LEBEDEV_35; 
      else if(N ==  590) this->algOrder_ = LEBEDEV_41; 
      else if(N ==  770) this->algOrder_ = LEBEDEV_47; 
      else if(N ==  974) this->algOrder_ = LEBEDEV_53; 
      else if(N == 1202) this->algOrder_ = LEBEDEV_59; 
      else if(N == 1454) this->algOrder_ = LEBEDEV_65; 
      else if(N == 1730) this->algOrder_ = LEBEDEV_71; 
      else if(N == 2030) this->algOrder_ = LEBEDEV_77; 
      else if(N == 2354) this->algOrder_ = LEBEDEV_83; 
      else if(N == 2702) this->algOrder_ = LEBEDEV_89; 
      else if(N == 3074) this->algOrder_ = LEBEDEV_95; 
      else if(N == 3470) this->algOrder_ = LEBEDEV_101;
      else if(N == 3890) this->algOrder_ = LEBEDEV_107;
      else if(N == 4334) this->algOrder_ = LEBEDEV_113;
      else if(N == 4802) this->algOrder_ = LEBEDEV_119;
      else if(N == 5294) this->algOrder_ = LEBEDEV_125;
      else if(N == 5810) this->algOrder_ = LEBEDEV_131;
      else
        CErr("Invalid Lebedev Grid Specification");

      this->loadAlgebraicPoints();
    };

    ~Lebedev(){ };
    IntegrationPoint operator[](size_t i){ 
      return IntegrationPoint(this->sphPoints_[i],this->sphWeights_[i]);
    };
};

class TwoDGrid2 : public Grid2 {
  protected:
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

    virtual inline IntegrationPoint operator[](size_t i) {
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

      IntegrationPoint totalIntPt(totalPoint,totalWeight);

      return totalIntPt;
    };

};

enum ATOMIC_PARTITION {
  BECKE,
  FRISCH
};

class AtomicGrid : public TwoDGrid2 {
  std::array<double,3> center_;
  double scalingFactor_;
  ATOMIC_PARTITION partitionScheme_;
  std::vector<std::array<double,3> > otherCenters_;

  public:
    AtomicGrid(size_t nPtsRad, size_t nPtsAng, 
        GRID_TYPE GTypeRad = EULERMAC, 
        GRID_TYPE GTypeAng = LEBEDEV, 
        std::array<double,3> center = {0.0,0.0,0.0},
        ATOMIC_PARTITION partitionScheme = BECKE, 
        double scalingFactor = 1.0) : 
      TwoDGrid2(nPtsRad,nPtsAng,GTypeRad,GTypeAng),
      center_(center),
      partitionScheme_(partitionScheme),
      scalingFactor_(scalingFactor) { };

    inline IntegrationPoint operator[](size_t i) {
      IntegrationPoint rawPoint = TwoDGrid2::operator[](i);

      // Re-center
      rawPoint.pt.set<0>(bg::get<0>(rawPoint.pt) + center_[0]);
      rawPoint.pt.set<1>(bg::get<1>(rawPoint.pt) + center_[1]);
      rawPoint.pt.set<2>(bg::get<2>(rawPoint.pt) + center_[2]);
      
      // Rescale radius
      rawPoint.pt.set<0>(bg::get<0>(rawPoint.pt) * scalingFactor_);
      rawPoint.pt.set<1>(bg::get<1>(rawPoint.pt) * scalingFactor_);
      rawPoint.pt.set<2>(bg::get<2>(rawPoint.pt) * scalingFactor_);

      // Rescale weight (w/o Partition Weights)
      rawPoint.weight *= scalingFactor_;

      return rawPoint;

    };
};

};
#endif
