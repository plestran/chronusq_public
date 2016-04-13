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

    virtual IntegrationPoint operator[](size_t) = 0;
    size_t npts(){return this->nPts_;};

    // Integrate Functions
    template<typename T>
    inline T integrate(std::function< T(cartGP) > func) {

      cout << "Integrate 1" << endl;
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

      cout << "Integrate 2" << endl;
      result = func((*this)[0]);
      for(auto iPt = 1; iPt < this->nPts_; iPt++)
        result += func((*this)[iPt]);

    };

    template <typename T>
    inline void integrate(std::function< void(IntegrationPoint,T&) > func,
        T& result) {
      cout << "Integrate 3" << endl;
      for(auto iPt = 0; iPt < this->nPts_; iPt++)
//        if((*this)[iPt].weight > 1e-6){
          func((*this)[iPt],result);
//        }
    };


    // Print function
    void printGrid(std::ostream&);
}; // class Grid2

class OneDGrid2 : public Grid2 {
  public:
    OneDGrid2(size_t npts = 0) : Grid2(npts){ };
    virtual ~OneDGrid2(){ };
    virtual IntegrationPoint operator[](size_t) = 0;
};

/**
 * Gauss-Chebyshev Quadrature of the first kind
 *
 * \int _{-1}^{1} \frac{f(x)}{\sqrt{1 - x^2}} \mathrm{d} x = \sum_n w_i f(x_i)
 *
 * x_i = \cos\left( \frac{2i - 1}{2n} \pi \right)
 * w_i = \frac{\pi}{n}
 *
 * where "n" is the number of integration points
 *
 */
class GaussChebFst : public OneDGrid2 {
  public:
    GaussChebFst(size_t npts = 0) : OneDGrid2(npts){ };
    ~GaussChebFst(){ };
    IntegrationPoint operator[](size_t i){ 
      // Generate points and weights for (-1,1)
      double pt  = std::cos( (2.0*(i+1)-1.0) / (2*this->nPts_) * math.pi );
      double wgt = (math.pi / this->nPts_);

      // Map (-1,1) -> (0,\inf)
        
      // As we're integrating f(x) not \frac{f(x)}{\sqrt{1 - x^2}}
      // we must factor the \sqrt{1-x^2} into the weights
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
      wgt *= 2.0 * std::sqrt(1 - pt * pt) / ( (1 - pt) * (1 - pt) );

      // Perform the coordinate shift
      //
      // r_i = \frac{(1+x_i)}{(1-x_i)}
      pt   = (1 + pt) / (1 - pt);

      return IntegrationPoint(pt,wgt);

    };
};

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
    GaussChebSnd(size_t npts = 0) : OneDGrid2(npts){ };
    ~GaussChebSnd(){ };
    IntegrationPoint operator[](size_t i){ 
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

class EulerMac : public OneDGrid2 {
  public:
    EulerMac(size_t npts) : OneDGrid2(npts){ };
    ~EulerMac(){ };
    IntegrationPoint operator[](size_t i){ 
      double pt  = (i + 1.0) * (i + 1.0) / 
        ((this->nPts_ - i) * (this->nPts_ - i));
      double wgt = 2.0 * (i + 1.0) * (this->nPts_ + 1.0) * pt * pt / 
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


class AtomicGrid2 : public TwoDGrid2 {
  double scalingFactor_;
  ATOMIC_PARTITION partitionScheme_;
  std::vector<std::array<double,3> > centers_;
  size_t centerIndx_;

  std::vector<double> partitionScratch_;

  double evalPartitionWeight(cartGP&);

  public:
    AtomicGrid2(size_t nPtsRad, size_t nPtsAng, 
        GRID_TYPE GTypeRad, GRID_TYPE GTypeAng, 
        ATOMIC_PARTITION partitionScheme, 
        std::vector<std::array<double,3> > centers,
        size_t centerIndx,
        double scalingFactor = 1.0) : 
      TwoDGrid2(nPtsRad,nPtsAng,GTypeRad,GTypeAng),
      partitionScheme_(partitionScheme),
      centers_(std::move(centers)),
      centerIndx_(centerIndx),
      scalingFactor_(scalingFactor) { 
        this->partitionScratch_.resize(this->centers_.size(),0.0);
      };

    inline IntegrationPoint operator[](size_t i) {
      IntegrationPoint rawPoint = TwoDGrid2::operator[](i);

      // Rescale radius
      rawPoint.pt.set<0>(bg::get<0>(rawPoint.pt) * scalingFactor_);
      rawPoint.pt.set<1>(bg::get<1>(rawPoint.pt) * scalingFactor_);
      rawPoint.pt.set<2>(bg::get<2>(rawPoint.pt) * scalingFactor_);
      double x = bg::get<0>(rawPoint.pt);
      double y = bg::get<1>(rawPoint.pt);
      double z = bg::get<2>(rawPoint.pt);
      double r = std::sqrt(x*x + y*y + z*z);
      // Re-center
      rawPoint.pt.set<0>(bg::get<0>(rawPoint.pt) + centers_[centerIndx_][0]);
      rawPoint.pt.set<1>(bg::get<1>(rawPoint.pt) + centers_[centerIndx_][1]);
      rawPoint.pt.set<2>(bg::get<2>(rawPoint.pt) + centers_[centerIndx_][2]);


      

      // Rescale weight (w/o Partition Weights)
      rawPoint.weight *= scalingFactor_ * r*r;

      //cout <<evalPartitionWeight(rawPoint.pt)<<endl;
//      if(rawPoint.weight > 1e-8)
      rawPoint.weight *= evalPartitionWeight(rawPoint.pt);
//      cout << " Weight Becke Done on center " << this->centerIndx_ <<endl;    
      return rawPoint;

    };
};

class Cube : public Grid2 {
  std::tuple<double,double,size_t> xRange_;
  std::tuple<double,double,size_t> yRange_;
  std::tuple<double,double,size_t> zRange_;

  public:
    Cube(std::tuple<double,double,size_t> xRange,
        std::tuple<double,double,size_t> yRange,
        std::tuple<double,double,size_t> zRange) :
        xRange_(xRange), yRange_(yRange), zRange_(zRange) { };

    inline IntegrationPoint operator[](size_t i) {
      // i -> (M,N,L)
      // M (x) = i mod nXPts
      // N (y) = (i div nXPts) mod nYPts
      // L (z) = i div (nXPts * nYPts)
      size_t zIndex = i % std::get<2>(this->zRange_);
      size_t yIndex = (i / std::get<2>(this->zRange_)) % 
        std::get<2>(this->yRange_);
      size_t xIndex = i / 
        (std::get<2>(this->zRange_) * std::get<2>(this->yRange_));

      double xPt = std::get<0>(this->xRange_) + 
        xIndex * (std::get<1>(this->xRange_) - std::get<0>(this->xRange_)) /
        (std::get<2>(this->xRange_) - 1);
      double yPt = std::get<0>(this->yRange_) + 
        yIndex * (std::get<1>(this->yRange_) - std::get<0>(this->yRange_)) /
        (std::get<2>(this->yRange_) - 1);
      double zPt = std::get<0>(this->zRange_) + 
        zIndex * (std::get<1>(this->zRange_) - std::get<0>(this->zRange_)) / 
        (std::get<2>(this->zRange_) - 1);

      cartGP pt(xPt,yPt,zPt);
      return IntegrationPoint(pt,1.0);

    };
};

};
#endif
