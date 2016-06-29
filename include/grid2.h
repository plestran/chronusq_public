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
  bool evalpt;
  IntegrationPoint(double pt_ = 0, double weight_ = 0) : 
    pt(pt_), weight(weight_), evalpt(true){ };
  IntegrationPoint(cartGP pt_ = cartGP(0.0,0.0,0.0), double weight_ = 0) : 
    pt(pt_), weight(weight_), evalpt(true){ };
};

// Classes
class Grid2 {
protected:
    size_t nPts_;  ///< number of grid points       
    std::vector<cartGP> gPoints_;
    std::vector<double> weights_;
    bool onTheFly_;
    bool haveGPs_;
public:
    Grid2(size_t npts = 0, bool onTheFly = true) : nPts_(npts),
      onTheFly_(onTheFly), haveGPs_(false){ 
      
    };

    virtual ~Grid2(){ };

    virtual IntegrationPoint operator[](size_t) = 0;
    size_t npts(){return this->nPts_;};
    virtual void generateGridPoints() {
      this->gPoints_.reserve(this->nPts_);
      this->weights_.reserve(this->nPts_);

      for(auto iPt = 0; iPt < nPts_; iPt++) {
        IntegrationPoint tmp = operator[](iPt);
        this->gPoints_.push_back(tmp.pt);
        this->weights_.push_back(tmp.weight);
      }
      this->haveGPs_ = true;
    };

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
      std::size_t NSkip(0);
      for(auto iPt = 0; iPt < this->nPts_; iPt++)
        if((*this)[iPt].evalpt){
          func((*this)[iPt],result);
        } else NSkip++;
      cout << "NSkip : " << NSkip << endl;
    };


    // Print function
    void printGrid(std::ostream&);
}; // class Grid2

class OneDGrid2 : public Grid2 {
  public:
    OneDGrid2(size_t npts = 0, bool onTheFly = true) : Grid2(npts,onTheFly){ };
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
    GaussChebFst(size_t npts = 0, bool onTheFly = true) : 
      OneDGrid2(npts,onTheFly){
        if(!onTheFly_) this->generateGridPoints();
      };

    ~GaussChebFst(){ };
    IntegrationPoint operator[](size_t i){ 
      if(!onTheFly_ && haveGPs_)
        return IntegrationPoint(gPoints_[i],weights_[i]);

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
    GaussChebSnd(size_t npts = 0, bool onTheFly = true) : 
      OneDGrid2(npts,onTheFly){
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
  void loadAlgebraicPoints();
  template <LEBEDEV_ALGEBRAIC_ORDER ORDER> void loadLebedev();

  public:
    Lebedev(size_t N) : OneDGrid2(N,false){
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
      this->haveGPs_ = true;
    };

    ~Lebedev(){ };
    IntegrationPoint operator[](size_t i){ 
      return IntegrationPoint(this->gPoints_[i],this->weights_[i]);
    };
};

class TwoDGrid2 : public Grid2 {
  protected:
    std::unique_ptr<OneDGrid2> GRad;
    std::unique_ptr<OneDGrid2> GAng;

  public:
    TwoDGrid2(size_t nPtsRad, size_t nPtsAng, GRID_TYPE GTypeRad, 
        GRID_TYPE GTypeAng, bool onTheFly = true) : 
      Grid2(nPtsRad*nPtsAng,true) {

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

    inline void generateGridPoints() { };

};

enum ATOMIC_PARTITION {
  BECKE,
  FRISCH
};


class AtomicGrid : public TwoDGrid2 {
  double scalingFactor_;
  ATOMIC_PARTITION partitionScheme_;
  std::vector<std::array<double,3> > centers_;
  size_t centerIndx_;
  RealMatrix *rIJ_;

  std::vector<double> partitionScratch_;

  double evalPartitionWeight(cartGP&);

  public:
    AtomicGrid(size_t nPtsRad, size_t nPtsAng, 
        GRID_TYPE GTypeRad, GRID_TYPE GTypeAng, 
        ATOMIC_PARTITION partitionScheme, 
        std::vector<std::array<double,3> > centers,
        RealMatrix *rIJ,
        size_t centerIndx,
        double scalingFactor = 1.0,
        bool onTheFly = true) : 
      TwoDGrid2(nPtsRad,nPtsAng,GTypeRad,GTypeAng,onTheFly),
      partitionScheme_(partitionScheme),
      centers_(std::move(centers)),
      rIJ_(rIJ),
      centerIndx_(centerIndx),
      scalingFactor_(scalingFactor) { 
        this->partitionScratch_.resize(this->centers_.size(),0.0);
    };

    inline IntegrationPoint operator[](size_t i) {
      // return the struc (Cart Pt, weight)
      IntegrationPoint rawPoint = TwoDGrid2::operator[](i);

      // Rescale radius
      rawPoint.pt.set<0>(bg::get<0>(rawPoint.pt) * scalingFactor_);
      rawPoint.pt.set<1>(bg::get<1>(rawPoint.pt) * scalingFactor_);
      rawPoint.pt.set<2>(bg::get<2>(rawPoint.pt) * scalingFactor_);

      // Re-center
      rawPoint.pt.set<0>(bg::get<0>(rawPoint.pt) + centers_[centerIndx_][0]);
      rawPoint.pt.set<1>(bg::get<1>(rawPoint.pt) + centers_[centerIndx_][1]);
      rawPoint.pt.set<2>(bg::get<2>(rawPoint.pt) + centers_[centerIndx_][2]);

      // Rescale weight (w/o Partition Weights)
      double r = bg::get<0>(GRad->operator[](i / GAng->npts()).pt)
        * scalingFactor_;
      rawPoint.weight *= scalingFactor_ * r*r;
      double partweight = 1;
//Screening no off APE
//      if(rawPoint.weight > 1e-8)
        partweight = evalPartitionWeight(rawPoint.pt);
        
      if(partweight < 1e-10) rawPoint.evalpt = false;
 
      rawPoint.weight *= partweight;
//Screening now off APE
//      if(partweight < 1e-6) rawPoint.evalpt = false;
      return rawPoint;

    };

    inline void setCenter(size_t i) { this->centerIndx_ = i; };
    inline void setScalingFactor(double x) { this->scalingFactor_ = x; };

    inline size_t& center(){ return this->centerIndx_; };
    inline double& scalingFactor(){ return this->scalingFactor_; };
};

class Cube : public Grid2 {
  std::tuple<double,double,size_t> xRange_;
  std::tuple<double,double,size_t> yRange_;
  std::tuple<double,double,size_t> zRange_;
  double xRes_;
  double yRes_;
  double zRes_;

  public:
    Cube(std::tuple<double,double,size_t> xRange,
        std::tuple<double,double,size_t> yRange,
        std::tuple<double,double,size_t> zRange) :
        xRange_(xRange), yRange_(yRange), zRange_(zRange), Grid2(0,false) { 
        
        xRes_ = (std::get<1>(this->xRange_) - std::get<0>(this->xRange_)) / 
          (std::get<2>(this->xRange_) - 1);
        yRes_ = (std::get<1>(this->yRange_) - std::get<0>(this->yRange_)) / 
          (std::get<2>(this->yRange_) - 1);
        zRes_ = (std::get<1>(this->zRange_) - std::get<0>(this->zRange_)) / 
          (std::get<2>(this->zRange_) - 1);

        nPts_ = std::get<2>(this->xRange_) * std::get<2>(this->yRange_)
          * std::get<2>(this->zRange_);
    };

    inline IntegrationPoint operator[](size_t i) {
      // i -> (M,N,L)
      // M (z) = i mod nZPts
      // N (y) = (i div nZPts) mod nYPts
      // L (x) = i div (nZPts * nYPts)
      size_t zIndex = i % std::get<2>(this->zRange_);
      size_t yIndex = (i / std::get<2>(this->zRange_)) % 
        std::get<2>(this->yRange_);
      size_t xIndex = i / 
        (std::get<2>(this->zRange_) * std::get<2>(this->yRange_));

      double xPt = std::get<0>(this->xRange_) + xIndex * this->xRes_;
      double yPt = std::get<0>(this->yRange_) + yIndex * this->yRes_;
      double zPt = std::get<0>(this->zRange_) + zIndex * this->zRes_;

      cartGP pt(xPt,yPt,zPt);
      return IntegrationPoint(pt,1.0);

    };
    inline void generateGridPoints() { };

    template<typename T, typename Mol>
    inline void genCubeFile(T func, std::string cubeFileName,
        Mol &molecule) {
      std::ofstream cubeFile(cubeFileName);
      // Print Cube File Header
      cubeFile << "ChronusQ CubeFile" << endl;
      cubeFile << "OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" << endl;
      cubeFile << std::left << std::fixed;
      cubeFile << std::setw(6) << molecule.nAtoms();
      cubeFile << std::setw(10) << std::get<0>(this->xRange_);
      cubeFile << std::setw(10) << std::get<0>(this->yRange_);
      cubeFile << std::setw(10) << std::get<0>(this->zRange_);
      cubeFile << endl;

      cubeFile << std::setw(6) << std::get<2>(this->xRange_);
      cubeFile << std::setw(10) << this->xRes_;
      cubeFile << std::setw(10) << 0.0;
      cubeFile << std::setw(10) << 0.0;
      cubeFile << endl;

      cubeFile << std::setw(6) << std::get<2>(this->yRange_);
      cubeFile << std::setw(10) << 0.0;
      cubeFile << std::setw(10) << this->yRes_;
      cubeFile << std::setw(10) << 0.0;
      cubeFile << endl;

      cubeFile << std::setw(6) << std::get<2>(this->zRange_);
      cubeFile << std::setw(10) << 0.0;
      cubeFile << std::setw(10) << 0.0;
      cubeFile << std::setw(10) << this->zRes_;
      cubeFile << endl;

      for(auto iAtm = 0; iAtm < molecule.nAtoms(); iAtm++){
        cubeFile << std::setw(6) << molecule.atomicZ(iAtm);
        cubeFile << std::setw(10) << 0.0;
        cubeFile << std::setw(10) << (*molecule.cart())(0,iAtm);
        cubeFile << std::setw(10) << (*molecule.cart())(1,iAtm);
        cubeFile << std::setw(10) << (*molecule.cart())(2,iAtm);
        cubeFile << endl;
      };

      for(auto iPt = 0; iPt < this->nPts_; iPt++){
        std::stringstream ss;
        ss << std::scientific << func(Cube::operator[](iPt).pt);
        std::string token(ss.str());
        std::replace(token.begin(),token.end(),'e','E');

        cubeFile << std::setw(15) << token;
        if( iPt % 6 == 5) cubeFile << endl;
//      else if(iPt != 0 && iPt % std::get<2>(this->zRange_) == 0) 
//        cubeFile << endl;
      };

    };
};

};
#endif
