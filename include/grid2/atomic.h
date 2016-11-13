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
#ifndef INCLUDED_ATOMIC_H
#define INCLUDED_ATOMIC_H
#include <grid2/twodgrid.h>

namespace ChronusQ {

enum ATOMIC_PARTITION {
  BECKE,
  FRISCH
};


class AtomicGrid : public TwoDGrid2 {
  double scalingFactor_;
  ATOMIC_PARTITION partitionScheme_;
  std::vector<std::array<double,3> > centers_;
  size_t centerIndx_;
  double nearestNeighbor_;
  double radCutOff_;
  RealMatrix *rIJ_;

  std::vector<std::vector<double>> partitionScratch_;

  double evalPartitionWeight(cartGP&);
  double hBecke(double x);
  double zFrisch(double x,double a);
  double gBecke(double x);
  double gFrisch(double x);
  inline double g(double x){
    if(     this->partitionScheme_ == BECKE)  return gBecke(x);
    else if(this->partitionScheme_ == FRISCH) return gFrisch(x);
  }

  public:
    void findNearestNeighbor();
    AtomicGrid(size_t nPtsRad, size_t nPtsAng, 
        GRID_TYPE GTypeRad, GRID_TYPE GTypeAng, 
        ATOMIC_PARTITION partitionScheme, 
        std::vector<std::array<double,3> > centers,
        RealMatrix *rIJ,
        size_t centerIndx,
        double screenTol = 0.0,
        double radCutOff = 1e6,
        double scalingFactor = 1.0,
        bool onTheFly = true) : 
      TwoDGrid2(nPtsRad,nPtsAng,screenTol,GTypeRad,GTypeAng,onTheFly),
      partitionScheme_(partitionScheme),
      centers_(std::move(centers)),
      rIJ_(rIJ),
      centerIndx_(centerIndx),
      radCutOff_(radCutOff),
      scalingFactor_(scalingFactor) { 
        this->partitionScratch_.resize(omp_get_max_threads());
        for(auto i = 0; i < omp_get_max_threads();i++)
          this->partitionScratch_[i].resize(this->centers_.size());
  
    };

    inline IntegrationPoint operator[](size_t i) {
      // return the struc (Cart Pt, weight)
      IntegrationPoint rawPoint = TwoDGrid2::operator[](i);

      // Rescale radius
      rawPoint.pt.set<0>(bg::get<0>(rawPoint.pt) * scalingFactor_);
      rawPoint.pt.set<1>(bg::get<1>(rawPoint.pt) * scalingFactor_);
      rawPoint.pt.set<2>(bg::get<2>(rawPoint.pt) * scalingFactor_);

      // Re-center
      rawPoint.pt.set<0>(bg::get<0>(rawPoint.pt)+centers_[centerIndx_][0]);
      rawPoint.pt.set<1>(bg::get<1>(rawPoint.pt)+centers_[centerIndx_][1]);
      rawPoint.pt.set<2>(bg::get<2>(rawPoint.pt)+centers_[centerIndx_][2]);

      // Rescale weight (w/o Partition Weights)
      double r = bg::get<0>(GRad->operator[](i / GAng->npts()).pt)
        * scalingFactor_;
      if(r > this->radCutOff_ && this->screenON_) {
//      if(r > this->radCutOff_) {
        rawPoint.evalpt = false;
//        cout << "Skipped " <<"r " << r << "rCut " << this->radCutOff_ <<endl;
        return rawPoint;
      }
      rawPoint.weight *= scalingFactor_ * r*r;
      double partweight = 1;
//Screening no off APE
//      if(rawPoint.weight > 1e-8)
//    if(std::abs(rawPoint.weight) < std::numeric_limits<double>::epsilon()) 
        partweight = evalPartitionWeight(rawPoint.pt);
//    else
//      partweight = 0.0;
        
//    if(partweight < 1e-10) rawPoint.evalpt = false;
 
      rawPoint.weight *= partweight;

      if((std::abs(rawPoint.weight) < this->screenTol_) && this->screenON_) 
//      if((std::abs(rawPoint.weight) < this->screenTol_) ) 
        rawPoint.evalpt = false;

//Screening now off APE
//      if(partweight < 1e-6) rawPoint.evalpt = false;
      return rawPoint;

    };

    inline void setCenter(size_t i) { this->centerIndx_ = i; };
    inline void setScalingFactor(double x) { this->scalingFactor_ = x; };
    inline void setRadCutOff(double x) {this->radCutOff_ = x; };

    inline size_t& center(){ return this->centerIndx_; };
    inline double& scalingFactor(){ return this->scalingFactor_; };
};

}
#endif
