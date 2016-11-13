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
#ifndef INCLUDED_CUBE_H
#define INCLUDED_CUBE_H
#include <grid2/grid2_def.h>

namespace ChronusQ{
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
        xRange_(xRange), yRange_(yRange), zRange_(zRange), Grid2(0,0.0,false) { 
        
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

}

#endif

