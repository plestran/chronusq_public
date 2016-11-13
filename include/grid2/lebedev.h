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
#ifndef INCLUDED_LEBEDEV_H
#define INCLUDED_LEBEDEV_H

#include <grid2/odgrid.h>

namespace ChronusQ {
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
    Lebedev(size_t N) : OneDGrid2(N,0.0,false){
/*
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
*/
      if(N ==  170)      this->algOrder_ = LEBEDEV_21; 
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
/*
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
*/
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
}

#endif
