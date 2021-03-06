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
#include <grid2/lebedev.h>

namespace ChronusQ { 
  void Lebedev::loadAlgebraicPoints() {
    this->gPoints_.reserve(this->nPts_);
    this->weights_.reserve(this->nPts_);
/*
    if(this->algOrder_ == LEBEDEV_3)         this->loadLebedev<LEBEDEV_3  >();
    else if(this->algOrder_ == LEBEDEV_5)    this->loadLebedev<LEBEDEV_5  >();
    else if(this->algOrder_ == LEBEDEV_7)    this->loadLebedev<LEBEDEV_7  >();
    else if(this->algOrder_ == LEBEDEV_9)    this->loadLebedev<LEBEDEV_9  >();
    else if(this->algOrder_ == LEBEDEV_11)   this->loadLebedev<LEBEDEV_11 >();
    else if(this->algOrder_ == LEBEDEV_13)   this->loadLebedev<LEBEDEV_13 >();
    else if(this->algOrder_ == LEBEDEV_15 )  this->loadLebedev<LEBEDEV_15 >();
    else if(this->algOrder_ == LEBEDEV_17 )  this->loadLebedev<LEBEDEV_17 >();
    else if(this->algOrder_ == LEBEDEV_19 )  this->loadLebedev<LEBEDEV_19 >();
    else if(this->algOrder_ == LEBEDEV_21 )  this->loadLebedev<LEBEDEV_21 >();
*/
    if(this->algOrder_ == LEBEDEV_21 )  this->loadLebedev<LEBEDEV_21 >();
    else if(this->algOrder_ == LEBEDEV_23 )  this->loadLebedev<LEBEDEV_23 >();
    else if(this->algOrder_ == LEBEDEV_25 )  this->loadLebedev<LEBEDEV_25 >();
    else if(this->algOrder_ == LEBEDEV_27 )  this->loadLebedev<LEBEDEV_27 >();
    else if(this->algOrder_ == LEBEDEV_29 )  this->loadLebedev<LEBEDEV_29 >();
    else if(this->algOrder_ == LEBEDEV_31 )  this->loadLebedev<LEBEDEV_31 >();
    else if(this->algOrder_ == LEBEDEV_35 )  this->loadLebedev<LEBEDEV_35 >();
    else if(this->algOrder_ == LEBEDEV_41 )  this->loadLebedev<LEBEDEV_41 >();
    else if(this->algOrder_ == LEBEDEV_47 )  this->loadLebedev<LEBEDEV_47 >();
    else if(this->algOrder_ == LEBEDEV_53 )  this->loadLebedev<LEBEDEV_53 >();
    else if(this->algOrder_ == LEBEDEV_59 )  this->loadLebedev<LEBEDEV_59 >();
    else if(this->algOrder_ == LEBEDEV_65 )  this->loadLebedev<LEBEDEV_65 >();
  //else if(this->algOrder_ == LEBEDEV_71 )  this->loadLebedev<LEBEDEV_71 >();
  //else if(this->algOrder_ == LEBEDEV_77 )  this->loadLebedev<LEBEDEV_77 >();
  //else if(this->algOrder_ == LEBEDEV_83 )  this->loadLebedev<LEBEDEV_83 >();
  //else if(this->algOrder_ == LEBEDEV_89 )  this->loadLebedev<LEBEDEV_89 >();
  //else if(this->algOrder_ == LEBEDEV_95 )  this->loadLebedev<LEBEDEV_95 >();
  //else if(this->algOrder_ == LEBEDEV_101)  this->loadLebedev<LEBEDEV_101>();
  //else if(this->algOrder_ == LEBEDEV_107)  this->loadLebedev<LEBEDEV_107>();
  //else if(this->algOrder_ == LEBEDEV_113)  this->loadLebedev<LEBEDEV_113>();
  //else if(this->algOrder_ == LEBEDEV_119)  this->loadLebedev<LEBEDEV_119>();
  //else if(this->algOrder_ == LEBEDEV_125)  this->loadLebedev<LEBEDEV_125>();
  //else if(this->algOrder_ == LEBEDEV_131)  this->loadLebedev<LEBEDEV_131>();
  };
}; // namespace ChronusQ
