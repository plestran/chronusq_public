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
#ifndef INCLUDED_ODGRID_H
#define INCLUDED_ODGRID_H

#include <grid2/grid2_def.h>

namespace ChronusQ {

class OneDGrid2 : public Grid2 {
  public:
    OneDGrid2(size_t npts = 0, double screenTol = 0.0, bool onTheFly = true) : Grid2(npts,screenTol,onTheFly){ };
    virtual ~OneDGrid2(){ };
    virtual IntegrationPoint operator[](size_t) = 0;
};
}
#endif
