/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explictly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
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
#ifndef  INCLUDED_WORKERS
#define  INCLUDED_WORKERS
//#include "mpi.h"
#include "global.h"
#include "fileio.h"
#include "atoms.h"
#include "molecule.h"
#include "matrix.h"
#include "basisset.h"
#include "aointegrals.h"
#include "mointegrals.h"
#include "sdresponse.h"
#include "singleslater.h"
#include "tools.h"
#include "classtools.h"

namespace ChronusQ {
int atlas(int, char**, GlobalMPI*);
int worker(GlobalMPI*);
} // namespace ChronusQ

#endif
