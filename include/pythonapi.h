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
#include <workers.h>
#include <global.h>
#include <memory.h>
#include <atoms.h>
#include <numdiff.h>
#include <response.h>

using namespace boost::python;
using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::AOIntegrals;
using ChronusQ::MOIntegrals;
using ChronusQ::RealTime;
using ChronusQ::CErr;
using ChronusQ::SDResponse;
using ChronusQ::CQMemManager;
using ChronusQ::NumericalDifferentiation;
using ChronusQ::RESPONSE_TYPE;

namespace ChronusQ{
  void Wrapper_readInput(FileIO&,Molecule&,BasisSet&,Controls&,BasisSet&);

  void Wrapper_CErr_Default(FileIO &);
  void Wrapper_CErr_Message(FileIO &, std::string );
  inline void CQSetNumThreads(int n) {
#ifdef _OPENMP
    omp_set_num_threads(n);
#endif
     int h = 0;
  };
  int getAtomicNumber(int);
  void Wrapper_initCQ(int,boost::python::list);
}
