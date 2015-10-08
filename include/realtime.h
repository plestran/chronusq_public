/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
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
#ifndef INCLUDED_REALTIME
#define INCLUDED_REALTIME
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <controls.h>
#include <aointegrals.h>
#include <singleslater.h>
#include <basisset.h>

namespace ChronusQ {
template<typename T>
class RealTime {
  FileIO *        fileio_;
  Controls *      controls_;
  BasisSet *	  basisset_;
  AOIntegrals *   aointegrals_;
  SingleSlater<T> *  groundState_;

  std::unique_ptr<SingleSlater<dcomplex>> ssPropagator_;

  int	nBasis_;
  int   isClosedShell_;
  int 	nOccA_;
  int 	nOccB_;
  int   nTCS_;
  int   Ref_;
  int	maxSteps_;	// Maximum number of steps
  int	nSkip_;		// Number of Steps to skpi printing
  int	initDensity_;  	// Initial density
  int   swapMOA_;	// Alpha MOs to swap, Format: MMMNNN, swap M with N
  int   swapMOB_;	// Beta MOs to swap, Format: MMMNNN, swap M with N
  int   typeOrtho_;   	// Type of orthonormal transformation
  int   methFormU_;	// Method of forming unitary transformation matrix

  double stepSize_;	// Input step size
  double deltaT_;	// Actual step size
  double currentTime_;	// Current time

  bool	frozenNuc_;     // Whether to freeze nuclei

  std::unique_ptr<std::array<double,3>> EDField_;

  std::unique_ptr<ComplexMatrix>  oTrans1_;
  std::unique_ptr<ComplexMatrix>  oTrans2_;
//  std::unique_ptr<ComplexMatrix>  PA_;
//  std::unique_ptr<ComplexMatrix>  PB_;
  std::unique_ptr<ComplexMatrix>  POA_;
  std::unique_ptr<ComplexMatrix>  POAsav_;
  std::unique_ptr<ComplexMatrix>  POB_;
  std::unique_ptr<ComplexMatrix>  POBsav_;
  std::unique_ptr<ComplexMatrix>  FOA_;
  std::unique_ptr<ComplexMatrix>  FOB_;
  std::unique_ptr<ComplexMatrix>  initMOA_;
  std::unique_ptr<ComplexMatrix>  initMOB_;
  std::unique_ptr<ComplexMatrix>  uTransA_;
  std::unique_ptr<ComplexMatrix>  uTransB_;
  std::unique_ptr<ComplexMatrix>  scratch_;

public:

  // constructor & destructor
  RealTime(){;};
  ~RealTime() {;};

  // pseudo-constructor
  void iniRealTime(Molecule *,BasisSet *,FileIO *,Controls *,AOIntegrals *,
                   SingleSlater<T> *);
  void iniDensity(); // initialize density
//  void formComplexFock();
  void formEDField();
  void printRT();
  void formUTrans();
  void doPropagation();
};

#include <realtime_alloc.h>
} // namespace ChronusQ
#endif
