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
  bool  isClosedShell_;

  std::array<double,3> EDField_;

  // FIXME: Need documentation for what these things actually are!
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
  
  inline void checkWorkers(){
    if(this->fileio_  == NULL) 
      CErr("Fatal: Must initialize RealTime with FileIO Object");
    if(this->controls_ == NULL) 
      CErr("Fatal: Must initialize RealTime with Controls Object",
           this->fileio_->out);
    if(this->aointegrals_== NULL)
      CErr("Fatal: Must initialize RealTime with AOIntegrals Object",
           this->fileio_->out);
    if(this->groundState_== NULL)
      CErr("Fatal: Must initialize RealTime with SingleSlater Reference Object",
           this->fileio_->out);
  }

  inline void checkMeta(){
    this->checkWorkers();
    if(this->nBasis_ == 0)
      CErr( "Fatal: RealTime Object Initialized with NBasis = 0",
        this->fileio_->out);
    if(this->Ref_ == SingleSlater<double>::_INVALID) 
      CErr("Fatal: RealTime reference not valid!",this->fileio_->out);
  }

public:

  // constructor & destructor
  RealTime(){
    this->nBasis_      = 0;
    this->Ref_         = 0;
    this->nOccA_       = 0;
    this->nOccB_       = 0;
    this->nSkip_       = 0;

    this->deltaT_      = 0.0;
    this->currentTime_ = 0.0;

    this->fileio_      = NULL;
    this->controls_    = NULL;
    this->aointegrals_ = NULL;
    this->groundState_ = NULL;

    this->oTrans1_ = nullptr;
    this->oTrans2_ = nullptr;
    this->POA_     = nullptr;
    this->POAsav_  = nullptr;
    this->POB_     = nullptr;
    this->POBsav_  = nullptr;
    this->FOA_     = nullptr;
    this->FOB_     = nullptr;
    this->initMOA_ = nullptr;
    this->initMOB_ = nullptr;
    this->uTransA_ = nullptr;
    this->uTransB_ = nullptr;
    this->scratch_ = nullptr;

    this->isClosedShell_ = false;

    // Standard Values
    this->frozenNuc_   = true;
    this->EDField_     = {0.0,0.0,0.0};
    this->nTCS_        = 1;
    this->maxSteps_    = 50;
    this->stepSize_    = 0.05;
    this->typeOrtho_   = 1;
    this->initDensity_ = 0;
    this->swapMOA_     = 0;
    this->swapMOB_     = 0;
    this->methFormU_   = 1;
  };
  ~RealTime() {;};

  inline void communicate(FileIO &fileio, Controls &cont, AOIntegrals &aoints, 
                SingleSlater<T> &groundState){
    this->fileio_      = &fileio;
    this->controls_    = &cont;
    this->aointegrals_ = &aoints;
    this->groundState_ = &groundState;
  }

  inline void initMeta(){
    this->nBasis_        = this->groundState_->nBasis();
    this->isClosedShell_ = this->groundState_->isClosedShell;
    this->Ref_           = this->groundState_->Ref();
    this->nTCS_          = this->groundState_->nTCS();
    this->nOccA_         = this->groundState_->nOccA();
    this->nOccB_         = this->groundState_->nOccB();
  }

  void alloc();

  // pseudo-constructor
  void iniRealTime(FileIO *,Controls *,AOIntegrals *,SingleSlater<T> *);
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
