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
  int   IEnvlp_;        // Type of envelope function
  double Ex_;           // Magnitude of electric field x-component
  double Ey_;           // Magnitude of electric field y-component
  double Ez_;           // Magnitude of electric field z-component
  double TOn_;          // Time field is turned on (fs)
  double TOff_;         // Time field is turned off (fs)
  double Freq_;         // Frequency of field (eV)
  double Phase_;        // Phase offset of field (rad)
  double Sigma_;        // Width of Gaussian envelope (eV)

  double stepSize_;	// Input step size
  double deltaT_;	// Actual step size
  double currentTime_;	// Current time

  bool	frozenNuc_;     // Whether to freeze nuclei
  bool  isClosedShell_;

  // Memory pointers
  dcomplex * SCR;         // Total Scratch space for RT
  dcomplex * oTrans1Mem_; // Orthogonalizing Transformation Matrix from Overlap
  dcomplex * oTrans2Mem_; // Orthogonalizing Transformation Matrix from Overlap Inverse
  dcomplex * POAMem_;     // Density Alpha in Orthonormal Basis
  dcomplex * POBMem_;     // Density Beta in Orthonormal Basis
  dcomplex * POAsavMem_;  // Saved copy of Density Alpha in Orthonormal Basis
  dcomplex * POBsavMem_;  // Saved copy of Density Beta in Orthonormal Basis
  dcomplex * FOAMem_;     // Fock Matrix Alpha in Orthonormal Basis 
  dcomplex * FOBMem_;     // Fock Matrix Beta in Orthonormal Basis
  dcomplex * initMOAMem_; // Ground State MO Alpha in orthonormal basis
  dcomplex * initMOBMem_; // Ground State MO Beta in orthonormal basis
  dcomplex * uTransAMem_; // Unitary Transformation Matrix [exp(-i*dt*F)] Alpha
  dcomplex * uTransBMem_; // Unitary Transformation Matrix [exp(-i*dt*F)] Beta
  dcomplex * scratchMem_; // NBas x NBas scratch Matrix

  double   * REAL_LAPACK_SCR;
  dcomplex * CMPLX_LAPACK_SCR;

  // Memory Lengths
  int lenScr_;
  int lenOTrans1_;
  int lenOTrans2_;
  int lenPOA_;
  int lenPOB_;
  int lenPOAsav_;
  int lenPOBsav_;
  int lenFOA_;
  int lenFOB_;
  int lenInitMOA_;
  int lenInitMOB_;
  int lenUTransA_;
  int lenUTransB_;
  int lenScratch_;

  int lenREAL_LAPACK_SCR;
  int lenCMPLX_LAPACK_SCR;
  int lWORK;

  std::array<double,3> EDField_;
  
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
    if(this->lenScr_ == 0)
      CErr("Fatal: RealTime given no scratch space",this->fileio_->out);
  }

  void initRTPtr();
  void initMemLen();
  void initMem();
  void initMaps();

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

    this->initRTPtr();

    this->isClosedShell_ = false;

    // Standard Values
    this->frozenNuc_   = true;
    this->nTCS_        = 1;
    this->maxSteps_    = 10;
    this->stepSize_    = 0.05;
    this->typeOrtho_   = Lowdin;
    this->initDensity_ = 0;
    this->swapMOA_     = 0;
    this->swapMOB_     = 0;
    this->methFormU_   = EigenDecomp;
    this->EDField_     = {0.0,0.0,0.0};
    this->Freq_        = 0.0;
    this->Phase_       = 0.0;
    this->Sigma_       = 0.0;
    this->TOn_         = 0.0;
    this->TOff_        = 1.0e4;
    this->IEnvlp_      = Constant;
  };
  ~RealTime() {;};

  enum ORTHO {
    Lowdin,
    Cholesky,
    Canonical 
  };

  enum FORM_U {
    EigenDecomp,
    Taylor
  };

  enum ENVELOPE {
    Constant, 
    LinRamp,
    Gaussian,
    Step,
    SinSq
  };

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

    this->initMemLen();
  }

  void alloc();

  // Getters

  // Setters
  inline void setMaxSteps(int i){ this->maxSteps_  = i;};
  inline void setStepSize(int i){ this->stepSize_ = i;};
  inline void setOrthoTyp(RealTime<T>::ORTHO i){ this->typeOrtho_ = i;};
  inline void setInitDen(int i){ this->initDensity_ = i;};
  inline void setSwapMOA(int i){ this->swapMOA_     = i;};
  inline void setSwapMOB(int i){ this->swapMOB_     = i;};
  inline void setFormU(RealTime<T>::FORM_U i){ this->methFormU_ = i;};
  inline void setEnvelope(RealTime<T>::ENVELOPE i){ this->IEnvlp_ = i;};
  inline void setFieldAmp(std::array<double,3> x){ 
    this->Ex_ = x[0];
    this->Ey_ = x[1];
    this->Ez_ = x[2];
  };
  inline void setTOn(double x){   this->TOn_   = x;};
  inline void setTOff(double x){  this->TOff_  = x;};
  inline void setFreq(double x){  this->Freq_  = x;};
  inline void setPhase(double x){ this->Phase_ = x;};
  inline void setSigma(double x){ this->Sigma_ = x;};
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
#include <realtime_print.h>
} // namespace ChronusQ
#endif
