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
#ifndef INCLUDED_CONTROLS
#define INCLUDED_CONTROLS
#include <global.h>
#include <cerr.h>
#include <tools.h>

/****************************/
/* Error Messages 6000-6999 */
/****************************/

namespace ChronusQ {
class Controls {

public:

  int   	printLevel;	// print level
  bool  	energyOnly;     // compute energy only
  bool  	optWaveFunction;// optimize wave function
  bool  	optGeometry;    // optimize geometry
  bool  	firstDer;      	// compute the first derivative
  bool  	secondDer; 	// compute the second derivative
  bool  	HF;             // use Hartree-Fock
  bool  	DFT;            // use density functional theory
  bool  	hybridDFT; 	// DFT is a hybrid functional
  bool          doTCS;          // Use Two-Component Spinor Basis
  bool  	restart;        // restart the calculation
  bool  	directTwoE;     // if direct two-electron will performed
  bool          buildn4eri;     // Build N^4 AO ERI tensor
  bool          doDF;           // Density fitting (RI) flag
  bool          doSDR;          // Perform Response calculation
  bool          doCUHF;         // To CUHF (Eigenfunction of S x S)
  bool          doComplex;      // Find complex solution
  bool          doUnit;         // Print consolidated infor for unit tests
  double 	thresholdS;
  double 	thresholdAB;
  double	thresholdSchawrtz;
  double        SCFdenTol_;
  double        SCFeneTol_;
  std::array<double,3> field_;
  int           SCFmaxIter_;
  int    	guess;         	// how to get the initial guess
  int           nthreads;       // Number of OpenMP threads
  int           SDMethod;
  bool          doRealTime;     // do RT-TDSCF?
  int           rtMaxSteps;     // Max number steps for RT-TDSCF
  double        rtTimeStep;     // Size of time step for RT-TDSCF
  int           rtTypeOrtho;    // Type of orthogonalization for RT-TDSCF.
  int           rtInitDensity;  // Initial Density for RT-TDSCF
  int           rtSwapMOA;      // which alpha MOs to swap for RT-TDSCF
  int           rtSwapMOB;      // which beta MOs to swap for RT-TDSCF
  int           rtMethFormU;    // Which method to form U propagator for RT-TDSCF
  std::array<double,3> rtField_;// Real Time Electric field magnitude Ex,Ey,Ez 
  double        rtFreq_;        // Frequency of electric field
  double        rtPhase_;       // Phase offset for electromagnetic field
  double        rtSigma_;       // Phase offset for electromagnetic field
  double        rtTOn_;         // Time field is applied
  double        rtTOff_;        // Time field is turned off 
  int           rtEnvelope_;    // Envelope function for RT EM field
  int           SDNSek;
  int           unitTest;       // Which class of unit test to perform
  std::string   gauFChkName;	// Gaussian formatted checkpoint filename
  std::string   gauMatElName;   // Gaussian raw matrix element file

  // Enum for Unit Tests
  enum{
    UnitSCF,
    UnitResp,
    UnitRT
  };

  Controls(){;};
  ~Controls(){;};
  void iniControls();
  void readSMP(int &);
  void readPSCF(std::fstream &,std::fstream &);
  void readDebug(std::string);
  void printSettings(std::ostream &out=cout);
//inline void printSettings(){this->printSettings(cout);};
};
} // namespace ChronusQ
#endif
