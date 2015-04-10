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
#include "global.h"

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
  bool  	restart;        // restart the calculation
  bool  	directTwoE;     // if direct two-electron will performed
  double 	thresholdS;
  double 	thresholdAB;
  double	thresholdSchawrtz;
  int    	guess;         	// how to get the initial guess
  char   	gauFChkName[MAXNAMELEN];	// Gaussian formatted checkpoint filename

  Controls(){;};
  ~Controls(){;};
  void iniControls();
};
} // namespace ChronusQ
#endif
