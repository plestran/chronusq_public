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
