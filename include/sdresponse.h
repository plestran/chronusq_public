#ifndef INCLUDED_SDRESPONSE
#define INCLUDED_SDRESPONSE
#include "global.h"
#include "matrix.h"
#include "molecule.h"
#include "controls.h"
#include "mointegrals.h"
#include "singleslater.h"
#include "basisset.h"

/****************************/
/* Error Messages 5000-5999 */
/****************************/

namespace ChronusQ {
class SDResponse {
  int       nBasis_;
  int       **R2Index_;
  int       nStates_;
  ChronusQ::BasisSet     	*basisSet_;
  ChronusQ::Molecule    	*molecule_;
  ChronusQ::FileIO       	*fileio_;
  ChronusQ::Controls     	*controls_;
  ChronusQ::MOIntegrals   *mointegrals_;
  SingleSlater  *singleSlater_;

public:
 
  // constructor & destructor
  SDResponse(){;};
  ~SDResponse() {;};
  // pseudo-constructor
  void iniSDResponse(ChronusQ::Molecule*,ChronusQ::BasisSet*,ChronusQ::MOIntegrals*,ChronusQ::FileIO*,ChronusQ::Controls*,SingleSlater*);

  void computeExcitedStates();         // compute the total electronic energy
  void printExcitedStateEnergies(); 
  void printInfo();

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagSDResponse);
  void mpiRecv(int,int tag=tagSDResponse);
};
}
#endif
