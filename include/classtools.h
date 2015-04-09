#ifndef  INCLUDED_CLASSTOOLS
#define  INCLUDED_CLASSTOOLS
#include "global.h"
#include "fileio.h"
#include "molecule.h"
#include "basisset.h"
#include "controls.h"

/*****************************/
/*Error Messages 15000-19999 */
/*****************************/

// read input files and initialize everything
void readInput(FileIO*,ChronusQ::Molecule*,BasisSet*,Controls*);

// trace of product of two symmetric matrices
double traceSymm(Matrix<double>*,Matrix<double>*);

#endif
