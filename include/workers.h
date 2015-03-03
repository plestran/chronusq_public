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

int atlas(int, char**, GlobalMPI*);
int worker(GlobalMPI*);

#endif
