#ifndef INCLUDED_GLOBAL
#define INCLUDED_GLOBAL

#define MAXNAMELEN 50
#define MAXANGULARMOMENTUM 6
#define MAXCONTRACTION 10
#define MAXATOMS 1000
//#include "oompi.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <sys/stat.h>
#include <pthread.h>
//#include <error.h>
#include <time.h>
//#include <new.h>
#include <cstring>
#include <complex>

//using namespace std;

using std::cout;
using std::endl;
using std::cin;
using std::fstream;
using std::ostream;
using std::ios;
using std::abs;
using std::nothrow;
using std::ifstream;

typedef std::complex<double> dcomplex;

//----------------//
//number constants//
//----------------//
struct Math {
  double zero, one, two, three, four, five, six, seven, eight, nine, ten, half, quarter;
  double sqrt2;
  double pi,pi32,sqrt2pi54; //pi, pi^{3/2} sqrt(2)*pi^{5/4}
  double small;
};
const Math math = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 0.5, 0.25,
		   1.4142135623731,
		   3.14159265358979,5.56832799683171,5.91496717279561,
		   1.0e-10};
//factorials n!
static double Factorial[8] ={
1,1,2, 6, 24, 120, 720, 5040 };

//double factorials n!!
static double dFactorial[21] = {
  1.00000000,
  1.00000000,
  3.00000000,
  15.00000000,
  105.00000000,
  945.00000000,
  10395.00000000,
  135135.00000000,
  2027025.00000000,
  34459425.00000000,
  654729075.00000000,
  13749310575.00000000,
  316234143225.00000000,
  7905853580625.00000000,
  213458046676875.00000000,
  6190283353629375.00000000,
  191898783962510624.00000000,
  6332659870762850304.00000000,
  221643095476699758592.00000000,
  8200794532637890838528.00000000,
  319830986772877752139776.00000000,
};

//------------------//
//Physical constants//
//------------------//
struct Phys {
  //Bohr radius per Angstrom
  double bohr;
  //number of cartesian AO's in a shell
};
const Phys phys = {0.5291772083000001};

//----------------//
//String constants//
//----------------//
const char bannerTop[100]="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";
const char bannerMid[100]="--------------------------------------------------------------------------------";
const char bannerEnd[100]="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";

//------------------//
// IO block numbers //
//------------------//
enum {blockControlFlags,blockMolecule,blockBasisSet, // 1,2,3
      blockSingleSlater,blockIntegrals};  // 4,5,6

//----------//
//MPI global//
//----------//

struct GlobalMPI {
  int  myid;
  int  size;
//  char nodeName[MPI_MAX_PROCESSOR_NAME];;
  int  nodeNameLen;
};
//---------//
//MPI Tags //
//---------//
enum {tagMolecule,tagMatrix,tagBasisSet,tagSingleSlater,tagIntegrals,tagSDResponse};


#endif
