#ifndef  INCLUDED_TOOLS
#define  INCLUDED_TOOLS
//#include <gsl/gsl_sf_erf.h>
#include "global.h"

/*****************************/
/*Error Messages 10000-14999 */
/*****************************/

// factorial function:  t! = 1*2*3...*t
double  factorial (int);

// double factorial:  (2t-1)!! = 1*3*5*..*(2t-1)
double  doubleFact (int);

// (base)^order:  0**0 == 1
double  powerInt (double,int);

// polynomial coeff: (x+a)^l= sum(coeff * x^i * a^(l-i))
double  polyCoeff (int, int);

// (x+pA)^a*(x+pB)^b = sum[coeff * x^k]; k<=a+b
double  kCoeff (int, int, int, double, double);

// make an uppercase copy of s:
void strupr (char*);

// make a lowercase copy of s:
void strlwr (char*);

// convert atomic symbol and mass to an index in the table of atom in atoms.h
int HashAtom (char*, int);

// convert shell symbol or angular momentum to number of AOs in the shell
int HashNAOs (int);
int HashNAOs (char*);

// convert shell symbol to angular momentum
int HashL (char*);

// convert angular momentum to AO index in a shell
int HashIAO(int,int*);

#endif
