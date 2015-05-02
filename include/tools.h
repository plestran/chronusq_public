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
#ifndef  INCLUDED_TOOLS
#define  INCLUDED_TOOLS
//#include <gsl/gsl_sf_erf.h>
#include <global.h>
#include <cerr.h>

/*****************************/
/*Error Messages 10000-14999 */
/*****************************/

namespace ChronusQ {
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

// make a lowercase of the string:
std::string stringlower(std::string);

//make a uppercaser of the string:
std::string stringupper(std::string);
// convert atomic symbol and mass to an index in the table of atom in atoms.h
int HashAtom (std::string, int);

// convert shell symbol or angular momentum to number of AOs in the shell
int HashNAOs (int);
int HashNAOs (char*);

// convert shell symbol to angular momentum
int HashL (std::string);

// convert angular momentum to AO index in a shell
int HashIAO(int,int*);


} // namespace ChronusQ
#endif
