/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explictly 
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
#include "classtools.h"
using ChronusQ::Molecule;
using ChronusQ::BasisSet;

/************************************************/
/* read input files and initialize everything   */
/* memory allocations are done here too         */
/************************************************/
namespace ChronusQ {
void readInput(FileIO *fileio, Molecule *mol, BasisSet *basis,Controls *controls) {
  int i, j, n, readInt;
  char readString[MAXNAMELEN];
  fileio->in >> readString;
  while(!(fileio->in.eof())) {
    strupr(readString);
    if(!strcmp(readString,"/*")) while(strcmp(readString,"*/")) fileio->in >> readString;
    else if(!strcmp(readString,"$CHARGE")) {
      fileio->in >> readInt;
      mol->readCharge(readInt);
    } else if(!strcmp(readString,"$SPIN")) {
      fileio->in >> readInt;
      mol->readSpin(readInt);
    } else if(!strcmp(readString,"$EXTRA")) {
      fileio->in>>readString;
      strupr(readString);
      if(!strcmp(readString,"RESTART")) controls->restart=true;
    } else if(!strcmp(readString,"$PRINT")) {
      fileio->in >> (controls->printLevel);
    } else if(!strcmp(readString,"$METHOD")) {
      fileio->in >> readString;
      strupr(readString);
      if(!strcmp(readString,"HF")) controls->HF=true;
      else {
	controls->HF=false;
	controls->DFT=true;
      };
    } else if(!strcmp(readString,"$GEOM")) {
      try{mol->readMolecule(fileio);}
      catch(int msg) {
	fileio->out<<"Reading molecular information failed! E#: "<<msg<<endl;
	exit(1);
      }; 
    } else if(!strcmp(readString,"$BASIS")) {
      basis->readBasisSet(fileio,mol);
    } else if(!strcmp(readString,"$GUESS")) {
      fileio->in>>readString;
      strupr(readString);
      if(!strcmp(readString,"INPUT")) {
        controls->guess = 1;
      } else if(!strcmp(readString,"GAUFCHK")) {
        controls->guess = 3;
        fileio->in>>controls->gauFChkName;
      };
    };
    fileio->in >> readString;
  };
};
/********************************/
/* trace of producto of two     */
/* symmetric Matrix<double>*ces           */
/********************************/
double traceSymm(Matrix<double>* a, Matrix<double>* b) {
  if(a->len()!=b->len()) throw 15001;
  double tmpVal = 0.0;
  int i;
  if(a->format()==0) for(i=0;i<a->len();i++) tmpVal+=a->data()[i]*b->data()[i];
  else if (a->format()==1){
    for(i=0;i<a->len();i++) {
      tmpVal+=double(2.0)*a->data()[i]*b->data()[i];
    };
    for(i=0;i<a->rows();i++) {
      tmpVal-=(*a)(i,i)*(*b)(i,i);
    };
  };
  return tmpVal;
};
} // namespace ChronusQ
