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
#include <classtools.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;

/************************************************/
/* read input files and initialize everything   */
/* memory allocations are done here too         */
/************************************************/
namespace ChronusQ {
void readInput(std::shared_ptr<FileIO> fileio, std::shared_ptr<Molecule> mol, 
               std::shared_ptr<BasisSet> basis,std::shared_ptr<Controls> controls) {
  int i, j, n, readInt;
  std::string readString;
  fileio->in >> readString;
  while(!(fileio->in.eof())) {
    readString=stringupper(readString);
    if(!readString.compare("/*")) while(readString.compare("*/")) fileio->in >> readString;
    else if(!readString.compare("$CHARGE")) {
      fileio->in >> readInt;
      mol->readCharge(readInt);
    } else if(!readString.compare("$SPIN")) {
      fileio->in >> readInt;
      mol->readSpin(readInt);
    } else if(!readString.compare("$EXTRA")) {
      fileio->in>>readString;
      readString=stringupper(readString);
      if(!readString.compare("RESTART")) controls->restart=true;
    } else if(!readString.compare("$PRINT")) {
      fileio->in >> (controls->printLevel);
    } else if(!readString.compare("$METHOD")) {
      fileio->in >> readString;
      readString=stringupper(readString);
      if(!readString.compare("HF")) controls->HF=true;
      else {
	controls->HF=false;
	controls->DFT=true;
      };
    } else if(!readString.compare("$GEOM")) {
      fileio->in >> readString;
      readString=stringupper(readString);  
      fstream *geomRead;
      if(!readString.compare("READ")) {
        geomRead = &fileio->in;
      } else if(!readString.compare("FILE")) {
        fileio->in >> readString;
        geomRead = new fstream(readString,ios::in);
        if(geomRead->fail()) CErr("Unable to open "+std::string(readString),fileio->out); 
        else fileio->out << "Reading geometry from " << readString << endl;
      } else {
        CErr("Unrecognized GEOM option: " + std::string(readString),fileio->out);
      }
      mol->readMolecule(fileio,*geomRead);
      if(!(geomRead==&fileio->in)){
//      fileio->out << "Closing " << readString << endl;
        geomRead->close();
        delete geomRead;
      }
    } else if(!readString.compare("$BASIS")) {
      basis->readBasisSet(fileio,mol);
    } else if(!readString.compare("$NSMP")) {
      fileio->in >> readInt;
      controls->readSMP(readInt);
    } else if(!readString.compare("$GUESS")) {
      fileio->in>>readString;
      readString=stringupper(readString);  
      if(!readString.compare("INPUT")) {
        controls->guess = 1;
      } else if(!readString.compare("GAUFCHK")) {
        controls->guess = 3;
        fileio->in>>controls->gauFChkName;
      };
    } else if(!readString.compare("$SCF")) {
      fileio->in>>readString;
      readString=stringupper(readString);
      if(!readString.compare("OFF"))
        controls->optWaveFunction = false;
      else if(!readString.compare("ON"))
        controls->optWaveFunction = true;
    } else if(!readString.compare("$DEBUG")) {
      fileio->in>>readString;
      readString=stringupper(readString);
      controls->readDebug(readString);
    };
    fileio->in >> readString;
  };
};
/********************************/
/* trace of producto of two     */
/* symmetric RealMatrices       */
/********************************/
double traceSymm(RealMatrix* a, RealMatrix* b) {
  if(a->size()!=b->size()) CErr("Only able to trace matricies of the same size"); // FIXME this is depreciated
  double tmpVal = 0.0;
  int i;
  for(i=0;i<a->size();i++) tmpVal+=a->data()[i]*b->data()[i];

  return tmpVal;
};
} // namespace ChronusQ
