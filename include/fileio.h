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
#ifndef INCLUDED_FILEIO
#define INCLUDED_FILEIO
#define MAXBLOCK 1000
#include <global.h>
#include <tools.h>

/****************************/
/* Error Messages 1000-1999 */
/****************************/
namespace ChronusQ {
class FileIO {

  char *name_in;                 // name of the input file
  char *name_out;                // name of the output file
  char *name_scr;                // name of the scratch file
  char *name_bin;                // name of the binary file
  int   sizeInt_;                // size of an integer in terms of sizeof(char)
  int   sizeLong_;               // size of a long integer in terms of sizeof(char)
  int   sizeFloat_;              // size of a float point in terms of sizeof(char)
  int   sizeDouble_;             // size of a double precision float point in terms of sizeof(char)
  long  block[MAXBLOCK];         // block pointers
  long  readPointer[2*MAXBLOCK]; // read pointers
  long  writePointer[2*MAXBLOCK];// write pointers

public:

  fstream in;                    // file handler of the input file
  fstream out;                   // file handler of the output file
  fstream scr;                   // file handler of the scratch file
  fstream bin;                   // file handler of the binary file

  // constructor and destructor
  FileIO(char *);
  void iniFileIO(bool);
  ~FileIO();

  // open files
  inline void openIn() {
    if(this->in.is_open()) this->in.close();
    this->in.open(this->name_in,ios::in);
    if(this->in.fail()) throw 1020;
  };

  // determine if a block is initialized
  inline bool isOpen(int blockNumber) {
    if(blockNumber<0||blockNumber>MAXBLOCK) throw 1008;
    if(block[blockNumber]<=0) return false;
    return true;
  };

  // access to private data
  int sizeInt()    { return this->sizeInt_;};
  int sizeLong()   { return this->sizeLong_;};
  int sizeFloat()  { return this->sizeFloat_;};
  int sizeDouble() { return this->sizeDouble_;};

  // read||write block pointers from||to binary file
  void iniBlock(int);
  void readBlock();
  void writeBlock();

  // read||write scratch and binary files
  void io(char*, int, char*, char*, int, int offset=-1,char*offsetType=NULL);
  void io(char*, int, char*, int*, int, int offset=-1,char*offsetType=NULL);
  void io(char*, int, char*, long*, int, int offset=-1,char*offsetType=NULL);
  void io(char*, int, char*, float*, int, int offset=-1,char*offsetType=NULL);
  void io(char*, int, char*, double*, int, int offset=-1,char*offsetType=NULL);
  void rw(char*, int, char*, char*, int, int offset=-1);
  int charOffset(int,char*);
};

} // namespace ChronusQ
#endif
