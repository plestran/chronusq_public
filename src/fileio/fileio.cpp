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
#include <fileio.h>
using ChronusQ::FileIO;
//-------------//
// constructor //
//-------------//
FileIO::FileIO(std::string nm_input) {
  char   testChar;
  int    testInt;
  long   testLong;
  float  testFloat;
  double testDouble;
  this->sizeInt_    = sizeof(testInt)   /sizeof(testChar);
  this->sizeLong_   = sizeof(testLong)  /sizeof(testChar);
  this->sizeFloat_  = sizeof(testFloat) /sizeof(testChar);
  this->sizeDouble_ = sizeof(testDouble)/sizeof(testChar);

  if(nm_input.empty()) {
    CErr("Fatal: Input File Required");
  } else {
    this->name_in = nm_input + ".inp";
    this->in.open(name_in,ios::in);
    if(this->in.fail()) CErr("Unable to open "+this->name_in);

    this->name_out = nm_input + ".out";
    this->out.open(name_out,ios::out);
    if(this->out.fail()) CErr("Unable to open "+this->name_out);

    // scratch file and binary files will be initialized in iniFileIO
    this->name_scr = nm_input + ".scr";
    this->name_bin = nm_input + ".bin";
  };
};
FileIO::FileIO(std::vector<std::string> nm_input) {
  char   testChar;
  int    testInt;
  long   testLong;
  float  testFloat;
  double testDouble;
  this->sizeInt_    = sizeof(testInt)   /sizeof(testChar);
  this->sizeLong_   = sizeof(testLong)  /sizeof(testChar);
  this->sizeFloat_  = sizeof(testFloat) /sizeof(testChar);
  this->sizeDouble_ = sizeof(testDouble)/sizeof(testChar);

  if(nm_input.size()==0) CErr("Fatal: Input File Required");
  cout << "HERE" << endl;

  std::string inputTag = "--inp=";
  std::string outputTag = "--out=";
  std::string scrTag = "--scr=";
  std::string binTag = "--bin=";
  for(auto i = 0; i < nm_input.size(); ++i) {
    if(!nm_input[i].compare(0,inputTag.length(),inputTag)){
      nm_input[i].erase(0,inputTag.length());
      this->name_in = nm_input[i];
    } else if(!nm_input[i].compare(0,outputTag.length(),outputTag)){
      nm_input[i].erase(0,outputTag.length());
      this->name_out = nm_input[i];
    } else if(!nm_input[i].compare(0,scrTag.length(),scrTag)){
      nm_input[i].erase(0,scrTag.length());
      this->name_scr = nm_input[i];
    } else if(!nm_input[i].compare(0,binTag.length(),binTag)){
      nm_input[i].erase(0,binTag.length());
      this->name_bin = nm_input[i];
    } else CErr("Input \""+nm_input[i]+"\" not recognized");
  }
  cout << "HERE" << endl;

  if(this->name_in.empty()) CErr("Fatal: Must specify an input file");
  if(this->name_out.empty()) this->name_out = this->name_in + ".out";
  if(this->name_scr.empty()) this->name_scr = this->name_in + ".scr";
  if(this->name_bin.empty()) this->name_bin = this->name_in + ".bin";

  this->in.open(name_in,ios::in);
  this->out.open(name_out,ios::out);

  if(this->in.fail()) CErr("Unable to open "+this->name_in);
  if(this->out.fail()) CErr("Unable to open "+this->name_out);
};
//------------//
// destructor //
//------------//
FileIO::~FileIO() {
  if(in.is_open()) in.close();
  if(scr.is_open()) scr.close();
  if(bin.is_open()) bin.close();
  if(scr.is_open()) scr.close();
  if(bin.is_open()) bin.close();
  if(out.is_open()) out.close();
};
/*
//-------------------//
// initialize FileIO //
//-------------------//
void FileIO::iniFileIO(bool restart) {
  int i;
  scr.open(name_scr,ios::out|ios::binary);
  scr.close();
  scr.open(name_scr,ios::in|ios::out|ios::binary);
  if(scr.fail()) throw 1004;

  if(restart) bin.open(name_bin,ios::in|ios::out|ios::binary);
  else {
    bin.open(name_bin,ios::out|ios::binary);
    bin.close();
    bin.open(name_bin,ios::in|ios::out|ios::binary);
  };
  if(bin.fail()) throw 1005;

  if(restart) {
    readBlock();
  } else {
    for(i=0;i<MAXBLOCK;i++) block[i] = -1;
    writeBlock();
  };
  for(i=0;i<MAXBLOCK;i++) {
    readPointer[i]  = 0;
    writePointer[i] = 0;
    readPointer[i+MAXBLOCK]  = 0;
    writePointer[i+MAXBLOCK] = 0;
  };
};
//-------------------------------------------------//
// initialize blocks                               //
//   only assign block pointer                     //
//   should be used with immediate write operation //
//-------------------------------------------------//
void FileIO::iniBlock(int blockNumber){
  if(blockNumber<0||blockNumber>MAXBLOCK) throw 1006;
  bin.seekp(0,ios::end);
  block[blockNumber]=bin.tellp();
  block[blockNumber]+=1;
};
//-------------------------------------------------//
// read||write block pointers from||to binary file //
//   block information is always stored at the     //
//   beginning of the binary file                  //
//-------------------------------------------------//
void FileIO::readBlock(){
  int charLen = MAXBLOCK*sizeLong_;
  char *charStorage;
  charStorage = (char *)block;
  bin.seekg(0,ios::beg);
  bin.read(charStorage,charLen);
};
void FileIO::writeBlock(){
  int charLen = MAXBLOCK*sizeLong_;
  char *charStorage;
  charStorage = (char *)block;
  bin.seekp(0,ios::beg);
  bin.write(charStorage,charLen);
};
//---------------------------------------------------------------------------------------------//
// read & write scratch and binary files                                                       //
//  io() is a wrapper to convert different types of variables to character storage requirement //
//  rw() is operates read|write on scratch and binary files                                    //
//---------------------------------------------------------------------------------------------//
void FileIO::io(char *op, int blockNumber, char *file, char *storage, int len, int offset, char *offsetType) {
  int cOffset = offset;
  if(offset>0&&offsetType!=NULL) cOffset = charOffset(offset,offsetType);
  try {this->rw(op, blockNumber, file, storage, len, cOffset);}
  catch(int msg) {
    this->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);
  };
};
void FileIO::io(char *op, int blockNumber, char *file, int *storage, int len, int offset, char *offsetType) {
  int charLen = len*sizeInt_;
  int cOffset = offset*sizeInt_;
  if(offset>0&&offsetType!=NULL) cOffset = charOffset(offset,offsetType);
  char *charStorage;
  charStorage = (char *)storage;
  try {this->rw(op, blockNumber, file, charStorage, charLen, cOffset);}
  catch(int msg) {
    this->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);
  };
};
void FileIO::io(char *op, int blockNumber, char *file, long *storage, int len, int offset, char *offsetType) {
  int charLen = len*sizeLong_;
  int cOffset = offset*sizeLong_;
  if(offset>0&&offsetType!=NULL) cOffset = charOffset(offset,offsetType);
  char *charStorage;
  charStorage = (char *)storage;
  try {this->rw(op, blockNumber, file, charStorage, charLen, cOffset);}
  catch(int msg) {
    this->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);
  };
};
void FileIO::io(char *op, int blockNumber, char *file, float *storage, int len, int offset, char *offsetType) {
  int charLen = len*sizeFloat_;
  int cOffset = offset*sizeFloat_;
  if(offset>0&&offsetType!=NULL) cOffset = charOffset(offset,offsetType);
  char *charStorage;
  charStorage = (char *)storage;
  try {this->rw(op, blockNumber, file, charStorage, charLen, cOffset);}
  catch(int msg) {
    this->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);
  };
};
void FileIO::io(char *op, int blockNumber, char *file, double *storage, int len, int offset, char *offsetType) {
  int charLen = len*sizeDouble_;
  int cOffset = offset*sizeDouble_;
  if(offset>0&&offsetType!=NULL) cOffset = charOffset(offset,offsetType);
  char *charStorage;
  charStorage = (char *)storage;
  try {this->rw(op, blockNumber, file, charStorage, charLen, cOffset);}
  catch(int msg) {
    this->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);
  };
};
void FileIO::rw(char *op, int blockNumber, char *file, char *charStorage, int charLen, int offset) {
  if(block[blockNumber]<0) {
    out<<"Error: block "<<blockNumber<<" is empty!"<<endl;
    throw 1007;
  };
  strupr(op);
  strupr(file);
  if(!strcmp(op,"R")||!strcmp(op,"READ")) {
    if(!strcmp(file,"SCR")) {
      if(offset>=0) readPointer[blockNumber+MAXBLOCK] = offset;
      this->scr.seekg(block[blockNumber]+readPointer[blockNumber]);
      this->scr.read(charStorage,charLen);
      readPointer[blockNumber+MAXBLOCK] += charLen;
    } else if(!strcmp(file,"BIN")) {
      if(offset>=0) readPointer[blockNumber] = offset;
      this->bin.seekg(block[blockNumber]+readPointer[blockNumber],ios::beg);
      this->bin.read(charStorage,charLen);
      readPointer[blockNumber] += charLen;
    } else throw 1017;
  } else if(!strcmp(op,"W")||!strcmp(op,"WRITE")) {
    if(!strcmp(file,"SCR")) {
      if(offset>=0) writePointer[blockNumber+MAXBLOCK] = offset;
      this->scr.seekp(block[blockNumber]+writePointer[blockNumber]);
      this->scr.write(charStorage,charLen);
      writePointer[blockNumber+MAXBLOCK] += charLen;
    } else if(!strcmp(file,"BIN")) {
      if(offset>=0) writePointer[blockNumber] = offset;
      this->bin.seekp(block[blockNumber]+writePointer[blockNumber],ios::beg);
      this->bin.write(charStorage,charLen);
      writePointer[blockNumber] += charLen;
    } else throw 1018;
  } else throw 1019;
};
//----------------------------------------//
// return offset in terms of sizeof(char) //
//----------------------------------------//
int FileIO::charOffset(int offset, char *offsetType) {
  strupr(offsetType);
  if(!strcmp(offsetType,"CHAR")) return offset;
  else if(!strcmp(offsetType,"INT")) return offset*sizeInt_;
  else if(!strcmp(offsetType,"LONG")) return offset*sizeInt_;
  else if(!strcmp(offsetType,"FLOAT")) return offset*sizeInt_;
  else if(!strcmp(offsetType,"DOUBLE")) return offset*sizeInt_;
  out<<"Unrecognized offset type! E#:"<<1020<<endl;
  exit(1);
};
*/
