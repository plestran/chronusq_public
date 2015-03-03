#include "fileio.h"
//-------------//
// constructor //
//-------------//
FileIO::FileIO(char *nm_input) {
  char   testChar;
  int    testInt;
  long   testLong;
  float  testFloat;
  double testDouble;
  sizeInt_    = sizeof(testInt)   /sizeof(testChar);
  sizeLong_   = sizeof(testLong)  /sizeof(testChar);
  sizeFloat_  = sizeof(testFloat) /sizeof(testChar);
  sizeDouble_ = sizeof(testDouble)/sizeof(testChar);

  name_in   = new char[MAXNAMELEN];
  name_out  = new char[MAXNAMELEN];
  name_scr  = new char[MAXNAMELEN];
  name_bin  = new char[MAXNAMELEN];

  if(nm_input==NULL) {
    this->out<<"Error: input filename required! "<<endl;
    throw 1001;
  } else {
    strcpy(name_in,nm_input);
    strcat(name_in,".inp");
    in.open(name_in,ios::in);
    if(in.fail()) throw 1002;

    strcpy(name_out,nm_input);
    strcat(name_out,".out");
    out.open(name_out,ios::out);
    if(out.fail()) throw 1003;

    // scratch file and binary files will be initialized in iniFileIO
    strcpy(name_scr,nm_input);
    strcat(name_scr,".scr");
    strcpy(name_bin,nm_input);
    strcat(name_bin,".bin");
  };
};
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
//------------//
// destructor //
//------------//
FileIO::~FileIO() {
  if(in.is_open()) in.close();
  if(scr.is_open()) scr.close();
  if(bin.is_open()) bin.close();
  if(scr.is_open()) scr.close();
  if(bin.is_open()) bin.close();
  if(remove(name_scr)!=0) out<<"Error when deleting scratch file!"<<endl;
  if(out.is_open()) out.close();
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
