#include "matrix.h"
using ChronusQ::FileIO;
using ChronusQ::Matrix;
//-----------------------//
// print matrix elements //
//-----------------------//
namespace ChronusQ {
template<>
void Matrix<double>::printAll(int list, ostream &output) {
  int i,j,k,n,end,endLT;
  output.precision(10);
  output.fill(' ');
  output.setf(ios::right,ios::scientific);
  output.setf(ios::fixed,ios::floatfield);
  output<<endl<<name_<<" :"<<endl;
  output<<bannerTop;

  if(this->format_==0) {
    if(trans_=='N'){
      for(i=0;i<cols_;i+=list) {
        output<<endl;
        end = list;
        output<<std::setw(5)<<" ";
        if((i+list)>=cols_) end = cols_ - i;
        for(k=i;k<i+end;k++) output<<std::setw(15)<<k+1;
        output<<endl;
        for(j=0;j<rows_;j++) {
          output<<std::setw(5)<<j+1;
          for(n=i;n<i+end;n++) output<<std::setw(15)<<data_[n*rows_+j]; // n in the column index (dbwy)
          output<<endl;
        };
      };
    } else if(trans_=='T'){
      for(i=0;i<rows_;i+=list) {
        output<<endl;
        end = list;
        output<<std::setw(5)<<" ";
        if((i+list)>=rows_) end = rows_ - i;
        for(k=i;k<i+end;k++) output<<std::setw(15)<<k+1;
        output<<endl;
        for(j=0;j<cols_;j++) {
          output<<std::setw(5)<<j+1;
          for(n=i;n<i+end;n++) output<<std::setw(15)<<data_[n+j*rows_]; // n in the column index (dbwy)
          output<<endl;
        };
      };
    } else { throw 3013;}
  }
  else if (this->format_==1) {
    for(i=0;i<cols_;i+=list) {
      output<<endl;
      end = i+list;
      output<<std::setw(5)<<" ";
      if((i+list)>=cols_) end = cols_;
      for(k=i;k<end;k++) output<<std::setw(15)<<k+1;
      output<<endl;
      for(j=i;j<rows_;j++) {
	output<<std::setw(5)<<j+1;
	endLT = j+1;
	if(j>=end) endLT = end;
	for(n=i;n<endLT;n++) output<<std::setw(15)<<data_[n*(this->rows_)-n*(n-1)/2+j-n];
	output<<endl;
      };
    };
  };
  output<<bannerEnd<<endl;
};
//--------------------------------//
// read from binary files //
//--------------------------------//
template<>
void Matrix<double>::ioRead(FileIO *fileio,int blockNumber,char *fileName,int charOffset) {
  const int nInteger = 6;
  int storage[nInteger];
  try { fileio->io("R",blockNumber,fileName,storage,nInteger,charOffset);}
  catch (int msg) {
    fileio->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);    
  };
  rows_      = storage[0];
  cols_      = storage[1];
  len_       = storage[2];
  format_    = storage[3];
  size_      = storage[4];
  haveEigen_ = storage[5];
  try { fileio->io("R",blockNumber,fileName,name_,MAXNAMELEN);}
  catch (int msg) {
    fileio->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);    
  };

  if(data_!=NULL) delete[] data_;
  data_ = new (nothrow) double[len_];
  if (data_==NULL) throw 3003;
  try { fileio->io("R",blockNumber,fileName,data_,len_);}
  catch (int msg) {
    fileio->out<<"Operation on file failed in Matrix<double>::ioWrite! E#:"<<msg<<endl;
    exit(1);    
  };

  if(haveEigen_) {
    if(eigenvector_!=NULL) delete[] eigenvector_;
    if(eigenvalue_ !=NULL) delete[] eigenvalue_;
    eigenvector_ = new (nothrow) double[rows_*rows_];
    eigenvalue_  = new (nothrow) double[rows_];
    if(eigenvector_==NULL) throw 3004;
    if(eigenvalue_ ==NULL) throw 3005;
    try { fileio->io("R",blockNumber,fileName,eigenvector_,rows_*rows_);}
    catch (int msg) {
      fileio->out<<"Operation on file failed in Matrix<double>::ioWrite! E#:"<<msg<<endl;
      exit(1);    
    };
    try { fileio->io("R",blockNumber,fileName,eigenvalue_,rows_);}
    catch (int msg) {
      fileio->out<<"Operation on file failed in Matrix<double>::ioWrite! E#:"<<msg<<endl;
      exit(1);    
    };
  };
};
//-----------------------//
// write to binary files //
//----------------------//
template<>
void Matrix<double>::ioWrite(FileIO *fileio,int blockNumber,char *fileName,int charOffset) {
  const int nInteger = 6;
  int storage[nInteger];
  storage[0] = rows_;
  storage[1] = cols_;
  storage[2] = len_;
  storage[3] = format_;
  storage[4] = size_;
  storage[5] = haveEigen_;
  try { fileio->io("W",blockNumber,fileName,storage,nInteger,charOffset);}
  catch (int msg) {
    fileio->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);    
  };
  try { fileio->io("W",blockNumber,fileName,name_,MAXNAMELEN);}
  catch (int msg) {
    fileio->out<<"Operation on file failed! E#:"<<msg<<endl;
    exit(1);    
  };
  try { fileio->io("W",blockNumber,fileName,data_,len_);}
  catch (int msg) {
    fileio->out<<"Operation on file failed in Matrix<double>::ioWrite! E#:"<<msg<<endl;
    exit(1);    
  };
  if(haveEigen_) {
    try { fileio->io("W",blockNumber,fileName,eigenvector_,rows_*rows_);}
    catch (int msg) {
      fileio->out<<"Operation on file failed in Matrix<double>::ioWrite! E#:"<<msg<<endl;
      exit(1);    
    };
    try { fileio->io("W",blockNumber,fileName,eigenvalue_,rows_);}
    catch (int msg) {
      fileio->out<<"Operation on file failed in Matrix<double>::ioWrite! E#:"<<msg<<endl;
      exit(1);    
    };
  };
};
} // namespace ChronusQ
