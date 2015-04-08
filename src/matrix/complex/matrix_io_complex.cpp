#include "matrix.h"
//-----------------------//
// print matrix elements //
//-----------------------//
template<>
void Matrix<dcomplex>::printAll(int list, ostream &output) {
  int i,j,k,n,end,endLT;
  output.precision(10);
  output.fill(' ');
  output.setf(ios::right,ios::scientific);
  output.setf(ios::fixed,ios::floatfield);
  output<<endl<<"Re["<<name_<<"] :"<<endl;
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
          for(n=i;n<i+end;n++) output<<std::setw(15)<<real(data_[n*rows_+j]); // n in the column index (dbwy)
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
          for(n=i;n<i+end;n++) output<<std::setw(15)<<real(data_[n+j*rows_]); // n in the column index (dbwy)
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
	for(n=i;n<endLT;n++) output<<std::setw(15)<<real(data_[n*(this->rows_)-n*(n-1)/2+j-n]);
	output<<endl;
      };
    };
  };

  output<<bannerEnd<<endl;
  output<<endl<<"Im["<<name_<<"] :"<<endl;
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
          for(n=i;n<i+end;n++) output<<std::setw(15)<<imag(data_[n*rows_+j]); // n in the column index (dbwy)
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
          for(n=i;n<i+end;n++) output<<std::setw(15)<<imag(data_[n+j*rows_]); // n in the column index (dbwy)
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
	for(n=i;n<endLT;n++) output<<std::setw(15)<<imag(data_[n*(this->rows_)-n*(n-1)/2+j-n]);
	output<<endl;
      };
    };
  };
  output<<bannerEnd<<endl;
};
