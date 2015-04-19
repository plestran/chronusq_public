#include "global.h"
#ifndef EIGINTER_INC
#define EIGINTER_INC
namespace Eigen {
  template<typename Derived>
  void prettyPrint(std::ostream & output, const Derived& m, std::string str){
  int list = 5;
  int i,j,k,n,end,endLT;
  output.precision(10);
  output.fill(' ');
  output.setf(std::ios::right,std::ios::scientific);
  output.setf(std::ios::fixed,std::ios::floatfield);
  output << std::endl << str + ": " << std::endl;
  output<<bannerTop;

  for(i=0;i<m.cols();i+=list) {
    output<<std::endl;
    end = list;
    output<<std::setw(5)<<" ";
    if((i+list)>=m.cols()) end = m.cols() - i;
    for(k=i;k<i+end;k++) output<<std::setw(15)<<k+1;
    output<<std::endl;
    for(j=0;j<m.rows();j++) {
      output<<std::setw(5)<<j+1;
      for(n=i;n<i+end;n++) output<<std::setw(15)<<m.data()[j*m.rows()+n]; // n in the column index (dbwy)
      output<<std::endl;
    };
  };
  output<<bannerEnd<<std::endl;
  }
}
#endif
