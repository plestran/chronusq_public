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
#ifndef EIGINTER_INC
#define EIGINTER_INC
//----------------//
//String constants//
//----------------//
const std::string bannerTop="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";
const std::string bannerMid="--------------------------------------------------------------------------------";
const std::string bannerEnd="xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx";


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
