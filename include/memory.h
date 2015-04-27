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
#ifndef INCLUDED_MEMORY
#define INCLUDED_MEMORY
#include <global.h>

namespace ChronusQ {
class Memory {
  unsigned long long iend_,len_;
  double *memory_;

public:

  //constructor
  Memory(long const, char *, ostream &output=cout);
  ~Memory() {delete[] memory_;};

  //allocation
  double *allocate(long const, ostream &output=cout);

  //check memory usage
  inline void check(ostream &output=cout, string *s=NULL) {
    if(s!=NULL) output<<s<<endl;
    cout<<iend_<<" used; "<<len_-iend_<<" available."<<endl;
  };

  inline bool intSij(int op, int i, int j) {
    ;
  };

  //rewind pointer to the beginning
  inline void release() {iend_ = 0;};

  //print elements
  void print(long const, long const, long col=5, ostream &output=cout, char *s=NULL);
};
} // namespace ChronusQ

#endif
