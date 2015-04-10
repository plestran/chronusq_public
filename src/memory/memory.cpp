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
#include "memory.h"
using ChronusQ::Memory;
//constructor
Memory::Memory(long n, char *unit, ostream output) {
  double test;
  int   i;
  unsigned long long const mw=8883608,gw=8589934592;
  unsigned long long const mb=1048576,gb=1073741824;
  i = sizeof(test);
  strlwr(unit);
  if(!strcmp(unit,"mw")) len_ = n*mw/i;
  if(!strcmp(unit,"mb")) len_ = n*mb/i;
  if(!strcmp(unit,"gw")) len_ = n*gb/i;
  if(!strcmp(unit,"gb")) len_ = n*gw/i;
  output<<endl<<n<<" "<<unit<<" of memory requested."<<endl;
  output<<endl<<len_*i<<" bytes for "<<len_<<" double precision numbers."<<endl;
  memory_ = new (nothrow) double[len_];
  if(memory_==NULL) throw 2001;
  iend_ = 0;
};

//print elements
void Memory::print(long const start, long const final, long col, ostream &output, char *s) {
  long i,j,end;
  if(s!=NULL) output<<s<<endl;
  for(i=start;i<=final;i+=col) {
    end = col;
    if((i+col)>final) end = final - i;
    for(j=0;j<=end;j++) output<<this->memory_[i+j]<<'/t';
    output<<'/n';
  };
};

