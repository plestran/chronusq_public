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

