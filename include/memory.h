#ifndef INCLUDED_MEMORY
#define INCLUDED_MEMORY
#include "global.h"

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

#endif
