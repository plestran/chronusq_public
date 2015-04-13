#include "basisset.h"

#ifdef USE_LIBINT

using ChronusQ::BasisSet;
typedef ChronusQ::Shell CShell;
typedef libint2::Shell LShell;

namespace ChronusQ {
std::vector<LShell> convShell() {
  std::vector<LShell> shells;
  shells.push_back(
     {
       {3.425250910, 0.623913730, 0.168855400},
       {
         {0, false, {0.15432897, 0.53532814, 0.44463454}},
         {1, false, {0.15432897, 0.53532814, 0.44463454}}
       },
     {{0,0,0}}
     }
  );

  cout << shells[0] << endl;
  return shells;

}
} // namespace ChronusQ
#endif
