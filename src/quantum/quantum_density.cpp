#include <quantum.h>

namespace ChronusQ {
  template<>
  void Quantum<double>::complexMyScale(){ };

  template<>
  void Quantum<dcomplex>::complexMyScale(){
    (*this->onePDMMy_) *= dcomplex(0.0,1.0);
  };
};
