#include <workers.h>
#include <global.h>

using ChronusQ::SDResponse;

namespace ChronusQ {
  template<>
  boost::python::list SDResponse<double>::Wrapper_excitationEnergies() {
    boost::python::list result;
    for(auto i = 0; i < this->nSek_; i++)
      result.append((*this->omega_)(i));
    return result;
  }

  template<>
  boost::python::list SDResponse<double>::Wrapper_oscStrengths() {
    boost::python::list result;
    for(auto i = 0; i < this->nSek_; i++)
      result.append((*this->oscStrength_)(0,i+1));
    return result;
  }

  template<>
  boost::python::list SDResponse<dcomplex>::Wrapper_excitationEnergies() {
    boost::python::list result;
    for(auto i = 0; i < this->nSek_; i++)
      result.append((*this->omega_)(i));
    return result;
  }

  template<>
  boost::python::list SDResponse<dcomplex>::Wrapper_oscStrengths() {
    boost::python::list result;
    for(auto i = 0; i < this->nSek_; i++)
      result.append((*this->oscStrength_)(0,i+1));
    return result;
  }
}; //namespace ChronusQ
