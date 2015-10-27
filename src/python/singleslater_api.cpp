#include <workers.h>
#include <global.h>

using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::AOIntegrals;

namespace ChronusQ {
  template<>
  void SingleSlater<double>::Wrapper_iniSingleSlater(Molecule &mol, BasisSet &basis,
         AOIntegrals &ints, FileIO &fileio, Controls &controls) {
     this->iniSingleSlater( &mol, &basis, &ints, &fileio, &controls); 
  }
  template<>
  void SingleSlater<dcomplex>::Wrapper_iniSingleSlater(Molecule &mol, BasisSet &basis,
         AOIntegrals &ints, FileIO &fileio, Controls &controls) {
     this->iniSingleSlater( &mol, &basis, &ints, &fileio, &controls); 
  }

  template<>
  boost::python::list SingleSlater<double>::Wrapper_dipole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++)
      result.append((*this->dipole_)(i));
    return result;
  }

  template<>
  boost::python::list SingleSlater<double>::Wrapper_quadrupole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++) {
      boost::python::list oneDim;
      for(auto j = 0; j < 3; j++) 
        oneDim.append((*this->quadpole_)(i,j));
      result.append(oneDim);
    }
    return result;
  }

  template<>
  boost::python::list SingleSlater<double>::Wrapper_octupole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++) {
      boost::python::list twoDim;
      for(auto j = 0; j < 3; j++) {
        boost::python::list oneDim;
        for(auto k = 0; k < 3; k++)
          oneDim.append((*this->octpole_)(i,j,k));
        twoDim.append(oneDim);
      }
      result.append(twoDim);
    }
    return result;
  }

  template<>
  boost::python::list SingleSlater<dcomplex>::Wrapper_dipole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++)
      result.append((*this->dipole_)(i));
    return result;
  }

  template<>
  boost::python::list SingleSlater<dcomplex>::Wrapper_quadrupole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++) {
      boost::python::list oneDim;
      for(auto j = 0; j < 3; j++) 
        oneDim.append((*this->quadpole_)(i,j));
      result.append(oneDim);
    }
    return result;
  }

  template<>
  boost::python::list SingleSlater<dcomplex>::Wrapper_octupole(){
    boost::python::list result;
    for(auto i = 0; i < 3; i++) {
      boost::python::list twoDim;
      for(auto j = 0; j < 3; j++) {
        boost::python::list oneDim;
        for(auto k = 0; k < 3; k++)
          oneDim.append((*this->octpole_)(i,j,k));
        twoDim.append(oneDim);
      }
      result.append(twoDim);
    }
    return result;
  }
};


