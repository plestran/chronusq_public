#include <workers.h>
#include <global.h>

using ChronusQ::RealTime;

namespace ChronusQ {

/*
  template<>
  boost::python::list RealTime<double>::Wrapper_recs(){
    boost::python::list result;
    for(auto i = this->propInfo.begin(); i != this->propInfo.end(); i++){
      Wrapper_PropInfo tmpRec(*i);
      result.append(tmpRec);
    }
    return result;
  }

  template<>
  boost::python::list RealTime<dcomplex>::Wrapper_recs(){
    boost::python::list result;
    for(auto i = this->propInfo.begin(); i != this->propInfo.end(); i++){
      Wrapper_PropInfo tmpRec(*i);
      result.append(tmpRec);
    }
    return result;
  }
*/

  template<>
  boost::python::list RealTime<double>::lastDipole(){
    boost::python::list result;
    for(auto i = 0; i < 4; i++)
      result.append(this->propInfo[this->propInfo.size()-1].dipole[i]);
    return result;
  }

  template<>
  double RealTime<double>::lastEnergy(){
    return this->propInfo[this->propInfo.size()-1].energy;
  }

  template<>
  double RealTime<double>::getTimeStep(){
    return this->propInfo[this->propInfo.size()-1].timeStep;
  }
  template<>
  boost::python::list RealTime<dcomplex>::lastDipole(){
    boost::python::list result;
    for(auto i = 0; i < 4; i++)
      result.append(this->propInfo[this->propInfo.size()-1].dipole[i]);
    return result;
  }

  template<>
  double RealTime<dcomplex>::lastEnergy(){
    return this->propInfo[this->propInfo.size()-1].energy;
  }

  template<>
  double RealTime<dcomplex>::getTimeStep(){
    return this->propInfo[this->propInfo.size()-1].timeStep;
  }
  



}; // namespace ChronusQ

