#include <boost/python.hpp>
#include <singleslater.h>
using namespace boost::python;
using ChronusQ::SingleSlater;

BOOST_PYTHON_MODULE(SingleSlater){
  class_<SingleSlater<double>,boost::noncopyable>("SingleSlater_double",init<>())
    .def("iniSingleSlater", &SingleSlater<double>::iniSingleSlater)
    .def("printInfo"      , &SingleSlater<double>::printInfo      )
  ;
};

