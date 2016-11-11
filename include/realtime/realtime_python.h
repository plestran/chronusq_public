template <typename T>
boost::python::list RealTime<T>::Python_propEnergy() {
  boost::python::list result;
  for(auto i : this->propInfo){
    result.append(i.energy);
  }
  return result;
};

template <typename T>
boost::python::list RealTime<T>::Python_propDipole() {
  boost::python::list result;
  for(auto i : this->propInfo){
    boost::python::list dipole;
    for(auto x : i.dipole)
      dipole.append(x);
    result.append(dipole);
  }
  return result;
};
