#include <grid2.h>
namespace ChronusQ {
  template<>
  void Lebedev::loadLebedev<LEBEDEV_3>(){
    this->gPoints_.push_back(
        cartGP( 1.000000000000000, 0.000000000000000, 0.000000000000000)
        );     
    this->gPoints_.push_back(
        cartGP(-1.000000000000000, 0.000000000000000, 0.000000000000000)
        );     
    this->gPoints_.push_back(
        cartGP( 0.000000000000000, 1.000000000000000, 0.000000000000000)
        );     
    this->gPoints_.push_back(
        cartGP( 0.000000000000000,-1.000000000000000, 0.000000000000000)
        );     
    this->gPoints_.push_back(
        cartGP( 0.000000000000000, 0.000000000000000, 1.000000000000000)
        );     
    this->gPoints_.push_back(
        cartGP( 0.000000000000000, 0.000000000000000,-1.000000000000000)
        );     
    this->weights_.push_back(0.166666666666667);
    this->weights_.push_back(0.166666666666667);
    this->weights_.push_back(0.166666666666667);
    this->weights_.push_back(0.166666666666667);
    this->weights_.push_back(0.166666666666667);
    this->weights_.push_back(0.166666666666667);
  };
};
