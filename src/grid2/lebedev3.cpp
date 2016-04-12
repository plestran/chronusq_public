#include <grid2.h>
namespace ChronusQ {
  template<>
  void Lebedev::loadLebedev<LEBEDEV_3>(){
    this->sphPoints_.push_back(
        cartGP( 1.000000000000000, 0.000000000000000, 0.000000000000000)
        );     
    this->sphPoints_.push_back(
        cartGP(-1.000000000000000, 0.000000000000000, 0.000000000000000)
        );     
    this->sphPoints_.push_back(
        cartGP( 0.000000000000000, 1.000000000000000, 0.000000000000000)
        );     
    this->sphPoints_.push_back(
        cartGP( 0.000000000000000,-1.000000000000000, 0.000000000000000)
        );     
    this->sphPoints_.push_back(
        cartGP( 0.000000000000000, 0.000000000000000, 1.000000000000000)
        );     
    this->sphPoints_.push_back(
        cartGP( 0.000000000000000, 0.000000000000000,-1.000000000000000)
        );     
    this->sphWeights_.push_back(0.166666666666667);
    this->sphWeights_.push_back(0.166666666666667);
    this->sphWeights_.push_back(0.166666666666667);
    this->sphWeights_.push_back(0.166666666666667);
    this->sphWeights_.push_back(0.166666666666667);
    this->sphWeights_.push_back(0.166666666666667);
  };
};
