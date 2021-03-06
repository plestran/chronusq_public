#include <grid2/lebedev.h>
namespace ChronusQ {
  template<>
  void Lebedev::loadLebedev<LEBEDEV_5>(){
    this->gPoints_.push_back(
        cartGP(1.000000000000000,0.000000000000000,0.000000000000000)
        );
    this->weights_.push_back(0.066666666666667);
    this->gPoints_.push_back(
        cartGP(-1.000000000000000,0.000000000000000,0.000000000000000)
        );     
    this->weights_.push_back(0.066666666666667);
    this->gPoints_.push_back(
        cartGP(0.000000000000000,1.000000000000000,0.000000000000000)
        );     
    this->weights_.push_back(0.066666666666667);
    this->gPoints_.push_back(
        cartGP(0.000000000000000,-1.000000000000000,0.000000000000000)
        );     
    this->weights_.push_back(0.066666666666667);
    this->gPoints_.push_back(
        cartGP(0.000000000000000,0.000000000000000,1.000000000000000)
        );     
    this->weights_.push_back(0.066666666666667);
    this->gPoints_.push_back(
        cartGP(0.000000000000000,0.000000000000000,-1.000000000000000)
        );     
    this->weights_.push_back(0.066666666666667);
    this->gPoints_.push_back(
        cartGP(0.577350269189626,0.577350269189626,0.577350269189626)
        );     
    this->weights_.push_back(0.075000000000000);
    this->gPoints_.push_back(
        cartGP(0.577350269189626,0.577350269189626,-0.577350269189626)
        );     
    this->weights_.push_back(0.075000000000000);
    this->gPoints_.push_back(
        cartGP(0.577350269189626,-0.577350269189626,0.577350269189626)
        );     
    this->weights_.push_back(0.075000000000000);
    this->gPoints_.push_back(
        cartGP(0.577350269189626,-0.577350269189626,-0.577350269189626)
        );     
    this->weights_.push_back(0.075000000000000);
    this->gPoints_.push_back(
        cartGP(-0.577350269189626,0.577350269189626,0.577350269189626)
        );     
    this->weights_.push_back(0.075000000000000);
    this->gPoints_.push_back(
        cartGP(-0.577350269189626,0.577350269189626,-0.577350269189626)
        );     
    this->weights_.push_back(0.075000000000000);
    this->gPoints_.push_back(
        cartGP(-0.577350269189626,-0.577350269189626,0.577350269189626)
        );     
    this->weights_.push_back(0.075000000000000);
    this->gPoints_.push_back(
        cartGP(-0.577350269189626,-0.577350269189626,-0.577350269189626)
        );     
    this->weights_.push_back(0.075000000000000);
  };
};
