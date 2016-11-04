#include <grid2/lebedev.h>
namespace ChronusQ {
  template<>
  void Lebedev::loadLebedev<LEBEDEV_11>(){
    this->gPoints_.push_back(cartGP(1.000000000000000,0.000000000000000,0.000000000000000));     
    this->gPoints_.push_back(cartGP(-1.000000000000000,0.000000000000000,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.000000000000000,1.000000000000000,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.000000000000000,-1.000000000000000,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.000000000000000,0.000000000000000,1.000000000000000));     
    this->gPoints_.push_back(cartGP(0.000000000000000,0.000000000000000,-1.000000000000000));     
    this->gPoints_.push_back(cartGP(0.000000000000000,0.707106781186547,0.707106781186548));     
    this->gPoints_.push_back(cartGP(0.000000000000000,0.707106781186548,-0.707106781186547));     
    this->gPoints_.push_back(cartGP(0.000000000000000,-0.707106781186547,0.707106781186548));     
    this->gPoints_.push_back(cartGP(0.000000000000000,-0.707106781186548,-0.707106781186547));     
    this->gPoints_.push_back(cartGP(0.707106781186547,0.000000000000000,0.707106781186548));     
    this->gPoints_.push_back(cartGP(0.707106781186548,0.000000000000000,-0.707106781186547));     
    this->gPoints_.push_back(cartGP(-0.707106781186547,0.000000000000000,0.707106781186548));     
    this->gPoints_.push_back(cartGP(-0.707106781186548,0.000000000000000,-0.707106781186547));     
    this->gPoints_.push_back(cartGP(0.707106781186548,0.707106781186547,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.707106781186548,-0.707106781186547,0.000000000000000));     
    this->gPoints_.push_back(cartGP(-0.707106781186547,0.707106781186548,0.000000000000000));     
    this->gPoints_.push_back(cartGP(-0.707106781186547,-0.707106781186548,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.577350269189626,0.577350269189626,0.577350269189626));     
    this->gPoints_.push_back(cartGP(0.577350269189626,0.577350269189626,-0.577350269189626));     
    this->gPoints_.push_back(cartGP(0.577350269189626,-0.577350269189626,0.577350269189626));     
    this->gPoints_.push_back(cartGP(0.577350269189626,-0.577350269189626,-0.577350269189626));     
    this->gPoints_.push_back(cartGP(-0.577350269189626,0.577350269189626,0.577350269189626));     
    this->gPoints_.push_back(cartGP(-0.577350269189626,0.577350269189626,-0.577350269189626));     
    this->gPoints_.push_back(cartGP(-0.577350269189626,-0.577350269189626,0.577350269189626));     
    this->gPoints_.push_back(cartGP(-0.577350269189626,-0.577350269189626,-0.577350269189626));     
    this->gPoints_.push_back(cartGP(0.301511344577764,0.301511344577764,0.904534033733291));     
    this->gPoints_.push_back(cartGP(0.301511344577764,0.301511344577764,-0.904534033733291));     
    this->gPoints_.push_back(cartGP(0.301511344577764,-0.301511344577764,0.904534033733291));     
    this->gPoints_.push_back(cartGP(0.301511344577764,-0.301511344577764,-0.904534033733291));     
    this->gPoints_.push_back(cartGP(-0.301511344577764,0.301511344577764,0.904534033733291));     
    this->gPoints_.push_back(cartGP(-0.301511344577764,0.301511344577764,-0.904534033733291));     
    this->gPoints_.push_back(cartGP(-0.301511344577764,-0.301511344577764,0.904534033733291));     
    this->gPoints_.push_back(cartGP(-0.301511344577764,-0.301511344577764,-0.904534033733291));     
    this->gPoints_.push_back(cartGP(0.301511344577764,0.904534033733291,0.301511344577764));     
    this->gPoints_.push_back(cartGP(0.301511344577764,-0.904534033733291,0.301511344577764));     
    this->gPoints_.push_back(cartGP(0.301511344577764,0.904534033733291,-0.301511344577764));     
    this->gPoints_.push_back(cartGP(0.301511344577764,-0.904534033733291,-0.301511344577764));     
    this->gPoints_.push_back(cartGP(-0.301511344577764,0.904534033733291,0.301511344577764));     
    this->gPoints_.push_back(cartGP(-0.301511344577764,-0.904534033733291,0.301511344577764));     
    this->gPoints_.push_back(cartGP(-0.301511344577764,0.904534033733291,-0.301511344577764));     
    this->gPoints_.push_back(cartGP(-0.301511344577764,-0.904534033733291,-0.301511344577764));     
    this->gPoints_.push_back(cartGP(0.904534033733291,0.301511344577764,0.301511344577764));     
    this->gPoints_.push_back(cartGP(-0.904534033733291,0.301511344577763,0.301511344577764));     
    this->gPoints_.push_back(cartGP(0.904534033733291,0.301511344577764,-0.301511344577764));     
    this->gPoints_.push_back(cartGP(-0.904534033733291,0.301511344577763,-0.301511344577764));     
    this->gPoints_.push_back(cartGP(0.904534033733291,-0.301511344577764,0.301511344577764));     
    this->gPoints_.push_back(cartGP(-0.904534033733291,-0.301511344577763,0.301511344577764));     
    this->gPoints_.push_back(cartGP(0.904534033733291,-0.301511344577764,-0.301511344577764));     
    this->gPoints_.push_back(cartGP(-0.904534033733291,-0.301511344577763,-0.301511344577764));     
    this->weights_.push_back(0.012698412698413);
    this->weights_.push_back(0.012698412698413);
    this->weights_.push_back(0.012698412698413);
    this->weights_.push_back(0.012698412698413);
    this->weights_.push_back(0.012698412698413);
    this->weights_.push_back(0.012698412698413);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.022574955908289);
    this->weights_.push_back(0.021093750000000);
    this->weights_.push_back(0.021093750000000);
    this->weights_.push_back(0.021093750000000);
    this->weights_.push_back(0.021093750000000);
    this->weights_.push_back(0.021093750000000);
    this->weights_.push_back(0.021093750000000);
    this->weights_.push_back(0.021093750000000);
    this->weights_.push_back(0.021093750000000);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
    this->weights_.push_back(0.020173335537919);
  };
};