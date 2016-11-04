#include <grid2/lebedev.h>
namespace ChronusQ {
  template<>
  void Lebedev::loadLebedev<LEBEDEV_21>(){
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
    this->gPoints_.push_back(cartGP(0.255125262111413,0.255125262111413,0.932642590312691));     
    this->gPoints_.push_back(cartGP(0.255125262111413,0.255125262111413,-0.932642590312691));     
    this->gPoints_.push_back(cartGP(0.255125262111413,-0.255125262111413,0.932642590312691));     
    this->gPoints_.push_back(cartGP(0.255125262111413,-0.255125262111413,-0.932642590312691));     
    this->gPoints_.push_back(cartGP(-0.255125262111413,0.255125262111413,0.932642590312691));     
    this->gPoints_.push_back(cartGP(-0.255125262111413,0.255125262111413,-0.932642590312691));     
    this->gPoints_.push_back(cartGP(-0.255125262111413,-0.255125262111413,0.932642590312691));     
    this->gPoints_.push_back(cartGP(-0.255125262111413,-0.255125262111413,-0.932642590312691));     
    this->gPoints_.push_back(cartGP(0.255125262111413,0.932642590312691,0.255125262111413));     
    this->gPoints_.push_back(cartGP(0.255125262111413,-0.932642590312691,0.255125262111413));     
    this->gPoints_.push_back(cartGP(0.255125262111413,0.932642590312691,-0.255125262111414));     
    this->gPoints_.push_back(cartGP(0.255125262111413,-0.932642590312691,-0.255125262111414));     
    this->gPoints_.push_back(cartGP(-0.255125262111414,0.932642590312691,0.255125262111413));     
    this->gPoints_.push_back(cartGP(-0.255125262111414,-0.932642590312691,0.255125262111413));     
    this->gPoints_.push_back(cartGP(-0.255125262111414,0.932642590312691,-0.255125262111414));     
    this->gPoints_.push_back(cartGP(-0.255125262111414,-0.932642590312691,-0.255125262111414));     
    this->gPoints_.push_back(cartGP(0.932642590312691,0.255125262111413,0.255125262111413));     
    this->gPoints_.push_back(cartGP(-0.932642590312691,0.255125262111413,0.255125262111413));     
    this->gPoints_.push_back(cartGP(0.932642590312691,0.255125262111413,-0.255125262111414));     
    this->gPoints_.push_back(cartGP(-0.932642590312691,0.255125262111413,-0.255125262111414));     
    this->gPoints_.push_back(cartGP(0.932642590312691,-0.255125262111413,0.255125262111413));     
    this->gPoints_.push_back(cartGP(-0.932642590312691,-0.255125262111413,0.255125262111413));     
    this->gPoints_.push_back(cartGP(0.932642590312691,-0.255125262111413,-0.255125262111414));     
    this->gPoints_.push_back(cartGP(-0.932642590312691,-0.255125262111413,-0.255125262111414));     
    this->gPoints_.push_back(cartGP(0.674360146036277,0.674360146036277,0.300793595137701));     
    this->gPoints_.push_back(cartGP(0.674360146036277,0.674360146036277,-0.300793595137702));     
    this->gPoints_.push_back(cartGP(0.674360146036277,-0.674360146036277,0.300793595137701));     
    this->gPoints_.push_back(cartGP(0.674360146036277,-0.674360146036277,-0.300793595137702));     
    this->gPoints_.push_back(cartGP(-0.674360146036277,0.674360146036277,0.300793595137701));     
    this->gPoints_.push_back(cartGP(-0.674360146036277,0.674360146036277,-0.300793595137702));     
    this->gPoints_.push_back(cartGP(-0.674360146036277,-0.674360146036277,0.300793595137701));     
    this->gPoints_.push_back(cartGP(-0.674360146036277,-0.674360146036277,-0.300793595137702));     
    this->gPoints_.push_back(cartGP(0.674360146036277,0.300793595137702,0.674360146036276));     
    this->gPoints_.push_back(cartGP(0.674360146036277,-0.300793595137702,0.674360146036276));     
    this->gPoints_.push_back(cartGP(0.674360146036276,0.300793595137701,-0.674360146036277));     
    this->gPoints_.push_back(cartGP(0.674360146036276,-0.300793595137701,-0.674360146036277));     
    this->gPoints_.push_back(cartGP(-0.674360146036277,0.300793595137702,0.674360146036276));     
    this->gPoints_.push_back(cartGP(-0.674360146036277,-0.300793595137702,0.674360146036276));     
    this->gPoints_.push_back(cartGP(-0.674360146036276,0.300793595137702,-0.674360146036277));     
    this->gPoints_.push_back(cartGP(-0.674360146036276,-0.300793595137702,-0.674360146036277));     
    this->gPoints_.push_back(cartGP(0.300793595137701,0.674360146036277,0.674360146036276));     
    this->gPoints_.push_back(cartGP(-0.300793595137701,0.674360146036277,0.674360146036276));     
    this->gPoints_.push_back(cartGP(0.300793595137701,0.674360146036276,-0.674360146036277));     
    this->gPoints_.push_back(cartGP(-0.300793595137701,0.674360146036276,-0.674360146036277));     
    this->gPoints_.push_back(cartGP(0.300793595137701,-0.674360146036277,0.674360146036276));     
    this->gPoints_.push_back(cartGP(-0.300793595137701,-0.674360146036277,0.674360146036276));     
    this->gPoints_.push_back(cartGP(0.300793595137701,-0.674360146036276,-0.674360146036277));     
    this->gPoints_.push_back(cartGP(-0.300793595137701,-0.674360146036276,-0.674360146036277));     
    this->gPoints_.push_back(cartGP(0.431891069671941,0.431891069671941,0.791795559393492));     
    this->gPoints_.push_back(cartGP(0.431891069671941,0.431891069671941,-0.791795559393492));     
    this->gPoints_.push_back(cartGP(0.431891069671941,-0.431891069671941,0.791795559393492));     
    this->gPoints_.push_back(cartGP(0.431891069671941,-0.431891069671941,-0.791795559393492));     
    this->gPoints_.push_back(cartGP(-0.431891069671941,0.431891069671941,0.791795559393492));     
    this->gPoints_.push_back(cartGP(-0.431891069671941,0.431891069671941,-0.791795559393492));     
    this->gPoints_.push_back(cartGP(-0.431891069671941,-0.431891069671941,0.791795559393492));     
    this->gPoints_.push_back(cartGP(-0.431891069671941,-0.431891069671941,-0.791795559393492));     
    this->gPoints_.push_back(cartGP(0.431891069671941,0.791795559393492,0.431891069671941));     
    this->gPoints_.push_back(cartGP(0.431891069671941,-0.791795559393492,0.431891069671941));     
    this->gPoints_.push_back(cartGP(0.431891069671941,0.791795559393492,-0.431891069671941));     
    this->gPoints_.push_back(cartGP(0.431891069671941,-0.791795559393492,-0.431891069671941));     
    this->gPoints_.push_back(cartGP(-0.431891069671941,0.791795559393492,0.431891069671941));     
    this->gPoints_.push_back(cartGP(-0.431891069671941,-0.791795559393492,0.431891069671941));     
    this->gPoints_.push_back(cartGP(-0.431891069671941,0.791795559393492,-0.431891069671941));     
    this->gPoints_.push_back(cartGP(-0.431891069671941,-0.791795559393492,-0.431891069671941));     
    this->gPoints_.push_back(cartGP(0.791795559393492,0.431891069671941,0.431891069671941));     
    this->gPoints_.push_back(cartGP(-0.791795559393492,0.431891069671941,0.431891069671941));     
    this->gPoints_.push_back(cartGP(0.791795559393492,0.431891069671941,-0.431891069671941));     
    this->gPoints_.push_back(cartGP(-0.791795559393492,0.431891069671941,-0.431891069671941));     
    this->gPoints_.push_back(cartGP(0.791795559393492,-0.431891069671941,0.431891069671941));     
    this->gPoints_.push_back(cartGP(-0.791795559393492,-0.431891069671941,0.431891069671941));     
    this->gPoints_.push_back(cartGP(0.791795559393492,-0.431891069671941,-0.431891069671941));     
    this->gPoints_.push_back(cartGP(-0.791795559393492,-0.431891069671941,-0.431891069671941));     
    this->gPoints_.push_back(cartGP(0.261393136033599,0.965232421976448,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.261393136033599,-0.965232421976448,0.000000000000000));     
    this->gPoints_.push_back(cartGP(-0.261393136033599,0.965232421976448,0.000000000000000));     
    this->gPoints_.push_back(cartGP(-0.261393136033599,-0.965232421976448,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.965232421976448,0.261393136033599,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.965232421976448,-0.261393136033599,0.000000000000000));     
    this->gPoints_.push_back(cartGP(-0.965232421976448,0.261393136033599,0.000000000000000));     
    this->gPoints_.push_back(cartGP(-0.965232421976448,-0.261393136033599,0.000000000000000));     
    this->gPoints_.push_back(cartGP(0.261393136033599,0.000000000000000,0.965232421976448));     
    this->gPoints_.push_back(cartGP(0.261393136033599,0.000000000000000,-0.965232421976448));     
    this->gPoints_.push_back(cartGP(-0.261393136033599,0.000000000000000,0.965232421976448));     
    this->gPoints_.push_back(cartGP(-0.261393136033599,0.000000000000000,-0.965232421976448));     
    this->gPoints_.push_back(cartGP(0.965232421976448,0.000000000000000,0.261393136033599));     
    this->gPoints_.push_back(cartGP(0.965232421976448,0.000000000000000,-0.261393136033599));     
    this->gPoints_.push_back(cartGP(-0.965232421976448,0.000000000000000,0.261393136033599));     
    this->gPoints_.push_back(cartGP(-0.965232421976448,0.000000000000000,-0.261393136033599));     
    this->gPoints_.push_back(cartGP(0.000000000000000,0.261393136033599,0.965232421976448));     
    this->gPoints_.push_back(cartGP(0.000000000000000,0.261393136033599,-0.965232421976448));     
    this->gPoints_.push_back(cartGP(0.000000000000000,-0.261393136033599,0.965232421976448));     
    this->gPoints_.push_back(cartGP(0.000000000000000,-0.261393136033599,-0.965232421976448));     
    this->gPoints_.push_back(cartGP(0.000000000000000,0.965232421976448,0.261393136033599));     
    this->gPoints_.push_back(cartGP(0.000000000000000,0.965232421976448,-0.261393136033599));     
    this->gPoints_.push_back(cartGP(0.000000000000000,-0.965232421976448,0.261393136033599));     
    this->gPoints_.push_back(cartGP(0.000000000000000,-0.965232421976448,-0.261393136033599));     
    this->gPoints_.push_back(cartGP(0.499045316179604,0.144663074432511,0.854415804684659));     
    this->gPoints_.push_back(cartGP(0.499045316179603,0.144663074432511,-0.854415804684659));     
    this->gPoints_.push_back(cartGP(0.499045316179604,-0.144663074432511,0.854415804684659));     
    this->gPoints_.push_back(cartGP(0.499045316179603,-0.144663074432511,-0.854415804684659));     
    this->gPoints_.push_back(cartGP(-0.499045316179604,0.144663074432512,0.854415804684659));     
    this->gPoints_.push_back(cartGP(-0.499045316179603,0.144663074432511,-0.854415804684659));     
    this->gPoints_.push_back(cartGP(-0.499045316179604,-0.144663074432512,0.854415804684659));     
    this->gPoints_.push_back(cartGP(-0.499045316179603,-0.144663074432511,-0.854415804684659));     
    this->gPoints_.push_back(cartGP(0.499045316179604,0.854415804684659,0.144663074432512));     
    this->gPoints_.push_back(cartGP(0.499045316179604,0.854415804684659,-0.144663074432512));     
    this->gPoints_.push_back(cartGP(0.499045316179604,-0.854415804684659,0.144663074432512));     
    this->gPoints_.push_back(cartGP(0.499045316179604,-0.854415804684659,-0.144663074432512));     
    this->gPoints_.push_back(cartGP(-0.499045316179604,0.854415804684659,0.144663074432512));     
    this->gPoints_.push_back(cartGP(-0.499045316179604,0.854415804684659,-0.144663074432512));     
    this->gPoints_.push_back(cartGP(-0.499045316179604,-0.854415804684659,0.144663074432512));     
    this->gPoints_.push_back(cartGP(-0.499045316179604,-0.854415804684659,-0.144663074432512));     
    this->gPoints_.push_back(cartGP(0.144663074432511,0.499045316179604,0.854415804684659));     
    this->gPoints_.push_back(cartGP(0.144663074432511,0.499045316179603,-0.854415804684659));     
    this->gPoints_.push_back(cartGP(0.144663074432511,-0.499045316179604,0.854415804684659));     
    this->gPoints_.push_back(cartGP(0.144663074432511,-0.499045316179603,-0.854415804684659));     
    this->gPoints_.push_back(cartGP(-0.144663074432512,0.499045316179604,0.854415804684659));     
    this->gPoints_.push_back(cartGP(-0.144663074432511,0.499045316179603,-0.854415804684659));     
    this->gPoints_.push_back(cartGP(-0.144663074432512,-0.499045316179604,0.854415804684659));     
    this->gPoints_.push_back(cartGP(-0.144663074432511,-0.499045316179603,-0.854415804684659));     
    this->gPoints_.push_back(cartGP(0.144663074432511,0.854415804684659,0.499045316179604));     
    this->gPoints_.push_back(cartGP(0.144663074432511,0.854415804684659,-0.499045316179604));     
    this->gPoints_.push_back(cartGP(0.144663074432511,-0.854415804684659,0.499045316179604));     
    this->gPoints_.push_back(cartGP(0.144663074432511,-0.854415804684659,-0.499045316179604));     
    this->gPoints_.push_back(cartGP(-0.144663074432511,0.854415804684659,0.499045316179604));     
    this->gPoints_.push_back(cartGP(-0.144663074432511,0.854415804684659,-0.499045316179604));     
    this->gPoints_.push_back(cartGP(-0.144663074432511,-0.854415804684659,0.499045316179604));     
    this->gPoints_.push_back(cartGP(-0.144663074432511,-0.854415804684659,-0.499045316179604));     
    this->gPoints_.push_back(cartGP(0.854415804684659,0.499045316179604,0.144663074432512));     
    this->gPoints_.push_back(cartGP(0.854415804684659,0.499045316179604,-0.144663074432512));     
    this->gPoints_.push_back(cartGP(0.854415804684659,-0.499045316179604,0.144663074432512));     
    this->gPoints_.push_back(cartGP(0.854415804684659,-0.499045316179604,-0.144663074432512));     
    this->gPoints_.push_back(cartGP(-0.854415804684659,0.499045316179603,0.144663074432512));     
    this->gPoints_.push_back(cartGP(-0.854415804684659,0.499045316179603,-0.144663074432512));     
    this->gPoints_.push_back(cartGP(-0.854415804684659,-0.499045316179603,0.144663074432512));     
    this->gPoints_.push_back(cartGP(-0.854415804684659,-0.499045316179603,-0.144663074432512));     
    this->gPoints_.push_back(cartGP(0.854415804684659,0.144663074432511,0.499045316179604));     
    this->gPoints_.push_back(cartGP(0.854415804684659,0.144663074432511,-0.499045316179604));     
    this->gPoints_.push_back(cartGP(0.854415804684659,-0.144663074432511,0.499045316179604));     
    this->gPoints_.push_back(cartGP(0.854415804684659,-0.144663074432511,-0.499045316179604));     
    this->gPoints_.push_back(cartGP(-0.854415804684659,0.144663074432511,0.499045316179604));     
    this->gPoints_.push_back(cartGP(-0.854415804684659,0.144663074432511,-0.499045316179604));     
    this->gPoints_.push_back(cartGP(-0.854415804684659,-0.144663074432511,0.499045316179604));     
    this->gPoints_.push_back(cartGP(-0.854415804684659,-0.144663074432511,-0.499045316179604));     
    this->weights_.push_back(0.005544842902037);
    this->weights_.push_back(0.005544842902037);
    this->weights_.push_back(0.005544842902037);
    this->weights_.push_back(0.005544842902037);
    this->weights_.push_back(0.005544842902037);
    this->weights_.push_back(0.005544842902037);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006071332770671);
    this->weights_.push_back(0.006383674773515);
    this->weights_.push_back(0.006383674773515);
    this->weights_.push_back(0.006383674773515);
    this->weights_.push_back(0.006383674773515);
    this->weights_.push_back(0.006383674773515);
    this->weights_.push_back(0.006383674773515);
    this->weights_.push_back(0.006383674773515);
    this->weights_.push_back(0.006383674773515);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.005183387587748);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006317929009814);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.006201670006589);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005477143385137);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
    this->weights_.push_back(0.005968383987681);
  };
};
