#include "controls.h"

/****************************/
/* Error Messages 6000-6999 */
/****************************/
//constructor and initialization
void Controls::iniControls(){
  this->printLevel = 0;
  this->energyOnly = false;
  this->optWaveFunction = true;
  this->optGeometry = false;
  this->firstDer = false;
  this->secondDer = false;
  this->HF = true;
  this->DFT = false;
  this->hybridDFT = false;
  this->restart = false;
  this->thresholdS = 1.0e-10;
  this->thresholdAB = 1.0e-6;
  this->thresholdSchawrtz = 1.0e-14;
  this->guess = 0;
  this->directTwoE = false;
};
