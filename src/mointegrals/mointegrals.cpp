#include "mointegrals.h"
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
//---------------------
// initialize MOIntegrals
//---------------------
void MOIntegrals::iniMOIntegrals(   Molecule *molecule, BasisSet *basisSet, FileIO *fileio, 
                                    Controls *controls, AOIntegrals *aointegrals, SingleSlater *singleSlater){
  this->molecule_ = molecule;
  this->basisSet_ = basisSet;
  this->fileio_   = fileio;
  this->controls_ = controls;
  this->aointegrals_ = aointegrals;
  this->singleSlater_ = singleSlater;

  this->haveMOiajb = false;
  this->haveMOijab = false;
  this->haveMOijka = false;
  this->haveMOijkl = false;
  this->haveMOiabc = false;
  this->haveMOabcd = false;
};

