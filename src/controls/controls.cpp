/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#include <controls.h>
using ChronusQ::Controls;
/****************************/
/* Error Messages 6000-6999 */
/****************************/
//constructor and initialization
void Controls::iniControls(){
  this->printLevel =        0;
  this->energyOnly =        false;
  this->optWaveFunction =   true;
  this->optGeometry =       false;
  this->firstDer =          false;
  this->secondDer =         false;
  this->HF =                true;
  this->DFT =               false;
  this->hybridDFT =         false;
  this->doTCS =             false;
  this->restart =           false;
  this->thresholdS =        1.0e-10;
  this->thresholdAB =       1.0e-6;
  this->thresholdSchawrtz = 1.0e-14;
  this->guess =             0;
  this->directTwoE =        true;
  this->buildn4eri =        false;
  this->doDF =              false;
  this->doSDR     =         false;
  this->doCUHF    =         false;
  this->doComplex =         false;
  this->doUnit    =         false;
  this->SDMethod  =         0;
  this->doRealTime =        false;
  this->rtMaxSteps =        10;
  this->rtTimeStep =        0.05;
  this->rtTypeOrtho =       0;
  this->rtInitDensity =     0;
  this->rtSwapMOA     =     0; 
  this->rtSwapMOB     =     0; 
  this->rtMethFormU   =     0;
  this->rtField_      =     {0.0,0.0,0.0};
  this->rtFreq_   =         0.0;
  this->rtPhase_   =        0.0;
  this->rtSigma_   =        0.0;
  this->rtTOn_ =            0.0;
  this->rtTOff_ =       10000.0;
  this->rtEnvelope_ =       0;
  this->SCFdenTol_ =        1e-10;
  this->SCFeneTol_ =        1e-12;
  this->SCFmaxIter_ =       128;
  this->unitTest    =       0;
  this->field_      =       {0.0,0.0,0.0};
  this->nthreads    = 1;
 
  this->SCFdenTol_ = 1e-10;
  this->SCFeneTol_ = 1e-12;
  this->SCFmaxIter_ = 128;
  this->unitTest    = 0;
  this->field_      = {0.0,0.0,0.0};
};

void Controls::readSMP(int &n) {
#ifdef _OPENMP
  this->nthreads = n;
  omp_set_num_threads(n);
#endif
}

void Controls::readPSCF(std::fstream &in, std::fstream &out){
  this->doSDR = true;
  std::string readString;
  in >> readString;
  readString = stringupper(readString);
  if(!readString.compare("CIS"))        this->SDMethod = 1;
  else if(!readString.compare("RPA"))   this->SDMethod = 2;
  else if(!readString.compare("PPRPA")) this->SDMethod = 3;
  else if(!readString.compare("PPATDA")) this->SDMethod = 4;
  else if(!readString.compare("PPCTDA")) this->SDMethod = 5;
  else if(!readString.compare("STAB")) this->SDMethod = 6;
  else CErr("Input PSCF Option Not Recgnized",out);

  if(this->SDMethod != 6) in >> this->SDNSek;
  else                    this->SDNSek = 6;
}

void Controls::readDebug(std::string str){
  if(!str.compare("AOERI")){
    this->directTwoE = false;
    this->buildn4eri = true;
  } else {
    cout << "Debug option \"" << str << "\" not recognized\n";
  }
}

void Controls::printSettings(ostream &out){
//out << bannerTop << endl;
  out << endl << "ChronusQ Control Settings:" << endl << bannerMid << endl;

  out << std::setw(35) << std::left << "  Print Level:" << this->printLevel << endl;
  out << std::setw(35) << std::left << "  Shared Memory Threads (OpenMP)" << this->nthreads << endl;
  out << std::setw(35) << std::left << "  Performing SCF:";
  if(this->optWaveFunction) {
    out << "YES" << endl;
    out << std::setw(35) << std::left << "    Schwartz Bound Threshold:" << std::scientific << this->thresholdSchawrtz << endl;
    out << std::setw(35) << std::left << "    Direct:";
    if(this->directTwoE) out << "YES" << endl;
    else out << "NO" <<endl;
  } else out << "NO" << endl;

  out << std::setw(35) << std::left << "  Building Full ERIs:";
  if(this->buildn4eri) out << "YES" << endl;
  else out << "NO" << endl;

  out << std::setw(35) << std::left << "  Density Fitting:";
  if(this->doDF) out << "ON" << endl;
  else out << "OFF" << endl;

  out << std::fixed;
//out << bannerEnd << endl;
}
