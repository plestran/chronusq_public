/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
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
  this->doDipole =          true;
  this->doQuadpole =        true;
  this->doOctpole =         true;
  this->doSDR     =         false;
  this->doCUHF    =         false;
  this->SDMethod  =         0;
  this->SCFdenTol_ = 1e-10;
  this->SCFeneTol_ = 1e-12;
  this->SCFmaxIter_ = 128;
#ifdef USE_LIBINT
  // Bootstrap Libint env
  libint2::init(); 
#endif
  this->nthreads = 1;
#ifdef USE_OMP
  // Set up Thread Pool
  this->nthreads = omp_get_max_threads();
  omp_set_num_threads(this->nthreads);
#endif
};

void Controls::readSMP(int &n) {
#ifdef USE_OMP
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
  else                    this->SDNSek = 3;
}

void Controls::readDebug(std::string str){
  if(!str.compare("AOERI")){
    this->directTwoE = false;
    this->buildn4eri = true;
  } else {
    cout << "Debug option \"" << str << "\" not recognized\n";
  }
}

void Controls::printSettings(fstream &out){
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

  out << std::setw(35) << std::left << "  Compute Dipole Moment:";
  if(this->doDipole) out << "YES" << endl;
  else out << "NO" << endl;
  out << std::setw(35) << std::left << "  Compute Quadrupole Moment:";
  if(this->doDipole) out << "YES" << endl;
  else out << "NO" << endl;
  out << std::setw(35) << std::left << "  Compute Octupole Moment:";
  if(this->doDipole) out << "YES" << endl;
  else out << "NO" << endl;
  out << std::fixed;
//out << bannerEnd << endl;
}
