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

void Controls::readDebug(std::string str){
  if(!str.compare("AOERI")){
    this->directTwoE = false;
    this->buildn4eri = true;
  } else {
    cout << "Debug option \"" << str << "\" not recognized\n";
  }
}
