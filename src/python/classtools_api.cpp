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
#include <workers.h>
#include <global.h>

using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::AOIntegrals;

namespace ChronusQ {
  void Wrapper_readInput(FileIO &fileio, Molecule &mol, BasisSet &basis, 
    Controls &controls, BasisSet &dfBasis) {
    readInput(&fileio,&mol,&basis,&controls,&dfBasis);
  }

  void Wrapper_initCQ(int argc, boost::python::list argv){
    char ** argv_ = new char *[argc];
    for(auto i = 0; i < argc; i++){
      std::string str = boost::python::extract<std::string>(argv[i]);
      auto len = str.length();
      argv_[i] = new char[len+1];
      strncpy(argv_[i],&str[0],len);
      argv_[i][len] = '\0'; // Termination character
    }
    initCQ(argc,argv_);
    for(auto i = 0; i < argc; i++){
      delete [] argv_[i];
    }
    delete [] argv_;
  };
};

