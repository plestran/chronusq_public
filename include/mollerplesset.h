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
#ifndef INCLUDED_MOLLERPLESSET
#define INCLUDED_MOLLERPLESSET
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <controls.h>
#include <mointegrals.h>
#include <singleslater.h>
#include <basisset.h>


namespace ChronusQ {
  class MollerPlesset {
    int nBasis_;
    int nTCS_;
    int Ref_;
    int order_;
    int nOA_;
    int nVA_;
    int nOB_;
    int nVB_;
    int nO_;         // NOA + NOB
    int nV_;         // NVA + NVB
    BasisSet *      basisSet_;
    Molecule *      molecule_;
    FileIO *        fileio_;
    Controls *      controls_;
    MOIntegrals *   mointegrals_;
    SingleSlater<double> *  singleSlater_;
  public:
    MollerPlesset(){;};
    ~MollerPlesset(){;};
    void iniMollerPlesset(Molecule *,BasisSet *,MOIntegrals *,FileIO *,Controls *,
                          SingleSlater<double> *);
    double MP2();
    double MP3();
  };
}

#endif
