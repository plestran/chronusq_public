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
#include <wavefunction.h>
namespace ChronusQ {
/*
template<>
template<>
WaveFunction<dcomplex>::WaveFunction(const WaveFunction<double> &other) :
  Quantum<dcomplex>::Quantum<dcomplex>(
    dynamic_cast<const Quantum<double>&>(other)),
  basisset_     (other.basisset()   ),               
  molecule_     (other.molecule()   ),               
  fileio_       (other.fileio()     ),                 
  aointegrals_  (other.aointegrals()),            
  nBasis_       (other.nBasis()     ),
  nShell_       (other.nShell()     ),
  nTT_          (other.nTT()        ),
  nO_           (other.nO()         ),
  nV_           (other.nV()         ),
  nOA_          (other.nOA()        ),
  nOB_          (other.nOB()        ),
  nVA_          (other.nVA()        ),
  nVB_          (other.nVB()        ),
  multip_       (other.multip()     ),
  energyNuclei_ (other.energyNuclei()),
  totalEnergy_  (other.totalEnergy()){

  this->alloc();
  this->moA_->real() = *other.moA();
  if(this->nTCS_ == 2 and !this->isClosedShell)
    this->moB_->real() = *other.moB();
}
*/

}
