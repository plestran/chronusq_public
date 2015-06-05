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
#ifndef INCLUDED_SDRESPONSE
#define INCLUDED_SDRESPONSE
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <controls.h>
#include <mointegrals.h>
#include <singleslater.h>
#include <basisset.h>

/****************************/
/* Error Messages 5000-5999 */
/****************************/

namespace ChronusQ {
class SDResponse {
  int       nBasis_;
  int       **R2Index_;
  int       nStates_;
  int       RHF_;
  friend class SingleSlater<double>;

  RealMatrix      XMO;
  RealMatrix      PDiag;
  std::unique_ptr<RealMatrix>    CISTransDen_;
  std::unique_ptr<RealMatrix>    CISEnergy_;
  std::unique_ptr<RealMatrix>    TransDipole_;
  BasisSet *      basisSet_;
  Molecule *      molecule_;
  FileIO *        fileio_;
  Controls *      controls_;
  MOIntegrals *   mointegrals_;
  SingleSlater<double> *  singleSlater_;
  RealTensor4d *  aoERI_;
  RealTensor3d *  elecDipole_;

public:
 
  // constructor & destructor
  SDResponse(){;};
  ~SDResponse() {;};
  // pseudo-constructor
  void iniSDResponse(Molecule *,BasisSet *,
                     MOIntegrals *,FileIO *,
                     Controls *,SingleSlater<double> *);

  inline int nOVA(){return this->singleSlater_->nOVA();};
  inline int nOVB(){return this->singleSlater_->nOVB();};
  void computeExcitedStates();         // compute the total electronic energy
  void printExcitedStateEnergies(); 
  void printInfo();
  void formRM();
  void DavidsonCIS();
  RealMatrix formRM2(RealMatrix &XMO);
  RealMatrix ReturnDiag();
  RealMatrix Guess(RealMatrix &PDiag);
  void TransDipole();
  double OscStrength();

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagSDResponse);
  void mpiRecv(int,int tag=tagSDResponse);
};
} // namespace ChronusQ
#endif
