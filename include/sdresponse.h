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
  int       nStates_;
  int       RHF_;
  int       nSek_;
  int       nGuess_;
  int       iMeth_;
  int       iPPRPA_;
  bool      haveDag_;
  int       nOA_;
  int       nVA_;
  int       nOB_;
  int       nVB_;
  int       nOAVA_;
  int       nOBVB_;
  int       nOAVB_;
  int       nOBVA_;
  int       nSingleDim_;
  double    rMu_;

  std::unique_ptr<RealCMMatrix>    transDen_;
  std::unique_ptr<RealMatrix>      oscStrength_;
  std::unique_ptr<VectorXd>        omega_;
  std::unique_ptr<RealTensor3d>    transDipole_;
  BasisSet *      basisSet_;
  Molecule *      molecule_;
  FileIO *        fileio_;
  Controls *      controls_;
  MOIntegrals *   mointegrals_;
  SingleSlater<double> *  singleSlater_;
  RealTensor4d *  aoERI_;
  RealTensor3d *  elecDipole_;

  std::unique_ptr<RealMatrix> tMO_;
  std::unique_ptr<RealCMMatrix> rmDiag_;
  std::unique_ptr<RealMatrix>  davGuess_;

public:
  enum{
    __invalid,
    CIS,
    RPA,
    PPRPA,
    PPATDA,
    PPCTDA,
    CCSD
  };
 
  // constructor & destructor
  SDResponse(){;};
  ~SDResponse() {;};
  // pseudo-constructor
  void iniSDResponse(Molecule *,BasisSet *,MOIntegrals *,FileIO *,Controls *,
                     SingleSlater<double> *);

  void computeExcitedStates();         // compute the total electronic energy
  void formRM();
  RealMatrix formRM2(RealMatrix &XMO);
  #include <sdresponse_getset.h>
  #include <sdresponse_qnrelated.h>
  #include <sdresponse_io.h>
  #include <sdresponse_misc.h>
  #include <sdresponse_prop.h>

};
} // namespace ChronusQ
#endif
