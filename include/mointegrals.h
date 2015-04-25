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
#ifndef  INCLUDED_MOINTEGRAL
#define  INCLUDED_MOINTEGRAL
//#include <gsl/gsl_sf_erf.h>
#include "global.h"
#include "basisset.h"
//#include "matrix.h"
#include "molecule.h"
#include "fileio.h"
#include "controls.h"
#include "tools.h"
#include "aointegrals.h"
#include "singleslater.h"

/****************************/
/* Error Messages 8000-8999 */
/****************************/
namespace ChronusQ {
class MOIntegrals{
  int       **iaIndex_;
  int       **ijIndex_;
  int       **abIndex_;

  ChronusQ::BasisSet     	*basisSet_;
  ChronusQ::Molecule    	*molecule_;
  ChronusQ::FileIO       	*fileio_;
  ChronusQ::Controls     	*controls_;
  ChronusQ::AOIntegrals   *aointegrals_;
  SingleSlater  *singleSlater_;

public:
  // these should be protected
  RealMatrix    *iajb_;
  RealMatrix    *ijab_;
  RealMatrix    *ijka_;
  RealMatrix    *ijkl_;
  RealMatrix    *iabc_;
  RealMatrix    *abcd_;

  bool      haveMOiajb;
  bool      haveMOijab;
  bool      haveMOijka;
  bool      haveMOijkl;
  bool      haveMOiabc;
  bool      haveMOabcd;
 
  MOIntegrals(){;};
  ~MOIntegrals(){
    if(      iajb_!=NULL) delete iajb_;
    if(      ijab_!=NULL) delete ijab_;
    if(      ijka_!=NULL) delete ijka_;
    if(      ijkl_!=NULL) delete ijkl_;
    if(      iabc_!=NULL) delete iabc_;
    if(      abcd_!=NULL) delete abcd_;
  };
  
  // initialization function
  void iniMOIntegrals(ChronusQ::Molecule*,ChronusQ::BasisSet*,ChronusQ::FileIO*,ChronusQ::Controls*,ChronusQ::AOIntegrals*,SingleSlater*);

  inline double &iajb(int i, int a, int j, int b){
    return (*iajb_)(this->iaIndex_[i][a],this->iaIndex_[j][b]);
  };
  inline double &iajb(int ia, int jb){
    return (*iajb_)(ia,jb);
  };

  void formiajb();
  void formijab();
  void formijka();
  void formijkl();
  void formiabc();
  void formabcd();
};
} // namespace ChronusQ

#endif
