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
#include <global.h>
#include <cerr.h>
#include <basisset.h>
#include <molecule.h>
#include <fileio.h>
#include <controls.h>
#include <tools.h>
#include <aointegrals.h>
#include <singleslater.h>

/****************************/
/* Error Messages 8000-8999 */
/****************************/
namespace ChronusQ {
class MOIntegrals{
  BasisSet *      basisSet_;
  Molecule *      molecule_;
  FileIO *        fileio_;
  Controls *      controls_;
  AOIntegrals *   aointegrals_;
  SingleSlater<double> *  singleSlater_;
/*
  std::unique_ptr<RealMatrix>    iajb_;
  std::unique_ptr<RealMatrix>    ijab_;
  std::unique_ptr<RealMatrix>    ijka_;
  std::unique_ptr<RealMatrix>    ijkl_;
  std::unique_ptr<RealMatrix>    iabc_;
  std::unique_ptr<RealMatrix>    abcd_;
*/
  std::unique_ptr<RealTensor4d>    iajb_;
  std::unique_ptr<RealTensor4d>    ijab_;
  std::unique_ptr<RealTensor4d>    ijka_;
  std::unique_ptr<RealTensor4d>    ijkl_;
  std::unique_ptr<RealTensor4d>    iabc_;
  std::unique_ptr<RealTensor4d>    abcd_;
  std::unique_ptr<RealTensor4d>    iajbAAAA_;
  std::unique_ptr<RealTensor4d>    ijabAAAA_;
  std::unique_ptr<RealTensor4d>    ijkaAAAA_;
  std::unique_ptr<RealTensor4d>    ijklAAAA_;
  std::unique_ptr<RealTensor4d>    iabcAAAA_;
  std::unique_ptr<RealTensor4d>    abcdAAAA_;
  std::unique_ptr<RealTensor4d>    iajbAABB_;
  std::unique_ptr<RealTensor4d>    ijabAABB_;
  std::unique_ptr<RealTensor4d>    ijkaAABB_;
  std::unique_ptr<RealTensor4d>    ijklAABB_;
  std::unique_ptr<RealTensor4d>    iabcAABB_;
  std::unique_ptr<RealTensor4d>    abcdAABB_;
  std::unique_ptr<RealTensor4d>    iajbBBBB_;
  std::unique_ptr<RealTensor4d>    ijabBBBB_;
  std::unique_ptr<RealTensor4d>    ijkaBBBB_;
  std::unique_ptr<RealTensor4d>    ijklBBBB_;
  std::unique_ptr<RealTensor4d>    iabcBBBB_;
  std::unique_ptr<RealTensor4d>    abcdBBBB_;

  std::unique_ptr<RealTensor2d>    locMOOcc_;
  std::unique_ptr<RealTensor2d>    locMOVir_;
  std::unique_ptr<RealTensor2d>    locMOAOcc_;
  std::unique_ptr<RealTensor2d>    locMOAVir_;
  std::unique_ptr<RealTensor2d>    locMOBOcc_;
  std::unique_ptr<RealTensor2d>    locMOBVir_;

  int nBasis_;
  int Ref_;
  int nTCS_;
  int nOA_;
  int nOB_;
  int nVA_;
  int nVB_;
  int nO_;
  int nV_;

  bool iajbIsDBar;
  bool ijabIsDBar;
  bool ijkaIsDBar;
  bool ijklIsDBar;
  bool iabcIsDBar;
  bool abcdIsDBar;

  bool haveLocMO;

  void getLocMO();

public:

  bool      haveMOiajb;
  bool      haveMOijab;
  bool      haveMOijka;
  bool      haveMOijkl;
  bool      haveMOiabc;
  bool      haveMOabcd;
 
  MOIntegrals(){;};
  ~MOIntegrals(){;};
  
  // initialization function
  void iniMOIntegrals(Molecule *,BasisSet *,
                      FileIO *,Controls *,
                      AOIntegrals *,SingleSlater<double> *);

  void formIAJB(bool);
  void formIJAB(bool);
  void formIJKA(bool);
  void formIJKL(bool);
  void formIABC(bool);
  void formABCD(bool);

  inline double IAJB(int i,int a,int j,int b,std::string spn="AAAA"){
    if(this->Ref_ == SingleSlater<double>::TCS){
      return (*this->iajb_)(i,a,j,b);
    } else {
      if(this->singleSlater_->isClosedShell){
        if(!spn.compare("AAAA") || !spn.compare("BBBB"))
          return (*this->iajbAAAA_)(i,a,j,b);
        else if(!spn.compare("AABB"))
          return (*this->iajbAABB_)(i,a,j,b);
        else CErr(spn+" is not a recognized spin order for IAJB",this->fileio_->out);
      } else {
        if(!spn.compare("AAAA"))
          return (*this->iajbAAAA_)(i,a,j,b);
        else if(!spn.compare("AABB"))
          return (*this->iajbAABB_)(i,a,j,b);
        else if(!spn.compare("BBBB"))
          return (*this->iajbBBBB_)(i,a,j,b);
        else CErr(spn+" is not a regocnized spin order for IAJB",this->fileio_->out);
      }
    }
  } 
};
} // namespace ChronusQ

#endif
