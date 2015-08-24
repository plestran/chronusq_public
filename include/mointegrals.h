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
template<typename T>
class MOIntegrals{
  typedef Tensor<T,Range4d> TTensor4d; 
  typedef Tensor<T,Range2d> TTensor2d; 
  BasisSet *      basisSet_;
  Molecule *      molecule_;
  FileIO *        fileio_;
  Controls *      controls_;
  AOIntegrals *   aointegrals_;
  SingleSlater<T> *  singleSlater_;

  std::unique_ptr<TTensor4d>    iajb_;
  std::unique_ptr<TTensor4d>    ijab_;
  std::unique_ptr<TTensor4d>    ijka_;
  std::unique_ptr<TTensor4d>    ijkl_;
  std::unique_ptr<TTensor4d>    iabc_;
  std::unique_ptr<TTensor4d>    abcd_;
  std::unique_ptr<TTensor4d>    iajbAAAA_;
  std::unique_ptr<TTensor4d>    ijabAAAA_;
  std::unique_ptr<TTensor4d>    ijkaAAAA_;
  std::unique_ptr<TTensor4d>    ijklAAAA_;
  std::unique_ptr<TTensor4d>    iabcAAAA_;
  std::unique_ptr<TTensor4d>    abcdAAAA_;
  std::unique_ptr<TTensor4d>    iajbAABB_;
  std::unique_ptr<TTensor4d>    ijabAABB_;
  std::unique_ptr<TTensor4d>    ijkaAABB_;
  std::unique_ptr<TTensor4d>    ijklAABB_;
  std::unique_ptr<TTensor4d>    iabcAABB_;
  std::unique_ptr<TTensor4d>    abcdAABB_;
  std::unique_ptr<TTensor4d>    iajbBBBB_;
  std::unique_ptr<TTensor4d>    ijabBBBB_;
  std::unique_ptr<TTensor4d>    ijkaBBBB_;
  std::unique_ptr<TTensor4d>    ijklBBBB_;
  std::unique_ptr<TTensor4d>    iabcBBBB_;
  std::unique_ptr<TTensor4d>    abcdBBBB_;

  std::unique_ptr<TTensor2d>    locMOOcc_;
  std::unique_ptr<TTensor2d>    locMOVir_;
  std::unique_ptr<TTensor2d>    locMOAOcc_;
  std::unique_ptr<TTensor2d>    locMOAVir_;
  std::unique_ptr<TTensor2d>    locMOBOcc_;
  std::unique_ptr<TTensor2d>    locMOBVir_;

  // Only for use with complex (artifact of BTAS)
  std::unique_ptr<RealTensor2d>    reLocMOOcc_;
  std::unique_ptr<RealTensor2d>    reLocMOVir_;
  std::unique_ptr<RealTensor2d>    reLocMOAOcc_;
  std::unique_ptr<RealTensor2d>    reLocMOAVir_;
  std::unique_ptr<RealTensor2d>    reLocMOBOcc_;
  std::unique_ptr<RealTensor2d>    reLocMOBVir_;
  std::unique_ptr<RealTensor2d>    imLocMOOcc_;
  std::unique_ptr<RealTensor2d>    imLocMOVir_;
  std::unique_ptr<RealTensor2d>    imLocMOAOcc_;
  std::unique_ptr<RealTensor2d>    imLocMOAVir_;
  std::unique_ptr<RealTensor2d>    imLocMOBOcc_;
  std::unique_ptr<RealTensor2d>    imLocMOBVir_;

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
  void testLocMO();

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
                      AOIntegrals *,SingleSlater<T> *);

  void formIAJB(bool);
  void formIJAB(bool);
  void formIJKA(bool);
  void formIJKL(bool);
  void formIABC(bool);
  void formABCD(bool);

  inline T IAJB(int i,int a,int j,int b,std::string spn="AAAA"){
    if(this->Ref_ == SingleSlater<T>::TCS){
      return (*this->iajb_)(i,a,j,b);
    } else {
      if(this->singleSlater_->isClosedShell){
//cout << "(" << i << " " << a + this->nOA_ << " | " << j << " " << this->nOA_+b << ") " <<(*this->iajbAABB_)(i,a,j,b) << endl;
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
  inline T IJAB(int i,int j,int a,int b,std::string spn="AAAA"){
    if(this->Ref_ == SingleSlater<T>::TCS){
      return (*this->ijab_)(i,j,a,b);
    } else {
      if(this->singleSlater_->isClosedShell){
//cout << "(" << i << " " << a + this->nOA_ << " | " << j << " " << this->nOA_+b << ") " <<(*this->ijabAABB_)(i,j,a,b) << endl;
        if(!spn.compare("AAAA") || !spn.compare("BBBB"))
          return (*this->ijabAAAA_)(i,j,a,b);
        else if(!spn.compare("AABB"))
          return (*this->ijabAABB_)(i,j,a,b);
        else CErr(spn+" is not a recognized spin order for ijab",this->fileio_->out);
      } else {
        if(!spn.compare("AAAA"))
          return (*this->ijabAAAA_)(i,j,a,b);
        else if(!spn.compare("AABB"))
          return (*this->ijabAABB_)(i,j,a,b);
        else if(!spn.compare("BBBB"))
          return (*this->ijabBBBB_)(i,j,a,b);
        else CErr(spn+" is not a regocnized spin order for ijab",this->fileio_->out);
      }
    }
  } 
  inline T IABJ(int i,int a,int b,int j,std::string spn="AAAA"){
    if(this->Ref_ == SingleSlater<T>::TCS){
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
  inline T ABCD(int a,int b,int c,int d,std::string spn="AAAA"){
    if(this->Ref_ == SingleSlater<T>::TCS){
//    cout << "DjM " << this->abcd_->size() << endl;
      return (*this->abcd_)(a,b,c,d);
    } else {
      if(this->singleSlater_->isClosedShell){
//cout << "(" << a + this->nOA_ << " " << b + this->nOA_ << " | " << c + this->nOA_ << " " << this->nOA_+d << ") " <<(*this->abcdAABB_)(a,b,c,d) << endl;
        if(!spn.compare("AAAA") || !spn.compare("BBBB"))
          return (*this->abcdAAAA_)(a,b,c,d);
        else if(!spn.compare("AABB"))
          return (*this->abcdAABB_)(a,b,c,d);
        else CErr(spn+" is not a recognized spin order for abcd",this->fileio_->out);
      } else {
        if(!spn.compare("AAAA"))
          return (*this->abcdAAAA_)(a,b,c,d);
        else if(!spn.compare("AABB"))
          return (*this->abcdAABB_)(a,b,c,d);
        else if(!spn.compare("BBBB"))
          return (*this->abcdBBBB_)(a,b,c,d);
        else CErr(spn+" is not a regocnized spin order for abcd",this->fileio_->out);
      }
    }
  } 
  inline T IJKL(int i,int j,int k,int l,std::string spn="AAAA"){
    if(this->Ref_ == SingleSlater<T>::TCS){
//    cout << "DjM " << this->ijkl_->size() << endl;
      return (*this->ijkl_)(i,j,k,l);
    } else {
      if(this->singleSlater_->isClosedShell){
//cout << "(" << a + this->nOA_ << " " << b + this->nOA_ << " | " << c + this->nOA_ << " " << this->nOA_+d << ") " <<(*this->ijklAABB_)(i,j,k,l) << endl;
        if(!spn.compare("AAAA") || !spn.compare("BBBB"))
          return (*this->ijklAAAA_)(i,j,k,l);
        else if(!spn.compare("AABB"))
          return (*this->ijklAABB_)(i,j,k,l);
        else CErr(spn+" is not a recognized spin order for ijkl",this->fileio_->out);
      } else {
        if(!spn.compare("AAAA"))
          return (*this->ijklAAAA_)(i,j,k,l);
        else if(!spn.compare("AABB"))
          return (*this->ijklAABB_)(i,j,k,l);
        else if(!spn.compare("BBBB"))
          return (*this->ijklBBBB_)(i,j,k,l);
        else CErr(spn+" is not a regocnized spin order for ijkl",this->fileio_->out);
      }
    }
  } 
}; //template Class MOIntegrals
} // namespace ChronusQ

#endif
