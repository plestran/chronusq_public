/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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
#include <global.h>
#include <molecule.h>
#include <basisset.h>
#include <aointegrals.h>
#include <singleslater.h>
#include <response.h>
#include <tools.h>
#include <classtools.h>


namespace ChronusQ {

  enum DiffType {
    TwoPointSymmetric
  };

  struct Derivatives {
    double          GS_ENERGY;
    double          GS_GRAD;
    Eigen::VectorXd ES_ENERGY;
    Eigen::VectorXd ES_GRAD;
    Eigen::VectorXd ES_GS_NACME;
    RealMatrix      ES_ES_NACME;
  };

  template<typename T>
  class NumericalDifferentiation {
    typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMatrix;

    DiffType          diffType_;
    Molecule        * molecule_undisplaced_;
    SingleSlater<T> * singleSlater_undisplaced_;
    Response<T>     * response_undisplaced_;

    std::vector<     Derivatives > dervData_;
    std::unique_ptr< Response<T> > generated_response_;

    bool generateESObjs_;
    RESPONSE_TYPE respType_;
    int responseDiffRoot_;
    int responseNRoots_;

    void checkPhase(SingleSlater<T>&,SingleSlater<T>&,TMatrix&);
    void checkPhase(SingleSlater<T>&,SingleSlater<T>&);
    void checkPhase(Response<T>&,Response<T>&,TMatrix&);
    void checkPhase(Response<T>&,Response<T>&);
    void checkPhase(TMatrix&, TMatrix&);
    void checkPhase(TMatrix&, TMatrix&, TMatrix&);

    void checkDegeneracies(SingleSlater<T>&);

  public:
    #include <numdiff_constructors.h>

    bool computeGSGradient;
    bool computeESGradient;
    bool computeES2GSNACME;
    bool computeES2ESNACME;
    bool doAllCartesianDOF;

    double step;
    
    // Setters
    inline void generateESObjs(){this->generateESObjs_ = true;};
    inline void setRespNRoots(int n){this->responseNRoots_ = n;};
    inline void setRespType(RESPONSE_TYPE type){this->respType_ = type;};
    inline void setRespRoot(int n){this->responseDiffRoot_ = n;};

    inline void setDiffType(DiffType diffType){this->diffType_ = diffType;}; 
    inline void setSingleSlater(SingleSlater<T> &ss){
      this->setMolecule(*ss.molecule());
      this->singleSlater_undisplaced_ = &ss;
      this->computeGSGradient         = true;
    };
    inline void setMolecule(Molecule &mol){
      this->molecule_undisplaced_ = &mol;
    };

    // Procedural functions
    inline void differentiate(){
      if(this->doAllCartesianDOF) this->cartesianDiff();
    }
    
    void cartesianDiff();
    void generateDispGeom(Molecule &, Molecule &, int, int);
    void computeGS(SingleSlater<T>&);
    void computeES(Response<T>&);
    double GSGradient(SingleSlater<T>&,SingleSlater<T>&);
    Eigen::VectorXd ESGradient(Response<T>&,Response<T>&);


    Eigen::VectorXd ES2GSNACME_CIS(SingleSlater<T>&,SingleSlater<T>&,
      Response<T>&,Response<T>&,TMatrix&,TMatrix&,TMatrix&,TMatrix&);
    Eigen::VectorXd ES2GSNACME_CIS(SingleSlater<T>&,SingleSlater<T>&,
      TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&);

    inline Eigen::VectorXd ES2GSNACME(SingleSlater<T> &ss_p1, 
      SingleSlater<T> &ss_m1, Response<T> &resp_p1, Response<T> &resp_m1, 
      TMatrix &SAO_0_p1, TMatrix &SAO_0_m1, TMatrix &SMO_0_p1, 
      TMatrix &SMO_0_m1){

      if(this->respType_ == RESPONSE_TYPE::CIS)
        return this->ES2GSNACME_CIS(ss_p1,ss_m1,resp_p1,resp_m1,SAO_0_p1,
          SAO_0_m1,SMO_0_p1,SMO_0_m1);
    }

    inline Eigen::VectorXd ES2GSNACME(SingleSlater<T> &ss_p1, 
      SingleSlater<T> &ss_m1, TMatrix &T_0, TMatrix &T_p1, TMatrix &T_m1, 
      TMatrix &SAO_0_p1, TMatrix &SAO_0_m1, TMatrix &SMO_0_p1, 
      TMatrix &SMO_0_m1){

      if(this->respType_ == RESPONSE_TYPE::CIS)
        return this->ES2GSNACME_CIS(ss_p1,ss_m1,T_0,T_p1,T_m1,SAO_0_p1,
          SAO_0_m1,SMO_0_p1,SMO_0_m1);
      else if(this->respType_ == RESPONSE_TYPE::PPTDA)
        return Eigen::VectorXd::Zero(this->responseNRoots_);
    }

    RealMatrix ES2ESNACME_CIS(SingleSlater<T>&,SingleSlater<T>&,
      Response<T>&,Response<T>&,TMatrix&,TMatrix&,TMatrix&,TMatrix&);
//  RealMatrix ES2ESNACME_PPTDA(SingleSlater<T>&,SingleSlater<T>&,
//    Response<T>&,Response<T>&,TMatrix&,TMatrix&,TMatrix&,TMatrix&){;};
    RealMatrix ES2ESNACME_CIS(SingleSlater<T>&,SingleSlater<T>&,
      TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&);
    RealMatrix ES2ESNACME_PPTDA(SingleSlater<T>&,SingleSlater<T>&,
      TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&,TMatrix&);

    inline RealMatrix ES2ESNACME(SingleSlater<T> &ss_p1, 
      SingleSlater<T> &ss_m1, Response<T> &resp_p1, Response<T> &resp_m1, 
      TMatrix &SAO_0_p1, TMatrix &SAO_0_m1, TMatrix &SMO_0_p1, 
      TMatrix &SMO_0_m1){
      if(this->respType_ == RESPONSE_TYPE::CIS)
        return this->ES2ESNACME_CIS(ss_p1,ss_m1,resp_p1,resp_m1,SAO_0_p1,
          SAO_0_m1,SMO_0_p1,SMO_0_m1);
      else if(this->respType_ == RESPONSE_TYPE::PPTDA)
        return this->ES2ESNACME_PPTDA(ss_p1,ss_m1,resp_p1,resp_m1,SAO_0_p1,
          SAO_0_m1,SMO_0_p1,SMO_0_m1);
    }

    inline RealMatrix ES2ESNACME(SingleSlater<T> &ss_p1, 
      SingleSlater<T> &ss_m1, TMatrix &T_0, TMatrix &T_p1, TMatrix &T_m1, 
      TMatrix &SAO_0_p1, TMatrix &SAO_0_m1, TMatrix &SMO_0_p1, 
      TMatrix &SMO_0_m1){
      if(this->respType_ == RESPONSE_TYPE::CIS)
        return this->ES2ESNACME_CIS(ss_p1,ss_m1,T_0,T_p1,T_m1,SAO_0_p1,
          SAO_0_m1,SMO_0_p1,SMO_0_m1);
      else if(this->respType_ == RESPONSE_TYPE::PPTDA)
        return this->ES2ESNACME_PPTDA(ss_p1,ss_m1,T_0,T_p1,T_m1,SAO_0_p1,
          SAO_0_m1,SMO_0_p1,SMO_0_m1);
    }


  }; // class NumericalDifferentiation
  #include <numdiff_procedural.h>
}; // namespace ChronusQ
