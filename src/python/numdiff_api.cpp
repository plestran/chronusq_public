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
#include <workers.h>
#include <numdiff.h>

namespace ChronusQ {
  template<>
  double NumericalDifferentiation<double>::Wrapper_GSEnergy(){
    return this->singleSlater_undisplaced_->totalEnergy;
  };

  template<>
  boost::python::list NumericalDifferentiation<double>::Wrapper_GSGrad(){
    boost::python::list result;
    for(auto IX = 0; IX < 3*this->molecule_undisplaced_->nAtoms(); IX++){
      result.append(this->dervData_[IX].GS_GRAD);
    }
    return result;
  };

  template<>
  boost::python::list NumericalDifferentiation<double>::Wrapper_ESEnergy(){
    boost::python::list result;
    for(auto ISt = 0; ISt < this->responseNRoots_; ISt++)
    if(this->respType_ == RESPONSE_TYPE::CIS){
      result.append(this->response_undisplaced_->
          template frequencies<SINGLETS>()(ISt));
    } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
      result.append( this->response_undisplaced_->
          template frequencies<A_PPTDA_SINGLETS>()(ISt)
      );
    }
    return result;
  };

  template<>
  boost::python::list NumericalDifferentiation<double>::Wrapper_ESGrad(){
    boost::python::list result;
    for(auto IX = 0; IX < 3*this->molecule_undisplaced_->nAtoms(); IX++){
      boost::python::list esGrad;
      for(auto ISt = 0; ISt < this->responseNRoots_; ISt++)
        esGrad.append(this->dervData_[IX].ES_GRAD(ISt));
      result.append(esGrad);
    }
    return result;
  };

  template<>
  boost::python::list NumericalDifferentiation<double>::Wrapper_ESGSNACME(){
    boost::python::list result;
    if(!this->computeES2GSNACME || this->respType_ == RESPONSE_TYPE::PPTDA)
      return result;

    for(auto IX = 0; IX < 3*this->molecule_undisplaced_->nAtoms(); IX++){
      boost::python::list esgsNACME;
      for(auto ISt = 0; ISt < this->responseNRoots_; ISt++)
        esgsNACME.append(this->dervData_[IX].ES_GS_NACME(ISt));
      result.append(esgsNACME);
    }
    return result;
  };


  template<>
  boost::python::list NumericalDifferentiation<double>::Wrapper_ESESNACME(){
    boost::python::list result;

    if(!this->computeES2ESNACME )
      return result;

    for(auto IX = 0; IX < 3*this->molecule_undisplaced_->nAtoms(); IX++){
      boost::python::list esesNACME;
      for(auto ISt = 0; ISt < this->responseNRoots_; ISt++){
	boost::python::list row;
        for(auto JSt = 0; JSt < this->responseNRoots_; JSt++)
          row.append(this->dervData_[IX].ES_ES_NACME(ISt,JSt));
	esesNACME.append(row);
      }
      result.append(esesNACME);
    }
    return result;
  };

}; // namespace ChronusQ
