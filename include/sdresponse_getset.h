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

  // Setters
  inline void setNSek(int n){   this->nSek_  = n; this->nGuess_ = 2*n; };
  inline void setMeth(int n){   this->iMeth_ = n; this->initMeth();    };
  inline void setNGuess(int n){ this->nGuess_ = n;                     };
  inline void setPPRPA(int n){ this->iPPRPA_ = n;                      };

  // Getters
  inline int           nGuess(){     return this->nGuess_;               };
  inline int           nSek(){       return this->nSek_;                 };
  inline int           iMeth(){      return this->iMeth_;                };
  inline int           nSingleDim(){ return this->nSingleDim_;           };
  inline int           nOVA(){       return this->singleSlater_->nOVA(); };
  inline int           nOVB(){       return this->singleSlater_->nOVB(); };
  inline VectorXd*     omega(){      return this->omega_.get();          };
  inline TCMMatrix*    transDen(){   return this->transDen_.get();       };
  inline RealCMMatrix* rmDiag(){     return this->rmDiag_.get();         };
  inline TMatrix*      davGuess(){   return this->davGuess_.get();       };
  inline FileIO*       fileio(){     return this->fileio_;               };
