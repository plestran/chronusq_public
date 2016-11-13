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
#include <aointegrals.h>
using ChronusQ::AOIntegrals;

//---------------------
// do DKH0
//---------------------
void AOIntegrals::DKH0(){ 
/*
    This whole routine is crap and needs to be restructured -- JJG
*/
    auto NSQ         = this->nBasis_ * this->nBasis_;
    double C         = phys.SPEED_OF_LIGHT;

    char JOBZ = 'V';
    char UPLO = 'L';
    int INFO;

    double *A     = new double[NSQ]; 
    double *W     = new double[this->nBasis_];  
    double *WORK  = new double[NSQ];
    double *SCR   = new double[NSQ];
    double *OVLP1 = new double[NSQ];
    double *OVLP2 = new double[NSQ];

    RealVecMap E(W,this->nBasis_);
    RealMap    V(A,this->nBasis_,this->nBasis_);
    RealMap    S(WORK,this->nBasis_,this->nBasis_); // Requires WORK to be NBSq
    RealMap    T(SCR,this->nBasis_,this->nBasis_); // Requires WORK to be NBSq
    RealMap    X(OVLP1,this->nBasis_,this->nBasis_); // get us the overlap for lowdin decomp
    RealMap    Xp(OVLP2,this->nBasis_,this->nBasis_); // get us the overlap for lowdin decomp

    E.setZero();
    V.setZero();
    S.setZero();

    std::memcpy(SCR,this->kinetic_->data(),
      this->nBasis_*this->nBasis_*sizeof(double));

    // put the kinetic energy in the orthonormal basis
    X = (*this->overlap_).pow(-0.5); // Make this more efficient... FIXME
    Xp = (*this->overlap_).pow(0.5); // Make this more efficient... FIXME
    V = X.transpose() * T * X;

    dsyev_(&JOBZ,&UPLO,&this->nBasis_,A,&this->nBasis_,W,WORK,&NSQ,&INFO);

//    V.transposeInPlace(); // BC Col major
    S.setZero(); // S will become our new T matrix
  
    for(auto i = 0; i < this->nBasis_; i++) {
      cout << E(i) << endl;
      S(i,i) = std::sqrt(2.0*E(i)*C*C + C*C*C*C) - C*C;
    }
    
    T = Xp.transpose() * V * S * V.adjoint() * Xp; 
 
    std::memcpy(this->kinetic_->data(),SCR,
      this->nBasis_*this->nBasis_*sizeof(double));
    prettyPrint(this->fileio_->out,T, "new kinetic");

    delete[] A;
    delete[] W;
    delete[] WORK;
    delete[] SCR;


};

