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
#include <sdresponse.h>
using ChronusQ::SDResponse;

namespace ChronusQ {
template<>
void SDResponse<dcomplex>::printPrinciple(int iSt){
  double printTol = 0.1;

  this->fileio_->out << "  Principle Transitions   ( tol = 0.1 )" << endl;
  for(auto ia = 0; ia < this->nOAVA_; ia++){
    auto xIA_Alpha = ia; auto yIA_Alpha = ia + this->nSingleDim_/2;

    auto alphaOccOrb = (xIA_Alpha % this->nOA_) + 1;
    auto alphaVirOrb = (xIA_Alpha / this->nOA_) + this->nOA_ + 1;

    double absXIA_Alpha, absYIA_Alpha;

    absXIA_Alpha = std::norm((*this->transDen_)(xIA_Alpha,iSt));
    if(this->iMeth_ == RPA)
      absYIA_Alpha = std::norm((*this->transDen_)(yIA_Alpha,iSt));

    if(absXIA_Alpha > printTol)
        this->fileio_->out << "    "
                           << alphaOccOrb << "A -> " << alphaVirOrb << "A   "
                           << std::fixed << std::setw(10) << std::right 
                           << std::norm((*this->transDen_)(xIA_Alpha,iSt)) 
                           << std::fixed << std::setw(10) << std::right 
                           << std::arg((*this->transDen_)(xIA_Alpha,iSt)) << endl;
    if(this->iMeth_ == RPA){
      if(absYIA_Alpha > printTol)
          this->fileio_->out << "    "
                             << alphaOccOrb << "A <- " << alphaVirOrb << "A   "
                             << std::fixed << std::setw(10) << std::right 
                             << std::norm((*this->transDen_)(yIA_Alpha,iSt))
                             << std::fixed << std::setw(10) << std::right 
                             << std::arg((*this->transDen_)(yIA_Alpha,iSt)) << endl;
    }
  }
  for(auto ia = this->nOAVA_; ia < this->nOAVA_ + this->nOBVB_; ia++){
    auto xIA_Beta = ia; auto yIA_Beta = ia + this->nSingleDim_/2;

    auto betaOccOrb = ((xIA_Beta - this->nOAVA_) % this->nOB_) + 1;
    auto betaVirOrb = ((xIA_Beta - this->nOAVA_) / this->nOB_) + this->nOB_ + 1;

    double absXIA_Beta, absYIA_Beta;

    absXIA_Beta = std::norm((*this->transDen_)(xIA_Beta,iSt));
    if(this->iMeth_ == RPA)
      absYIA_Beta = std::norm((*this->transDen_)(yIA_Beta,iSt));

    if(absXIA_Beta > printTol)
        this->fileio_->out << "    "
                           << betaOccOrb << "B -> " << betaVirOrb << "B   "
                           << std::fixed << std::setw(10) << std::right 
                           << std::norm((*this->transDen_)(xIA_Beta,iSt))
                           << std::fixed << std::setw(10) << std::right 
                           << std::arg((*this->transDen_)(xIA_Beta,iSt)) << endl;
    if(this->iMeth_ == RPA){
      if(absYIA_Beta > printTol)
          this->fileio_->out << "    "
                             << betaOccOrb << "B <- " << betaVirOrb << "B   "
                             << std::fixed << std::setw(10) << std::right
                             << std::norm((*this->transDen_)(yIA_Beta,iSt))
                             << std::fixed << std::setw(10) << std::right 
                             << std::arg((*this->transDen_)(yIA_Beta,iSt)) << endl;
    }
  }
  this->fileio_->out << bannerMid << endl << endl;
}

template<>
void SDResponse<dcomplex>::printInfo() {
};//printInfo

template<>
void SDResponse<dcomplex>::printExcitedStateEnergies(){
  this->fileio_->out << bannerTop << endl;
  if(this->iMeth_ == CIS)
    this->fileio_->out << "CIS";
  else if(this->iMeth_ == RPA)
    this->fileio_->out << "RPA";
  else if(this->iMeth_ == STAB)
    this->fileio_->out << "Stability";
  this->fileio_->out << " Diagonalization for lowest " << this->nSek_ << " eigenstates" << endl;
  this->fileio_->out << bannerMid << endl << endl;
  
  for(auto iSt = 0; iSt < this->nSek_; iSt++){
    double Omega = (*this->omega_)(iSt);
    this->fileio_->out << "Excited State " << iSt+1 << ":" << endl;
    this->fileio_->out << "  \u03C9 = " << std::setw(10) << std::setprecision(7) << std::fixed << Omega                   << " Eh   ";
    this->fileio_->out << "  \u03C9 = " << std::setw(10) << std::setprecision(7) << std::fixed << Omega*phys.eVPerHartree << " eV   ";
    this->fileio_->out << "  \u03C9 = " << std::setw(10) << std::setprecision(7) << std::fixed << Omega*phys.nmPerHartree << " nm   " << endl;
    this->fileio_->out << "  f(" << 0 << "," << iSt+1 << ") = " << (*this->oscStrength_)(0,iSt+1) << endl;
    this->printPrinciple(iSt);
  }
} //printExcitedStateEnergies()
} // namespace ChronusQ
