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
#include <basisset.h>

namespace ChronusQ{

  /**
   *  Print the header for BasisSet::printInfo
   */ 
  void BasisSet::printHeader(){
    this->fileio_->out << endl << "Basis Set Information:" << endl;;
    this->fileio_->out << bannerTop << endl;
  } // BasisSet::printHeader

  /**
   *  Print BasisSet metadata:
   *    nBasis
   *    nPrimitive
   *    maxPrim
   *    nShell
   *    nLShell
   *    nShellPair
   */ 
  void BasisSet::printMeta(){
    
    this->fileio_->out << "  Basis Set path: " << this->basisPath_ << endl;
    this->fileio_->out << endl;

    this->fileio_->out << "  nBasis      =  " <<  this->nBasis_     << endl; 
    this->fileio_->out << "  nPrimitive  =  " <<  this->nPrimitive_ << endl; 
    if(this->printLevel_ > 1) {
      this->fileio_->out << "  maxPrim     =  " <<  this->maxPrim_    << endl; 
      this->fileio_->out << "  maxL        =  " <<  this->maxL_       << endl; 
    }
    this->fileio_->out << "  nShell      =  " <<  this->nShell_     << endl; 
    if(this->printLevel_ > 2) { 
      for(auto i = 0; i < this->nLShell_.size(); i++){
        this->fileio_->out << "  n" << HashS(i) << "Shell     =  "
                           << this->nLShell_[i] << endl;
      } // loop i
     
      this->fileio_->out << "  nShellPair  =  " <<  this->nShellPair_ << endl; 
    }
    this->fileio_->out << "  Using Cartesian Functions?: ";
    if(this->forceCart_) this->fileio_->out << "Yes";
    else                 this->fileio_->out << "No";
    this->fileio_->out << endl;

    this->fileio_->out << bannerEnd<< endl << endl;
  } // BasisSet::printMeta


  /**
   *  Print Basis Set definition
   */ 
  void BasisSet::printBasis(){
    this->fileio_->out.precision(6);
    this->fileio_->out.fill(' ');
    this->fileio_->out.setf(ios::right,ios::adjustfield);
    this->fileio_->out.setf(ios::fixed,ios::floatfield);
    this->fileio_->out << endl << "Cartesian Basis Functions:" << endl;
    this->fileio_->out << endl;
    this->fileio_->out << std::setw(16) << "        Shell Type" << "    Center" << 
      std::setw(8) << "L" << std::setw(26) << "exponent" << std::setw(18) << 
      "coefficient" << endl;
    this->fileio_->out << bannerMid << endl;

    auto iShell = 0;
    for(auto shell : this->shells_){
      this->fileio_->out << std::setw(5)  << iShell + 1;
      this->fileio_->out << std::setw(8)  << HashS(shell.contr[0].l);
      this->fileio_->out << std::setw(13) << this->mapSh2Cen_[iShell] << "  "; 
      this->fileio_->out << std::setw(8)  << shell.contr[0].l;
 
      for(auto iPrim = 0; iPrim < shell.alpha.size(); iPrim++){
        if(iPrim != 0) this->fileio_->out << std::setw(36) << " ";

        this->fileio_->out << std::setw(26) << shell.alpha[iPrim];
        this->fileio_->out << std::setw(18) << shell.contr[0].coeff[iPrim];
        this->fileio_->out << endl;
      } // loop iPrim

      this->fileio_->out << endl;
      iShell++;
    } // loop shell
  } // BasisSet::printBasis

  /**
   *  Print important information about a BasisSet object
   */ 
  void BasisSet::printInfo(){
    this->printHeader();
    this->printMeta();
    if(this->printLevel_ > 3) this->printBasis();
  } // BasisSet::printInfo
} // namespace ChronusQ

