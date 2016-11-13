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
#include <gauinterface.h>

void GauJob::genInput(){
  this->gauInput_.open(this->gauFName_+".com", std::ios::out);
  if(!this->gauInput_){
    std::cout << "Error opening " << this->gauFName_+".com " << " for writing" << std::endl;
    exit(1);
  }

  std::string keywords = "#p ";
  keywords += "hf/"+this->basisName_+" ";
  if(this->doOpt_) keywords += "opt ";
  keywords += "output=RawMatrixElement";
  
  this->gauInput_ << keywords << std::endl;
  this->gauInput_ << std::endl << this->gauFName_ << std::endl;
  this->gauInput_ << std::endl << this->charge_ << " " << this->multip_<< std::endl;
  for(auto iAtm = 0; iAtm < this->atoms_.size(); iAtm++){
    this->gauInput_ << std::setw(5) << std::left << this->atoms_[iAtm];
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
      this->gauInput_ << std::setw(12) << std::setprecision(8) << std::fixed << this->cart_[iXYZ + iAtm*3];
    }
    this->gauInput_ << std::endl;
  }
  this->gauInput_ << std::endl << this->gauFName_ + ".matel" << std::endl;
  this->gauInput_ << "\n";
}
void GauJob::run(){
  system(("gdv "+this->gauFName_+".com").c_str());  
  this->matEl_ = std::unique_ptr<GauMatEl>(new GauMatEl(this->gauFName_+".matel"));

}
