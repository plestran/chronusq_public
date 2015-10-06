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
#include <molecule.h>
using ChronusQ::Molecule;
//-------------//
// initializer //
//-------------//
void Molecule::iniMolecule(int nAtoms, std::ostream &out) {
  if(nAtoms < 1) CErr("No Atoms given to build Molecule",out);
  this->nAtoms_= nAtoms;
  this->cart_  = std::unique_ptr<RealMatrix>(new RealMatrix(3, this->nAtoms_));
  this->index_ = new int[this->nAtoms_];
};
//--------------------------------------------//
// Read molecular information from input file //
//--------------------------------------------//
void Molecule::readMolecule(FileIO * fileio, std::istream &geomRead){
  int i, j, n, readInt;
  std::string readString;
  geomRead >> readInt;
  iniMolecule(readInt,fileio->out);
  nTotalE_ = 0;
  for(i=0;i<nAtoms_;i++) {
    geomRead >> readString;
    geomRead >> readInt;
    if((n=HashAtom(readString,readInt))!=-1) index_[i] = n;
    else {
      CErr("Error: invalid atomic symbol or mass number!",fileio->out);
    };
    nTotalE_ += atom[n].atomicNumber;
    geomRead >> (*cart_)(0,i);
    (*cart_)(0,i) = (*cart_)(0,i)/phys.bohr;
    geomRead >> (*cart_)(1,i);
    (*cart_)(1,i) = (*cart_)(1,i)/phys.bohr;
    geomRead >> (*cart_)(2,i);
    (*cart_)(2,i) = (*cart_)(2,i)/phys.bohr;
  };
  double sqrAB;
  energyNuclei_ = 0.0;
  for(i=0;i<nAtoms_;i++) 
    for(j=i+1;j<nAtoms_;j++) {
      sqrAB = 0.0;
      for(n=0;n<3;n++) sqrAB += ((*cart_)(n,i)-(*cart_)(n,j))*((*cart_)(n,i)-(*cart_)(n,j));
      energyNuclei_ += atom[index_[i]].atomicNumber*atom[index_[j]].atomicNumber/sqrt(sqrAB);
  };
  this->computeRij();
  this->toCOM(0);
  this->computeI();

};
//---------------------------------------------------//
// Print out molecular carteisan coordinates in bohr //
//---------------------------------------------------//
void Molecule::printInfo(std::ostream &out) {
  out.precision(8);
  out.fill(' ');
  out.setf(ios::right,ios::adjustfield);
  out.setf(ios::fixed,ios::floatfield);
  out<<"\nMolecular Information:"<<endl;
  out<<std::setw(15)<<"nAtoms ="<<std::setw(8)<<nAtoms_<<std::setw(5)<<" "
	     <<std::setw(20)<<"Charge ="<<std::setw(8)<<charge_<<endl;
  out<<std::setw(15)<<"nElectrons ="<<std::setw(8)<<nTotalE_<<endl;
  out<<"\nCartesian coordinates (bohr):"<<endl;
  out<<bannerTop<<endl;
  out<<std::setw(24)<<" "<<std::setw(15)<<"X"<<std::setw(15)<<"Y"<<std::setw(15)<<"Z"<<endl;
  out<<bannerMid<<endl;
  for(int i=0;i<nAtoms_;i++)
    out<<std::setw(8)<<i+1<<std::setw(8)<<atom[index_[i]].symbol<<std::setw(8)<<atom[index_[i]].atomicNumber
	       <<std::setw(15)<<(*cart_)(0,i)<<std::setw(15)<<(*cart_)(1,i)<<std::setw(15)<<(*cart_)(2,i)<<endl;
  out<<bannerMid<<endl;
  prettyPrint(out,*this->momentOfInertia_,"Moment of Inertia Tensor (AMU-bohr\u00B2)");
  prettyPrint(out,*this->rIJ_,"Interatomic Distance Matrix (bohr)");
};

void Molecule::toCOM(int Iop){
  int nAtoms;
  this->COM_ = std::unique_ptr<VectorXd>(new VectorXd(3)); 
  double TotW = 0;
  if(Iop == 0){
    for(auto iA = 0; iA < this->nAtoms_; iA++){
       TotW += elements[this->index_[iA]].mass;
       (*this->COM_) += (this->cart_->col(iA)) * 
                        elements[this->index_[iA]].mass;
    }
  }
  else if(Iop == 1){
    for(auto iA=0; iA < this->nAtoms_; iA++){
       TotW += elements[this->index_[iA]].atomicNumber;
       (*this->COM_)+=(this->cart_->col(iA)) * 
                      elements[this->index_[iA]].atomicNumber;
     }
  }
}

void Molecule::computeI(){
  this->momentOfInertia_ = std::unique_ptr<RealMatrix>(new RealMatrix(3,3));

  // Assuming X,Y,Z unit vectors as intertial frame
  RealMatrix E = RealMatrix::Identity(3,3); 

  for(auto iAtm = 0; iAtm < this->nAtoms_; iAtm++){
    *this->momentOfInertia_ += elements[this->index_[iAtm]].mass*(
                     this->cart_->col(iAtm).dot(this->cart_->col(iAtm))*E -
                     this->cart_->col(iAtm)*this->cart_->col(iAtm).transpose());
  }
  
}
void Molecule::computeRij(){
  this->rIJ_ = 
    std::unique_ptr<RealMatrix>(new RealMatrix(this->nAtoms_,this->nAtoms_));
  for(auto iAtm = 0; iAtm < this->nAtoms_; iAtm++)
  for(auto jAtm = 0; jAtm < iAtm;          jAtm++){
    (*this->rIJ_)(iAtm,jAtm) = (cart_->col(iAtm) - cart_->col(jAtm)).norm();
  }
  (*this->rIJ_) = this->rIJ_->selfadjointView<Lower>();
}
