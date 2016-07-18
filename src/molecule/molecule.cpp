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
#include <molecule.h>
using ChronusQ::Molecule;
//-------------//
// initializer //
//-------------//
void Molecule::alloc(std::ostream &out) {
  if(this->nAtoms_ < 1) CErr("No Atoms given to build Molecule",out);
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
  this->nAtoms_ = readInt;
  this->alloc(fileio->out);
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
  this->computeNucRep();
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
  if(this->printLevel_ > 0) {
    out << "Molecular Information:" << endl;
    out << bannerTop<< endl;
 
    out << std::setw(15) << "nAtoms =" << std::setw(8) << this->nAtoms_
    		<< std::setw(30) << "nElectrons   =" << std::setw(8) << this->nTotalE_
        << endl;
 
    out << std::setw(15) << "Charge =" << std::setw(8) << this->charge_ 
     		<< std::setw(30) << "Multiplicity ="  << std::setw(8) << this->multip_
        << endl;

    out << std::setw(26) << "Nuclear Repulsion =" << std::setw(24) 
        << std::setprecision(5) << std::scientific 
        << this->energyNuclei_ << " Eh"
        << endl;
    out << bannerEnd << endl;

    out << endl << "Cartesian Coordinates (bohr):" << endl;
    out << bannerTop<< endl;
    out << std::setw(18) << "Atom" << std::setw(21) << "X" 
				<< std::setw(15) << "Y" 	 << std::setw(15) << "Z" << endl;
    out << bannerMid << endl;
 
    for(auto i = 0; i < this->nAtoms_; i++)
      out << std::setw(8)  << i+1 << std::setw(8) << atom[index_[i]].symbol
          << std::setw(8)  << atom[index_[i]].atomicNumber
          << std::setw(15) << (*cart_)(0,i) 
          << std::setw(15) << (*cart_)(1,i)
          << std::setw(15) << (*cart_)(2,i)
          << endl;
    out << bannerMid << endl;
  }

  if(this->printLevel_ > 3){
    prettyPrint(out,*this->momentOfInertia_,
                "Moment of Inertia Tensor (AMU-bohr\u00B2)");
 
    prettyPrint(out,*this->rIJ_,"Interatomic Distance Matrix (bohr)");
  }
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

void Molecule::computeNucRep(){

  this->energyNuclei_ = 0.0;
  for(auto iAtm = 0       ; iAtm < this->nAtoms_; iAtm++)
  for(auto jAtm = iAtm + 1; jAtm < this->nAtoms_; jAtm++){
    auto rDist = this->cart_->col(iAtm) - this->cart_->col(jAtm);
    auto sqrDistAB = rDist.dot(rDist);
    this->energyNuclei_ += atom[this->index_[iAtm]].atomicNumber *
                           atom[this->index_[jAtm]].atomicNumber /
                           std::sqrt(sqrDistAB);
  }
}

void Molecule::generateFiniteWidthNuclei(){
  //
  // Generate libint2::Shell structs to hold a gaussian charge distribution
  // for each of the nuclei. The scheme is taken from Dyall, et al.
  //
  // L. Visscher and K.G. Dyall; Atomic Data and Nuclear Data Tables; 67,
  //   207-224 (1997)
  //
  //
  for(auto iAtm = 0; iAtm < this->nAtoms_; iAtm++){
    double varience = 
      0.836 * std::pow(elements[index_[iAtm]].massNumber,1.0/3.0)
      + 0.570; // fm
    varience *= 1e-5; // Ang
    varience /= phys.bohr; // Bohr

    varience *= varience;

    std::vector<double> zeta;
    zeta.push_back(3.0 / (2.0*varience));

    std::vector<double> cont;
    cont.push_back(elements[index_[iAtm]].atomicNumber);
    this->finiteWidthNuclei_.push_back(
        libint2::Shell{ zeta, {{0,false,cont}},
        {{(*this->cart())(0,iAtm),
           (*this->cart())(1,iAtm),
           (*this->cart())(2,iAtm)}}
        }
      );

    // Handle the fact that libint likes to make things square normalized
    // not normalized
    finiteWidthNuclei_.back().contr[0].coeff[0] = 
      elements[index_[iAtm]].atomicNumber * std::pow(zeta[0] / math.pi,1.5);

  };
};
