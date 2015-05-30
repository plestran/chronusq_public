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
/******************************
 * Print Energy Contributions *
 ******************************/
template<typename T>
void SingleSlater<T>::printEnergy(){
  this->fileio_->out<<"\nEnergy Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"E(one electron) = "<<std::setw(15)<<this->energyOneE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"E(two electron) = "<<std::setw(15)<<this->energyTwoE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"E(nuclear repulsion) = "<<std::setw(15)<<this->energyNuclei<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<"E(total) = "<<std::setw(15)<<this->totalEnergy<<std::setw(5)<<" Eh "<<endl;
};

/**********************************
 * Print Wavefunction Information *
 **********************************/
template<typename T>
void SingleSlater<T>::printInfo() {
  this->fileio_->out<<"\nSingle Slater Determinent Wave Function Information:"<<endl;
  this->fileio_->out<<std::setw(15)<<"nOccA ="<<std::setw(8)<<this->nOccA_<<std::setw(5)<<" "<<std::setw(20)<<"nVirA ="<<std::setw(8)<<this->nVirA_<<endl;
  this->fileio_->out<<std::setw(15)<<"nOccB ="<<std::setw(8)<<this->nOccB_<<std::setw(5)<<" "<<std::setw(20)<<"nVirB ="<<std::setw(8)<<this->nVirB_<<endl;
  this->fileio_->out<<std::setw(15)<<"Multiplicity ="<<std::setw(8)<<this->spin_<<endl;
};
/***********************************************
 * Print the Multipole Moments (Electric Only) *
 ***********************************************/
template<typename T>
void SingleSlater<T>::printMultipole(){
  this->fileio_->out << bannerTop << endl;
  this->fileio_->out << std::setw(50) << std::left <<"Electric Dipole Moment"
                        << "(Debye)" << endl;
  this->fileio_->out << std::left << std::setw(5) <<"X=" 
                     << std::fixed << std::right << std::setw(20) 
                     << (*this->dipole_)(0,0)/phys.debye << endl;
  this->fileio_->out << std::left << std::setw(5) <<"Y=" 
                     << std::fixed << std::right << std::setw(20) 
                     << (*this->dipole_)(1,0)/phys.debye << endl;
  this->fileio_->out << std::left << std::setw(5) <<"Z=" 
                     << std::fixed << std::right << std::setw(20) 
                     << (*this->dipole_)(2,0)/phys.debye << endl;
  if(this->controls_->doQuadpole) {
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left << "Electric Quadrupole Moment" 
                       <<  "(Debye-\u212B)" << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(0,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(0,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(0,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(1,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(1,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(1,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(2,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(2,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->quadpole_)(2,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left << "Electric Quadrupole Moment (Traceless)" 
                       <<  "(Debye-\u212B)" << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(0,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(0,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(0,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(1,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(1,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(1,2)*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(2,0)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(2,1)*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->tracelessQuadpole_)(2,2)*phys.bohr/phys.debye << endl;
  }
  if(this->controls_->doOctpole) {
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left << "Electric Octupole Moment" 
                       << "(Debye-\u212B\u00B2)" << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,0,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,0,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,0,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,1,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,1,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,1,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,2,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,2,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(0,2,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,0,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,0,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,0,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,1,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,1,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,1,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,2,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,2,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(1,2,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,0,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,0,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,0,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,1,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,1,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,1,2)*phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,2,0)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,2,1)*phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << (*this->octpole_)(2,2,2)*phys.bohr*phys.bohr/phys.debye << endl;
  }
  this->fileio_->out << endl << bannerEnd << endl;
}
