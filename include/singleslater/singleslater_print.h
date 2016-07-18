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
/******************************
 * Print Energy Contributions *
 ******************************/
template<typename T>
void SingleSlater<T>::printEnergy(){
  this->fileio_->out<<"\nEnergy Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<std::fixed<<"E(one electron) = "<<std::setw(15)<<this->energyOneE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<std::fixed<<"E(two electron) = "<<std::setw(15)<<this->energyTwoE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<std::fixed<<"E(nuclear repulsion) = "<<std::setw(15)<<this->energyNuclei<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<std::fixed<<"E(total) = "<<std::setw(15)<<this->totalEnergy<<std::setw(5)<<" Eh "<<endl;
};
/******************************
 * Print Energy Contributions *
 ******************************/

/**********************************
 * Print Wavefunction Information *
 **********************************/
template<typename T>
void SingleSlater<T>::printInfo() {
  this->fileio_->out<<"\nSingle Slater Determinent Wave Function Information:"<<endl;
  this->fileio_->out<<std::setw(15)<<"nOccA ="<<std::setw(8)<<this->nOccA_<<std::setw(5)<<" "<<std::setw(20)<<"nVirA ="<<std::setw(8)<<this->nVirA_<<endl;
  this->fileio_->out<<std::setw(15)<<"nOccB ="<<std::setw(8)<<this->nOccB_<<std::setw(5)<<" "<<std::setw(20)<<"nVirB ="<<std::setw(8)<<this->nVirB_<<endl;
  this->fileio_->out<<std::setw(15)<<"Multiplicity ="<<std::setw(8)<<this->multip_<<endl;
};
/***********************************************
 * Print the Multipole Moments (Electric Only) *
 ***********************************************/
template<typename T>
void SingleSlater<T>::printMultipole(){

  double xSmall = 1e-10;
  if(this->maxMultipole_ >= 1) {
    for(auto i = 0; i < 3; i++){
      if(std::abs(this->elecDipole_[i]) < xSmall)
        this->elecDipole_[i] = std::abs(this->elecDipole_[i]);
      if(this->maxMultipole_ >= 2) {
        for(auto j = 0; j < 3; j++){
          if(std::abs(this->elecQuadpole_[i][j]) < xSmall){
            this->elecQuadpole_[i][j] = std::abs(this->elecQuadpole_[i][j]);
            this->elecTracelessQuadpole_[i][j] = 
              std::abs(this->elecTracelessQuadpole_[i][j]);
          }
          for(auto k = 0; k < 3; k++){
            if(std::abs(this->elecOctpole_[i][j][k]) < xSmall){
              this->elecOctpole_[i][j][k] = 
                std::abs(this->elecOctpole_[i][j][k]);
            }
          }
        }
      }
    }
  }

  if(this->maxMultipole_ >= 1) {
    this->fileio_->out << "\nMultipole Information:" << endl;
    this->fileio_->out << bannerTop << endl;
    this->fileio_->out << std::setw(50) << std::left <<"Electric Dipole Moment"
                          << "(Debye)" << endl;
    this->fileio_->out << std::left << std::setw(5) <<"X=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecDipole_[0]/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" Y=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecDipole_[1]/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" Z=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecDipole_[2]/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"Tot=" 
                       << std::fixed << std::right << std::setw(20) 
                       << std::sqrt(
                             this->elecDipole_[2]*this->elecDipole_[2] + 
                             this->elecDipole_[1]*this->elecDipole_[1] + 
                             this->elecDipole_[0]*this->elecDipole_[0]  
                           )/phys.debye << endl;
   }

  if(this->maxMultipole_ >= 2) {
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left << "Electric Quadrupole Moment" 
                       <<  "(Debye-\u212B)" << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[0][0]*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[0][1]*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[0][2]*phys.bohr/phys.debye 
                       << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[1][0]*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[1][1]*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[1][2]*phys.bohr/phys.debye 
                       << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[2][0]*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[2][1]*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecQuadpole_[2][2]*phys.bohr/phys.debye 
                       << endl;
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left 
                       << "Electric Quadrupole Moment (Traceless)" <<  "(Debye-\u212B)" 
                       << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[0][0] * 
                         phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[0][1] * 
                         phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[0][2] * 
                         phys.bohr/phys.debye 
                       << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[1][0] * 
                         phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[1][1] * 
                         phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[1][2] * 
                         phys.bohr/phys.debye 
                       << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[2][0] * 
                         phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[2][1] * 
                         phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecTracelessQuadpole_[2][2] * 
                         phys.bohr/phys.debye 
                       << endl;
  }

  if(this->maxMultipole_ >= 3) {
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(50) << std::left 
                       << "Electric Octupole Moment" 
                       << "(Debye-\u212B\u00B2)" << endl;

    this->fileio_->out << std::left << std::setw(5) <<"XXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][0][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][0][1] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][0][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][1][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][1][1] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][1][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"XZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][2][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][2][1] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" XZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[0][2][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][0][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][0][1] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][0][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][1][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][1][1] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][1][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"YZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][2][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][2][1] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" YZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[1][2][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZXX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][0][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZXY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][0][1] *
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZXZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][0][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZYX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][1][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZYY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][1][1] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZYZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][1][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
    this->fileio_->out << std::left << std::setw(5) <<"ZZX=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][2][0] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZY=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][2][1] * 
                         phys.bohr*phys.bohr/phys.debye;
    this->fileio_->out << std::left << std::setw(5) <<" ZZZ=" 
                       << std::fixed << std::right << std::setw(20) 
                       << this->elecOctpole_[2][2][2] * 
                         phys.bohr*phys.bohr/phys.debye << endl;
  }
  this->fileio_->out << bannerEnd << endl << endl;
}

template<typename T>
void SingleSlater<T>::printSExpect(){
  this->fileio_->out << "Spin Information:" << endl;
  this->fileio_->out << "  <Sx> = " << this->Sx_ << endl;
  this->fileio_->out << "  <Sy> = " << this->Sy_ << endl;
  this->fileio_->out << "  <Sz> = " << this->Sz_ << endl;
  this->fileio_->out << "  <S\u00B2> = " << this->Ssq_ << endl;
};

template<typename T>
void SingleSlater<T>::printSCFHeader(ostream &output){
  output << bannerTop << endl;
  output << "Self Consistent Field (SCF) Settings:" << endl << endl;
//cout << std::setprecision(6);

  output << std::setw(38) << std::left << "  SCF Type:" << this->SCFType_ 
         << endl;

  output << std::setw(38)   << std::left << "  Density Convergence Tolerence:" 
         << std::scientific << std::setprecision(6) << this->denTol_ << endl;

  output << std::setw(38)   << std::left << "  Energy Convergence Tolerence:" 
         << std::scientific << std::setprecision(6) << this->eneTol_ << endl;

  output << std::setw(38) << std::left << "  Maximum Number of SCF Cycles:" 
         << this->maxSCFIter_ << endl;

  output << std::setw(38) << std::left << "  Integral Contraction Algorithm:";
  if(this->aointegrals_->integralAlgorithm == AOIntegrals::DIRECT)
    output << "Direct";
  else if (this->aointegrals_->integralAlgorithm == AOIntegrals::DENFIT)
    output << "Density-Fitting (BTAS)";
  else if (this->aointegrals_->integralAlgorithm == AOIntegrals::INCORE)
    output << "In-Core (BTAS)";
  output << endl;

  output << std::setw(38) << std::left << "  Initial Guess:";
  if(this->guess_ == SAD)
    output << "Superposition of Atomic Densities";
  else if(this->guess_ == CORE)
    output << "Core Hamiltonian";
  else if(this->guess_ == READ)
    output << "Read";
  output << endl;
  output << std::setw(38) << std::left << "  DIIS Extrapolation Algorithm:";
  if(this->doDIIS) output << "CDIIS";
  else             output << "No DIIS Extrapolation";
  output << endl;


  if(this->isDFT){
    output << std::setw(38) << std::left << "  Density Functional:";

    if(this->DFTKernel_ == USERDEFINED)
      output << "User Defined";
    else if(this->DFTKernel_ == LSDA)
      output << "LSDA";

    output << endl;

    
    if(this->DFTKernel_ == USERDEFINED){
      output << std::setw(38) << std::left << "    Exchange Kernel:";
      output << this->dftFunctionals_[0]->name;
      output << endl;

      output << std::setw(38) << std::left << "    Correlation Kernel:";
      output << this->dftFunctionals_[1]->name;
      output << endl;
   }

   
    output << std::setw(38) << std::left << "    Radial Grid:";
    if(this->dftGrid_ == EULERMACL)
      output << "Euler-Maclaurin";
    else if(this->dftGrid_ == GAUSSCHEB)
      output << "Gauss-Chebyshev (1st Kind)";
    output << "  (" << this->nRadDFTGridPts_ << ")";
    output << endl;

    output << std::setw(38) << std::left << "    Angular Grid:";
    output << "Lebedev";
    output << "  (" << this->nAngDFTGridPts_ << ")" << endl;

   
    output << std::setw(38) << std::left << "    Quadrature Weight Scheme:";
    if(this->weightScheme_ == BECKE)
      output << "Becke";
    else if(this->weightScheme_ == FRISCH)
      output << "Gaussian (as in the company...)";
    output << endl;

    if(this->screenVxc){
       output << std::setw(38) << std::left 
              << "    Quadrature Screen Tolerence:"
              << std::scientific << std::setprecision(6) << this->epsScreen 
              << endl;
    }
    
  }

  std::array<double,3> null = {{0,0,0}};
  if(this->elecField_ != null) {
    output << endl;
    output << std::setw(38) << std::left << "  Static Electric Field (Dipole):"
           << "{" << this->elecField_[0] << ", " << this->elecField_[1] << ", "
           << this->elecField_[2] << "}" << endl;
  }


  output << endl << bannerMid << endl;
  output << std::setw(16) << "SCF Iteration";
  output << std::setw(18) << "Energy (Eh)";
  output << std::setw(18) << "\u0394E (Eh)";
  if(this->Ref_ == TCS)
    output << std::setw(18) << "|\u0394P|";
  else {
    output << std::setw(18) << "|\u0394P(\u03B1)|";
    if(!this->isClosedShell)
      output << std::setw(18) << "|\u0394P(\u03B2)|";
  }
  output << endl;
  output << std::setw(16) << "-------------";
  output << std::setw(18) << "-----------";
  output << std::setw(18) << "-------";
  if(this->Ref_ == TCS)
    output << std::setw(18) << "----";
  else {
    output << std::setw(18) << "-------";
    if(!this->isClosedShell)
      output << std::setw(18) << "-------";
  }
  output << endl;
}

template<typename T>
void SingleSlater<T>::printSCFIter(int iter, double EDel,double PARMS,double PBRMS){
  this->fileio_->out << std::setw(16) << std::left 
                     << "  SCFIt: " + std::to_string(iter+1);
  this->fileio_->out << std::setw(18) << std::fixed << std::setprecision(10)
                     << std::left << this->totalEnergy;
  this->fileio_->out << std::setw(14) << std::scientific << std::right 
                     << std::setprecision(7) << EDel;
  this->fileio_->out << "   ";
  this->fileio_->out << std::setw(13) << std::scientific << std::right 
                     << std::setprecision(7) << PARMS;
  if(!this->isClosedShell && this->Ref_ != TCS) {
    this->fileio_->out << "   ";
    this->fileio_->out << std::setw(13) << std::scientific << std::right 
                       << std::setprecision(7) << PBRMS;
  }
  

  this->fileio_->out << endl;
}      
