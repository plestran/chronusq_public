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
/******************************
 * Print Energy Contributions *
 ******************************/
template<typename T>
void SingleSlater<T>::printEnergy(){
  this->fileio_->out<<"\nEnergy Information:"<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<std::fixed<<"E(one electron) = "<<std::setw(15)<<this->energyOneE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<std::fixed<<"E(two electron) = "<<std::setw(15)<<this->energyTwoE<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<std::fixed<<"E(nuclear repulsion) = "<<std::setw(15)<<this->energyNuclei_<<std::setw(5)<<" Eh "<<endl;
  this->fileio_->out<<std::right<<std::setw(30)<<std::fixed<<"E(total) = "<<std::setw(15)<<this->totalEnergy_<<std::setw(5)<<" Eh "<<endl;
};
/******************************
 * Print Energy Contributions *
 ******************************/

/**********************************
 * Print Wavefunction Information *
 **********************************/
template<typename T>
void SingleSlater<T>::printInfo() {
  this->fileio_->out << "\nSingle Slater Determinent Wave Function Information:" << endl;
  this->fileio_->out << std::setw(15) << "nOA = " << std::setw(8) << this->nOA_
                     << std::setw(5)  << " "      
                     << std::setw(20) << "nVA = " << std::setw(8) << this->nVA_ << endl;
  this->fileio_->out << std::setw(15) << "nOB = " << std::setw(8) << this->nOB_
                     << std::setw(5)  << " "      
                     << std::setw(20) << "nVB = " << std::setw(8) << this->nVB_ << endl;
  this->fileio_->out << std::setw(15) << "Multiplicity = " << std::setw(8) << this->multip_ << endl;
};
/***********************************************
 * Print the Multipole Moments (Electric Only) *
 ***********************************************/
template<typename T>
void SingleSlater<T>::printMultipole(){

  this->fileio_->out << std::left << std::setprecision(10);
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
    this->fileio_->out << bannerTop << endl << endl;;
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
    this->fileio_->out << endl << endl;
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
    this->fileio_->out << endl << endl;
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
    this->fileio_->out << endl << endl;
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
  this->fileio_->out << bannerTop << endl << endl;
  this->fileio_->out << std::setprecision(5) << std::fixed;
  this->fileio_->out << "  <Sx> = " << std::setw(10) << std::right 
                     << this->Sx_ << endl;
  this->fileio_->out << "  <Sy> = " << std::setw(10) << std::right 
                     << this->Sy_ << endl;
  this->fileio_->out << "  <Sz> = " << std::setw(10) << std::right 
                     << this->Sz_ << endl;
  this->fileio_->out << "  <S\u00B2> = " << std::setw(10) << std::right 
                     << this->Ssq_ << endl;
};

template<typename T>
void SingleSlater<T>::printSCFHeader(ostream &output){
  output << BannerTop << endl;
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
  else if(this->guess_ == RANDOM)
    output << "Random";
  output << endl;
  output << std::setw(38) << std::left << "  DIIS Extrapolation Algorithm:";
  if(this->doDIIS) output << "CDIIS";
  else             output << "No DIIS Extrapolation";
  output << endl;
/*
  output << std::setw(38) << std::left << "  Imaginary Time Propagation:";
  if(this->doITP) output << "True, with dt = " << this->dt;
  else            output << "False"; 
  output << endl;
*/
  if(this->doDamp) {
    output << std::setw(38) << std::left << "  Damping Parameter:";
    output << std::scientific << std::setprecision(4) << this->dampParam;
    output << endl;
  }

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

      if(dftFunctionals_.size() > 1){
        output << std::setw(38) << std::left << "    Correlation Kernel:";
        output << this->dftFunctionals_[1]->name;
        output << endl;
      }
   }

   
    output << std::setw(38) << std::left << "    Radial Grid:";
    if(this->dftGrid_ == EULERMAC)
      output << "Euler-Maclaurin";
    else if(this->dftGrid_ == GAUSSCHEBFST)
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
  if(this->nTCS_ == 2)
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
  if(this->nTCS_ == 2)
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
                     << std::left << this->totalEnergy_;
  this->fileio_->out << std::setw(14) << std::scientific << std::right 
                     << std::setprecision(7) << EDel;
  this->fileio_->out << "   ";
  this->fileio_->out << std::setw(13) << std::scientific << std::right 
                     << std::setprecision(7) << PARMS;
  if(!this->isClosedShell && this->nTCS_ == 1) {
    this->fileio_->out << "   ";
    this->fileio_->out << std::setw(13) << std::scientific << std::right 
                       << std::setprecision(7) << PBRMS;
  }
  

  this->fileio_->out << endl;
}      

template <typename T>
void SingleSlater<T>::printDensity() {
  prettyPrintSmart(this->fileio_->out,(*this->onePDMScalar_),
    "Density (Scalar)");
  if(this->nTCS_ == 2 or !this->isClosedShell)
    prettyPrintSmart(this->fileio_->out,(*this->onePDMMz_),
      "Density (Mz)");
  if(this->nTCS_ == 2) {
    prettyPrintSmart(this->fileio_->out,(*this->onePDMMx_),
      "Density (Mx)");
    prettyPrintSmart(this->fileio_->out,(*this->onePDMMy_),
      "Density (My)");
  }
};

template <typename T>
void SingleSlater<T>::printFock() {
  prettyPrintSmart(this->fileio_->out,(*this->fockScalar_),
    "Fock (Scalar)");
  if(this->nTCS_ == 2 or !this->isClosedShell)
    prettyPrintSmart(this->fileio_->out,(*this->fockMz_),
      "Fock (Mz)");
  if(this->nTCS_ == 2) {
    prettyPrintSmart(this->fileio_->out,(*this->fockMx_),
      "Fock (Mx)");
    prettyPrintSmart(this->fileio_->out,(*this->fockMy_),
      "Fock (My)");
  }
};

template <typename T>
void SingleSlater<T>::printPT() {
  prettyPrintSmart(this->fileio_->out,(*this->PTScalar_),
    "PT (Scalar)");
  if(this->nTCS_ == 2 or !this->isClosedShell)
    prettyPrintSmart(this->fileio_->out,(*this->PTMz_),
      "PT (Mz)");
  if(this->nTCS_ == 2) {
    prettyPrintSmart(this->fileio_->out,(*this->PTMx_),
      "PT (Mx)");
    prettyPrintSmart(this->fileio_->out,(*this->PTMy_),
      "PT (My)");
  }
};

template <typename T>
void SingleSlater<T>::printSCFResults(){
  this->fileio_->out << endl << endl;
  this->fileio_->out << "SCF Results:" << endl;
  this->fileio_->out << BannerTop << endl << endl;;

  this->fileio_->out << "Orbital Eigenvalues";
  if(this->nTCS_ == 1)
    this->fileio_->out << " (Alpha)";
  this->fileio_->out << " / Eh :" << endl << bannerTop << endl;

  int no = this->nTCS_ == 2 ? this->nO_ : this->nOA_;
  for(auto e = 0 ; e !=  this->epsA_->size(); ++e){
    if(e == 0) this->fileio_->out << "Occupied:" << endl;
    else if(e == no) this->fileio_->out << endl << "Unoccupied:" << endl;
    this->fileio_->out << std::setw(13) << std::scientific 
                       << std::setprecision(4) <<
                        (*this->epsA_)(e);
    if(e != 0 and (e+1) % 5 == 0) this->fileio_->out << endl;
  }
  this->fileio_->out << endl;

  if(this->nTCS_ == 1 and !this->isClosedShell) {
    this->fileio_->out << endl;
    this->fileio_->out << "Orbital Eigenvalues (Beta) / Eh :";
    this->fileio_->out << endl << bannerTop << endl;

    for(auto e = 0 ; e !=  this->epsB_->size(); ++e){
      if(e == 0) this->fileio_->out << "Occupied:" << endl;
      else if(e == no) this->fileio_->out << endl << "Unoccupied:" << endl;
      this->fileio_->out << std::setw(13) << std::scientific 
                         << std::setprecision(4) <<
                          (*this->epsB_)(e);
      if(e != 0 and (e+1) % 5 == 0) this->fileio_->out << endl;
    }
  }
}

template <typename T>
void SingleSlater<T>::printTimings() {

    this->fileio_->out << endl << "SCF Timings: "<<endl << BannerTop 
      << endl;
    this->fileio_->out << endl << "Average timing over SCF iterations" << endl 
                       << bannerMid << endl;


    this->fileio_->out << std::left << std::setw(60) <<
      "Average Wall time for Fock build:";
    this->fileio_->out << std::left << std::setw(15) 
      << this->avgFockD_.count() << " sec" << endl;

    this->fileio_->out << std::left << std::setw(60) <<
      "  Average Wall time for G[P]:";
    this->fileio_->out << std::left << std::setw(15) 
      << this->avgPTD_.count() << " sec" << endl;

    if(this->isDFT) {
      this->fileio_->out << std::left << std::setw(60) <<
        "  Average Wall time for EXC / VXC:";
      this->fileio_->out << std::left << std::setw(15) 
        << this->avgVXCD_.count() << " sec" << endl;
    }

    this->fileio_->out << std::left << std::setw(60) <<
      "Average Wall time for Fock diagonalization:";
    this->fileio_->out << std::left << std::setw(15) 
      << this->avgDiagD_.count() << " sec" << endl;

    this->fileio_->out << std::left << std::setw(60) <<
      "Average Wall time spent on orthonormal transformation:";
    this->fileio_->out << std::left << std::setw(15) 
      << this->avgOrthoD_.count() << " sec" << endl;

    if(this->doDIIS){
      this->fileio_->out << std::left << std::setw(60) <<
        "Average Wall time for [F,P]:";
      this->fileio_->out << std::left << std::setw(15) 
        << this->avgCommD_.count() << " sec" << endl;
    }

    this->fileio_->out << std::left << std::setw(60) << " "
      << std::left << std::setw(15) << "---------------" 
      << "----" << endl;

    this->fileio_->out << std::left << std::setw(60) 
      << "Total wall time for SCF procedure:" 
      << std::left << std::setw(15) << this->SCFD_.count() << " sec" 
      << endl;

    this->fileio_->out << BannerEnd << endl;
};
