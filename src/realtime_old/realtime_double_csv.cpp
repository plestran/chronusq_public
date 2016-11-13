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
#include <realtime.h>
#include <aointegrals.h>
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;
using ChronusQ::RealTime;

namespace ChronusQ {

template<>
void RealTime<double>::writeDipoleCSV(PropInfo & rec, long int & iStep){
    if (iStep == 0) {
      *csvs[0] << "Time Step (a.u.), Energy (Eh), Dipole X debye), Dipole Y (debye),"
               << "Dipole Z (debye), Dipole Tot (debye)" << endl;
    }
    *csvs[0] << std::fixed << std::setprecision(10) << rec.timeStep << ", " << rec.energy << ", " << rec.dipole[0] 
             << ", " << rec.dipole[1] << ", " << rec.dipole[2] << ", " << rec.dipole[3] << endl;
};

template<>
void RealTime<double>::writeAppliedFieldCSV(PropInfo & rec, long int & iStep){
    if (iStep == 0) {
      *csvs[1] << "Time Step (a.u.), Ex (a.u), Ey (a.u.), Ez (a.u.), || E total || (a.u.)" << endl;
    }
    *csvs[1] << std::fixed << std::setprecision(10) << rec.timeStep << ", " << rec.appliedfield[0] 
             << ", " << rec.appliedfield[1] << ", " << rec.appliedfield[2] << ", " << rec.appliedfield[3] << endl;
};

template<>
void RealTime<double>::writeMullikenCSV(PropInfo & rec, long int & iStep){
  if (iStep == 0) {   
    *csvs[2] << std::setw(14) << "Atom number";
    for(auto iAtm = 0; iAtm < this->ssPropagator_->molecule()->nAtoms(); iAtm++) {
      *csvs[2] << std::setw(14) << iAtm;
    }
    *csvs[2] << endl;
    *csvs[2] << std::setw(14) << "Atom symbol";
    for(auto iAtm = 0; iAtm < this->ssPropagator_->molecule()->nAtoms(); iAtm++) {
      *csvs[2] << std::setw(14) << elements[this->ssPropagator_->molecule()->index(iAtm)].symbol;
    }
    *csvs[2] << endl;
    *csvs[2] << std::setw(14) << "Time (a.u.)" << endl;
  }
    *csvs[2] << std::fixed << std::setw(14) << std::setprecision(10) << rec.timeStep;
    for(auto iAtm = 0; iAtm < this->ssPropagator_->molecule()->nAtoms(); iAtm++) {
      *csvs[2] << ", " << std::setw(14) << rec.mullPop[iAtm];
    }
    *csvs[2] << endl;
};

template<>
void RealTime<double>::writeOrbitalCSV(PropInfo & ref, long int & iStep){
  if(iStep == 0) {
    *csvs[3] << std::setw(10) << "Time (au)";
    for(auto idx = 0; idx != this->nBasis_*this->nTCS_; idx++){
      *csvs[3] << ", " << std::setw(10) << idx + 1;
    }
    *csvs[3] << endl;
  }
  *csvs[3] << std::fixed << std::setw(10) << std::setprecision(6) << ref.timeStep; 
  for(auto idx = 0; idx != this->nBasis_*this->nTCS_; idx++) {
    *csvs[3] << ", " << std::setw(10) << std::setprecision(6) << ref.orbitalOccA[idx];
  }
  *csvs[3] << endl;
  if(!this->isClosedShell_ && this->Ref_ != SingleSlater<double>::TCS){
    if(iStep == 0) {
      for(auto idx = 0; idx != this->nBasis_*this->nTCS_; idx++){
        *csvs[4] << ", " << std::setw(10) << idx + 1;
      }
      *csvs[4] << endl;
    }
    *csvs[4] << std::fixed << std::setw(10) << std::setprecision(6) << ref.timeStep; 
    for(auto idx = 0; idx != this->nBasis_*this->nTCS_; idx++) {
      *csvs[4] << ", " << std::setw(10) << std::setprecision(6) << ref.orbitalOccB[idx];
    }
    *csvs[4] << endl;
  }
};

template<>
void RealTime<double>::tarCSVFiles(){
#ifdef UNIX
  this->fileio_->out << "Tarring CSV files for Real-Time Propagation..." 
                     << endl;

  if(!system(NULL)) CErr("Cannot find available processor",this->fileio_->out);
  std::string command = "tar -cf " + this->fileio_->fileName() + ".tar ";
  for(auto i = this->csvFiles.begin(); i != this->csvFiles.end(); i++)
    command += i->second + " ";
  auto info = system(command.c_str());
  command = "rm "; 
  for(auto i = this->csvFiles.begin(); i != this->csvFiles.end(); i++)
    command += i->second + " ";
  info = system(command.c_str());
#endif
};

} //namespace ChronusQ
