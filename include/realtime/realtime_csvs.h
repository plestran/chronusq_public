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
template <typename T>
void RealTime<T>::initCSV() {
  // Create/open CSVs for printing results
  csvs_.push_back(
    new std::ofstream(fileio_->fileName() + "_RealTime_Dipole.csv")
  );
  csvFiles_[csvs_[0]] = fileio_->fileName() + "_RealTime_Dipole.csv";

  csvs_.push_back(
    new std::ofstream(fileio_->fileName() + "_RealTime_AppliedField.csv")
  );
  csvFiles_[csvs_[1]] = fileio_->fileName() + "_RealTime_AppliedField.csv";

  csvs_.push_back(
    new std::ofstream(fileio_->fileName() + "_RealTime_Mulliken.csv")
  );
  csvFiles_[csvs_[2]] = fileio_->fileName() + "_RealTime_Mulliken.csv";

  csvs_.push_back(
    new std::ofstream(fileio_->fileName() + "_RealTime_OrbOcc_Alpha.csv")
  );
  csvFiles_[csvs_[3]] = fileio_->fileName() + "_RealTime_OrbOcc_Alpha.csv";

  if(!ssPropagator_->isClosedShell && ssPropagator_->nTCS() != 2){
    csvs_.push_back(
      new std::ofstream(fileio_->fileName() + "_RealTime_OrbOcc_Beta.csv")
    );
    csvFiles_[csvs_[4]] = fileio_->fileName() + "_RealTime_OrbOcc_Beta.csv";
  }


  // Print CSV Headers
  // DIPOLE
  *csvs_[0] << "Time Step (a.u.), Energy (Eh), Dipole X debye),"
           << "Dipole Y (debye), Dipole Z (debye), Dipole Tot (debye)" 
           << endl;

  // APPLIED FIELD
  *csvs_[1] << "Time Step (a.u.), Ex (a.u), Ey (a.u.), Ez (a.u.),"
           << "|| E total || (a.u.)" << endl;

  // MULLIKEN
  *csvs_[2] << std::setw(14) << "Atom number";
  for(auto iAtm = 0; iAtm < ssPropagator_->molecule()->nAtoms(); iAtm++) 
    *csvs_[2] << std::setw(14) << iAtm;
  *csvs_[2] << endl;
  *csvs_[2] << std::setw(14) << "Atom symbol";
  for(auto iAtm = 0; iAtm < ssPropagator_->molecule()->nAtoms(); iAtm++) 
    *csvs_[2] << std::setw(14) 
             << elements[ssPropagator_->molecule()->index(iAtm)].symbol;
  *csvs_[2] << endl;
  *csvs_[2] << std::setw(14) << "Time (a.u.)" << endl;

  // ORBITAL OCCUPATION
  auto NBT = ssPropagator_->basisset()->nBasis() * ssPropagator_->nTCS();
  *csvs_[3] << std::setw(10) << "Time (au)";
  for(auto idx = 0; idx != NBT; idx++)
    *csvs_[3] << ", " << std::setw(10) << idx + 1;
  *csvs_[3] << endl;
  if(!ssPropagator_->isClosedShell && ssPropagator_->nTCS() != 2){
    *csvs_[4] << std::setw(10) << "Time (au)";
    for(auto idx = 0; idx != NBT; idx++)
      *csvs_[4] << ", " << std::setw(10) << idx + 1;
    *csvs_[4] << endl;
  }
};


template <typename T>
void RealTime<T>::writeDipoleCSV(){
    *csvs_[0] << std::fixed << std::setprecision(10) 
             << propInfo.back().timeStep << ", " 
             << propInfo.back().energy << ", " 
             << propInfo.back().dipole[0] << ", " 
             << propInfo.back().dipole[1] << ", " 
             << propInfo.back().dipole[2] << ", " 
             << propInfo.back().dipole[3] << endl;
};

template <typename T>
void RealTime<T>::writeAppliedFieldCSV(){
    *csvs_[1] << std::fixed << std::setprecision(10) 
             << propInfo.back().timeStep << ", " 
             << propInfo.back().appliedfield[0] << ", " 
             << propInfo.back().appliedfield[1] << ", " 
             << propInfo.back().appliedfield[2] << ", " 
             << propInfo.back().appliedfield[3] << endl;
};

template <typename T>
void RealTime<T>::writeMullikenCSV(){
    *csvs_[2] << std::fixed << std::setw(14) << std::setprecision(10) 
             << propInfo.back().timeStep;
    for(auto iAtm = 0; iAtm < ssPropagator_->molecule()->nAtoms(); iAtm++) 
      *csvs_[2] << ", " << std::setw(14) << propInfo.back().mullPop[iAtm];

    *csvs_[2] << endl;
};

template <typename T>
void RealTime<T>::writeOrbitalCSV(){
  auto NBT = ssPropagator_->basisset()->nBasis() * ssPropagator_->nTCS();
  *csvs_[3] << std::fixed << std::setw(10) << std::setprecision(6) 
           << propInfo.back().timeStep; 
  for(auto idx = 0; idx != NBT; idx++) 
    *csvs_[3] << ", " << std::setw(10) << std::setprecision(6) 
             << propInfo.back().orbitalOccA[idx];
  *csvs_[3] << endl;

  if(!ssPropagator_->isClosedShell && ssPropagator_->nTCS() != 2){
    *csvs_[4] << std::fixed << std::setw(10) << std::setprecision(6) 
             << propInfo.back().timeStep; 
    for(auto idx = 0; idx != NBT; idx++) 
      *csvs_[4] << ", " << std::setw(10) << std::setprecision(6) 
               << propInfo.back().orbitalOccB[idx];
    *csvs_[4] << endl;
  }
};

template <typename T>
void RealTime<T>::tarCSVFiles(){
#ifdef UNIX
  fileio_->out << "Tarring CSV files for Real-Time Propagation..." << endl;
  if(!system(NULL)) CErr("Cannot find available processor",fileio_->out);
  std::string command = "tar -cf " + fileio_->fileName() + ".tar ";
  for(auto i = csvFiles_.begin(); i != csvFiles_.end(); i++)
    command += i->second + " ";
  auto info = system(command.c_str());
  command = "rm "; 
  for(auto i = csvFiles_.begin(); i != csvFiles_.end(); i++)
    command += i->second + " ";
  info = system(command.c_str());
#endif
};

