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
#include <workers.h>
#include <davidson.h>
using namespace ChronusQ;

int ChronusQ::atlas(int argc, char *argv[], GlobalMPI *globalMPI) {
  time_t currentTime;

  // Pointers for important storage 
  auto molecule     	= std::unique_ptr<Molecule>(new Molecule());
  auto basisset     	= std::unique_ptr<BasisSet>(new BasisSet());
  auto dfBasisset     	= std::unique_ptr<BasisSet>(new BasisSet());
  auto controls     	= std::unique_ptr<Controls>(new Controls());
  auto aointegrals	= std::unique_ptr<AOIntegrals>(new AOIntegrals());
  auto hartreeFock	= std::unique_ptr<SingleSlater<double>>(new SingleSlater<double>());
  std::unique_ptr<FileIO> fileIO;

  // Initialize the FileIO object
  std::vector<std::string> argv_string;
  for(auto i = 1; i < argc; ++i) if(argv[i][0]=='-') argv_string.push_back(argv[i]);
  if(argv_string.size()==0) fileIO = std::unique_ptr<FileIO>(new FileIO(argv[1]));
  else fileIO = std::unique_ptr<FileIO>(new FileIO(argv_string));


  // print out the starting time of the job
  time(&currentTime);
  fileIO->out<<"Job started: "<<ctime(&currentTime)<<endl;
  //fileIO->out<<"Central control process is on "<<globalMPI->nodeName<<endl;

  // Initialize default settings and read input
  controls->iniControls();
  readInput(fileIO.get(),molecule.get(),basisset.get(),controls.get(),dfBasisset.get());
//  fileIO->iniFileIO(controls->restart);

  // print out molecular and basis set information
  molecule->printInfo(fileIO.get(),controls.get());
  basisset->printInfo_libint(fileIO.get(),controls.get());

  dfBasisset->printInfo_libint(fileIO.get(),controls.get());

  aointegrals->iniAOIntegrals(molecule.get(),basisset.get(),fileIO.get(),controls.get(),dfBasisset.get());
  cout << "HERE" << endl;
  hartreeFock->iniSingleSlater(molecule.get(),basisset.get(),aointegrals.get(),fileIO.get(),controls.get());
  cout << "HERE" << endl;
  hartreeFock->printInfo();
  cout << "HERE" << endl;
  if(controls->guess==0) hartreeFock->formGuess();
  else if(controls->guess==1) hartreeFock->readGuessIO();
  else if(controls->guess==2) {
    GauMatEl matEl(controls->gauMatElName);
    hartreeFock->readGuessGauMatEl(matEl);
  }
  else if(controls->guess==3) hartreeFock->readGuessGauFChk(controls->gauFChkName);
  cout << "HERE" << endl;
  hartreeFock->formFock();
  aointegrals->printTimings();
  hartreeFock->computeEnergy();
  if(controls->optWaveFunction) hartreeFock->SCF();
  else fileIO->out << "**Skipping SCF Optimization**" << endl; 
  hartreeFock->computeMultipole();

//if(controls->doDF) aointegrals->compareRI();
/*
  MOIntegrals *moIntegrals = new MOIntegrals();
  moIntegrals->iniMOIntegrals(molecule,basisset,fileIO,controls,aointegrals,hartreeFock);

  SDResponse *sdResponse = new SDResponse();
  sdResponse->iniSDResponse(molecule,basisset,moIntegrals,fileIO,controls,hartreeFock);

  sdResponse->computeExcitedStates();
*/
  time(&currentTime);
  fileIO->out<<"\nJob finished: "<<ctime(&currentTime)<<endl;
  SingleSlater<dcomplex> newSS(hartreeFock.get());
  newSS.printInfo();
  prettyPrint(cout,*newSS.densityA(),"New D");
#ifdef USE_LIBINT
  libint2::cleanup();
#endif


return  1;
};


