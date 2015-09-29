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
#include <mollerplesset.h>
#include <grid.h>
using namespace ChronusQ;

int ChronusQ::atlas(int argc, char *argv[], GlobalMPI *globalMPI) {
  time_t currentTime;

  // Pointers for important storage and job control
  auto molecule           = std::unique_ptr<Molecule>(new Molecule());
  auto basisset           = std::unique_ptr<BasisSet>(new BasisSet());
  auto dfBasisset         = std::unique_ptr<BasisSet>(new BasisSet());
  auto controls           = std::unique_ptr<Controls>(new Controls());
  auto aointegrals        = std::unique_ptr<AOIntegrals>(new AOIntegrals());
  auto mointegralsReal    = std::unique_ptr<MOIntegrals<double>>(new MOIntegrals<double>());
  auto mointegralsComplex = std::unique_ptr<MOIntegrals<dcomplex>>(new MOIntegrals<dcomplex>());
  auto hartreeFockReal	  = std::unique_ptr<SingleSlater<double>>(  new SingleSlater<double>()  );
  auto hartreeFockComplex = std::unique_ptr<SingleSlater<dcomplex>>(new SingleSlater<dcomplex>());
  auto realtime		  = std::unique_ptr<RealTime>(new RealTime());
  auto sdResponseReal     = std::unique_ptr<SDResponse<double>>(new SDResponse<double>());
  auto sdResponseComplex  = std::unique_ptr<SDResponse<dcomplex>>(new SDResponse<dcomplex>());
//auto twoDGrid     	= std::unique_ptr<TwoDGrid>(new TwoDGrid());
  std::unique_ptr<FileIO> fileIO;
  std::unique_ptr<GauJob> gauJob;

  // Initialize the FileIO object
  std::vector<std::string> argv_string;
  for(auto i = 1; i < argc; ++i) if(argv[i][0]=='-') argv_string.push_back(argv[i]);
  if(argv_string.size()==0) fileIO = std::unique_ptr<FileIO>(new FileIO(argv[1]));
  else fileIO = std::unique_ptr<FileIO>(new FileIO(argv_string));


  // print out the starting time of the job
  time(&currentTime);
  fileIO->out<<"Job started: "<<ctime(&currentTime)<<endl;

  // Initialize default settings and read input
  controls->iniControls();
  readInput(fileIO.get(),molecule.get(),basisset.get(),controls.get(),dfBasisset.get());

  // print out molecular and basis set information
  controls->printSettings(fileIO->out);
  molecule->printInfo(fileIO.get(),controls.get());
  basisset->printInfo();


  // Initialize memory for AO integral storage
  aointegrals->iniAOIntegrals(molecule.get(),basisset.get(),fileIO.get(),controls.get(),
    dfBasisset.get());

  // Initialize memory for wave function related quantities
  // Logic for Real / Complex
  if(!controls->doComplex){
    hartreeFockReal->iniSingleSlater(molecule.get(),basisset.get(),aointegrals.get(),
      fileIO.get(),controls.get());
    hartreeFockReal->printInfo();
  } else {
    hartreeFockComplex->iniSingleSlater(molecule.get(),basisset.get(),aointegrals.get(),
      fileIO.get(),controls.get());
    hartreeFockComplex->printInfo();
  }

  // Form the initial guess for the wave function. Either form a new
  // guess from scratch or read off input / gaussian
  //
  // ** Note that guess from Gaussian is buggy **
  if(!controls->doComplex){
    if(controls->guess==0) hartreeFockReal->formGuess();
    else if(controls->guess==1) hartreeFockReal->readGuessIO();
    else if(controls->guess==2) {
      GauMatEl matEl(controls->gauMatElName);
      hartreeFockReal->readGuessGauMatEl(matEl);
    }
    else if(controls->guess==3) hartreeFockReal->readGuessGauFChk(controls->gauFChkName);
  } else {
    if(controls->guess==0) hartreeFockComplex->formGuess();
    else CErr("Cannot Read Guess for Complex Wavefunctions (NYI)",fileIO->out);
  }

  // Optimize wave function (?)
  if(!controls->doComplex){
    // Form initial (primer) Fock matrix
    hartreeFockReal->formFock();
    aointegrals->printTimings();

    // Compute initial energy
    hartreeFockReal->computeEnergy();

    // Optionally optimize the wavefunction through SCF (unless told to skip)
    if(controls->optWaveFunction)  hartreeFockReal->SCF();
    else fileIO->out << "**Skipping SCF Optimization**" << endl; 

    // Compute the Electric Multipole Moments
    hartreeFockReal->computeMultipole();
  } else {
    // Form initial (primer) Fock matrix
    hartreeFockComplex->formFock();
    aointegrals->printTimings();

    // Compute initial energy
    hartreeFockComplex->computeEnergy();

    // Optionally optimize the wavefunction through SCF (unless told to skip)
    if(controls->optWaveFunction)  hartreeFockComplex->SCF();
    else fileIO->out << "**Skipping SCF Optimization**" << endl; 

    // Compute the Electric Multipole Moments
    hartreeFockComplex->computeMultipole();
  }

  // Run PSCF Calculations
  if(controls->doSDR) {
    if(!controls->doComplex){
      mointegralsReal->iniMOIntegrals(molecule.get(),basisset.get(),fileIO.get(),controls.get(),aointegrals.get(),hartreeFockReal.get());
      sdResponseReal->setPPRPA(1);
      sdResponseReal->iniSDResponse(molecule.get(),basisset.get(),mointegralsReal.get(),fileIO.get(),
                                controls.get(),hartreeFockReal.get());
      
      sdResponseReal->IterativeRPA();
    } else {
      mointegralsComplex->iniMOIntegrals(molecule.get(),basisset.get(),fileIO.get(),controls.get(),aointegrals.get(),hartreeFockComplex.get());
      sdResponseComplex->setPPRPA(1);
      sdResponseComplex->iniSDResponse(molecule.get(),basisset.get(),mointegralsComplex.get(),fileIO.get(),
                                controls.get(),hartreeFockComplex.get());
      
      sdResponseComplex->IterativeRPA();
    }
  }
  if(controls->doUnit) printUnitInfo(controls.get(),hartreeFockReal.get(),sdResponseReal.get());

/*
////// APS ////
  twoDGrid->iniTwoDGrid(fileIO.get(),molecule.get(),basisset.get(),aointegrals.get(),hartreeFockReal.get(),100,194);
////// APE ////
*/
/*
//fds
  realtime->iniRealTime(molecule.get(),basisset.get(),fileIO.get(),controls.get(),aointegrals.get(),hartreeFockReal.get());
  fileIO->out<<"\niniRealTime Done: "<<ctime(&currentTime)<<endl;
  realtime->iniDensity();
  fileIO->out<<"\niniDensity Done: "<<ctime(&currentTime)<<endl;
  realtime->doPropagation();
  fileIO->out<<"\ndoPropagation Done: "<<ctime(&currentTime)<<endl;
//fds
*/

  // Cleanup Libint env
#ifdef USE_LIBINT
  libint2::cleanup();
#endif


  return  1;
};


