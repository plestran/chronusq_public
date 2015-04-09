#include "workers.h"
using namespace ChronusQ;

int atlas(int argc, char *argv[], GlobalMPI *globalMPI) {
  int i,j,k,l;
  time_t currentTime;
  Molecule    	*molecule     	= new Molecule();
  BasisSet     	*basisset     	= new BasisSet();
  Controls     	*controls     	= new Controls();
  AOIntegrals	*aointegrals	= new AOIntegrals();
  SingleSlater	*hartreeFock	= new SingleSlater();
  FileIO       	*fileIO;

  try { fileIO=new FileIO(argv[1]);}
  catch(int msg) {
    cout<<"Unable to open file! E#:"<<msg<<endl;
    exit(1);
  };

  // print out the starting time of the job
  time(&currentTime);
  fileIO->out<<"Job started: "<<ctime(&currentTime)<<endl;
  //fileIO->out<<"Central control process is on "<<globalMPI->nodeName<<endl;

  // read input
  controls->iniControls();
  readInput(fileIO,molecule,basisset,controls);
  fileIO->iniFileIO(controls->restart);

  // print out molecular and basis set information
  molecule->printInfo(fileIO,controls);
  basisset->printInfo(fileIO,controls);
  aointegrals->iniAOIntegrals(molecule,basisset,fileIO,controls);
  hartreeFock->iniSingleSlater(molecule,basisset,aointegrals,fileIO,controls);
  hartreeFock->printInfo();
  if(controls->guess==0) hartreeFock->formGuess();
  else if(controls->guess==1) hartreeFock->readGuessIO();
  else if(controls->guess==2) ;
  else if(controls->guess==3) hartreeFock->readGuessGauFChk(controls->gauFChkName);
  hartreeFock->formFock();
  hartreeFock->computeEnergy();
  hartreeFock->SCF();
  MOIntegrals *moIntegrals = new MOIntegrals();
  moIntegrals->iniMOIntegrals(molecule,basisset,fileIO,controls,aointegrals,hartreeFock);

  SDResponse *sdResponse = new SDResponse();
  sdResponse->iniSDResponse(molecule,basisset,moIntegrals,fileIO,controls,hartreeFock);

  sdResponse->computeExcitedStates();

  time(&currentTime);
  fileIO->out<<"\nJob finished: "<<ctime(&currentTime)<<endl;
  delete  molecule;
  delete  basisset;
  delete  fileIO;
  delete  aointegrals;
  delete  controls;
  return  1;
};


