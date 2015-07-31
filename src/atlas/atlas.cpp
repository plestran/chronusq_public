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
#include <grid.h>
using namespace ChronusQ;

  //  Test Function for One-dimensional grid
void func(Grid *g){
  // Reference numeric integration computed with Mathematica
  double ref = 0.18947234582049224;
   (*g).genGrid();
//   (*g).printGrid(); 
   cout << "Test Integral value= "<< (*g).integrate() << endl;
//   cout << "Test Integral err  = "<< std::scientific <<std::abs((*g).integrate() - ref)/ref << endl;
   cout << "Test Integral err  = "<<std::abs((*g).integrate() - ref)/ref << endl;
}

void func2D(OneDGrid *g){
   (*g).genGrid();
   (*g).printGrid(); 
   cout << "Test Integral value= "<< (*g).integrate() << endl;
}

//void func2D_plus(OneDGrid *gr, OneDGrid *gs){
//   (*gr).genGrid();
//   (*gs).genGrid();
//     TwoDGrid(*gr, *gs) grs2;
//   TwoDGrid *G3(*gr,*gs);
//   cout << "Test Integral value= "<< (g2d).integrate() << endl;
//}

int ChronusQ::atlas(int argc, char *argv[], GlobalMPI *globalMPI) {
  time_t currentTime;

  // Pointers for important storage 
  auto molecule     	= std::unique_ptr<Molecule>(new Molecule());
  auto basisset     	= std::unique_ptr<BasisSet>(new BasisSet());
  auto dfBasisset     	= std::unique_ptr<BasisSet>(new BasisSet());
  auto controls     	= std::unique_ptr<Controls>(new Controls());
  auto aointegrals	= std::unique_ptr<AOIntegrals>(new AOIntegrals());
  auto hartreeFock	= std::unique_ptr<SingleSlater<double>>(new SingleSlater<double>());
  auto sdResponse       = std::unique_ptr<SDResponse>(new SDResponse());
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
  //fileIO->out<<"Central control process is on "<<globalMPI->nodeName<<endl;

  // Initialize default settings and read input
  controls->iniControls();
  readInput(fileIO.get(),molecule.get(),basisset.get(),controls.get(),dfBasisset.get());
//  fileIO->iniFileIO(controls->restart);

  // print out molecular and basis set information
  controls->printSettings(fileIO->out);
  molecule->printInfo(fileIO.get(),controls.get());
  basisset->printInfo_libint(fileIO.get(),controls.get());

  dfBasisset->printInfo_libint(fileIO.get(),controls.get());

  aointegrals->iniAOIntegrals(molecule.get(),basisset.get(),fileIO.get(),controls.get(),dfBasisset.get());
  hartreeFock->iniSingleSlater(molecule.get(),basisset.get(),aointegrals.get(),fileIO.get(),controls.get());
  hartreeFock->printInfo();
  if(controls->guess==0) hartreeFock->formGuess();
  else if(controls->guess==1) hartreeFock->readGuessIO();
  else if(controls->guess==2) {
    GauMatEl matEl(controls->gauMatElName);
    hartreeFock->readGuessGauMatEl(matEl);
  }
  else if(controls->guess==3) hartreeFock->readGuessGauFChk(controls->gauFChkName);
//APS I have MO Please check in which controls call the following function
//hartreeFock->matchord();
//APE
  hartreeFock->formFock();
  aointegrals->printTimings();
  hartreeFock->computeEnergy();
  std::shared_ptr<MOIntegrals> moIntegrals = std::make_shared<MOIntegrals>();
  if(controls->optWaveFunction) {
    hartreeFock->doCUHF = false;
    hartreeFock->SCF();
  }
  //MOIntegrals *moIntegrals = new MOIntegrals();
  //moIntegrals->iniMOIntegrals(molecule,basisset,fileIO,controls,aointegrals,hartreeFock);
  else fileIO->out << "**Skipping SCF Optimization**" << endl; 
  hartreeFock->computeMultipole();
  if(controls->doSDR) {
    sdResponse->iniSDResponse(molecule.get(),basisset.get(),moIntegrals.get(),fileIO.get(),
                              controls.get(),hartreeFock.get());
    sdResponse->IterativeRPA();
  }

//if(controls->doDF) aointegrals->compareRI();
/*
  MOIntegrals *moIntegrals = new MOIntegrals();
  moIntegrals->iniMOIntegrals(molecule,basisset,fileIO,controls,aointegrals,hartreeFock);


 int Iop=0;
 molecule->toCOM(Iop);  // call object molecule pointing to function toCOM-Iop=0 Center of Mass
 Iop=1;
 molecule->toCOM(Iop);  // call object molecule pointing to function toCOM-Iop=1 Center of Nuclear Charges
*/
//Test one dimensional grid
  fileIO->out << "**AP One dimensional grid test**" << endl;
  int Ngridr =   500;
  int NLeb    = 38;
  double radius = 1.0;
   GaussChebyshev1stGrid Rad(Ngridr,0.0,radius);
   LebedevGrid GridLeb(NLeb);
   
//   func(&Rad);
//   func2D(&GridLeb);
//
   Rad.genGrid();
   GridLeb.genGrid();
   TwoDGrid G3(&Rad,&GridLeb);
   cout << "Sphere Int = " << G3.integrate() <<endl;
   cout << "Sphere Err = " << std::abs(G3.integrate()-(4.0*radius*radius*radius*math.pi/3.0)) <<endl;
//   G3.printGrid();
  time(&currentTime);
  fileIO->out<<"\nJob finished: "<<ctime(&currentTime)<<endl;
/*
  SingleSlater<dcomplex> newSS(hartreeFock.get());
  newSS.formFock();
*/
/*
  double *tmp = new double[3*2];
  for(auto i =0; i < 6; i++) tmp[i] = 0.0;
  tmp[0] = 0.9;
  std::vector<int> atm;
  atm.push_back(1);
  atm.push_back(1);
  GauJob job(false,"sto-3g",tmp,atm,0,1);
  job.run();
*/
/*
  Eigen::SelfAdjointEigenSolver<RealMatrix> ES;
  ES.compute((*hartreeFock->densityA())+(*hartreeFock->densityB())/2);
  cout << ES.eigenvalues() << endl;
  cout << endl <<ES.eigenvectors()*(*hartreeFock->densityA())*ES.eigenvectors().transpose() << endl;
  cout << endl << (hartreeFock->moB()->transpose())*(*aointegrals->overlap_)*(*hartreeFock->densityB())*(*aointegrals->overlap_)*(*hartreeFock->moB()) << endl;
  RealMatrix X = aointegrals->overlap_->pow(-0.5);
  ES.compute(X*((*hartreeFock->densityA())+(*hartreeFock->densityB()))*X.transpose()/2);
  cout << endl << ES.eigenvalues() << endl;
  cout << endl << ES.eigenvectors().transpose()*(*hartreeFock->densityA())*ES.eigenvectors() << endl;
*/
#ifdef USE_LIBINT
  libint2::cleanup();
#endif


return  1;
};


