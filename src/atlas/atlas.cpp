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

RealMatrix AX(const RealMatrix &A, const RealMatrix &B) {return A*B;};

int ChronusQ::atlas(int argc, char *argv[], GlobalMPI *globalMPI) {
  int i,j,k,l;
  time_t currentTime;
  auto molecule     	= std::unique_ptr<Molecule>(new Molecule());
  auto basisset     	= std::unique_ptr<BasisSet>(new BasisSet());
  auto dfBasisset     	= std::unique_ptr<BasisSet>(new BasisSet());
  auto controls     	= std::unique_ptr<Controls>(new Controls());
  auto aointegrals	= std::unique_ptr<AOIntegrals>(new AOIntegrals());
  auto hartreeFock	= std::unique_ptr<SingleSlater>(new SingleSlater());
  std::unique_ptr<FileIO> fileIO;
  cout << dfBasisset.get() << endl;

  std::vector<std::string> argv_string;
  for(auto i = 1; i < argc; ++i) if(argv[i][0]=='-') argv_string.push_back(argv[i]);
  if(argv_string.size()==0) fileIO = std::unique_ptr<FileIO>(new FileIO(argv[1]));
  else fileIO = std::unique_ptr<FileIO>(new FileIO(argv_string));


  // print out the starting time of the job
  time(&currentTime);
  fileIO->out<<"Job started: "<<ctime(&currentTime)<<endl;
  //fileIO->out<<"Central control process is on "<<globalMPI->nodeName<<endl;
  // read input
  controls->iniControls();
  readInput(fileIO.get(),molecule.get(),basisset.get(),controls.get(),dfBasisset.get());
//  fileIO->iniFileIO(controls->restart);

  // print out molecular and basis set information
  molecule->printInfo(fileIO.get(),controls.get());
  basisset->printInfo_libint(fileIO.get(),controls.get());
  dfBasisset->printInfo_libint(fileIO.get(),controls.get());
  aointegrals->iniAOIntegrals(molecule.get(),basisset.get(),fileIO.get(),controls.get());
  hartreeFock->iniSingleSlater(molecule.get(),basisset.get(),aointegrals.get(),fileIO.get(),controls.get());
  hartreeFock->printInfo();
#ifdef USE_LIBINT
  aointegrals->computeSchwartz();
  if(controls->buildn4eri) aointegrals->computeAOTwoE();
#endif
  if(controls->guess==0) hartreeFock->formGuess();
  else if(controls->guess==1) hartreeFock->readGuessIO();
  else if(controls->guess==2) {
    GauMatEl matEl(controls->gauMatElName);
    hartreeFock->readGuessGauMatEl(matEl);
  }
  else if(controls->guess==3) hartreeFock->readGuessGauFChk(controls->gauFChkName);
//APS I have MO
  fileIO->out<<"APS Print 0"<<endl;
  hartreeFock->matchord();
//APE
  hartreeFock->formFock();
  aointegrals->printTimings();
  hartreeFock->computeEnergy();
  if(controls->optWaveFunction) hartreeFock->SCF();
  else fileIO->out << "**Skipping SCF Optimization**" << endl; 
  hartreeFock->computeMultipole();

/*
  MOIntegrals *moIntegrals = new MOIntegrals();
  moIntegrals->iniMOIntegrals(molecule,basisset,fileIO,controls,aointegrals,hartreeFock);

  SDResponse *sdResponse = new SDResponse();
  sdResponse->iniSDResponse(molecule,basisset,moIntegrals,fileIO,controls,hartreeFock);

  sdResponse->computeExcitedStates();
*/
//APS
  int Iop=0;
  molecule->toCOM(Iop);  // call object molecule pointing to function toCOM-Iop=0 Center of Mass
  Iop=1;
  molecule->toCOM(Iop);  // call object molecule pointing to function toCOM-Iop=1 Center of Nuclear Charges
//APE
  time(&currentTime);
  fileIO->out<<"\nJob finished: "<<ctime(&currentTime)<<endl;
  int N = 500;
  int NSek = 15;
  std::shared_ptr<RealMatrix> A = std::make_shared<RealMatrix>(N,N);
  for(auto i = 0; i < N; i++) (*A)(i,i) = i+1;
  (*A) = (*A) + RealMatrix::Random(N,N);
  (*A) = A->selfadjointView<Eigen::Lower>();
  Eigen::EigenSolver<RealMatrix> ES;
  for(int i = 0; i < N; i++)
  for(int j = 0; j < N; j++) {
    (*A)(i,j) = std::abs((*A)(i,j));
  }
/*
  ES.compute(*A);
  Eigen::VectorXd E = ES.eigenvalues().real();
  RealMatrix VR = ES.eigenvectors().real();
  cout << endl << ES.eigenvectors() << endl;
  ES.compute(A->transpose());
  RealMatrix VL = ES.eigenvectors().real();
  cout << endl << ES.eigenvectors() << endl;
//VR.normCol();
//VL.normCol();

  Eigen::FullPivLU<RealMatrix> lu(VL.transpose()*VR);
  cout << endl << ES.eigenvectors() << endl;

  RealMatrix L = RealMatrix::Identity(N,N);
  L.triangularView<Eigen::StrictlyLower>() = lu.matrixLU();
  RealMatrix U = lu.matrixLU().triangularView<Eigen::Upper>();
  for(auto i = 0; i < N; i++){
    cout << endl << (*A)*VR.col(i) - E(i)*VR.col(i) << endl;
  }

  cout << endl << endl << L.inverse()*lu.permutationP()*VL.transpose()*VR*lu.permutationQ()*U.inverse() << endl;

  VR = VR*lu.permutationQ()*U.inverse();
  VL = VL*lu.permutationP().transpose()*L.inverse().transpose();
  cout << endl << VL.transpose()*VR << endl;


  for(auto i = 0; i < N; i++){
    cout << endl << (*A)*VR.col(i) - E(i)*VR.col(i) << endl;
  }

  cout << endl << VL.transpose()*VR << endl;
  biOrth(VL,VR);
  cout << endl << VL.transpose()*VR << endl;
//cout << *A << endl;
*/
//Davidson<double> dav(&AX,A,NSek,N);
//dav.run(fileIO->out);
//GauMatEl gau("file.out");
//double * x = NULL;
//gau.readRec(GauMatEl::overlap,x);
//cout << "OUTSIDE" << endl;
//cout << x[23] << endl << endl;
//gau.readRec(GauMatEl::dipole,x);
//cout << "OUTSIDE" << endl;
//cout << x[23] << endl << endl;
  

#ifdef USE_LIBINT
  libint2::cleanup();
#endif


return  1;
};


