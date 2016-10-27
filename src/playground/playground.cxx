//#define EIGEN_RUNTIME_NO_MALLOC
//#include <response.h>
//#include <workers.h>
//#include <pythonapi.h>
#include <global.h>
#include <pythonapi.h>
#include <singleslater.h>
#include <mointegrals.h>
#include <response.h>
#include <realtime.h>

struct MyStruct {
  RealMatrix VXCA;
  RealMatrix VXCB;
  double Energy;

  MyStruct(size_t N) : VXCA(N,N), VXCB(N,N), Energy(0.0){ 
    VXCA.setZero(); VXCB.setZero();
  };
};
#include <grid2.h>

using namespace ChronusQ;

enum MOLECULE_PRESETS {
  WATER,Methanol,H,HE,SO,OxMolecule,Li
};

template<MOLECULE_PRESETS T>
void loadPresets(Molecule&);

template<>
void loadPresets<WATER>(Molecule &mol) {
  mol.setNAtoms(3);
  mol.setCharge(0);
  mol.setNTotalE(10);
  mol.setMultip(1);
  mol.alloc();
  mol.setIndex(0,HashAtom("O",0));
  mol.setIndex(1,HashAtom("H",0));
  mol.setIndex(2,HashAtom("H",0));
  mol.setCart(0,0.000000000 ,-0.07579184359, 0.0);
  mol.setCart(1,0.866811829 ,0.6014357793  ,0.0);
  mol.setCart(2,-0.866811829, 0.6014357793 ,0.0);
};
template<>
void loadPresets<Methanol>(Molecule &mol) {
  mol.setCharge(0);
  mol.setNTotalE(18);
  mol.setMultip(1);
  mol.setNAtoms(6);
  mol.alloc(); //allocates all memory for the class

  mol.setIndex(0,HashAtom("C",0));
  mol.setIndex(1,HashAtom("H",0));
  mol.setIndex(2,HashAtom("O",0));
  mol.setIndex(3,HashAtom("H",0));
  mol.setIndex(4,HashAtom("H",0));
  mol.setIndex(5,HashAtom("H",0));

  mol.setCart(0,-1.013487,1.725956,1.257405);
  mol.setCart(1,-0.069872,1.679306,0.755096);
  mol.setCart(2,-1.488002,3.074931,1.256099);
  mol.setCart(3,-1.168237,3.527873,2.039804);
  mol.setCart(4,-1.718167,1.099499,0.751562);
  mol.setCart(5,-0.897365,1.389691,2.266534);
};
template<>
void loadPresets<H>(Molecule &mol){
  mol.setNAtoms(1);
  mol.setCharge(0);
  mol.setNTotalE(1);
  mol.setMultip(2);
  mol.alloc();
  mol.setIndex(0,HashAtom("H",0));
  mol.setCart(0,0.000000000 ,0.00000000000, 0.0);
};
template<>
void loadPresets<HE>(Molecule &mol){
  mol.setNAtoms(1);
  mol.setCharge(0);
  mol.setNTotalE(2);
  mol.setMultip(1);
  mol.alloc();
  mol.setIndex(0,HashAtom("He",0));
//  mol.setCart(0,0.000000000 ,-0.00000000000, 0.0);
  mol.setCart(0,0.000000000 ,0.00000000000, 0.0);
};
template<>
void loadPresets<SO>(Molecule &mol){
  mol.setNAtoms(1);
  mol.setCharge(0);
  mol.setNTotalE(8);
  mol.setMultip(1);
  mol.alloc();
  mol.setIndex(0,HashAtom("O",0));
//  mol.setCart(0,0.000000000 ,-0.00000000000, 0.0);
  mol.setCart(0,0.000000000 ,0.00000000000, 0.0);
};
template<>
void loadPresets<OxMolecule>(Molecule &mol){
  mol.setCharge(0);
  mol.setNTotalE(16);
  mol.setMultip(3);
  mol.setNAtoms(2);
  mol.alloc(); //allocates all memory for the class
  mol.setIndex(0,HashAtom("O",0));
  mol.setIndex(1,HashAtom("O",0));
  mol.setCart(0,-0.608586,0.0,0.0); //In Angstroms!
  mol.setCart(1,0.608586,0.0,0.0);
};
template<>
void loadPresets<Li>(Molecule &mol){
  mol.setNAtoms(1);
  mol.setCharge(0);
  mol.setNTotalE(3);
  mol.setMultip(2);
  mol.alloc();
  mol.setIndex(0,HashAtom("Li",0));
  mol.setCart(0,0.000000000 ,0.00000000000, 0.0);
};


int main(int argc, char **argv){

  CQMemManager memManager;
  Molecule moleculeLi;
  Molecule moleculeOxy;
  Molecule moleculeWat;
  BasisSet basisOxy;
  BasisSet basisLi;
  BasisSet basisWat;
  AOIntegrals aointsOxy;
  AOIntegrals aointsLi;
  AOIntegrals aointsWat;
  SingleSlater<double> singleSlater109;
  SingleSlater<double> singleSlater112;
  SingleSlater<double> singleSlater115;
  SingleSlater<double> singleSlater118;
  SingleSlater<double> singleSlater121;
  SingleSlater<double> singleSlater124;
  SingleSlater<double> singleSlater127;
  SingleSlater<double> singleSlater130;
  SingleSlater<double> singleSlater133;
  SingleSlater<double> singleSlater137;
  SingleSlater<double> singleSlater136;
  SingleSlater<double> singleSlater150;
  SingleSlater<double> singleSlater151;
  SingleSlater<double> singleSlater160;
  SingleSlater<double> singleSlater161;
  SingleSlater<double> singleSlater162;
  SingleSlater<double> singleSlater163;
  FileIO fileio("test.inp","test.out");

  memManager.setTotalMem(256e6);
  initCQ(argc,argv);
  CQSetNumThreads(1);
 
//////////////////////////////////////////////////////
  loadPresets<OxMolecule>(moleculeOxy);
  loadPresets<Li>(moleculeLi);
  loadPresets<WATER>(moleculeWat);
// Molecule
  moleculeWat.convBohr();
  moleculeWat.computeNucRep();
  moleculeWat.computeRij();
  moleculeWat.computeI();
  cout << "Mol Wat " <<endl;
  moleculeLi.convBohr();
  moleculeLi.computeNucRep();
  moleculeLi.computeRij();
  moleculeLi.computeI();
  cout << "Mol Li " <<endl;

  moleculeOxy.convBohr();
  moleculeOxy.computeNucRep();
  moleculeOxy.computeRij();
  moleculeOxy.computeI();

  cout << "Mol Oxy " <<endl;
  fileio.iniH5Files();

//  basisWat.findBasisFile("cc-pvtz");
  basisWat.findBasisFile("sto-3g");
  basisWat.communicate(fileio);
//  basisWat.forceCart();
  basisWat.parseGlobal();
  basisWat.constructLocal(&moleculeWat);
  basisWat.makeMaps(&moleculeWat);
  basisWat.renormShells();
  cout << "Bas Wat " <<endl;
  aointsWat.communicate(moleculeWat,basisWat,fileio,memManager);
  aointsWat.initMeta();
  aointsWat.alloc();
  cout << "AoInt Wat " <<endl;

  basisLi.findBasisFile("sto-3g");
  basisLi.communicate(fileio);
//  basisLi.forceCart();
  basisLi.parseGlobal();
  basisLi.constructLocal(&moleculeLi);
  basisLi.makeMaps(&moleculeLi);
  basisLi.renormShells();
  cout << "Bas Li " <<endl;
  aointsLi.communicate(moleculeLi,basisLi,fileio,memManager);
  aointsLi.initMeta();
  aointsLi.alloc();
  cout << "AoInt Li " <<endl;

//  basisOxy.findBasisFile("cc-pvtz");
  basisOxy.findBasisFile("sto-3g");
  basisOxy.communicate(fileio);
//  basisOxy.forceCart();
  basisOxy.parseGlobal();
  basisOxy.constructLocal(&moleculeOxy);
  basisOxy.makeMaps(&moleculeOxy);
  basisOxy.renormShells();
  cout << "Bas Oxy " <<endl;
  aointsOxy.communicate(moleculeOxy,basisOxy,fileio,memManager);
  aointsOxy.initMeta();
  aointsOxy.alloc();
  cout << "AoInt Oxy " <<endl;

/*
  singleSlater109.setRef("RSLATER");
//  singleSlater109.isDFT = true;
//  singleSlater109.isHF  = false;
  singleSlater109.setSCFEneTol(1e-12);
  singleSlater109.setSCFMaxIter(10000);
  singleSlater109.doDIIS = true;
  singleSlater109.setGuess(CORE);
//  singleSlater109.addSlater();
//  singleSlater109.setxHF(0.0);
  singleSlater109.communicate(moleculeWat,basisWat,aointsWat,fileio,memManager);
  singleSlater109.initMeta();
//  singleSlater109.genMethString();
  singleSlater109.alloc();
  singleSlater109.formGuess();
  singleSlater109.computeProperties();
  singleSlater109.printProperties();
  singleSlater109.SCF3();
  singleSlater109.computeProperties();
  singleSlater109.printProperties();
  singleSlater109.mullikenPop();

  cout << "Energy " << singleSlater109.totalEnergy() <<endl;
  cout << "PASSED 109? " << bool(std::abs(singleSlater109.totalEnergy() - (-74.0665552706)) <= 1e-8) << endl;
//  singleSlater121.setRef("UBLYP");
  singleSlater121.setRef("USLATER");
//  singleSlater121.isDFT = true;
//  singleSlater121.isHF  = false;
  singleSlater121.setSCFEneTol(1e-12);
  singleSlater121.setSCFMaxIter(10000);
  singleSlater121.doDIIS = true;
  singleSlater121.setGuess(CORE);
//  singleSlater121.addSlater();
//  singleSlater121.setxHF(0.0);
  singleSlater121.communicate(moleculeOxy,basisOxy,aointsOxy,fileio,memManager);
  singleSlater121.initMeta();
//  singleSlater121.genMethString();
  singleSlater121.alloc();
  singleSlater121.formGuess();
  singleSlater121.computeProperties();
  singleSlater121.printProperties();
  singleSlater121.SCF3();
  singleSlater121.computeProperties();
  singleSlater121.printProperties();
  singleSlater121.mullikenPop();

  cout << "Energy " << singleSlater121.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater121.totalEnergy() - (-146.0745567304)) <<endl;
  cout << "PASSED 121 - Oxy Slater STO-3G? " << bool(std::abs(singleSlater121.totalEnergy() - (-146.0745567304)) <= 1e-8) << endl;
///
  singleSlater127.setRef("ULSDA");
//  singleSlater127.isDFT = true;
//  singleSlater127.isHF  = false;
  singleSlater127.setSCFEneTol(1e-12);
  singleSlater127.setSCFMaxIter(10000);
  singleSlater127.doDIIS = true;
  singleSlater127.setGuess(CORE);
 // singleSlater127.createLSDA();
//  singleSlater127.setxHF(0.0);
  singleSlater127.communicate(moleculeOxy,basisOxy,aointsOxy,fileio,memManager);
  singleSlater127.initMeta();
//  singleSlater127.genMethString();
  singleSlater127.alloc();
  singleSlater127.formGuess();
  singleSlater127.computeProperties();
  singleSlater127.printProperties();
  singleSlater127.SCF3();
  singleSlater127.computeProperties();
  singleSlater127.printProperties();
  singleSlater127.mullikenPop();

  cout << "Energy " << singleSlater127.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater127.totalEnergy() - (-147.5135588885)) <<endl;
  cout << "PASSED 127 - Oxy LSDA STO-3G? " << bool(std::abs(singleSlater127.totalEnergy() - (-147.5135588885)) <= 1e-8) << endl;
///
  singleSlater133.setRef("USVWN5");
//  singleSlater133.isDFT = true;
//  singleSlater133.isHF  = false;
  singleSlater133.setSCFEneTol(1e-12);
  singleSlater133.setSCFMaxIter(10000);
  singleSlater133.doDIIS = true;
  singleSlater133.setGuess(CORE);
//  singleSlater133.addSlater();
//  singleSlater133.addVWN5();
//  singleSlater133.setxHF(0.0);
  singleSlater133.communicate(moleculeOxy,basisOxy,aointsOxy,fileio,memManager);
  singleSlater133.initMeta();
//  singleSlater133.genMethString();
  singleSlater133.alloc();
  singleSlater133.formGuess();
  singleSlater133.computeProperties();
  singleSlater133.printProperties();
  singleSlater133.SCF3();
  singleSlater133.computeProperties();
  singleSlater133.printProperties();
  singleSlater133.mullikenPop();

  cout << "Energy " << singleSlater133.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater133.totalEnergy() - (-147.1965659202)) <<endl;
  cout << "PASSED 133 - Oxy LSDA STO-3G? " << bool(std::abs(singleSlater133.totalEnergy() - (-147.1965659202)) <= 1e-8) << endl;
///
  singleSlater137.setRef("UB88");
//  singleSlater137.isDFT = true;
//  singleSlater137.isHF  = false;
  singleSlater137.setSCFEneTol(1e-12);
  singleSlater137.setSCFMaxIter(10000);
  singleSlater137.doDIIS = true;
  singleSlater137.setGuess(CORE);
//  singleSlater137.createB88();
//  singleSlater137.setxHF(0.0);
  singleSlater137.communicate(moleculeOxy,basisOxy,aointsOxy,fileio,memManager);
  singleSlater137.initMeta();
//  singleSlater137.genMethString();
  singleSlater137.alloc();
  singleSlater137.formGuess();
  singleSlater137.computeProperties();
  singleSlater137.printProperties();
  singleSlater137.SCF3();
  singleSlater137.computeProperties();
  singleSlater137.printProperties();
  singleSlater137.mullikenPop();

  cout << "Energy " << singleSlater137.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater137.totalEnergy() - (-147.6742271106)) <<endl;
  cout << "PASSED 137 - Oxy B88 STO-3G? " << bool(std::abs(singleSlater137.totalEnergy() - (-147.6742271106)) <= 1e-8) << endl;

  singleSlater151.setRef("UBLYP");
//  singleSlater151.isDFT = true;
//  singleSlater151.isHF  = false;
  singleSlater151.setSCFEneTol(1e-12);
  singleSlater151.setSCFMaxIter(10000);
  singleSlater151.doDIIS = true;
  singleSlater151.setGuess(CORE);
//  singleSlater151.addSlater();
//  singleSlater151.addB88();
//  singleSlater151.addLYP();
//  singleSlater151.setxHF(0.0);
  singleSlater151.communicate(moleculeOxy,basisOxy,aointsOxy,fileio,memManager);
  singleSlater151.initMeta();
//  singleSlater151.genMethString();
  singleSlater151.alloc();
  singleSlater151.formGuess();
  singleSlater151.computeProperties();
  singleSlater151.printProperties();
  singleSlater151.SCF3();
  singleSlater151.computeProperties();
  singleSlater151.printProperties();
  singleSlater151.mullikenPop();

  cout << "Energy " << singleSlater151.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater151.totalEnergy() - (-148.2438136080)) <<endl;
  cout << "PASSED 151 - Oxy BLYP STO-3G? " << bool(std::abs(singleSlater151.totalEnergy() - (-148.2438136080)) <= 1e-8) << endl;

  singleSlater161.setRef("UB3LYP");
//  singleSlater161.isDFT = true;
//  singleSlater161.isHF  = false;
  singleSlater161.setSCFEneTol(1e-12);
  singleSlater161.setSCFMaxIter(10000);
  singleSlater161.doDIIS = true;
  singleSlater161.setGuess(CORE);
//  singleSlater161.createB3LYP();
//  singleSlater161.setxHF(0.0);
  singleSlater161.communicate(moleculeOxy,basisOxy,aointsOxy,fileio,memManager);
  singleSlater161.initMeta();
//  singleSlater161.genMethString();
  singleSlater161.alloc();
  singleSlater161.formGuess();
  singleSlater161.computeProperties();
  singleSlater161.printProperties();
  singleSlater161.SCF3();
  singleSlater161.computeProperties();
  singleSlater161.printProperties();
  singleSlater161.mullikenPop();

  cout << "Energy " << singleSlater161.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater161.totalEnergy() - (-148.2725627255)) <<endl;
  cout << "PASSED 161 - Oxy B3LYP STO-3G? " << bool(std::abs(singleSlater151.totalEnergy() - (-148.2725627255)) <= 1e-8) << endl;
 
  singleSlater162.setRef("UBHandH");
//  singleSlater162.isDFT = true;
//  singleSlater162.isHF  = false;
  singleSlater162.setSCFEneTol(1e-12);
  singleSlater162.setSCFMaxIter(10000);
  singleSlater162.doDIIS = true;
  singleSlater162.setGuess(CORE);
//  singleSlater162.createBHandH();
//  singleSlater162.setxHF(0.0);
  singleSlater162.communicate(moleculeOxy,basisOxy,aointsOxy,fileio,memManager);
  singleSlater162.initMeta();
//  singleSlater162.genMethString();
  singleSlater162.alloc();
  singleSlater162.formGuess();
  singleSlater162.computeProperties();
  singleSlater162.printProperties();
  singleSlater162.SCF3();
  singleSlater162.computeProperties();
  singleSlater162.printProperties();
  singleSlater162.mullikenPop();

  cout << "Energy " << singleSlater162.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater162.totalEnergy() - (-148.2725627255)) <<endl;
  cout << "PASSED 162 - Oxy BHandH STO-3G? " << bool(std::abs(singleSlater162.totalEnergy() - (-148.2725627255)) <= 1e-8) << endl;

// 112
  singleSlater112.setRef("RLSDA");
//  singleSlater112.isDFT = true;
//  singleSlater112.isHF  = false;
  singleSlater112.setSCFEneTol(1e-12);
  singleSlater112.setSCFMaxIter(10000);
  singleSlater112.doDIIS = true;
  singleSlater112.setGuess(CORE);
//  singleSlater112.createLSDA();
  singleSlater112.communicate(moleculeWat,basisWat,aointsWat,fileio,memManager);
  singleSlater112.initMeta();
//  singleSlater112.genMethString();
  singleSlater112.alloc();
  singleSlater112.formGuess();
  singleSlater112.computeProperties();
  singleSlater112.printProperties();
  singleSlater112.SCF3();
  singleSlater112.computeProperties();
  singleSlater112.printProperties();
  singleSlater112.mullikenPop();

  cout << "Energy " << singleSlater112.totalEnergy() <<endl;
  cout << "PASSED 112? " << bool(std::abs(singleSlater112.totalEnergy() - (-74.9289919100)) <= 1e-8) << endl;
// TEST 115
  singleSlater115.setRef("RSVWN5");
//  singleSlater115.isDFT = true;
//  singleSlater115.isHF  = false;
  singleSlater115.setSCFEneTol(1e-12);
  singleSlater115.setSCFMaxIter(10000);
  singleSlater115.doDIIS = true;
  singleSlater115.setGuess(CORE);
//  singleSlater115.addSlater();
//  singleSlater115.addVWN5();
//  singleSlater115.setxHF(0.0);
  singleSlater115.communicate(moleculeWat,basisWat,aointsWat,fileio,memManager);
  singleSlater115.initMeta();
//  singleSlater115.genMethString();
  singleSlater115.alloc();
  singleSlater115.formGuess();
  singleSlater115.computeProperties();
  singleSlater115.printProperties();
  singleSlater115.SCF3();
  singleSlater115.computeProperties();
  singleSlater115.printProperties();
  singleSlater115.mullikenPop();

  cout << "Energy " << singleSlater115.totalEnergy() <<endl;
  cout << "PASSED 115? " << bool(std::abs(singleSlater115.totalEnergy() - (-74.7332742682)) <= 1e-8) << endl;
// TEST 136
  singleSlater136.setRef("RB88");
//  singleSlater136.isDFT = true;
//  singleSlater136.isHF  = false;
  singleSlater136.setSCFEneTol(1e-12);
  singleSlater136.setSCFMaxIter(10000);
  singleSlater136.doDIIS = true;
  singleSlater136.setGuess(CORE);
//  singleSlater136.createB88();
  singleSlater136.communicate(moleculeWat,basisWat,aointsWat,fileio,memManager);
  singleSlater136.initMeta();
//  singleSlater136.genMethString();
  singleSlater136.alloc();
  singleSlater136.formGuess();
  singleSlater136.computeProperties();
  singleSlater136.printProperties();
  singleSlater136.SCF3();
  singleSlater136.computeProperties();
  singleSlater136.printProperties();
  singleSlater136.mullikenPop();

  cout << "Energy " << singleSlater136.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater136.totalEnergy() - (-74.9512108608)) <<endl;
  cout << "PASSED 136? " << bool(std::abs(singleSlater136.totalEnergy() - (-74.9512108608)) <= 1e-8) << endl;
// TEST 150
  singleSlater150.setRef("RBLYP");
//  singleSlater150.isDFT = true;
//  singleSlater150.isHF  = false;
  singleSlater150.setSCFEneTol(1e-12);
  singleSlater150.setSCFMaxIter(10000);
  singleSlater150.doDIIS = true;
  singleSlater150.setGuess(CORE);
//  singleSlater150.addSlater();
//  singleSlater150.addB88();
//  singleSlater150.addLYP();
//  singleSlater150.setxHF(0.0);
  singleSlater150.communicate(moleculeWat,basisWat,aointsWat,fileio,memManager);
  singleSlater150.initMeta();
//  singleSlater150.genMethString();
  singleSlater150.alloc();
  singleSlater150.formGuess();
  singleSlater150.computeProperties();
  singleSlater150.printProperties();
  singleSlater150.SCF3();
  singleSlater150.computeProperties();
  singleSlater150.printProperties();
  singleSlater150.mullikenPop();

  cout << "Energy " << singleSlater150.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater150.totalEnergy() - (-75.2846829249)) <<endl;
  cout << "PASSED 150? " << bool(std::abs(singleSlater150.totalEnergy() - (-75.2846829249)) <= 1e-8) << endl;
// TEST B3LYP WATER 
  singleSlater160.setRef("RB3LYP");
//  singleSlater160.isDFT = true;
//  singleSlater160.isHF  = false;
  singleSlater160.setSCFEneTol(1e-12);
  singleSlater160.setSCFMaxIter(10000);
  singleSlater160.doDIIS = true;
  singleSlater160.setGuess(CORE);
//  singleSlater160.createB3LYP();
  singleSlater160.communicate(moleculeWat,basisWat,aointsWat,fileio,memManager);
  singleSlater160.initMeta();
//  singleSlater160.genMethString();
  singleSlater160.alloc();
  singleSlater160.formGuess();
  singleSlater160.computeProperties();
  singleSlater160.printProperties();
  singleSlater160.SCF3();
  singleSlater160.computeProperties();
  singleSlater160.printProperties();
  singleSlater160.mullikenPop();

  cout << "Energy " << singleSlater160.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater160.totalEnergy() - (-75.3122931451)) <<endl;
  cout << "PASSED 160 - WATER B3LYP STO3G? " << bool(std::abs(singleSlater160.totalEnergy() - (-75.3122931451)) <= 1e-8) << endl;

// TEST B3LYP WATER 
  singleSlater163.setRef("RBHandH");
//  singleSlater163.isDFT = true;
//  singleSlater163.isHF  = false;
  singleSlater163.setSCFEneTol(1e-12);
  singleSlater163.setSCFMaxIter(10000);
  singleSlater163.doDIIS = true;
  singleSlater163.setGuess(CORE);
//  singleSlater163.createBHandH();
  singleSlater163.communicate(moleculeWat,basisWat,aointsWat,fileio,memManager);
  singleSlater163.initMeta();
//  singleSlater163.genMethString();
  singleSlater163.alloc();
  singleSlater163.formGuess();
  singleSlater163.computeProperties();
  singleSlater163.printProperties();
  singleSlater163.SCF3();
  singleSlater163.computeProperties();
  singleSlater163.printProperties();
  singleSlater163.mullikenPop();

  cout << "Energy " << singleSlater163.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater163.totalEnergy() - (-75.3122931451)) <<endl;
  cout << "PASSED 163 - WATER BHandH STO3G? " << bool(std::abs(singleSlater163.totalEnergy() - (-75.3122931451)) <= 1e-8) << endl;

//  Lit UKS
  singleSlater118.setRef("USLATER");
//  singleSlater118.isDFT = true;
//  singleSlater118.isHF  = false;
  singleSlater118.setSCFEneTol(1e-12);
  singleSlater118.setSCFMaxIter(10000);
  singleSlater118.doDIIS = true;
  singleSlater118.setGuess(CORE);
//  singleSlater118.addSlater();
//  singleSlater118.setxHF(0.0);
  singleSlater118.communicate(moleculeLi,basisLi,aointsLi,fileio,memManager);
  singleSlater118.initMeta();
//  singleSlater118.genMethString();
  singleSlater118.alloc();
  singleSlater118.formGuess();
  singleSlater118.computeProperties();
  singleSlater118.printProperties();
  singleSlater118.SCF3();
  singleSlater118.computeProperties();
  singleSlater118.printProperties();
  singleSlater118.mullikenPop();

  cout << "Energy " << singleSlater118.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater118.totalEnergy() - (-7.0664115457)) <<endl;
  cout << "PASSED 118 - Li Slater-STO-3G? " << bool(std::abs(singleSlater118.totalEnergy() - (-7.0664115457)) <= 1e-8) << endl;
////
  singleSlater124.setRef("ULSDA");
//  singleSlater124.isDFT = true;
//  singleSlater124.isHF  = false;
  singleSlater124.setSCFEneTol(1e-12);
  singleSlater124.setSCFMaxIter(10000);
  singleSlater124.doDIIS = true;
  singleSlater124.setGuess(CORE);
//  singleSlater124.createLSDA();
  singleSlater124.communicate(moleculeLi,basisLi,aointsLi,fileio,memManager);
  singleSlater124.initMeta();
//  singleSlater124.genMethString();
  singleSlater124.alloc();
  singleSlater124.formGuess();
  singleSlater124.computeProperties();
  singleSlater124.printProperties();
  singleSlater124.SCF3();
  singleSlater124.computeProperties();
  singleSlater124.printProperties();
  singleSlater124.mullikenPop();

  cout << "Energy " << singleSlater124.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater124.totalEnergy() - (-7.2764020089)) <<endl;
  cout << "PASSED 124 - Li LSDA-STO-3G? " << bool(std::abs(singleSlater124.totalEnergy() - (-7.2764020089)) <= 1e-8) << endl;
////
  singleSlater130.setRef("USVWN5");
//  singleSlater130.isDFT = true;
//  singleSlater130.isHF  = false;
  singleSlater130.setSCFEneTol(1e-12);
  singleSlater130.setSCFMaxIter(10000);
  singleSlater130.doDIIS = true;
  singleSlater130.setGuess(CORE);
//  singleSlater130.addSlater();
//  singleSlater130.addVWN5();
//  singleSlater130.setxHF(0.0);
  singleSlater130.communicate(moleculeLi,basisLi,aointsLi,fileio,memManager);
  singleSlater130.initMeta();
//  singleSlater130.genMethString();
  singleSlater130.alloc();
  singleSlater130.formGuess();
  singleSlater130.computeProperties();
  singleSlater130.printProperties();
  singleSlater130.SCF3();
  singleSlater130.computeProperties();
  singleSlater130.printProperties();
  singleSlater130.mullikenPop();

  cout << "Energy " << singleSlater130.totalEnergy() <<endl;
  cout << "Energy " << std::abs(singleSlater130.totalEnergy() - (-7.2213016498)) <<endl;
  cout << "PASSED 130 - Li SVWN5-STO-3G? " << bool(std::abs(singleSlater130.totalEnergy() - (-7.2213016498)) <= 1e-8) << endl;
*/
// OXy
//TES
//T - 112 Real RKS-SCF-LSDA (SVWN3) - STO-3G
/*
  molecule2.convBohr();
  molecule2.computeNucRep();
  molecule2.computeRij();
  molecule2.computeI();

  singleSlater2.setRef(RHF);
  singleSlater2.isDFT = true;
  singleSlater2.isHF  = false;
  singleSlater2.setSCFEneTol(1e-12);
  singleSlater2.setSCFMaxIter(10000);
  singleSlater2.doDIIS = true;

  singleSlater2.setGuess(CORE);
//  singleSlater2.setGuess(READ);
//  fileio.doRestart = true;


  singleSlater2.isDFT = true;
  singleSlater2.isHF = false;
//singleSlater2.setExchKernel(SingleSlater<double>::EXCH::B88);
  //singleSlater2.setExchKernel(SingleSlater<double>::EXCH::NOEXCH);
//  singleSlater2.setCorrKernel(SingleSlater<double>::CORR::NOCORR);
//  singleSlater2.setCorrKernel(SingleSlater<double>::CORR::VWN5);
//  singleSlater2.addB88();
//  singleSlater2.addLYP();
//  singleSlater2.createB88();
//    singleSlater2.createB3LYP();
  singleSlater2.createLSDA();
//singleSlater2.addVWN5();
//singleSlater2.setPrintLevel(5);

//singleSlater2.setxHF(0.0);
//basis2.findBasisFile("sto-3g");
//basis2.findBasisFile("3-21g");
  basis2.findBasisFile("STO-3G");
//basis2.findBasisFile("cc-pVTZ");
  basis2.communicate(fileio);
  basis2.parseGlobal();
  basis2.constructLocal(&molecule2);
  basis2.makeMaps(&molecule2);
  basis2.renormShells();


  aoints2.communicate(molecule2,basis2,fileio,memManager);
  singleSlater2.communicate(molecule2,basis2,aoints2,fileio,memManager);
//moints.communicate(molecule2,basis2,fileio,aoints2,singleSlater2);

  aoints2.initMeta();
  singleSlater2.initMeta();
  singleSlater2.genMethString();

  aoints2.alloc();
  singleSlater2.alloc();

  singleSlater2.formGuess();
  singleSlater2.computeProperties();
  singleSlater2.printProperties();
//singleSlater2.formFock();
//singleSlater2.computeEnergy();
  singleSlater2.SCF3();
  singleSlater2.computeProperties();
  singleSlater2.printProperties();

  singleSlater2.mullikenPop();

  cout << "Energy " << singleSlater2.totalEnergy() <<endl;
  cout << "PASSED 112? " << bool(std::abs(singleSlater2.totalEnergy() - (-74.9289919100)) <= 1e-8) << endl;
///////////////////////////////////////////////////////////////////////////
//TEST - 115 Real RKS-SCF-SVWN5 - STO-3G
//loadPresets<H>(molecule);
//loadPresets<OxMolecule>(molecule);
  Molecule molecule3;
  BasisSet basis3;
  AOIntegrals aoints3;
  SingleSlater<double> singleSlater3;
  loadPresets<WATER>(molecule3);
//loadPresets<Methanol>(molecule3);
//loadPresets<HE>(molecule3);
//loadPresets<SO>(molecule3);
//loadPresets<Li>(molecule3);
  molecule3.convBohr();
  molecule3.computeNucRep();
  molecule3.computeRij();
  molecule3.computeI();

  singleSlater3.setRef(RHF);
  singleSlater3.isDFT = true;
  singleSlater3.isHF  = false;
  singleSlater3.setSCFEneTol(1e-12);
  singleSlater3.setSCFMaxIter(10000);
  singleSlater3.doDIIS = true;

  singleSlater3.setGuess(CORE);
//  singleSlater3.setGuess(READ);
//  fileio.doRestart = true;


  singleSlater3.isDFT = true;
  singleSlater3.isHF = false;
//singleSlater3.setExchKernel(SingleSlater<double>::EXCH::B88);
  //singleSlater3.setExchKernel(SingleSlater<double>::EXCH::NOEXCH);
//  singleSlater3.setCorrKernel(SingleSlater<double>::CORR::NOCORR);
//  singleSlater3.setCorrKernel(SingleSlater<double>::CORR::VWN5);
  singleSlater3.addSlater();
  singleSlater3.addVWN5();
//  singleSlater3.addLYP();
//  singleSlater3.createB88();
//    singleSlater3.createB3LYP();
//  singleSlater3.createLSDA();
//singleSlater3.addVWN5();
//singleSlater3.setPrintLevel(5);

  singleSlater3.setxHF(0.0);
//basis3.findBasisFile("sto-3g");
//basis3.findBasisFile("3-21g");
  basis3.findBasisFile("STO-3G");
//basis3.findBasisFile("cc-pVTZ");
  basis3.communicate(fileio);
  basis3.parseGlobal();
  basis3.constructLocal(&molecule3);
  basis3.makeMaps(&molecule3);
  basis3.renormShells();


  aoints3.communicate(molecule3,basis3,fileio,memManager);
  singleSlater3.communicate(molecule3,basis3,aoints3,fileio,memManager);
//moints.communicate(molecule3,basis3,fileio,aoints3,singleSlater3);

  aoints3.initMeta();
  singleSlater3.initMeta();
  singleSlater3.genMethString();

  aoints3.alloc();
  singleSlater3.alloc();

  singleSlater3.formGuess();
  singleSlater3.computeProperties();
  singleSlater3.printProperties();
//singleSlater3.formFock();
//singleSlater3.computeEnergy();
  singleSlater3.SCF3();
  singleSlater3.computeProperties();
  singleSlater3.printProperties();

  singleSlater3.mullikenPop();
  cout << "Energy " << singleSlater3.totalEnergy() <<endl;
  cout << "PASSED 115? " << bool(std::abs(singleSlater3.totalEnergy() - (-74.7332742682)) <= 1e-8) << endl;
//TEST - 118 Real UKS-SCF-SLATER - STO-3G

  Molecule molecule4;
  BasisSet basis4;
  AOIntegrals aoints4;
  SingleSlater<double> singleSlater4;
  loadPresets<Li>(molecule4);
  molecule4.convBohr();
  molecule4.computeNucRep();
  molecule4.computeRij();
  molecule4.computeI();

  singleSlater4.setRef(UHF);
  singleSlater4.isDFT = true;
  singleSlater4.isHF  = false;
  singleSlater4.setSCFEneTol(1e-12);
  singleSlater4.setSCFMaxIter(10000);
  singleSlater4.doDIIS = true;

  singleSlater4.setGuess(CORE);

  singleSlater4.isDFT = true;
  singleSlater4.isHF = false;
  singleSlater4.addSlater();
  singleSlater4.setxHF(0.0);
  basis4.findBasisFile("STO-3G");
  basis4.communicate(fileio);
  basis4.parseGlobal();
  basis4.constructLocal(&molecule4);
  basis4.makeMaps(&molecule4);
  basis4.renormShells();


  aoints4.communicate(molecule4,basis4,fileio,memManager);
  singleSlater4.communicate(molecule4,basis4,aoints4,fileio,memManager);

  aoints4.initMeta();
  singleSlater4.initMeta();
  singleSlater4.genMethString();

  aoints4.alloc();
  singleSlater4.alloc();

  singleSlater4.formGuess();
  singleSlater4.computeProperties();
  singleSlater4.printProperties();
  singleSlater4.SCF3();
  singleSlater4.computeProperties();
  singleSlater4.printProperties();

  singleSlater4.mullikenPop();
  cout << "Energy " << singleSlater4.totalEnergy() <<endl;
  cout << "PASSED 118? " << bool(std::abs(singleSlater4.totalEnergy() - (-7.0664115457)) <= 1e-8) << endl;

//TEST - 121 Real UKS-SCF-SLATER - STO-3G
  Molecule molecule5;
  BasisSet basis5;
  AOIntegrals aoints5;
  SingleSlater<double> singleSlater5;
  loadPresets<OxMolecule>(molecule5);
  molecule5.computeNucRep();
  molecule5.computeRij();
  molecule5.computeI();

  singleSlater5.setRef(UHF);
  singleSlater5.isDFT = true;
  singleSlater5.isHF  = false;
  singleSlater5.setSCFEneTol(1e-12);
  singleSlater5.setSCFMaxIter(10000);
  singleSlater5.doDIIS = true;

  singleSlater5.setGuess(CORE);

  singleSlater5.isDFT = true;
  singleSlater5.isHF = false;
  singleSlater5.addSlater();
  singleSlater5.setxHF(0.0);
  basis5.findBasisFile("STO-3G");
  basis5.communicate(fileio);
  basis5.parseGlobal();
  basis5.constructLocal(&molecule5);
  basis5.makeMaps(&molecule5);
  basis5.renormShells();


  aoints5.communicate(molecule5,basis5,fileio,memManager);
  singleSlater5.communicate(molecule5,basis5,aoints5,fileio,memManager);
  aoints5.initMeta();
  singleSlater5.initMeta();
  singleSlater5.genMethString();

  aoints5.alloc();
  singleSlater5.alloc();

  singleSlater5.formGuess();
  singleSlater5.computeProperties();
  singleSlater5.printProperties();
  singleSlater5.SCF3();
  singleSlater5.computeProperties();
  singleSlater5.printProperties();

  singleSlater5.mullikenPop();
  cout << "Energy " << singleSlater5.totalEnergy() <<endl;
  cout << "PASSED 121? " << bool(std::abs(singleSlater5.totalEnergy() - (-146.0745567304)) <= 1e-8) << endl;
*/
/*


  singleSlater.onePDMOrtho()[1]->swap(*singleSlater.onePDMOrtho()[3]);
  singleSlater.fockOrtho()[1]->swap(*singleSlater.fockOrtho()[3]);
*/
/*
  singleSlater.onePDMOrtho()[1]->swap(*singleSlater.onePDMOrtho()[2]);
  singleSlater.fockOrtho()[1]->swap(*singleSlater.fockOrtho()[2]);
*/

/*
  // Rotate 90 deg around y (should yield x)
  singleSlater.rotateDensities({1,0,0},math.pi/4);
  singleSlater.computeProperties();
  singleSlater.printProperties();


  ComplexMatrix P(2*singleSlater.onePDMMz()->rows(),2*singleSlater.onePDMMz()->rows());
  ComplexMatrix F(2*singleSlater.onePDMMz()->rows(),2*singleSlater.onePDMMz()->rows());

  std::vector<std::reference_wrapper<ComplexMap>> scatteredFock_;
  std::vector<std::reference_wrapper<ComplexMap>> scatteredDen_;

  ComplexMap PMap(P.data(),P.rows(),P.cols());
  ComplexMap FMap(F.data(),F.rows(),F.cols());

  ComplexMatrix TMP(basis.nBasis()*singleSlater.nTCS(),basis.nBasis()*singleSlater.nTCS());
  TMP.setZero();

  Quantum<dcomplex>::spinGather(FMap,scatteredFock_);  
  Quantum<dcomplex>::spinGather(PMap,scatteredDen_);  


  prettyPrintSmart(cout,*singleSlater.moA(),"MO");
  for(auto i = 0; i < 2*basis.nBasis(); i++) {
  for(auto j = 0; j < 2*basis.nBasis(); j+=2) {
    dcomplex a = (*singleSlater.moA())(j,i); 
    dcomplex b = (*singleSlater.moA())(j+1,i); 
    // Y Rotation ...

    (*singleSlater.moA())(j,i) 
       = std::cos(-math.pi / 4) * a  +  std::sin(-math.pi / 4) * b; 
    (*singleSlater.moA())(j+1,i) 
       = -std::sin(-math.pi / 4) * a  +  std::cos(-math.pi / 4) * b; 


    // X Rotation ...
    (*singleSlater.moA())(j,i) 
       = std::cos(math.pi / 4) * a  +  math.ii * std::sin(math.pi / 4) * b; 
    (*singleSlater.moA())(j+1,i) 
       = math.ii * std::sin(math.pi / 4) * a  +  std::cos(math.pi / 4) * b; 
  }
  }
  prettyPrintSmart(cout,*singleSlater.moA(),"MO");
  singleSlater.formDensity();
  singleSlater.computeProperties();
  singleSlater.printProperties();

//singleSlater.printFock();
//singleSlater.printPT();


  ComplexMatrix FPScalar2(FPScalar);
  FPScalar2.setZero();
  for(auto i = 0; i < singleSlater.fockOrtho().size(); i++) {
    FPScalar2 += 
      (*singleSlater.fockOrtho()[i]) * (*singleSlater.onePDMOrtho()[i]);
  } 
  FPScalar2 *= 0.5;

  ComplexMatrix FPMz2(FPMz);
  FPMz2  = (*singleSlater.fockOrtho()[0]) * (*singleSlater.onePDMOrtho()[1]);
  FPMz2 += (*singleSlater.fockOrtho()[1]) * (*singleSlater.onePDMOrtho()[0]);
  FPMz2 += math.ii * 
    ((*singleSlater.fockOrtho()[3]) * (*singleSlater.onePDMOrtho()[2]));
  FPMz2 -= math.ii * 
    ((*singleSlater.fockOrtho()[2]) * (*singleSlater.onePDMOrtho()[3]));
  FPMz2 *= 0.5;

  ComplexMatrix FPMy2(FPMy);
  FPMy2  = (*singleSlater.fockOrtho()[0]) * (*singleSlater.onePDMOrtho()[2]);
  FPMy2 += (*singleSlater.fockOrtho()[2]) * (*singleSlater.onePDMOrtho()[0]);
  FPMy2 += math.ii * 
    ((*singleSlater.fockOrtho()[1]) * (*singleSlater.onePDMOrtho()[3]));
  FPMy2 -= math.ii * 
    ((*singleSlater.fockOrtho()[3]) * (*singleSlater.onePDMOrtho()[1]));
  FPMy2 *= 0.5;

  ComplexMatrix FPMx2(FPMx);
  FPMx2  = (*singleSlater.fockOrtho()[0]) * (*singleSlater.onePDMOrtho()[3]);
  FPMx2 += (*singleSlater.fockOrtho()[3]) * (*singleSlater.onePDMOrtho()[0]);
  FPMx2 += math.ii * 
    ((*singleSlater.fockOrtho()[2]) * (*singleSlater.onePDMOrtho()[1]));
  FPMx2 -= math.ii * 
    ((*singleSlater.fockOrtho()[1]) * (*singleSlater.onePDMOrtho()[2]));
  FPMx2 *= 0.5;

//prettyPrintSmart(cout,FPScalar,"Scalar1");
//prettyPrintSmart(cout,FPScalar2,"Scalar2");
  prettyPrintSmart(cout,FPScalar - FPScalar2,"Diff Scalar");
  prettyPrintSmart(cout,FPMz - FPMz2,"Diff Mz");
  prettyPrintSmart(cout,FPMy - FPMy2,"Diff My");
  prettyPrintSmart(cout,FPMx - FPMx2,"Diff Mx");
*/
/*
  rt.communicate(singleSlater);
  rt.alloc();
//rt.setMaxSteps(827000); // roughly 1 ps of dynamics
  rt.setMaxSteps(1000); 
//rt.setTOff(0.10);
//rt.setEDFieldAmp({0.001,0.0,0.0});
  rt.setIEnvlp(Step);
  rt.doPropagation();
*/
/*
  cout << endl;
  MOIntegrals<double> moints;
  moints.communicate(singleSlater,memManager);
  moints.initMeta();
//moints.testMOInts();
  FOPPA<double> resp(DIAGONALIZATION,SPIN_SEPARATED,false,false);
  resp.communicate(singleSlater,memManager);
//resp.doFull();
  resp.setNSek(3);
  resp.setNGuess(10);
  resp.initMeta();
  resp.alloc();
  resp.runResponse();
*/
  finalizeCQ();
  return 0;
};

