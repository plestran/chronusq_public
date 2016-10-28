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

  finalizeCQ();
  return 0;
};

