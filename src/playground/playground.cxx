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
  WATER,Methanol,H,HE,SO,OxMolecule,Li,MnAcAc
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

template<>
void loadPresets<MnAcAc>(Molecule &molecule){

    molecule.setCharge(0);
    molecule.setNTotalE(87);
    molecule.setMultip(2);
    molecule.setNAtoms(15);
    molecule.alloc();

    molecule.setIndex(0,HashAtom("C",0));
    molecule.setIndex(1,HashAtom("O",0));
    molecule.setIndex(2,HashAtom("H",0));
    molecule.setIndex(3,HashAtom("C",0));
    molecule.setIndex(4,HashAtom("O",0));
    molecule.setIndex(5,HashAtom("H",0));
    molecule.setIndex(6,HashAtom("H",0));
    molecule.setIndex(7,HashAtom("C",0));
    molecule.setIndex(8,HashAtom("O",0));
    molecule.setIndex(9,HashAtom("H",0));
    molecule.setIndex(10,HashAtom("C",0));
    molecule.setIndex(11,HashAtom("O",0));
    molecule.setIndex(12,HashAtom("H",0));
    molecule.setIndex(13,HashAtom("H",0));
    molecule.setIndex(14,HashAtom("Mn",0));

    molecule.setCart(0,  0.000000,  1.230231, 1.676654);
    molecule.setCart(1,  0.000000,  0.836540, 0.491685);
    molecule.setCart(2,  0.000000,  2.268652, 1.888047);
    molecule.setCart(3, -0.000000, -1.230231, 1.676654);
    molecule.setCart(4, -0.000000, -0.836540, 0.491685);
    molecule.setCart(5, -0.000000, -2.268652, 1.888047);
    molecule.setCart(6,  0.000000,  0.000000, 2.514988);
    molecule.setCart(7,  1.230231, -0.000000,-1.676654);
    molecule.setCart(8,  0.836540, -0.000000,-0.491685);
    molecule.setCart(9,  2.268652, -0.000000,-1.888047);
    molecule.setCart(10,-1.230231,  0.000000,-1.676654);
    molecule.setCart(11,-0.836540,  0.000000,-0.491685);
    molecule.setCart(12,-2.268652,  0.000000,-1.888047);
    molecule.setCart(13, 0.000000,  0.000000,-2.514988);
    molecule.setCart(14, 0.000000,  0.000000, 0.000000);
}

int main(int argc, char **argv){

  CQMemManager memManager;
  Molecule molecule;
  BasisSet basis;
  AOIntegrals aoints;
  SingleSlater<dcomplex> singleSlater;
  RealTime<dcomplex> rt;
  FileIO fileio("test.inp","test.out");

  memManager.setTotalMem(256e6);
  initCQ(argc,argv);
  CQSetNumThreads(4);
  
//////////////////////////////////////////////////////
//loadPresets<H>(molecule);
//loadPresets<OxMolecule>(molecule);
//loadPresets<WATER>(molecule);
//loadPresets<Methanol>(molecule);
//loadPresets<HE>(molecule);
//loadPresets<SO>(molecule);
  loadPresets<Li>(molecule);
//loadPresets<MnAcAc>(molecule);
  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef("GHF");
//singleSlater.setSCFEneTol(1e-12);
  singleSlater.setSCFMaxIter(10000);
  singleSlater.doDIIS = true;
  singleSlater.doDamp = false;
//singleSlater.dampParam = 0.2;

  singleSlater.setGuess(CORE);

  fileio.iniH5Files();


//basis.forceCart();
  basis.findBasisFile("sto-3g");
//basis.findBasisFile("3-21g");
//basis.findBasisFile("6-31G");
//basis.findBasisFile("cc-pVTZ");
//basis.findBasisFile("cc-pVDZ");
  basis.communicate(fileio);
  basis.parseGlobal();
  basis.constructLocal(&molecule);
  basis.makeMaps(&molecule);
//basis.renormShells();


  aoints.setAlgorithm(AOIntegrals::INTEGRAL_ALGORITHM::INCORE);
  aoints.communicate(molecule,basis,fileio,memManager);
  singleSlater.communicate(molecule,basis,aoints,fileio,memManager);
//moints.communicate(molecule,basis,fileio,aoints,singleSlater);

  aoints.initMeta();
  singleSlater.initMeta();

  aoints.alloc();
  singleSlater.alloc();

//singleSlater.setPrintLevel(4);
  singleSlater.formGuess();
//singleSlater.printDensity();
  singleSlater.SCF3();
  singleSlater.computeProperties();
  singleSlater.printProperties();

  singleSlater.rotateDensities({1.0,0.0,0.0},math.pi/2);
  singleSlater.SCF3();
  singleSlater.computeProperties();
  singleSlater.printProperties();

  prettyPrintSmart(cout,*singleSlater.moA(),"MO");

/*
  rt.communicate(singleSlater);
  rt.alloc();
//rt.setMaxSteps(827000); // roughly 1 ps of dynamics
  rt.setMaxSteps(10); 
  rt.setTOff(0.0000001);
  rt.setEDFieldAmp({0.0005,0.0,0.0});
  rt.setIEnvlp(Step);
  rt.doPropagation();
*/
  MOIntegrals<dcomplex> moints;
  moints.communicate(singleSlater,memManager);
  moints.initMeta();
  moints.testMOInts();
/*
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

