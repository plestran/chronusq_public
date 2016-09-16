//#define EIGEN_RUNTIME_NO_MALLOC
//#include <response.h>
//#include <workers.h>
//#include <pythonapi.h>
#include <global.h>
#include <pythonapi.h>
#include <singleslater.h>
#include <mointegrals.h>
#include <response.h>

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
  WATER, HE,SO,Li
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
  Molecule molecule;
  BasisSet basis;
  AOIntegrals aoints;
  SingleSlater<double> singleSlater;
  FileIO fileio("test.inp","test.out");

  memManager.setTotalMem(256e6);
  initCQ(argc,argv);
  fileio.iniH5Files();
  CQSetNumThreads(1);
  
  loadPresets<WATER>(molecule);
//loadPresets<HE>(molecule);
//loadPresets<SO>(molecule);
//loadPresets<Li>(molecule);
  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(SingleSlater<double>::RHF);
  singleSlater.isClosedShell = true;
//singleSlater.setRef(SingleSlater<double>::UHF);
//singleSlater.isClosedShell = false;

/*
  singleSlater.isDFT = true;
  singleSlater.isHF = false;
//singleSlater.setExchKernel(SingleSlater<double>::EXCH::B88);
  //singleSlater.setExchKernel(SingleSlater<double>::EXCH::NOEXCH);
//  singleSlater.setCorrKernel(SingleSlater<double>::CORR::NOCORR);
//  singleSlater.setCorrKernel(SingleSlater<double>::CORR::VWN5);
  singleSlater.addB88();
  singleSlater.addLYP();
*/
//singleSlater.addSlater();
//singleSlater.addVWN5();
//singleSlater.setPrintLevel(5);

//basis.findBasisFile("sto-3g");
//basis.findBasisFile("3-21g");
  basis.findBasisFile("6-31G");
  basis.communicate(fileio);
  basis.parseGlobal();
  basis.constructLocal(&molecule);
  basis.makeMaps(&molecule);
  basis.renormShells();


  aoints.communicate(molecule,basis,fileio,memManager);
  singleSlater.communicate(molecule,basis,aoints,fileio,memManager);
//moints.communicate(molecule,basis,fileio,aoints,singleSlater);

  aoints.initMeta();
  aoints.integralAlgorithm = AOIntegrals::INCORE;
  singleSlater.initMeta();
  singleSlater.genMethString();
  singleSlater.setSCFEneTol(1e-12);

  aoints.alloc();
  singleSlater.alloc();

  singleSlater.formGuess();
//singleSlater.formFock();
//singleSlater.computeEnergy();
  singleSlater.SCF3();
  singleSlater.computeProperties();
  singleSlater.printProperties();


  cout << endl;
  MOIntegrals<double> moints;
  moints.communicate(singleSlater,memManager);
  moints.initMeta();
//moints.testMOInts();
  FOPPA<double> resp(DIAGONALIZATION,SPIN_ADAPTED,true);
  resp.communicate(singleSlater,memManager);
//resp.doFull();
  resp.setNSek(2);
  resp.setNGuess(8);
  resp.initMeta();
  resp.alloc();
  resp.runResponse();
  finalizeCQ();
  return 0;
};

