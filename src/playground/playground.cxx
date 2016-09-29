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
  Molecule molecule;
  BasisSet basis;
  AOIntegrals aoints;
  SingleSlater<dcomplex> singleSlater;
//RealTime<double> rt;
  FileIO fileio("test.inp","test.out");

  memManager.setTotalMem(256e6);
  initCQ(argc,argv);
  fileio.iniH5Files();
  CQSetNumThreads(1);
  
//loadPresets<H>(molecule);
  loadPresets<OxMolecule>(molecule);
//loadPresets<WATER>(molecule);
//loadPresets<Methanol>(molecule);
//loadPresets<HE>(molecule);
//loadPresets<SO>(molecule);
//loadPresets<Li>(molecule);
  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(X2C);
  singleSlater.setGuess(CORE);
  singleSlater.isClosedShell = false;
  singleSlater.setNTCS(2);
  singleSlater.setSCFEneTol(1e-12);
  singleSlater.setSCFMaxIter(10000);
  singleSlater.doDIIS = true;

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
//basis.findBasisFile("cc-pVTZ");
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
  aoints.doX2C = true;
  aoints.useFiniteWidthNuclei = true;
  aoints.setPrintLevel(3);
  aoints.twoEFudge = 2; //default is 1
  singleSlater.initMeta();
  singleSlater.genMethString();

  aoints.alloc();
  singleSlater.alloc();

  singleSlater.formGuess();
  singleSlater.computeProperties();
  singleSlater.printProperties();
//singleSlater.formFock();
//singleSlater.computeEnergy();
  singleSlater.SCF3();
  singleSlater.computeProperties();
  singleSlater.printProperties();

  singleSlater.onePDMMz()->swap(*singleSlater.onePDMMx());
  singleSlater.computeProperties();
  singleSlater.printProperties();
  singleSlater.doDIIS = false; 
  singleSlater.SCF3();
  singleSlater.computeProperties();
  singleSlater.printProperties();

/*
  singleSlater.setGuess(SingleSlater<dcomplex>::ONLY);
  singleSlater.onePDMMz()->swap(*singleSlater.onePDMMx());
  singleSlater.PTMz()->swap(*singleSlater.PTMx());
  singleSlater.computeProperties();
  singleSlater.printProperties();

  ComplexMatrix TMP(basis.nBasis()*singleSlater.nTCS(),basis.nBasis()*singleSlater.nTCS());
  TMP.setZero();

  ComplexMap TMPMap(TMP.data(),TMP.rows(),TMP.cols());
  std::vector<std::reference_wrapper<ComplexMap>> scattered;
  scattered.emplace_back(*singleSlater.onePDMScalar());
  scattered.emplace_back(*singleSlater.onePDMMz());
  scattered.emplace_back(*singleSlater.onePDMMy());
  scattered.emplace_back(*singleSlater.onePDMMx());
  Quantum<dcomplex>::spinGather(TMPMap,scattered);
*/

/*
//for(auto OPDM : singleSlater.onePDM()) *OPDM *= 0.5;
  singleSlater.printDensity();
  prettyPrintSmart(fileio.out,TMPMap,"Gathered");
  Quantum<dcomplex>::spinScatter(TMPMap,scattered);
  singleSlater.printDensity();
*/

/*
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

  prettyPrintSmart(cout,*singleSlater.moA(),"MO");
  singleSlater.setPrintLevel(4);
  singleSlater.doDIIS = false;
//singleSlater.formGuess();
  singleSlater.SCF3();
  singleSlater.computeProperties();
  singleSlater.printProperties();
  prettyPrintSmart(cout,*singleSlater.moA(),"MO");
*/

/*
  rt.communicate(singleSlater);
  rt.alloc();
//rt.setMaxSteps(827000); // roughly 1 ps of dynamics
  rt.setMaxSteps(1000); 
  rt.setTOff(0.10);
  rt.setEDFieldAmp({0.001,0.0,0.0});
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
  finalizeCQ();
*/
  return 0;
};

