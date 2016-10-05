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
  SingleSlater<double> singleSlater;
  RealTime<double> rt;
  FileIO fileio("test.inp","test.out");

  memManager.setTotalMem(256e6);
  initCQ(argc,argv);
  CQSetNumThreads(1);
  
//////////////////////////////////////////////////////
//loadPresets<H>(molecule);
//loadPresets<OxMolecule>(molecule);
  loadPresets<WATER>(molecule);
//loadPresets<Methanol>(molecule);
//loadPresets<HE>(molecule);
//loadPresets<SO>(molecule);
//loadPresets<Li>(molecule);
  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(RHF);
  singleSlater.isDFT = true;
  singleSlater.isHF  = false;
  singleSlater.setSCFEneTol(1e-12);
  singleSlater.setSCFMaxIter(10000);
  singleSlater.doDIIS = true;

  singleSlater.setGuess(CORE);
//  singleSlater.setGuess(READ);
//  fileio.doRestart = true;

  fileio.iniH5Files();

  singleSlater.isDFT = true;
  singleSlater.isHF = false;
//singleSlater.setExchKernel(SingleSlater<double>::EXCH::B88);
  //singleSlater.setExchKernel(SingleSlater<double>::EXCH::NOEXCH);
//  singleSlater.setCorrKernel(SingleSlater<double>::CORR::NOCORR);
//  singleSlater.setCorrKernel(SingleSlater<double>::CORR::VWN5);
//  singleSlater.addB88();
//  singleSlater.addLYP();
//  singleSlater.createB88();
    singleSlater.createB3LYP();
//singleSlater.addSlater();
//singleSlater.addVWN5();
//singleSlater.setPrintLevel(5);

//singleSlater.setxHF(0.0);
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

  singleSlater.mullikenPop();

  cout << "PASSED? " << bool((singleSlater.totalEnergy() - (-76.3671171302)) <= 1e-8) << endl;
///////////////////////////////////////////////////////////////////////////
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

