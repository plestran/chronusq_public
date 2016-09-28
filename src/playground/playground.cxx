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
  SingleSlater<dcomplex> singleSlater;
  RealTime<dcomplex> rt;
  FileIO fileio("test.inp","test.out");

  memManager.setTotalMem(256e6);
  initCQ(argc,argv);
  fileio.iniH5Files();
  CQSetNumThreads(1);
  
//loadPresets<WATER>(molecule);
//loadPresets<HE>(molecule);
//loadPresets<SO>(molecule);
  loadPresets<Li>(molecule);
  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

//singleSlater.setRef(RHF);
//singleSlater.isClosedShell = true;
//singleSlater.setRef(UHF);
  singleSlater.setRef(TCS);
  singleSlater.setNTCS(2);
  singleSlater.isClosedShell = false;

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

  basis.findBasisFile("sto-3g");
//basis.findBasisFile("3-21g");
//basis.findBasisFile("6-31G");
  basis.communicate(fileio);
  basis.parseGlobal();
  basis.constructLocal(&molecule);
  basis.makeMaps(&molecule);
  basis.renormShells();


  aoints.communicate(molecule,basis,fileio,memManager);
  singleSlater.communicate(molecule,basis,aoints,fileio,memManager);
//moints.communicate(molecule,basis,fileio,aoints,singleSlater);

  aoints.initMeta();
//aoints.integralAlgorithm = AOIntegrals::INCORE;
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

  for(auto iF : singleSlater.fockOrtho()) scatteredFock_.emplace_back(*iF);
  for(auto iF : singleSlater.onePDMOrtho()) scatteredDen_.emplace_back(*iF);

  Quantum<dcomplex>::spinGather(FMap,scatteredFock_);  
  Quantum<dcomplex>::spinGather(PMap,scatteredDen_);  


  ComplexMatrix FP = F*P;
  ComplexMap    FPMap(FP.data(),FP.rows(),FP.cols());
  ComplexMatrix FPScalar(FP.rows()/2,FP.cols()/2);
  ComplexMatrix FPMz(FP.rows()/2,FP.cols()/2);
  ComplexMatrix FPMy(FP.rows()/2,FP.cols()/2);
  ComplexMatrix FPMx(FP.rows()/2,FP.cols()/2);

  ComplexMap FPScalarMap(FPScalar.data(),FPScalar.rows(),FPScalar.cols());
  ComplexMap FPMzMap(FPMz.data(),FPMz.rows(),FPMz.cols());
  ComplexMap FPMyMap(FPMy.data(),FPMy.rows(),FPMy.cols());
  ComplexMap FPMxMap(FPMx.data(),FPMx.rows(),FPMx.cols());

  std::vector<std::reference_wrapper<ComplexMap>> scatteredFP_;
  scatteredFP_.emplace_back(FPScalarMap);
  scatteredFP_.emplace_back(FPMzMap);
  scatteredFP_.emplace_back(FPMyMap);
  scatteredFP_.emplace_back(FPMxMap);

  Quantum<dcomplex>::spinScatter(FPMap,scatteredFP_);


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
  rt.communicate(singleSlater);
  rt.alloc();
//rt.setMaxSteps(827000); // roughly 1 ps of dynamics
  rt.setMaxSteps(1000); 
  rt.setTOff(0.10);
  rt.setEDFieldAmp({0.001,0.0,0.0});
  rt.setIEnvlp(Step);
  rt.doPropagation();

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

