#include <response.h>
#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;
int main(){
  Molecule molecule;
  BasisSet basis;
  Controls controls;
  AOIntegrals aoints;
  MOIntegrals<double> moints;
  SingleSlater<double> singleSlater;
  Response<double> resp;
  FileIO fileio("test.inp","test.out");

  initCQ();
  controls.iniControls();
  fileio.iniH5Files();
  fileio.iniStdGroups();
  CQSetNumThreads(1);

  molecule.setCharge(0);
  molecule.setMultip(1);
  //molecule.setNAtoms(3);
  molecule.setNAtoms(2);
  molecule.alloc(fileio.out);

/*
  // Molecule Specification for Water
  molecule.setIndex(0,HashAtom("O",0));
  molecule.setIndex(1,HashAtom("H",0));
  molecule.setIndex(2,HashAtom("H",0));
  molecule.setCart(0,0.000000000 ,-0.07579184359, 0.0);
  molecule.setCart(1,0.866811829 ,0.6014357793  ,0.0);
  molecule.setCart(2,-0.866811829, 0.6014357793 ,0.0);
  molecule.setNTotalE(10);
*/

  // Molecule Specification ofr BH
  molecule.setIndex(0,HashAtom("B",0));
  molecule.setIndex(1,HashAtom("H",0));
  molecule.setCart(0,0.0,0.0,0.0);
  molecule.setCart(1,0.0,0.0,1.232);
  molecule.setNTotalE(6);

  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(SingleSlater<double>::RHF);
  singleSlater.isClosedShell = true;

  basis.communicate(fileio);
  //basis.findBasisFile("sto3g");
  basis.findBasisFile("3-21G");
  basis.parseGlobal();
  basis.constructLocal(&molecule);
  basis.makeMaps(1,&molecule);
  basis.renormShells();

  aoints.communicate(molecule,basis,fileio,controls);
  singleSlater.communicate(molecule,basis,aoints,fileio,controls);
  moints.communicate(molecule,basis,fileio,controls,aoints,singleSlater);

  aoints.initMeta();
  aoints.integralAlgorithm = AOIntegrals::INCORE;
  singleSlater.initMeta();
  singleSlater.genMethString();

  aoints.alloc();
  singleSlater.alloc();

  singleSlater.formGuess();
  singleSlater.formFock();
  singleSlater.computeEnergy();
  singleSlater.SCF();
  singleSlater.computeProperties();
  singleSlater.printProperties();
  
  moints.communicate(molecule,basis,fileio,controls,aoints,singleSlater);
  moints.initMeta();
  resp.communicate(singleSlater,moints,fileio);  
  resp.setMeth(RESPONSE_TYPE::PPTDA);
  resp.doSA();
  resp.doResponse();
  resp.setNSek(3);
  finalizeCQ(); 
  return 0;
};
