#include <response.h>
#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;

enum MOLECULE_PRESETS {
  WATER
};

template<MOLECULE_PRESETS T>
void loadPresets(Molecule&);

template<>
void loadPresets<WATER>(Molecule &mol) {
  mol.setNAtoms(3);
  mol.setCharge(2);
  mol.setNTotalE(8);
  mol.setMultip(1);
  mol.alloc();
  mol.setIndex(0,HashAtom("O",0));
  mol.setIndex(1,HashAtom("H",0));
  mol.setIndex(2,HashAtom("H",0));
  mol.setCart(0,0.000000000 ,-0.07579184359, 0.0);
  mol.setCart(1,0.866811829 ,0.6014357793  ,0.0);
  mol.setCart(2,-0.866811829, 0.6014357793 ,0.0);
};

int main(int argc, char **argv){
  Molecule molecule;
  BasisSet basis;
  Controls controls;
  AOIntegrals aoints;
  MOIntegrals<double> moints;
  SingleSlater<double> singleSlater;
  Response<double> resp;
  FileIO fileio("test.inp","test.out");

  initCQ(argc,argv);
  controls.iniControls();
  fileio.iniH5Files();
  fileio.iniStdGroups();
  CQSetNumThreads(1);
  
  loadPresets<WATER>(molecule);
  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(SingleSlater<double>::RHF);
  singleSlater.isClosedShell = true;

  basis.findBasisFile("sto3g");
  basis.communicate(fileio);
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
  //resp.doSA();
  resp.setNSek(3);
  resp.doFull();
  resp.doResponse();

  RealMatrix T = resp.transDen<AAB_PPTDA>();

  RealMatrix TAO(singleSlater.nBasis(),singleSlater.nBasis());
  RealMatrix TAO2(singleSlater.nBasis(),singleSlater.nBasis());
  RealMatrix TMO(singleSlater.nBasis(),singleSlater.nBasis());

  auto delta = [](int i, int j) -> int {
    if(i == j) return 1;
    else       return 0;
  };

  for(auto iSt = 0; iSt < 8 ; iSt++){
     // This builds the total difference density in the AO basis
     for(auto mu = 0; mu < singleSlater.nBasis(); mu++)
     for(auto nu = 0; nu < singleSlater.nBasis(); nu++)
     for(auto a = 0, ab = 0; a < resp.nVA(); a++      )
     for(auto b = 0        ; b < resp.nVB(); b++, ab++)
     for(auto c = 0, cd = 0; c < resp.nVA(); c++      )
     for(auto d = 0        ; d < resp.nVB(); d++, cd++){
       TAO2(mu,nu) += T.col(iSt)(ab) * T.col(iSt)(cd) * ( 
          (*singleSlater.moA())(mu,resp.nOA() + a) * 
          (*singleSlater.moA())(nu,resp.nOA() + c) * delta(b,d)
         + 
          (*singleSlater.moA())(mu,resp.nOA() + b) * 
          (*singleSlater.moA())(nu,resp.nOA() + d) * delta(a,c)
       ); 
     };
   
     // This builds the (alpha) difference density in the MO basis
     for(auto a = 0; a < resp.nVA(); a++      )
     for(auto b = 0        ; b < resp.nVB(); b++)
     for(auto c = 0, cd = 0; c < resp.nVA(); c++      )
     for(auto d = 0        ; d < resp.nVB(); d++, cd++)
     for(auto e = 0, ef = 0; e < resp.nVA(); e++      )
     for(auto f = 0        ; f < resp.nVB(); f++, ef++){
       TMO(resp.nOA() + a, resp.nOB() + b) += 
         T.col(iSt)(cd) * T.col(iSt)(ef) *
         delta(b,f) * delta(c,e) * delta(d,a);
     }
   
     // This transforms the (alpha) difference density into the AO basis
     TAO = (*singleSlater.moA()) * TMO * singleSlater.moA()->adjoint();

     cout <<TAO2.cwiseQuotient(TAO) << endl;
   
   //prettyPrint(cout,TMO,"TMO");
   //prettyPrint(cout,TAO,"TAO");
   //prettyPrint(cout,TAO2,"TAO2");
   //cout << T.col(iSt).dot(T.col(iSt)) << endl;
   //cout << (TAO * (*aoints.overlap_)).trace() << endl;
   //cout << (TAO2 * (*aoints.overlap_)).trace() << endl;
   
     RealMatrix XAO(singleSlater.nBasis(),singleSlater.nBasis());
     for(auto mu = 0; mu < singleSlater.nBasis(); mu++)
     for(auto nu = 0; nu < singleSlater.nBasis(); nu++)
     for(auto a = 0, ab = 0; a < resp.nVA(); a++      )
     for(auto b = 0        ; b < resp.nVB(); b++, ab++){
       XAO(mu,nu) += T.col(iSt)(ab) * 
         (*singleSlater.moA())(mu,resp.nOA() + a) *  
         (*singleSlater.moA())(nu,resp.nOA() + b); 
     }
   
   //prettyPrint(cout,T.col(iSt),"X");
   //prettyPrint(cout,XAO,"XAO");
     /*
     double TPDM = (TAO2 * (*aoints.overlap_)).trace() *
     0.25*((*singleSlater.densityA()) * (*aoints.overlap_)).trace()
     + 
     0.25*(TAO2 * (*aoints.overlap_) * (*singleSlater.densityA()) * (*aoints.overlap_))
        .trace() 
     +
     0.0*(XAO * (*aoints.overlap_) * XAO.adjoint() * (*aoints.overlap_)).trace();
     */
   
    /*
     double TPDM = -2.0* (TAO * (*aoints.overlap_) * (*singleSlater.densityA()) * (*aoints.overlap_)).trace() 
       + 2.0 * (XAO * (*aoints.overlap_) * XAO.adjoint() * (*aoints.overlap_)).trace();
       */
     double TPDM = -0.5* (TAO2 * (*aoints.overlap_) * (*singleSlater.densityA()) * (*aoints.overlap_)).trace() 
       - 2.0 * (XAO * (*aoints.overlap_) * XAO.adjoint() * (*aoints.overlap_)).trace();
   
     double OPDM = 2;
   
     cout << (3.0*OPDM + TPDM)/4 << endl;;
   //cout << OPDM << endl;;
   //cout << TPDM << endl;;
  }

   
  finalizeCQ(); 
  return 0;
};

