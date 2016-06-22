#define EIGEN_RUNTIME_NO_MALLOC
#include <response.h>
#include <workers.h>
#include <pythonapi.h>
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
  auto gaussian = [&](cartGP pt) -> double {
    double r = bg::get<0>(pt);
    return std::exp(-r*r);
  };
  auto sphGaussian = [&](cartGP pt) -> double {
    double x = bg::get<0>(pt);
    double y = bg::get<1>(pt);
    double z = bg::get<2>(pt);
    double r = std::sqrt(x*x + y*y + z*z);
    return std::exp(-r*r);
  };

  auto mat = [&](cartGP pt) -> RealMatrix {
    RealMatrix A(2,2);
    A(0,0) = gaussian(pt);
    A(1,1) = sphGaussian(pt);
    return A;
  };

  auto sphere = [&](cartGP pt) -> double {
    double x = bg::get<0>(pt);
    double y = bg::get<1>(pt);
    double z = bg::get<2>(pt);

    return x*x + y*y + z*z;
  };

  auto sphHarmonic = [&](cartGP pt) -> double {
    double x = bg::get<0>(pt);
    double y = bg::get<1>(pt);
    double z = bg::get<2>(pt);
    double r = std::sqrt(x*x + y*y + z*z);
    double ct = z / r;
//    return (3.0 * z*z) / (4 * math.pi * r*r);
//    cout << ct << endl;
    return 3*ct*ct;
  };

//  EulerMac G(300000);
//  double f = G.integrate<double>(gaussian); // exp(-x^2)
//  double g = G.integrate<double>(sphGaussian); // exp(-r^2) * r^2 over r
//  auto X = G.integrate<RealMatrix>(mat);
//
//  Lebedev L(302);
//  cout << L.npts() << endl;;
//
//  double h = L.integrate<double>(sphHarmonic);
//  cout << h << endl;
//
//  cout.precision(10);
//  std::cout <<  f << endl;
//  std::cout <<  g << endl;
//  cout << endl;
//  std::cout <<  X << endl;
//
//  TwoDGrid2 Sphere(75,302,EULERMAC,LEBEDEV);
//
//  cout << 4*math.pi*Sphere.integrate<double>(sphGaussian) << endl;
//
//  auto sphGaussian2 = [&](IntegrationPoint pt) -> double {
//    double x = bg::get<0>(pt.pt);
//    double y = bg::get<1>(pt.pt);
//    double z = bg::get<2>(pt.pt);
//    double r = std::sqrt(x*x + y*y + z*z);
//    return pt.weight *  std::exp(-r*r);
//  };
//
//  auto sphGaussian3 = [&](IntegrationPoint pt, double &result){
//    double x = bg::get<0>(pt.pt);
//    double y = bg::get<1>(pt.pt);
//    double z = bg::get<2>(pt.pt);
//    double r = std::sqrt(x*x + y*y + z*z);
//    result +=  pt.weight *  std::exp(-r*r);
//  };
//
//  double res, res2(1.0);
//  Sphere.integrate<double>(sphGaussian2,res);
//  Sphere.integrate<double>(sphGaussian3,res2);
//  cout << 4 * math.pi * res << endl;
//  cout << 4 * math.pi * (res2 - 1.0) << endl;
//
//  RealMatrix scr(2,2);
//  auto noallocMatIntegrate = [&](IntegrationPoint pt, RealMatrix &result){
//    double x = bg::get<0>(pt.pt);
//    double y = bg::get<1>(pt.pt);
//    double z = bg::get<2>(pt.pt);
//    double r = std::sqrt(x*x + y*y + z*z);
//    scr.setZero();
//    scr(0,0) = pt.weight * r * r * std::exp(-r*r);
//    scr(1,1) = pt.weight * r * r * std::exp(-r*r);
//
//    result += scr;
//  };
//
//  RealMatrix resMat(2,2);
//  Eigen::internal::set_is_malloc_allowed(false);
//  Sphere.integrate<RealMatrix>(noallocMatIntegrate,resMat);
//  Eigen::internal::set_is_malloc_allowed(true);
//  cout << 4 * math.pi * resMat << endl;
//
//  std::vector<std::array<double,3>> empty;
//
////AtomicGrid SphereAtom(75,302,EULERMAC,LEBEDEV,{0.0,0.0,0.0},BECKE,
////    empty);
////cout << 4 * math.pi * SphereAtom.integrate<double>(sphGaussian) << endl;
//
////AtomicGrid SphereAtom2(75,302,EULERMAC,LEBEDEV,{2.0,0.0,0.0},BECKE,
////    empty);
//
////SphereAtom.printGrid(cout);
////SphereAtom2.printGrid(cout);
//
//
//  // Test Molecular integration
//
//
//  std::array<double,3> center1 = {0.0,0.0,0.0};
//  std::array<double,3> center2 = {1.0,0.0,0.0};
//  std::array<double,3> center3 = {0.0,2.0,0.0};
//
//  std::vector<std::array<double,3>> other1;
//  std::vector<std::array<double,3>> other2;
//  std::vector<std::array<double,3>> other3;
//
//  other1.push_back(center2);
//  other1.push_back(center3);
//  other2.push_back(center1);
//  other2.push_back(center3);
//  other3.push_back(center1);
//  other3.push_back(center2);
//
////AtomicGrid Center1(75,590,EULERMAC,LEBEDEV,{0.0,0.0,0.0},BECKE,
////    other1);
////AtomicGrid Center2(75,590,EULERMAC,LEBEDEV,{1.0,0.0,0.0},BECKE,
////    other2);
////AtomicGrid Center3(75,590,EULERMAC,LEBEDEV,{0.0,2.0,0.0},BECKE,
////    other3);
//
////cout << "HERE" << endl;
////cout << 4 * math.pi * Center1.integrate<double>(sphGaussian) << endl;
////cout << 4 * math.pi * Center2.integrate<double>(sphGaussian) << endl;
////cout << 4 * math.pi * Center3.integrate<double>(sphGaussian) << endl;
////cout << 4 * math.pi * (Center1.integrate<double>(sphGaussian) +
////   Center2.integrate<double>(sphGaussian) +
////   Center3.integrate<double>(sphGaussian) ) << endl;
//
//  std::vector<std::array<double,3>> centers;
//  centers.push_back(center1);
//  centers.push_back(center2);
//  centers.push_back(center3);
//
//  AtomicGrid NCenter1(75,590,EULERMAC,LEBEDEV,BECKE,centers,0);
//  AtomicGrid NCenter2(75,590,EULERMAC,LEBEDEV,BECKE,centers,1);
//  AtomicGrid NCenter3(75,590,EULERMAC,LEBEDEV,BECKE,centers,2);
//
//  cout << "HERE" << endl;
//  cout << 4 * math.pi * NCenter1.integrate<double>(sphGaussian) << endl;
//  cout << 4 * math.pi * NCenter2.integrate<double>(sphGaussian) << endl;
//  cout << 4 * math.pi * NCenter3.integrate<double>(sphGaussian) << endl;
//  cout << 4 * math.pi * (NCenter1.integrate<double>(sphGaussian) +
//     NCenter2.integrate<double>(sphGaussian) +
//     NCenter3.integrate<double>(sphGaussian) ) << endl;



  CQMemManager memManager;
  Molecule molecule;
  BasisSet basis;
  AOIntegrals aoints;
  MOIntegrals<double> moints;
  SingleSlater<double> singleSlater;
  Response<double> resp;
  FileIO fileio("test.inp","test.out");

  memManager.setTotalMem(256e6);
  initCQ(argc,argv);
  fileio.iniH5Files();
  fileio.iniStdGroups();
  CQSetNumThreads(1);
  
  loadPresets<Li>(molecule);
  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(SingleSlater<double>::UHF);
  singleSlater.isClosedShell = false;

  basis.findBasisFile("sto3g");
//  basis.findBasisFile("6-31g");
  basis.communicate(fileio);
  basis.parseGlobal();
  basis.constructLocal(&molecule);
  basis.makeMaps(&molecule);
  basis.renormShells();


  aoints.communicate(molecule,basis,fileio,memManager);
  singleSlater.communicate(molecule,basis,aoints,fileio,memManager);
  moints.communicate(molecule,basis,fileio,aoints,singleSlater);

  aoints.initMeta();
  aoints.integralAlgorithm = AOIntegrals::INCORE;
  singleSlater.initMeta();
  singleSlater.genMethString();

  aoints.alloc();
  singleSlater.alloc();

  singleSlater.formGuess();
  singleSlater.formFock();
  singleSlater.computeEnergy();
  singleSlater.SCF2();
  singleSlater.computeProperties();
  singleSlater.printProperties();

  SingleSlater<double> newSS;
  newSS.setRef(SingleSlater<double>::TCS);
  newSS.isClosedShell = false;
  newSS.setNTCS(2);
  newSS.communicate(molecule,basis,aoints,fileio,memManager);
  newSS.initMeta();
  newSS.genMethString();
  newSS.alloc();

  Eigen::Map<RealMatrix,0,Eigen::Stride<Dynamic,Dynamic> >
    CA(newSS.moA()->data(),basis.nBasis(),2*basis.nBasis(),
        Eigen::Stride<Dynamic,Dynamic>(2*basis.nBasis(),2));
  Eigen::Map<RealMatrix,0,Eigen::Stride<Dynamic,Dynamic> >
    CB(newSS.moA()->data()+1,basis.nBasis(),2*basis.nBasis(),
        Eigen::Stride<Dynamic,Dynamic>(2*basis.nBasis(),2));

  (*newSS.onePDMScalar()) = (*singleSlater.onePDMScalar());
  (*newSS.onePDMMz()) = (*singleSlater.onePDMMz());

//(*newSS.onePDMMx()) = (*newSS.onePDMMz());
//newSS.onePDMMz()->setZero();

  newSS.gatherDensity();
//prettyPrintTCS(cout,*newSS.onePDMA(),"PA");
//prettyPrint(cout,*singleSlater.onePDMA(),"PA");
//prettyPrint(cout,*singleSlater.onePDMB(),"PB");

  newSS.haveMO = true;
  newSS.haveDensity = true;
  newSS.formFock();
//prettyPrintTCS(cout,*newSS.fockA(),"FA");
//prettyPrint(cout,*singleSlater.fockA(),"FA");
//prettyPrint(cout,*singleSlater.fockB(),"FB");
  
  newSS.computeEnergy();
  cout << "E = " << newSS.totalEnergy << endl;

  RealMatrix TMPAA =  
    (*singleSlater.fockA()) * (*aoints.overlap_) * (*singleSlater.onePDMA());
  RealMatrix TMPBB =  
    (*singleSlater.fockB()) * (*aoints.overlap_) * (*singleSlater.onePDMB());

  RealMatrix TMPSS = 
    (*newSS.fockScalar()) * (*aoints.overlap_) * (*newSS.onePDMScalar());
  RealMatrix TMPZZ = 
    (*newSS.fockMz()) * (*aoints.overlap_) * (*newSS.onePDMMz());
  RealMatrix TMPYY = // Need -?
    (*newSS.fockMy()) * (*aoints.overlap_) * (*newSS.onePDMMy());
  RealMatrix TMPXX = 
    (*newSS.fockMx()) * (*aoints.overlap_) * (*newSS.onePDMMx());

  RealMatrix TMPSZ = 
    (*newSS.fockScalar()) * (*aoints.overlap_) * (*newSS.onePDMMz());
  RealMatrix TMPSY = 
    (*newSS.fockScalar()) * (*aoints.overlap_) * (*newSS.onePDMMy());
  RealMatrix TMPSX = 
    (*newSS.fockScalar()) * (*aoints.overlap_) * (*newSS.onePDMMx());

  RealMatrix TMPZS = 
    (*newSS.fockMz()) * (*aoints.overlap_) * (*newSS.onePDMScalar());
  RealMatrix TMPYS = 
    (*newSS.fockMy()) * (*aoints.overlap_) * (*newSS.onePDMScalar());
  RealMatrix TMPXS = 
    (*newSS.fockMx()) * (*aoints.overlap_) * (*newSS.onePDMScalar());

  RealMatrix TMPXZ = 
    (*newSS.fockMx()) * (*aoints.overlap_) * (*newSS.onePDMMz());
  RealMatrix TMPZX = 
    (*newSS.fockMz()) * (*aoints.overlap_) * (*newSS.onePDMMx());

  RealMatrix TMPYZ = 
    (*newSS.fockMy()) * (*aoints.overlap_) * (*newSS.onePDMMz());
  RealMatrix TMPZY = 
    (*newSS.fockMz()) * (*aoints.overlap_) * (*newSS.onePDMMy());

  RealMatrix TMPXY = 
    (*newSS.fockMx()) * (*aoints.overlap_) * (*newSS.onePDMMy());
  RealMatrix TMPYX = 
    (*newSS.fockMy()) * (*aoints.overlap_) * (*newSS.onePDMMx());

  RealMatrix CommI = TMPSS + TMPZZ + TMPXX + TMPYY;
  RealMatrix CommZ = TMPSZ + TMPZS - TMPXY + TMPYX;
  RealMatrix CommX = TMPSX + TMPXS - TMPYZ + TMPZY;
  RealMatrix CommY = TMPSY + TMPYS + TMPZX - TMPXZ;

  cout << " TMPSS   " <<  TMPSS.norm() << endl; 
  cout << " TMPZZ   " <<  TMPZZ.norm() << endl; 
  cout << " TMPYY   " <<  TMPYY.norm() << endl; 
  cout << " TMPXX   " <<  TMPXX.norm() << endl; 
  cout << " TMPSZ   " <<  TMPSZ.norm() << endl; 
  cout << " TMPSY   " <<  TMPSY.norm() << endl; 
  cout << " TMPSX   " <<  TMPSX.norm() << endl; 
  cout << " TMPZS   " <<  TMPZS.norm() << endl; 
  cout << " TMPYS   " <<  TMPYS.norm() << endl; 
  cout << " TMPXS   " <<  TMPXS.norm() << endl; 
  cout << " TMPXZ   " <<  TMPXZ.norm() << endl; 
  cout << " TMPZX   " <<  TMPZX.norm() << endl; 
  cout << " TMPYZ   " <<  TMPYZ.norm() << endl; 
  cout << " TMPZY   " <<  TMPZY.norm() << endl; 
  cout << " TMPXY   " <<  TMPXY.norm() << endl; 
  cout << " TMPYX   " <<  TMPYX.norm() << endl; 

  RealMatrix CommTotal(*newSS.onePDMA());
  CommTotal.setZero();

//CommZ = CommX;
//CommX.setZero();
  std::vector<std::reference_wrapper<RealMatrix>> comComp;
  comComp.emplace_back(CommI);
  comComp.emplace_back(CommZ);
  comComp.emplace_back(CommY);
  comComp.emplace_back(CommX);

  Quantum<double>::spinGather(CommTotal,comComp);

//prettyPrint(cout,TMPAA,"COMMAA");
//prettyPrint(cout,TMPBB,"COMMBB");
//prettyPrintTCS(cout,0.5*CommTotal,"Comm");
//prettyPrint(cout,CommI,"I");
//prettyPrint(cout,CommX,"X");
//prettyPrint(cout,CommY,"Y");
//prettyPrint(cout,CommZ,"Z");

  finalizeCQ();
  return 0;
};

