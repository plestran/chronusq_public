#define EIGEN_RUNTIME_NO_MALLOC
#include <response.h>
#include <workers.h>
#include <pythonapi.h>
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
  WATER, HE,SO
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
//  auto sphGaussian2 = [&](ChronusQ::IntegrationPoint pt) -> double {
//    double x = bg::get<0>(pt.pt);
//    double y = bg::get<1>(pt.pt);
//    double z = bg::get<2>(pt.pt);
//    double r = std::sqrt(x*x + y*y + z*z);
//    return pt.weight *  std::exp(-r*r);
//  };
//
//  auto sphGaussian3 = [&](ChronusQ::IntegrationPoint pt, double &result){
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
//  auto noallocMatIntegrate = [&](ChronusQ::IntegrationPoint pt, RealMatrix &result){
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
////ChronusQ::AtomicGrid SphereAtom(75,302,EULERMAC,LEBEDEV,{0.0,0.0,0.0},BECKE,
////    empty);
////cout << 4 * math.pi * SphereAtom.integrate<double>(sphGaussian) << endl;
//
////ChronusQ::AtomicGrid SphereAtom2(75,302,EULERMAC,LEBEDEV,{2.0,0.0,0.0},BECKE,
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
////ChronusQ::AtomicGrid Center1(75,590,EULERMAC,LEBEDEV,{0.0,0.0,0.0},BECKE,
////    other1);
////ChronusQ::AtomicGrid Center2(75,590,EULERMAC,LEBEDEV,{1.0,0.0,0.0},BECKE,
////    other2);
////ChronusQ::AtomicGrid Center3(75,590,EULERMAC,LEBEDEV,{0.0,2.0,0.0},BECKE,
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
//  ChronusQ::AtomicGrid NCenter1(75,590,EULERMAC,LEBEDEV,BECKE,centers,0);
//  ChronusQ::AtomicGrid NCenter2(75,590,EULERMAC,LEBEDEV,BECKE,centers,1);
//  ChronusQ::AtomicGrid NCenter3(75,590,EULERMAC,LEBEDEV,BECKE,centers,2);
//
//  cout << "HERE" << endl;
//  cout << 4 * math.pi * NCenter1.integrate<double>(sphGaussian) << endl;
//  cout << 4 * math.pi * NCenter2.integrate<double>(sphGaussian) << endl;
//  cout << 4 * math.pi * NCenter3.integrate<double>(sphGaussian) << endl;
//  cout << 4 * math.pi * (NCenter1.integrate<double>(sphGaussian) +
//     NCenter2.integrate<double>(sphGaussian) +
//     NCenter3.integrate<double>(sphGaussian) ) << endl;



  Molecule molecule;
  BasisSet basis;
  AOIntegrals aoints;
  MOIntegrals<double> moints;
  SingleSlater<double> singleSlater;
  Response<double> resp;
  CQMemManager memManager;
  FileIO fileio("test.inp","test.out");

  memManager.setTotalMem(256e6);
  initCQ(argc,argv);
  fileio.iniH5Files();
  fileio.iniStdGroups();
  CQSetNumThreads(1);
  
  loadPresets<WATER>(molecule);
//loadPresets<HE>(molecule);
//  loadPresets<SO>(molecule);
  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(SingleSlater<double>::RHF);
  singleSlater.isClosedShell = true;
  singleSlater.isDFT = true;
  singleSlater.isHF = false;
  singleSlater.setExchKernel(SingleSlater<double>::EXCH::SLATER);
  //singleSlater.setExchKernel(SingleSlater<double>::EXCH::NOEXCH);
//  singleSlater.setCorrKernel(SingleSlater<double>::CORR::NOCORR);
  singleSlater.setCorrKernel(SingleSlater<double>::CORR::VWN3);
  singleSlater.setPrintLevel(5);

  basis.findBasisFile("sto3g");
//  basis.findBasisFile("3-21g");
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
  singleSlater.SCF();
  singleSlater.computeProperties();
  singleSlater.printProperties();

  RealMatrix SCRATCH2(singleSlater.nBasis(),singleSlater.nBasis());
  VectorXd   SCRATCH1(singleSlater.nBasis());
  std::chrono::duration<double> T1; 
  std::chrono::duration<double> T2; 
  std::chrono::duration<double> T3; 

  auto valVxc = [&](ChronusQ::IntegrationPoint pt, MyStruct &result) {
    // Evaluate the basis product in SCRATCH
//    if (pt.weight < singleSlater.dftFunctionals_[0]->epsScreen){
//    return 0.0; }
    SCRATCH1.setZero();
    cartGP GP = pt.pt;
    double rhoA;
    double rhoB;
    auto shMap = basis.MapGridBasis(GP); 
    if(shMap[0]) { 
//       cout << "Skip all pts " <<endl;
       return 0.0;}
    for(auto iShell = 0; iShell < basis.nShell(); iShell++){
      if(!shMap[iShell+1]) {
//       cout << basis.getradCutSh(iShell) << endl;
        continue;}
      int b_s = basis.mapSh2Bf(iShell);
      int size= basis.shells(iShell).size();

      libint2::Shell shTmp = basis.shells(iShell);
      double * buff = basis.basisDEval(0,shTmp,&pt.pt);
      RealMap bMap(buff,size,1);
      SCRATCH1.block(b_s,0,size,1) = bMap;

      delete [] buff;
    };

    if(SCRATCH1.norm() < 1e-8) return 0.0;
    SCRATCH2 = SCRATCH1 * SCRATCH1.transpose();
    rhoA = singleSlater.computeProperty<double,ALPHA>(SCRATCH2);
    rhoB = singleSlater.computeProperty<double,BETA>(SCRATCH2);
    for(auto i = 0; i < singleSlater.dftFunctionals_.size(); i++){
      DFTFunctional::DFTInfo kernelXC = 
        singleSlater.dftFunctionals_[i]->eval(rhoA, rhoB);
      result.VXCA   += pt.weight * SCRATCH2 * kernelXC.ddrhoA; 
      result.VXCB   += pt.weight * SCRATCH2 * kernelXC.ddrhoB; 
      result.Energy += pt.weight * (rhoA+rhoB) * kernelXC.eps;
    }
  };

  std::vector<std::array<double,3>> atomicCenters;

  for(auto iAtm = 0; iAtm < molecule.nAtoms(); iAtm++){
    atomicCenters.push_back(
        {(*molecule.cart())(0,iAtm),
         (*molecule.cart())(1,iAtm),
         (*molecule.cart())(2,iAtm)}
    );
  };


  ChronusQ::AtomicGrid AGrid(100,302,GAUSSCHEBFST,LEBEDEV,BECKE,atomicCenters,0,1.0,
      false);

  MyStruct res(singleSlater.nBasis());
  basis.radcut(1.0e-10, 50, 1.0e-7);
  for(auto iAtm = 0; iAtm < molecule.nAtoms(); iAtm++){
    AGrid.center() = iAtm;
    AGrid.scalingFactor()=0.5*elements[molecule.index(iAtm)].sradius/phys.bohr;
    std::chrono::duration<double> TVEX; 
    auto t1s = std::chrono::high_resolution_clock::now();
    AGrid.integrate<MyStruct>(valVxc,res);
    auto t1f = std::chrono::high_resolution_clock::now();
    TVEX += t1f - t1s;
    cout << "VEX time " << TVEX.count() << endl;
  };
//  double Cx = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));
  prettyPrint(cout,4*math.pi*res.VXCA,"A");
  prettyPrint(cout,4*math.pi*res.VXCB,"B");
  cout << "ENERGY " << 4*math.pi*res.Energy  << endl;


  finalizeCQ();
  return 0;
};

