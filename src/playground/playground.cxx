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
  singleSlater.setCorrKernel(SingleSlater<double>::CORR::NOCORR);
  singleSlater.setPrintLevel(5);

//  basis.findBasisFile("sto3g");
  basis.findBasisFile("6-31g");
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

  RealMatrix SCRATCH2(singleSlater.nBasis(),singleSlater.nBasis());
  VectorXd   SCRATCH1(singleSlater.nBasis());
  std::chrono::duration<double> T1; 
  std::chrono::duration<double> T2; 
  std::chrono::duration<double> T3; 
//  basis.radcut(1e-10,50,1e-7);
  auto density = [&](IntegrationPoint pt, double &result) {
    // Evaluate the basis product in SCRATCH
    SCRATCH1.setZero();
    auto t1s = std::chrono::high_resolution_clock::now();
    cartGP GP = pt.pt;
//    auto shMap = basis.MapGridBasis(GP); 
//    if(shMap[0]) { 
//       cout << "Skip all pts " <<endl;
//       return 0.0;}
    for(auto iShell = 0; iShell < basis.nShell(); iShell++){
//      if(!shMap[iShell+1]) {
//       cout << basis.getradCutSh(iShell) << endl;
//        continue;}
      int b_s = basis.mapSh2Bf(iShell);
      int size= basis.shells(iShell).size();

      libint2::Shell shTmp = basis.shells(iShell);
      double * buff = basis.basisDEval(0,shTmp,&pt.pt);
      RealMap bMap(buff,size,1);
      SCRATCH1.block(b_s,0,size,1) = bMap;

      delete [] buff;
    };
    auto t2s = std::chrono::high_resolution_clock::now();

    if(SCRATCH1.norm() < 1e-8) return 0.0;
    SCRATCH2 = SCRATCH1 * SCRATCH1.transpose();
   
    auto t3s = std::chrono::high_resolution_clock::now();
    result += pt.weight * singleSlater.computeProperty<double,TOTAL>(SCRATCH2); 
    auto t3f = std::chrono::high_resolution_clock::now();

    T1 += t2s - t1s;
    T2 += t3s - t2s;
    T3 += t3f - t3s;
  };


  auto valVxc = [&](IntegrationPoint pt, MyStruct &result) {
    // Evaluate the basis product in SCRATCH
    SCRATCH1.setZero();
    cartGP GP = pt.pt;
    double rhoA;
    double rhoB;
    DFTFunctional::DFTInfo kernelXC;
//    SlaterExchange * dftFun;
//    auto shMap = basis.MapGridBasis(GP); 
//    if(shMap[0]) { 
//       cout << "Skip all pts " <<endl;
//       return 0.0;}
    for(auto iShell = 0; iShell < basis.nShell(); iShell++){
//      if(!shMap[iShell+1]) {
//       cout << basis.getradCutSh(iShell) << endl;
//        continue;}
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
    kernelXC = singleSlater.dftFunctionals_[0]->eval(rhoA, rhoB);
    result.VXCA   += pt.weight * SCRATCH2 * kernelXC.ddrhoA; 
    result.VXCB   += pt.weight * SCRATCH2 * kernelXC.ddrhoB; 
    result.Energy += pt.weight * (rhoA+rhoB) * kernelXC.eps;
  };

  auto numOverlap = [&](IntegrationPoint pt, RealMatrix &result) {
    // Evaluate the basis product in SCRATCH
    for(auto iShell = 0; iShell < basis.nShell(); iShell++){
      int b_s = basis.mapSh2Bf(iShell);
      int size= basis.shells(iShell).size();

      libint2::Shell shTmp = basis.shells(iShell);
      double * buff = basis.basisDEval(0,shTmp,
          &pt.pt);
      RealMap bMap(buff,size,1);
      SCRATCH1.block(b_s,0,size,1) = bMap;

      delete [] buff;
    };

    SCRATCH2 = SCRATCH1 * SCRATCH1.transpose();
    result += pt.weight * SCRATCH2 ;
  };

  std::vector<std::array<double,3>> atomicCenters;

  for(auto iAtm = 0; iAtm < molecule.nAtoms(); iAtm++){
    atomicCenters.push_back(
        {(*molecule.cart())(0,iAtm),
         (*molecule.cart())(1,iAtm),
         (*molecule.cart())(2,iAtm)}
    );
  //cout << atomicCenters[iAtm][0] << "\t";
  //cout << atomicCenters[iAtm][1] << "\t";
  //cout << atomicCenters[iAtm][2] << "\t";
  };


  double rho = 0;
  RealMatrix NS(singleSlater.nBasis(),singleSlater.nBasis());
  NS.setZero();
//  coeff = singleSlater.dftFunctionals_[0]->getCxVx();
  auto t4s = std::chrono::high_resolution_clock::now();
  AtomicGrid AGrid(100,590,GAUSSCHEBFST,LEBEDEV,BECKE,atomicCenters,0,1.0,
      false);

  for(auto iAtm = 0; iAtm < molecule.nAtoms(); iAtm++){
    AGrid.center() = iAtm;
    AGrid.scalingFactor()=0.5*elements[molecule.index(iAtm)].sradius/phys.bohr;
    AGrid.integrate<double>(density,rho);
  };
  auto t4f = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> T4 = t4f - t4s;

  //prettyPrint(cout,4*math.pi*NS,"NS");
  //prettyPrint(cout,*aoints.overlap_,"S");
  cout << "RHO " << 4*math.pi*rho << endl;
  cout << "T1 " << T1.count() << endl;
  cout << "T2 " << T2.count() << endl;
  cout << "T3 " << T3.count() << endl;
  cout << "T4 " << T4.count() << endl;

  MyStruct res(singleSlater.nBasis());
  for(auto iAtm = 0; iAtm < molecule.nAtoms(); iAtm++){
    AGrid.center() = iAtm;
    AGrid.scalingFactor()=0.5*elements[molecule.index(iAtm)].sradius/phys.bohr;
    AGrid.integrate<MyStruct>(valVxc,res);
  };

  prettyPrint(cout,4*math.pi*res.VXCA,"A");
  prettyPrint(cout,4*math.pi*res.VXCB,"B");
  cout << "ENERGY " << 4*math.pi*res.Energy  << endl;


  finalizeCQ();
  return 0;
};

