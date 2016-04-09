#define EIGEN_RUNTIME_NO_MALLOC
#include <response.h>
#include <workers.h>
#include <pythonapi.h>
#include <grid2.h>

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

  EulerMac G(300000);
  double f = G.integrate<double>(gaussian); // exp(-x^2)
  double g = G.integrate<double>(sphGaussian); // exp(-r^2) * r^2 over r
  auto X = G.integrate<RealMatrix>(mat);

  Lebedev L(302);
  cout << L.npts() << endl;;

  double h = L.integrate<double>(sphHarmonic);
  cout << h << endl;

  cout.precision(10);
  std::cout <<  f << endl;
  std::cout <<  g << endl;
  cout << endl;
  std::cout <<  X << endl;

  TwoDGrid2 Sphere(75,302,EULERMAC,LEBEDEV);

  cout << 4*math.pi*Sphere.integrate<double>(sphGaussian) << endl;

  auto sphGaussian2 = [&](IntegrationPoint pt) -> double {
    double x = bg::get<0>(pt.pt);
    double y = bg::get<1>(pt.pt);
    double z = bg::get<2>(pt.pt);
    double r = std::sqrt(x*x + y*y + z*z);
    return pt.weight *  std::exp(-r*r);
  };

  auto sphGaussian3 = [&](IntegrationPoint pt, double &result){
    double x = bg::get<0>(pt.pt);
    double y = bg::get<1>(pt.pt);
    double z = bg::get<2>(pt.pt);
    double r = std::sqrt(x*x + y*y + z*z);
    result +=  pt.weight *  std::exp(-r*r);
  };

  double res, res2(1.0);
  Sphere.integrate<double>(sphGaussian2,res);
  Sphere.integrate<double>(sphGaussian3,res2);
  cout << 4 * math.pi * res << endl;
  cout << 4 * math.pi * (res2 - 1.0) << endl;

  RealMatrix scr(2,2);
  auto noallocMatIntegrate = [&](IntegrationPoint pt, RealMatrix &result){
    double x = bg::get<0>(pt.pt);
    double y = bg::get<1>(pt.pt);
    double z = bg::get<2>(pt.pt);
    double r = std::sqrt(x*x + y*y + z*z);
    scr.setZero();
    scr(0,0) = pt.weight * r * r * std::exp(-r*r);
    scr(1,1) = pt.weight * r * r * std::exp(-r*r);

    result += scr;
  };

  RealMatrix resMat(2,2);
  Eigen::internal::set_is_malloc_allowed(false);
  Sphere.integrate<RealMatrix>(noallocMatIntegrate,resMat);
  Eigen::internal::set_is_malloc_allowed(true);
  cout << 4 * math.pi * resMat << endl;

  std::vector<std::array<double,3>> empty;

  AtomicGrid SphereAtom(75,302,EULERMAC,LEBEDEV,{0.0,0.0,0.0},BECKE,
      empty);
  cout << 4 * math.pi * SphereAtom.integrate<double>(sphGaussian) << endl;

  AtomicGrid SphereAtom2(75,302,EULERMAC,LEBEDEV,{2.0,0.0,0.0},BECKE,
      empty);

  SphereAtom.printGrid(cout);
  SphereAtom2.printGrid(cout);


  // Test Molecular integration


  std::array<double,3> center1 = {0.0,0.0,0.0};
  std::array<double,3> center2 = {1.0,0.0,0.0};
  std::array<double,3> center3 = {0.0,2.0,0.0};

  std::vector<std::array<double,3>> other1;
  std::vector<std::array<double,3>> other2;
  std::vector<std::array<double,3>> other3;

  other1.push_back(center2);
  other1.push_back(center3);
  other2.push_back(center1);
  other2.push_back(center3);
  other3.push_back(center1);
  other3.push_back(center2);

  AtomicGrid Center1(75,302,EULERMAC,LEBEDEV,{0.0,0.0,0.0},BECKE,
      other1);
  AtomicGrid Center2(75,302,EULERMAC,LEBEDEV,{1.0,0.0,0.0},BECKE,
      other2);
  AtomicGrid Center3(75,302,EULERMAC,LEBEDEV,{0.0,2.0,0.0},BECKE,
      other3);

  cout << "HERE" << endl;
  cout << 4 * math.pi * Center1.integrate<double>(sphGaussian) << endl;
  cout << 4 * math.pi * Center2.integrate<double>(sphGaussian) << endl;
  cout << 4 * math.pi * Center3.integrate<double>(sphGaussian) << endl;
  cout << 4 * math.pi * (Center1.integrate<double>(sphGaussian) +
     Center2.integrate<double>(sphGaussian) +
     Center3.integrate<double>(sphGaussian) ) << endl;
  return 0;
};

