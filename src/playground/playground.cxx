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
    return r * r * std::exp(-r*r);
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
  double f = G.integrate<double>(gaussian);
  double g = G.integrate<double>(sphGaussian);
  auto X = G.integrate<RealMatrix>(mat);

  Lebedev L(302);
  cout << L.npts() << endl;;

  double one(0.0);
  for(auto i = 0; i < 302; i++){
    cout << std::setw(22) << std::setprecision(10) << bg::get<0>(L[i].pt);
    cout << std::setw(22) << std::setprecision(10) << bg::get<1>(L[i].pt);
    cout << std::setw(22) << std::setprecision(10) << bg::get<2>(L[i].pt);
    cout << std::setw(22) << std::setprecision(10) << L[i].weight;
    cout << endl;
    one += L[i].weight;
  };
  cout << one << endl;


  double h = L.integrate<double>(sphHarmonic);
  cout << h << endl;

  cout.precision(10);
  std::cout <<  f << endl;
  std::cout <<  g << endl;
  cout << endl;
  std::cout <<  X << endl;

  TwoDGrid2 Sphere(75,302,EULERMAC,LEBEDEV);

  cout << 4*math.pi*Sphere.integrate<double>(sphGaussian) << endl;
  return 0;
};

