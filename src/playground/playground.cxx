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
    double r = bg::get<0>(pt);
    return r * r * std::exp(-r*r);
  };

  auto mat = [&](cartGP pt) -> RealMatrix {
    RealMatrix A(2,2);
    A(0,0) = gaussian(pt);
    A(1,1) = sphGaussian(pt);
    return A;
  };

  EulerMac G(10000);
  double f = G.integrate<double>(gaussian);
  double g = G.integrate<double>(sphGaussian);
  auto X = G.integrate<RealMatrix>(mat);


  cout.precision(10);
  std::cout <<  f << endl;
  std::cout <<  g << endl;
  cout << endl;
  std::cout <<  X << endl;
  return 0;
};

