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
  GaussChebFst G(10000);
  double f = 0;
  for(auto i = 0; i < 10000; i++){
    IntegrationPoint pt = G[i];
    double r = bg::get<0>(pt.pt);

    std::cout << r << '\t';
    std::cout << pt.weight << std::endl;
    f += pt.weight * std::exp(-r*r);
  }
  cout.precision(10);
  std::cout <<  f << endl;
  return 0;
};

