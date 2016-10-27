#include <wavefunction.h>
namespace ChronusQ {
/*
template<>
template<>
WaveFunction<dcomplex>::WaveFunction(const WaveFunction<double> &other) :
  Quantum<dcomplex>::Quantum<dcomplex>(
    dynamic_cast<const Quantum<double>&>(other)),
  basisset_     (other.basisset()   ),               
  molecule_     (other.molecule()   ),               
  fileio_       (other.fileio()     ),                 
  aointegrals_  (other.aointegrals()),            
  nBasis_       (other.nBasis()     ),
  nShell_       (other.nShell()     ),
  nTT_          (other.nTT()        ),
  nO_           (other.nO()         ),
  nV_           (other.nV()         ),
  nOA_          (other.nOA()        ),
  nOB_          (other.nOB()        ),
  nVA_          (other.nVA()        ),
  nVB_          (other.nVB()        ),
  multip_       (other.multip()     ),
  energyNuclei_ (other.energyNuclei()),
  totalEnergy_  (other.totalEnergy()){

  this->alloc();
  this->moA_->real() = *other.moA();
  if(this->nTCS_ == 2 and !this->isClosedShell)
    this->moB_->real() = *other.moB();
}
*/

}
