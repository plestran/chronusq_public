#include <quantum.h>
namespace ChronusQ {
  template<>
  template<>
  Quantum<double>::Quantum(const Quantum<dcomplex> &other){
      this->elecDipole_            = 
        const_cast<Quantum<dcomplex>&>(other).elecDipole();
      this->elecQuadpole_          = 
        const_cast<Quantum<dcomplex>&>(other).elecQuadpole();
      this->elecTracelessQuadpole_ = 
        const_cast<Quantum<dcomplex>&>(other).elecTracelessQuadpole();
      this->elecOctpole_           = 
        const_cast<Quantum<dcomplex>&>(other).elecOctpole();
      this->nTCS_                  = 
        const_cast<Quantum<dcomplex>&>(other).nTCS();
      this->isClosedShell          = 
        const_cast<Quantum<dcomplex>&>(other).isClosedShell;
      this->maxMultipole_          = 
        const_cast<Quantum<dcomplex>&>(other).maxMultipole(); 

      this->densityA_ = std::unique_ptr<RealMatrix>(
          new RealMatrix(
            const_cast<Quantum<dcomplex>&>(other).densityA()->rows(),
            const_cast<Quantum<dcomplex>&>(other).densityA()->cols()
          )
        );
      (*this->densityA_) = 
        const_cast<Quantum<dcomplex>&>(other).densityA()->real();

      if(!this->isClosedShell && this->nTCS_ != 2){
        this->densityB_ = std::unique_ptr<RealMatrix>(
            new RealMatrix(
              const_cast<Quantum<dcomplex>&>(other).densityB()->rows(),
              const_cast<Quantum<dcomplex>&>(other).densityB()->cols()
            )
          );
        (*this->densityB_) = 
          const_cast<Quantum<dcomplex>&>(other).densityB()->real();
      }
  };
  template<>
  template<>
  Quantum<dcomplex>::Quantum(const Quantum<double> &other){
      this->elecDipole_            = 
        const_cast<Quantum<double>&>(other).elecDipole();
      this->elecQuadpole_          = 
        const_cast<Quantum<double>&>(other).elecQuadpole();
      this->elecTracelessQuadpole_ = 
        const_cast<Quantum<double>&>(other).elecTracelessQuadpole();
      this->elecOctpole_           = 
        const_cast<Quantum<double>&>(other).elecOctpole();
      this->nTCS_                  = 
        const_cast<Quantum<double>&>(other).nTCS();
      this->isClosedShell          = 
        const_cast<Quantum<double>&>(other).isClosedShell;
      this->maxMultipole_          = 
        const_cast<Quantum<double>&>(other).maxMultipole(); 

      this->densityA_ = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(
            const_cast<Quantum<double>&>(other).densityA()->rows(),
            const_cast<Quantum<double>&>(other).densityA()->cols()
          )
        );
      this->densityA_->real() = 
        (*const_cast<Quantum<double>&>(other).densityA());

      if(!this->isClosedShell && this->nTCS_ != 2){
        this->densityB_ = std::unique_ptr<ComplexMatrix>(
            new ComplexMatrix(
              const_cast<Quantum<double>&>(other).densityB()->rows(),
              const_cast<Quantum<double>&>(other).densityB()->cols()
            )
          );
        this->densityB_->real() = (*const_cast<Quantum<double>&>(other).densityB());
      }
  };
};
