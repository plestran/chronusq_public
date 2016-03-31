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

      this->onePDMA_ = std::unique_ptr<RealMatrix>(
          new RealMatrix(
            const_cast<Quantum<dcomplex>&>(other).onePDMA()->rows(),
            const_cast<Quantum<dcomplex>&>(other).onePDMA()->cols()
          )
        );
      (*this->onePDMA_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMA()->real();

      if(!this->isClosedShell && this->nTCS_ != 2){
        this->onePDMB_ = std::unique_ptr<RealMatrix>(
            new RealMatrix(
              const_cast<Quantum<dcomplex>&>(other).onePDMB()->rows(),
              const_cast<Quantum<dcomplex>&>(other).onePDMB()->cols()
            )
          );
        (*this->onePDMB_) = 
          const_cast<Quantum<dcomplex>&>(other).onePDMB()->real();
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

      this->onePDMA_ = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(
            const_cast<Quantum<double>&>(other).onePDMA()->rows(),
            const_cast<Quantum<double>&>(other).onePDMA()->cols()
          )
        );
      this->onePDMA_->real() = 
        (*const_cast<Quantum<double>&>(other).onePDMA());

      if(!this->isClosedShell && this->nTCS_ != 2){
        this->onePDMB_ = std::unique_ptr<ComplexMatrix>(
            new ComplexMatrix(
              const_cast<Quantum<double>&>(other).onePDMB()->rows(),
              const_cast<Quantum<double>&>(other).onePDMB()->cols()
            )
          );
        this->onePDMB_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMB());
      }
  };
};
