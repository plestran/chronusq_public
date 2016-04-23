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

    this->isScattered_ = 
      const_cast<Quantum<dcomplex>&>(other).isScattered();

    if(!this->isScattered_) {
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
    } else {
      this->onePDMScalar_ = std::unique_ptr<RealMatrix>(
          new RealMatrix(
            const_cast<Quantum<dcomplex>&>(other).onePDMScalar()->rows(),
            const_cast<Quantum<dcomplex>&>(other).onePDMScalar()->cols()
          )
        );
      (*this->onePDMScalar_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMScalar()->real();

      if(!this->isClosedShell || this->nTCS_ == 2) {
        this->onePDMMz_ = std::unique_ptr<RealMatrix>(
            new RealMatrix(
              const_cast<Quantum<dcomplex>&>(other).onePDMMz()->rows(),
              const_cast<Quantum<dcomplex>&>(other).onePDMMz()->cols()
            )
          );
        (*this->onePDMMz_) = 
          const_cast<Quantum<dcomplex>&>(other).onePDMMz()->real();
      }

      if(this->nTCS_ == 2) {
        this->onePDMMx_ = std::unique_ptr<RealMatrix>(
            new RealMatrix(
              const_cast<Quantum<dcomplex>&>(other).onePDMMx()->rows(),
              const_cast<Quantum<dcomplex>&>(other).onePDMMx()->cols()
            )
          );

        this->onePDMMy_ = std::unique_ptr<RealMatrix>(
            new RealMatrix(
              const_cast<Quantum<dcomplex>&>(other).onePDMMy()->rows(),
              const_cast<Quantum<dcomplex>&>(other).onePDMMy()->cols()
            )
          );


        (*this->onePDMMx_) = 
          const_cast<Quantum<dcomplex>&>(other).onePDMMx()->real();
        (*this->onePDMMy_) = 
          const_cast<Quantum<dcomplex>&>(other).onePDMMy()->real();
      };
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

    this->isScattered_ = 
      const_cast<Quantum<double>&>(other).isScattered();

    if(!this->isScattered_) {
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
    } else {
      this->onePDMScalar_ = std::unique_ptr<ComplexMatrix>(
          new ComplexMatrix(
            const_cast<Quantum<double>&>(other).onePDMScalar()->rows(),
            const_cast<Quantum<double>&>(other).onePDMScalar()->cols()
          )
        );
      this->onePDMScalar_->real() = 
        (*const_cast<Quantum<double>&>(other).onePDMScalar());

      if(!this->isClosedShell || this->nTCS_ == 2) {
        this->onePDMMz_ = std::unique_ptr<ComplexMatrix>(
            new ComplexMatrix(
              const_cast<Quantum<double>&>(other).onePDMMz()->rows(),
              const_cast<Quantum<double>&>(other).onePDMMz()->cols()
            )
          );
        this->onePDMMz_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMMz());
      }

      if(this->nTCS_ == 2) {
        this->onePDMMx_ = std::unique_ptr<ComplexMatrix>(
            new ComplexMatrix(
              const_cast<Quantum<double>&>(other).onePDMMx()->rows(),
              const_cast<Quantum<double>&>(other).onePDMMx()->cols()
            )
          );

        this->onePDMMy_ = std::unique_ptr<ComplexMatrix>(
            new ComplexMatrix(
              const_cast<Quantum<double>&>(other).onePDMMy()->rows(),
              const_cast<Quantum<double>&>(other).onePDMMy()->cols()
            )
          );


        this->onePDMMx_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMMx());
        this->onePDMMy_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMMy());
      };
    }
  };
};
