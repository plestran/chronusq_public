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

    auto NBT = const_cast<Quantum<dcomplex>&>(other).onePDMA()->rows(); 
    auto NB  = const_cast<Quantum<dcomplex>&>(other).onePDMScalar()->rows(); 
    auto NBTSq = NBT*NBT;
    auto NBSq = NB*NB; 

    this->onePDMA_ = std::unique_ptr<RealMap>(
        new RealMap( this->memManager_->malloc<double>(NBTSq),NBT,NBT)
      );
    (*this->onePDMA_) = 
      const_cast<Quantum<dcomplex>&>(other).onePDMA()->real();
    
    if(!this->isClosedShell && this->nTCS_ != 2){
      this->onePDMB_ = std::unique_ptr<RealMap>(
          new RealMap( this->memManager_->malloc<double>(NBTSq),NBT,NBT)
        );
      (*this->onePDMB_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMB()->real();
    }

    this->onePDMScalar_ = std::unique_ptr<RealMap>(
        new RealMap( this->memManager_->malloc<double>(NBSq),NB,NB)
      );
    (*this->onePDMScalar_) = 
      const_cast<Quantum<dcomplex>&>(other).onePDMScalar()->real();

    if(!this->isClosedShell || this->nTCS_ == 2) {
      this->onePDMMz_ = std::unique_ptr<RealMap>(
          new RealMap( this->memManager_->malloc<double>(NBSq),NB,NB)
        );
      (*this->onePDMMz_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMMz()->real();
    }

    if(this->nTCS_ == 2) {
      this->onePDMMx_ = std::unique_ptr<RealMap>(
          new RealMap( this->memManager_->malloc<double>(NBSq),NB,NB)
        );

      this->onePDMMy_ = std::unique_ptr<RealMap>(
          new RealMap( this->memManager_->malloc<double>(NBSq),NB,NB)
        );


      (*this->onePDMMx_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMMx()->real();
      (*this->onePDMMy_) = 
        const_cast<Quantum<dcomplex>&>(other).onePDMMy()->real();
    };
  };

  template<>
  template<>
  Quantum<dcomplex>::Quantum(const Quantum<double> &other){
    cout << "HERE 3" << endl;
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
    this->memManager_ =
      const_cast<Quantum<double>&>(other).memManager();

    this->isScattered_ = 
      const_cast<Quantum<double>&>(other).isScattered();

    auto NBT = const_cast<Quantum<double>&>(other).onePDMA()->rows(); 
    auto NBTSq = NBT*NBT;

    this->onePDMA_ = std::unique_ptr<ComplexMap>(
        new ComplexMap( this->memManager_->malloc<dcomplex>(NBTSq),NBT,NBT)
      );
    this->onePDMA_->setZero();

    this->onePDMA_->real() = (*const_cast<Quantum<double>&>(other).onePDMA());
    prettyPrintComplex(cout,*this->onePDMA_,"PA");

    if(!this->isClosedShell || this->nTCS_ == 2){
      auto NB  = const_cast<Quantum<double>&>(other).onePDMScalar()->rows(); 
      auto NBSq = NB*NB; 

      this->onePDMB_ = std::unique_ptr<ComplexMap>(
          new ComplexMap( this->memManager_->malloc<dcomplex>(NBTSq),NBT,NBT)
        );
      this->onePDMScalar_ = std::unique_ptr<ComplexMap>(
          new ComplexMap( this->memManager_->malloc<dcomplex>(NBSq),NB,NB)
        );
      this->onePDMMz_ = std::unique_ptr<ComplexMap>(
          new ComplexMap( this->memManager_->malloc<dcomplex>(NBSq),NB,NB)
        );

      this->onePDMB_->setZero();
      this->onePDMScalar_->setZero();
      this->onePDMMz_->setZero();

      this->onePDMScalar_->real() = 
        (*const_cast<Quantum<double>&>(other).onePDMScalar());
      this->onePDMMz_->real() = 
        (*const_cast<Quantum<double>&>(other).onePDMMz());

      if(this->nTCS_ == 1)
        this->onePDMB_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMB());
      else {
        this->onePDMMx_ = std::unique_ptr<ComplexMap>(
            new ComplexMap( this->memManager_->malloc<dcomplex>(NBSq),NB,NB)
          );
        this->onePDMMy_ = std::unique_ptr<ComplexMap>(
            new ComplexMap( this->memManager_->malloc<dcomplex>(NBSq),NB,NB)
          );

        this->onePDMMx_->setZero();
        this->onePDMMy_->setZero();

        this->onePDMMx_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMMx());
        this->onePDMMy_->real() = 
          (*const_cast<Quantum<double>&>(other).onePDMMy());

      } // NTCS check
    } // if not RHF/KS
  }; // COPY
};
