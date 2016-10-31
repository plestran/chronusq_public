#ifndef INCLUDED_WFN
#define INCLUDED_WFN
#include <global.h>
#include <cerr.h>
#include <molecule.h>
#include <aointegrals.h>
#include <quantum.h>

namespace ChronusQ {
template<typename T>
class WaveFunction : public Quantum<T> {
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMatrix;
  typedef Eigen::Map<TMatrix> TMap;
  typedef Eigen::Matrix<T,Dynamic,1,ColMajor> TVec;

protected:

  BasisSet *    basisset_;    ///< Basis Set
  Molecule *    molecule_;    ///< Molecular specificiations
  FileIO *      fileio_;      ///< Access to output file
  AOIntegrals * aointegrals_; ///< Molecular Integrals over GTOs (AO basis)

  // Dimensions to be copied or calculated from BasisSet
  int nBasis_;
  int nShell_;
  int nTT_;

  // Dimensions to be copied or calculated from Molecule
  int nO_;
  int nV_;
  int nOA_;
  int nOB_;
  int nVA_;
  int nVB_;
  int multip_;

  // Energies
  double energyNuclei_; ///< N-N Repulsion Energy
  double totalEnergy_;  ///< Sum of all energetic contributions

  // Matricies
  std::unique_ptr<TMap>  moA_; ///< Alpha or Full MO Coefficients
  std::unique_ptr<TMap>  moB_; ///< Beta MO Coefficient Matrix
  std::unique_ptr<RealMap>  epsA_;       ///< Alpha or Full Eigenenergies
  std::unique_ptr<RealMap>  epsB_;       ///< Beta Fock Eigenenergie

  inline void checkWorkers(){
    if(this->fileio_  == NULL) 
      CErr("Fatal: Must initialize WaveFunction with FileIO Object");
    if(this->basisset_ == NULL) 
      CErr("Fatal: Must initialize WaveFunction with BasisSet Object",
           this->fileio_->out);
    if(this->molecule_ == NULL) 
      CErr("Fatal: Must initialize WaveFunction with Molecule Object",
           this->fileio_->out);
    if(this->aointegrals_== NULL)
      CErr("Fatal: Must initialize WaveFunction with AOIntegrals Object",
           this->fileio_->out);
    if(this->memManager_ == NULL)
      CErr("Fatal: Must initialize WaveFunction with CQMemManager Object",
           this->fileio_->out);
    
  }

  inline void checkMeta(){
    this->checkWorkers();
    if(this->nBasis_ == 0 || this->nShell_ == 0)
      CErr(
        "Fatal: WaveFunction Object Initialized with NBasis = 0 or NShell = 0",
        this->fileio_->out);

    if((this->molecule_->nTotalE() % 2 == 0 && this->multip_ % 2 == 0) ||
       (this->molecule_->nTotalE() % 2 != 0 && this->multip_ % 2 != 0))
      CErr(std::string("Fatal: The specified multiplicity is impossible within")
           +std::string(" the number of electrons given"),this->fileio_->out);

    if(this->multip_ != 1 and this->isClosedShell)
      CErr("Closed Shell Reference and Spin-Multiplicity != 0 Incompatible",
        this->fileio_->out);
  }
public:

  WaveFunction() : Quantum<T>(),
    basisset_     (NULL),               
    molecule_     (NULL),               
    fileio_       (NULL),                 
    aointegrals_  (NULL),            
    nBasis_       (0),
    nShell_       (0),
    nTT_          (0),
    nO_           (0),
    nV_           (0),
    nOA_          (0),
    nOB_          (0),
    nVA_          (0),
    nVB_          (0),
    multip_       (0),
    energyNuclei_ (0.0),
    totalEnergy_  (0.0),
    moA_          (nullptr),        
    moB_          (nullptr) { ; }; 

  WaveFunction(WaveFunction &other) :
    Quantum<T>(dynamic_cast<Quantum<T>&>(other)),
    basisset_     (other.basisset_   ),               
    molecule_     (other.molecule_   ),               
    fileio_       (other.fileio_     ),                 
    aointegrals_  (other.aointegrals_),            
    nBasis_       (other.nBasis_     ),
    nShell_       (other.nShell_     ),
    nTT_          (other.nTT_        ),
    nO_           (other.nO_         ),
    nV_           (other.nV_         ),
    nOA_          (other.nOA_        ),
    nOB_          (other.nOB_        ),
    nVA_          (other.nVA_        ),
    nVB_          (other.nVB_        ),
    multip_       (other.multip_     ),
    energyNuclei_ (other.energyNuclei_),
    totalEnergy_  (other.totalEnergy_){

    this->alloc();
    *this->moA_ = *other.moA_;
    if(this->nTCS_ == 1 and !this->isClosedShell)
      *this->moB_ = *other.moB_;
  };

  template <typename U>
  WaveFunction(const U&);

  // Link up to all of the other worker classes
  inline void communicate(Molecule &mol, BasisSet&basis, AOIntegrals &aoints, 
    FileIO &fileio, CQMemManager &memManager){

    Quantum<T>::communicate(memManager);

    this->molecule_    = &mol;
    this->basisset_    = &basis;
    this->fileio_      = &fileio;
    this->aointegrals_ = &aoints;
  }

  inline void initMeta(){
    this->checkWorkers();

    this->nBasis_       = this->basisset_->nBasis();
    this->nTT_          = this->nBasis_ * (this->nBasis_ + 1) / 2;
    this->multip_       = this->molecule_->multip();
    this->nShell_       = this->basisset_->nShell();
    this->energyNuclei_ = this->molecule_->energyNuclei();

    int nTotalE  = this->molecule_->nTotalE();
    int nSingleE = this->multip_ - 1;
    this->nOB_ = (nTotalE - nSingleE)/2;
    this->nVB_ = this->nBasis_ - this->nOB_;
    this->nOA_ = this->nOB_ + nSingleE;
    this->nVA_ = this->nBasis_ - this->nOA_;

    this->nO_ = this->nOA_ + this->nOB_;
    this->nV_ = this->nVA_ + this->nVB_;
  }

  void alloc();

  // Polymorphism...
  virtual void computeEnergy() = 0;
    
  // Quantum dependencies
  virtual void formDensity() = 0;
  virtual void computeSSq()  = 0;

  // access to private data
  inline int nBasis()  const  { return this->nBasis_; };
  inline int nTT()     const  { return this->nTT_;    };
  inline int nShell()  const  { return this->nShell_; };
  inline int nOA()     const  { return this->nOA_;    };
  inline int nOB()     const  { return this->nOB_;    };     
  inline int nVA()     const  { return this->nVA_;    };
  inline int nVB()     const  { return this->nVB_;    };
  inline int nO()      const  { return this->nO_;     };     
  inline int nV()      const  { return this->nV_;     };
  inline int multip()  const  { return this->multip_; };

  inline const double &energyNuclei() const { return this->energyNuclei_; };
  inline const double &totalEnergy()  const { return this->totalEnergy_;  };

  inline TMap* moA()   const  { return this->moA_.get(); };
  inline TMap* moB()   const  { return this->moB_.get(); };
  inline RealMap* epsA()              { return this->epsA_.get();     };
  inline RealMap* epsB()              { return this->epsB_.get();     };
  inline BasisSet     * basisset()    const   { return this->basisset_;    };
  inline Molecule     * molecule()    const   { return this->molecule_;    };
  inline FileIO       * fileio()      const   { return this->fileio_;      };
  inline AOIntegrals  * aointegrals() const   { return this->aointegrals_; };
};
};
#endif
