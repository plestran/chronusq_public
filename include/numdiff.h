#include <global.h>
#include <molecule.h>
#include <basisset.h>
#include <aointegrals.h>
#include <singleslater.h>
#include <response.h>
#include <tools.h>
#include <classtools.h>


namespace ChronusQ {

  enum DiffType {
    TwoPointSymmetric
  };

  struct Derivatives {
    double          GS_ENERGY;
    double          GS_GRAD;
    Eigen::VectorXd ES_ENERGY;
    Eigen::VectorXd ES_GRAD;
    Eigen::VectorXd ES_GS_NACME;
    RealMatrix      ES_ES_NACME;
  };

  template<typename T>
  class NumericalDifferentiation {
    Molecule        * molecule_undisplaced_;
    SingleSlater<T> * singleSlater_undisplaced_;
    Response<T>     * response_undisplaced_;

    DiffType diffType_;
    std::vector<Derivatives> dervData_;
  public:
    #include <numdiff_constructors.h>

    bool computeGSGradient;
    bool computeESGradient;
    bool computeES2GSNACME;
    bool computeES2ESNACME;
    bool doAllCartesianDOF;

    double step;
    
    // Setters
    inline void setDiffType(DiffType diffType){this->diffType_ = diffType;}; 
    inline void setSingleSlater(SingleSlater<T> &ss){
      this->setMolecule(*ss.molecule());
      this->singleSlater_undisplaced_ = &ss;
      this->computeGSGradient         = true;
    };
    inline void setMolecule(Molecule &mol){
      this->molecule_undisplaced_ = &mol;
    };





    // Procedural functions
    inline void differentiate(){
      if(this->doAllCartesianDOF) this->cartesianDiff();
    }
    
    void cartesianDiff();
    void computeGS(BasisSet&,BasisSet&,SingleSlater<T>&,SingleSlater<T>&);
    void computeES();
    void GSGradient();
    void ESGradient();
    void ES2GSNACME();
    void ES2ESNACME();


  }; // class NumericalDifferentiation
  #include <numdiff_procedural.h>
}; // namespace ChronusQ
