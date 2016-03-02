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
    void setDiffType(DiffType diffType){this->diffType_ = diffType;}; 




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
