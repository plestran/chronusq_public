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
    DiffType          diffType_;
    Molecule        * molecule_undisplaced_;
    SingleSlater<T> * singleSlater_undisplaced_;
    Response<T>     * response_undisplaced_;

    std::vector<     Derivatives > dervData_;
    std::unique_ptr< Response<T> > generated_response_;

    bool generateESObjs_;
    RESPONSE_TYPE respType_;
    int responseDiffRoot_;
    int responseNRoots_;

  public:
    #include <numdiff_constructors.h>

    bool computeGSGradient;
    bool computeESGradient;
    bool computeES2GSNACME;
    bool computeES2ESNACME;
    bool doAllCartesianDOF;

    double step;
    
    // Setters
    inline void generateESObjs(){this->generateESObjs_ = true;};
    inline void setRespNRoots(int n){this->responseNRoots_ = n;};
    inline void setRespType(RESPONSE_TYPE type){this->respType_ = type;};
    inline void setRespRoot(int n){this->responseDiffRoot_ = n;};

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
    void generateDispGeom(Molecule &, Molecule &, int, int);
    void computeGS(SingleSlater<T>&);
    void computeES(Response<T>&);
    double GSGradient(SingleSlater<T>&,SingleSlater<T>&);
    Eigen::VectorXd ESGradient(Response<T>&,Response<T>&);
    void ES2GSNACME();
    void ES2ESNACME();


  }; // class NumericalDifferentiation
  #include <numdiff_procedural.h>
}; // namespace ChronusQ
