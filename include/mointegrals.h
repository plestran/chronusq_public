#ifndef  INCLUDED_MOINTEGRAL
#define  INCLUDED_MOINTEGRAL
//#include <gsl/gsl_sf_erf.h>
#include "global.h"
#include "basisset.h"
#include "matrix.h"
#include "molecule.h"
#include "fileio.h"
#include "controls.h"
#include "tools.h"
#include "aointegrals.h"
#include "singleslater.h"

/****************************/
/* Error Messages 8000-8999 */
/****************************/
class MOIntegrals{
  int       **iaIndex_;
  int       **ijIndex_;
  int       **abIndex_;

  ChronusQ::BasisSet     	*basisSet_;
  ChronusQ::Molecule    	*molecule_;
  ChronusQ::FileIO       	*fileio_;
  ChronusQ::Controls     	*controls_;
  ChronusQ::AOIntegrals   *aointegrals_;
  SingleSlater  *singleSlater_;

public:
  // these should be protected
  ChronusQ::Matrix<double>    *iajb_;
  ChronusQ::Matrix<double>    *ijab_;
  ChronusQ::Matrix<double>    *ijka_;
  ChronusQ::Matrix<double>    *ijkl_;
  ChronusQ::Matrix<double>    *iabc_;
  ChronusQ::Matrix<double>    *abcd_;

  bool      haveMOiajb;
  bool      haveMOijab;
  bool      haveMOijka;
  bool      haveMOijkl;
  bool      haveMOiabc;
  bool      haveMOabcd;
 
  MOIntegrals(){;};
  ~MOIntegrals(){
    if(      iajb_!=NULL) delete iajb_;
    if(      ijab_!=NULL) delete ijab_;
    if(      ijka_!=NULL) delete ijka_;
    if(      ijkl_!=NULL) delete ijkl_;
    if(      iabc_!=NULL) delete iabc_;
    if(      abcd_!=NULL) delete abcd_;
  };
  
  // initialization function
  void iniMOIntegrals(ChronusQ::Molecule*,ChronusQ::BasisSet*,ChronusQ::FileIO*,ChronusQ::Controls*,ChronusQ::AOIntegrals*,SingleSlater*);

  inline double &iajb(int i, int a, int j, int b){
    return (*iajb_)(this->iaIndex_[i][a],this->iaIndex_[j][b]);
  };
  inline double &iajb(int ia, int jb){
    return (*iajb_)(ia,jb);
  };

  void formiajb();
  void formijab();
  void formijka();
  void formijkl();
  void formiabc();
  void formabcd();
};

#endif
