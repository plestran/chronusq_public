/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#ifndef INCLUDED_RESP_DECL
#define INCLUDED_RESP_DECL

#include <pscf.h>
#include <qn.h>

namespace ChronusQ {

enum RESPONSE_MATRIX_PARTITION {
  FULL,
  SPIN_SEPARATED,
  SPIN_ADAPTED
};

struct ResponseSettings {
  RESPONSE_MATRIX_PARTITION part;
  bool doSinglets;
  bool doTriplets;
  bool doTDA;
};

inline ResponseSettings DefaultResponseSettings() {
  return ResponseSettings{FULL,true,false,true};
};

template <typename T>
class ResponseMatrix : public QNCallable<T> {
  // Useful Eigen Typedefs
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TMat> TMap;
  typedef Eigen::Map<TVec> TVecMap;

protected:
  CQMemManager *memManager_;
  PostSCF<T>   *pscf_;


/*
  T * fullMatMem_;
  T * solutionVecRMem_;
  T * solutionVecLMem_;
*/

  ResponseSettings sett_;

public:

  ResponseMatrix(PostSCF<T> *pscf, ResponseSettings sett): QNCallable<T>(),
    pscf_(pscf), sett_(sett) { 
    this->memManager_ = pscf_->memManager(); 
    this->checkMeta();
  };

  ResponseMatrix(PostSCF<T> *pscf): 
    ResponseMatrix(pscf,DefaultResponseSettings()){ };

  ResponseMatrix() : ResponseMatrix(NULL,DefaultResponseSettings()){ };


  virtual void formFull()  = 0;
  virtual void initMeta()  = 0;
  virtual void postSolve() = 0;

  inline void alloc() {
    cout << "In ResponseMatrix alloc" << endl;
    this->solutionVecR_ = 
      this->memManager_->template malloc<T>(this->nSek_*this->nSingleDim_);
    this->omega_ = this->memManager_->template malloc<double>(this->nSek_);
    this->diag_  = this->memManager_->template malloc<double>(this->nSingleDim_);
  } 

  inline void checkMeta() {
    if(this->pscf_->nTCS() == 2 and this->sett_.part != FULL){
      cout << "Switching Matrix partitioning to FULL for 2C method" << endl;
      this->sett_.part = FULL;
    } else if(this->pscf_->nTCS() == 1 and 
              !this->pscf_->isClosedShell and 
              this->sett_.part == SPIN_ADAPTED) {
      cout << "Switching Matrix partitioning to SPIN_SEPARATED for ";
      cout << "Unresticted method" << endl;
      this->sett_.part = SPIN_SEPARATED;
    }
  }
};

template <typename T>
class Response : public PostSCF<T> {
  // Useful Eigen Typedefs
  typedef Eigen::Matrix<T,Dynamic,Dynamic,ColMajor> TMat;
  typedef Eigen::Matrix<T,Dynamic,1> TVec;
  typedef Eigen::Map<TMat> TMap;
  typedef Eigen::Map<TVec> TVecMap;

protected:

  /** Job Control for Response calculation **/
  bool useIncoreInts_;   ///< Use Incore Integrals
  bool doFull_;          ///< Do full matrix problem (in-core)
  bool debugIter_;       ///< Diagonalize full incore matrix iteratively
  bool doTDA_;           ///< Invoke TDA

  size_t nSek_;
  size_t nGuess_;

  QNProblemType iJob_;   ///< Response Job Type
  RESPONSE_MATRIX_PARTITION part_;

  std::vector<ResponseMatrix<T>*> iMat_;

public:
  
  Response(QNProblemType typ, RESPONSE_MATRIX_PARTITION part, bool doTDA): 
    PostSCF<T>(),
    iJob_(typ), doTDA_(doTDA), part_(part),
    useIncoreInts_(false),
    doFull_       (false),
    debugIter_    (false),
    nSek_         (3),
    nGuess_       (6){ };

  Response() : Response(DIAGONALIZATION,FULL,false){ };

  Response(Response &other) :
    PostSCF<T>(dynamic_cast<PostSCF<T>&>(other))
    { ; };

  // Quantum compliant
  virtual void formDensity() = 0;
  virtual void computeSSq() = 0;


  inline void communicate(WaveFunction<T> &wfn, CQMemManager &memManager) {
    PostSCF<T>::communicate(wfn,memManager);
    cout << "In Response Communicate" << endl;
  };

  inline void initMeta(){
    PostSCF<T>::initMeta();
    cout << "In Response initMeta" << endl;
  }

  inline void runResponse(){
    std::function<H5::DataSet*(const H5::CompType&,std::string&,
      std::vector<hsize_t>&)> fileFactory = 
        std::bind(&FileIO::createScratchPartition,this->fileio_,
        std::placeholders::_1,std::placeholders::_2,std::placeholders::_3);

    this->formGuess();
    for(auto MAT : this->iMat_) {
      QuasiNewton2<T> qn(MAT,this->memManager_,fileFactory);
      if(this->doFull_) {
        MAT->formFull();
        MAT->setAlgorithm(FULL_SOLVE);
      } else 
        MAT->formFull(); // This is for debugging...?
    //qn.run();
    //MAT->postSolve();
    //TMap VR(MAT->solutionVecR(),MAT->nSingleDim(),this->nSek_);
    //prettyPrintSmart(cout,VR,"Solution R");
    //RealMap W(MAT->omega(),this->nSek_,1);
    //prettyPrintSmart(cout,phys.eVPerHartree*W,"E");
    }
  };

  inline void alloc() {
    for(auto MAT : this->iMat_) {
      MAT->initMeta();
      if(this->doFull_) MAT->setNSek(MAT->nSingleDim());
      else {             
        MAT->setNSek(this->nSek_);
        MAT->setNGuess(this->nGuess_);
      }
      MAT->alloc();
    }
  }


  inline void formGuess() {
    for(auto MAT : this->iMat_) {
      MAT->formDiag();
      MAT->formGuess();
    }
  }
  void doTDA()  { this->doTDA_  = true; };
  void doFull() { this->doFull_ = true; };
  void setNSek(size_t N) { this->nSek_ = N; };
  void setNGuess(size_t N) { this->nGuess_ = N; };


};
};

#endif

