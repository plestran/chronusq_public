/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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

template<typename Dervd>
inline double diffNorm( const Dervd &A, const Dervd &B ){
  return (A-B).norm();
}

template<typename Dervd>
inline double diffNormI( const Dervd &A ){
  return (A - Dervd::Identity(A.rows(),A.cols())).norm();
}

template<typename Dervd>
inline double selfInner( const Dervd &A){
  return A.dot(A);
}
template<typename T>
void NumericalDifferentiation<T>::cartesianDiff(){
  FileIO fileioTmp("NumDiffScr.inp",
//   this->singleSlater_undisplaced_->fileio()->fileNameOut(),
     "NumDiffScr.out","NumDiffScr.rst");

  fileioTmp.iniH5Files();
  fileioTmp.iniStdGroups();

  int nAtoms = this->molecule_undisplaced_->nAtoms(); 
  int multip = this->molecule_undisplaced_->multip(); 
  int charge = this->molecule_undisplaced_->charge();
  int nElec  = this->molecule_undisplaced_->nTotalE();
  std::string basisName = 
    this->singleSlater_undisplaced_->basisset()->basisPath();

  Molecule mol_p1, mol_m1;

  
  mol_p1.setNAtoms(nAtoms);
  mol_m1.setNAtoms(nAtoms);
  mol_p1.setMultip(multip);
  mol_m1.setMultip(multip);
  mol_p1.setCharge(charge);
  mol_m1.setCharge(charge);
  mol_p1.setNTotalE(nElec);
  mol_m1.setNTotalE(nElec);

  mol_p1.alloc(fileioTmp.out);
  mol_m1.alloc(fileioTmp.out);


  for(auto iAtm = 0; iAtm < nAtoms; iAtm++){
    mol_p1.index(iAtm) = this->molecule_undisplaced_->index(iAtm);
    mol_m1.index(iAtm) = this->molecule_undisplaced_->index(iAtm);
  }


  if(this->generateESObjs_){
    // Allocate Response object locally
    this->generated_response_ = 
      std::unique_ptr<Response<T>>(new Response<T>());

    this->response_undisplaced_ = this->generated_response_.get();

    // Temporary MOints for the undisplaced geometry
    MOIntegrals<T> moints;
    moints.communicate(*this->molecule_undisplaced_,
      *this->singleSlater_undisplaced_->basisset(),fileioTmp,
      *this->singleSlater_undisplaced_->aointegrals()->controls(),
      *this->singleSlater_undisplaced_->aointegrals(),
      *this->singleSlater_undisplaced_);

    this->response_undisplaced_->communicate(*this->singleSlater_undisplaced_,
      moints,fileioTmp); 


    moints.initMeta();
    this->singleSlater_undisplaced_->fileio()->out 
      << "Performing Response Calculation at Undisplaced Geometry" << endl;

    this->computeES(*this->response_undisplaced_);
  }


  // Copy T
  TMatrix T_0;
  if(this->computeESGradient){
    if(this->respType_ == RESPONSE_TYPE::CIS) {
      T_0 = this->response_undisplaced_->
        template transDen<SINGLETS>().block(0,0,
        this->response_undisplaced_->
          template nMatDim<SINGLETS>(),this->responseNRoots_);
    } else if(this->respType_ == RESPONSE_TYPE::PPTDA) {
      T_0 = this->response_undisplaced_->
        template transDen<A_PPTDA_SINGLETS>().block(0,0,
        this->response_undisplaced_->template nMatDim<A_PPTDA_SINGLETS>(),
        this->responseNRoots_);
    }
  }


  this->singleSlater_undisplaced_->fileio()->out
    << "  Summary of Results for Undisplaced Geometry:" << endl; 
  this->singleSlater_undisplaced_->fileio()->out
    << "    SCF Energy:" << endl;

  this->singleSlater_undisplaced_->fileio()->out
    << "    E(0) = "
    << std::setprecision(10) << this->singleSlater_undisplaced_->totalEnergy 
    << endl;

  if(this->computeESGradient){
    this->singleSlater_undisplaced_->fileio()->out 
      << "    ES Energies:" << endl;
    for(auto iRt = 0; iRt < this->responseNRoots_; iRt++){
      this->singleSlater_undisplaced_->fileio()->out 
        << "      W(0," << iRt << ") = "
        << std::setprecision(10); 

      if(this->respType_ == RESPONSE_TYPE::CIS){
        this->singleSlater_undisplaced_->fileio()->out 
          << this->response_undisplaced_->template frequencies<SINGLETS>()(iRt) 
          << endl;
      } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
        this->singleSlater_undisplaced_->fileio()->out 
          << this->response_undisplaced_->
            template frequencies<A_PPTDA_SINGLETS>()(iRt) 
          << endl;
      }
    }
  }

  this->singleSlater_undisplaced_->fileio()->out << endl;

/*
  BasisSet basis_p1, basis_m1;
  basis_p1.communicate(fileioTmp);
  basis_m1.communicate(fileioTmp);
  basis_p1.findBasisFile(this->singleSlater_undisplaced_->basisset()->basisPath());
  basis_m1.findBasisFile(this->singleSlater_undisplaced_->basisset()->basisPath());
  basis_p1.parseGlobal();
  basis_m1.parseGlobal();
  // so AOIntegrals / SS know about nBasis
  basis_p1.constructLocal(this->molecule_undisplaced_);
  basis_m1.constructLocal(this->molecule_undisplaced_);

  AOIntegrals aoints_p1, aoints_m1;
  aoints_p1.communicate(mol_p1,basis_p1,fileioTmp,
    *this->singleSlater_undisplaced_->aointegrals()->controls());
  aoints_m1.communicate(mol_m1,basis_m1,fileioTmp,
    *this->singleSlater_undisplaced_->aointegrals()->controls());

  aoints_p1.integralAlgorithm = 
    this->singleSlater_undisplaced_->aointegrals()->integralAlgorithm;
  aoints_m1.integralAlgorithm = 
    this->singleSlater_undisplaced_->aointegrals()->integralAlgorithm;

  aoints_p1.initMeta();
  aoints_m1.initMeta();
  aoints_p1.alloc();
  aoints_m1.alloc();


  SingleSlater<T> ss_p1, ss_m1;
  ss_p1.communicate(mol_p1,basis_p1,aoints_p1,fileioTmp,
    *this->singleSlater_undisplaced_->aointegrals()->controls());
  ss_m1.communicate(mol_m1,basis_m1,aoints_m1,fileioTmp,
    *this->singleSlater_undisplaced_->aointegrals()->controls());

  ss_p1.setRef(this->singleSlater_undisplaced_->Ref());
  ss_m1.setRef(this->singleSlater_undisplaced_->Ref());

  ss_p1.isClosedShell = this->singleSlater_undisplaced_->isClosedShell;
  ss_m1.isClosedShell = this->singleSlater_undisplaced_->isClosedShell;

  ss_p1.genMethString();
  ss_m1.genMethString();
  ss_p1.initMeta();
  ss_m1.initMeta();
  ss_p1.alloc();
  ss_m1.alloc();
*/
  
  

  std::ofstream outGSGrad("geom_gdv.gsgrad");
  std::ofstream outESGrad("geom_gdv.esgrad");
  for(auto iAtm = 0, IX = 0; iAtm < nAtoms; iAtm++)
  for(auto iXYZ = 0; iXYZ < 3     ; iXYZ++, IX++) {
    this->singleSlater_undisplaced_->fileio()->out <<
      "Numerically Evaluating First Derivative IX = " << IX + 1 <<
      " (IATM = " << iAtm + 1 << " , IXYZ = " << iXYZ + 1 <<
      ")" << endl;

    this->generateDispGeom(mol_p1,mol_m1,iXYZ,iAtm);

    mol_p1.printInfo(this->singleSlater_undisplaced_->fileio()->out);
    mol_m1.printInfo(this->singleSlater_undisplaced_->fileio()->out);

    BasisSet        basis_p1;  BasisSet        basis_m1;
    AOIntegrals     aoints_p1; AOIntegrals     aoints_m1;
    SingleSlater<T> ss_p1;     SingleSlater<T> ss_m1;
    Response<T>     resp_p1;   Response<T>     resp_m1;
    MOIntegrals<T>  moints_p1; MOIntegrals<T>  moints_m1;

    basis_p1.communicate(fileioTmp);
    basis_m1.communicate(fileioTmp);

    basis_p1.findBasisFile(basisName);
    basis_m1.findBasisFile(basisName);
 
 
    basis_p1.parseGlobal();
    basis_m1.parseGlobal();
 
    basis_p1.constructLocal(&mol_p1);
    basis_m1.constructLocal(&mol_m1);
    
    basis_p1.makeMaps(1,&mol_p1);
    basis_m1.makeMaps(1,&mol_m1);
    basis_p1.renormShells();
    basis_m1.renormShells();


    aoints_p1.communicate(mol_p1,basis_p1,fileioTmp,
      *this->singleSlater_undisplaced_->aointegrals()->controls());
    aoints_m1.communicate(mol_m1,basis_m1,fileioTmp,
      *this->singleSlater_undisplaced_->aointegrals()->controls());

    ss_p1.communicate(mol_p1,basis_p1,aoints_p1,fileioTmp,
      *this->singleSlater_undisplaced_->aointegrals()->controls());
    ss_m1.communicate(mol_m1,basis_m1,aoints_m1,fileioTmp,
      *this->singleSlater_undisplaced_->aointegrals()->controls());

    ss_p1.setRef(SingleSlater<T>::RHF);
    ss_m1.setRef(SingleSlater<T>::RHF);
 
    // For now, Response requires In-Core integrals ... FIXME
    if(this->computeESGradient){
      aoints_p1.integralAlgorithm = AOIntegrals::INTEGRAL_ALGORITHM::INCORE;
      aoints_m1.integralAlgorithm = AOIntegrals::INTEGRAL_ALGORITHM::INCORE;
    }

    aoints_p1.initMeta();
    aoints_m1.initMeta();

    ss_p1.isClosedShell = true;
    ss_m1.isClosedShell = true;
 
    ss_p1.initMeta();
    ss_m1.initMeta();
 
    ss_p1.genMethString();
    ss_m1.genMethString();
 
    aoints_p1.alloc();
    aoints_m1.alloc();
 
    ss_p1.alloc();
    ss_m1.alloc();
 
/*
    basis_p1.resetAll(); basis_m1.resetAll();

    basis_p1.constructLocal(&mol_p1);
    basis_m1.constructLocal(&mol_m1);
    basis_p1.makeMaps(1,&mol_p1);
    basis_m1.makeMaps(1,&mol_m1);

    ss_p1.initMeta();
    ss_m1.initMeta();
*/

/*
    ss_p1.moA()->setZero();
    ss_p1.densityA()->setZero();
    ss_p1.haveMO = false;
    ss_p1.haveDensity = false;
    ss_p1.havePT = false;
    aoints_p1.haveAOOneE = false;
    aoints_p1.haveAOTwoE = false;
*/
/*
    ss_m1.moA()->setZero();
    ss_m1.densityA()->setZero();
*/
    this->singleSlater_undisplaced_->fileio()->out 
      << "  Performing GS SCF Calculation at + Displaced Geometry" << endl;
    this->computeGS(ss_p1);
    this->singleSlater_undisplaced_->fileio()->out 
      << "  Performing GS SCF Calculation at - Displaced Geometry" << endl;
    this->computeGS(ss_m1);


    // If we're doing NACME, phase check the MO's so that the response
    // is done properly
      
    RealMatrix SAO_0_p1, SAO_0_m1;
    RealMatrix SMO_0_p1, SMO_0_m1;
    if(this->computeES2GSNACME || this->computeES2ESNACME){
      this->singleSlater_undisplaced_->fileio()->out << 
        "  Performing Phase Check on Displaced MOs" << endl;

      SAO_0_p1 = AOIntegrals::genSpx(
        *this->singleSlater_undisplaced_->basisset(),basis_p1
      );
      SAO_0_m1 = AOIntegrals::genSpx(
        *this->singleSlater_undisplaced_->basisset(),basis_m1
      );

      SMO_0_p1 =
        this->singleSlater_undisplaced_->moA()->adjoint() *
        SAO_0_p1 * (*ss_p1.moA());

      SMO_0_m1 =
        this->singleSlater_undisplaced_->moA()->adjoint() *
        SAO_0_m1 * (*ss_m1.moA());

      cout << endl;
      cout << "  Checking | C - C' | Before Phase Check:" << endl;
      
      cout << "  | C(X,Y) - C(X+DX,Y) | = " 
           << diffNorm((*this->singleSlater_undisplaced_->moA()),
               (*ss_p1.moA())) 
           << endl;  
      cout << "  | C(X,Y) - C(X-DX,Y) | = " 
           << diffNorm((*this->singleSlater_undisplaced_->moA()),
               (*ss_m1.moA())) 
           << endl;  


      this->checkPhase((*this->singleSlater_undisplaced_),ss_p1,SMO_0_p1);
      this->checkPhase((*this->singleSlater_undisplaced_),ss_m1,SMO_0_m1);
//    this->checkPhase((*this->singleSlater_undisplaced_->moA()),
//      (*ss_p1.moA()));
//    this->checkPhase((*this->singleSlater_undisplaced_->moA()),
//      (*ss_m1.moA()));

      cout << endl;
      cout << "  Checking | C - C' | After Phase Check:" << endl;
      
      cout << "  | C(X,Y) - C(X+DX,Y) | = " 
           << diffNorm((*this->singleSlater_undisplaced_->moA()),
               (*ss_p1.moA())) 
           << endl;  
      cout << "  | C(X,Y) - C(X-DX,Y) | = " 
           << diffNorm((*this->singleSlater_undisplaced_->moA()),
               (*ss_m1.moA())) 
           << endl;  


      // Recompute the MO Overlaps at phase corrected MOs
      SMO_0_p1 =
        this->singleSlater_undisplaced_->moA()->adjoint() *
        SAO_0_p1 * (*ss_p1.moA());

      SMO_0_m1 =
        this->singleSlater_undisplaced_->moA()->adjoint() *
        SAO_0_m1 * (*ss_m1.moA());
    }


    TMatrix T_p1, T_m1;
    if(this->computeESGradient) {
      moints_p1.communicate(mol_p1,basis_p1,fileioTmp,
        *this->singleSlater_undisplaced_->aointegrals()->controls(),
        aoints_p1,ss_p1);
      moints_m1.communicate(mol_m1,basis_m1,fileioTmp,
        *this->singleSlater_undisplaced_->aointegrals()->controls(),
        aoints_m1,ss_m1);

      resp_p1.communicate(ss_p1,moints_p1,fileioTmp);
      resp_m1.communicate(ss_m1,moints_m1,fileioTmp);

      moints_p1.initMeta();
      moints_m1.initMeta();

      this->singleSlater_undisplaced_->fileio()->out 
        << "  Performing Response Calculation at + Displaced Geometry" << endl;
      this->computeES(resp_p1);

      this->singleSlater_undisplaced_->fileio()->out 
        << "  Performing Response Calculation at - Displaced Geometry" << endl;
      this->computeES(resp_m1);

      // Copy T's
      if(this->respType_ == RESPONSE_TYPE::CIS){
        T_p1 = resp_p1.template transDen<SINGLETS>().block(0,0,
            resp_p1.template nMatDim<SINGLETS>(),this->responseNRoots_);
        T_m1 = resp_m1.template transDen<SINGLETS>().block(0,0,
            resp_m1.template nMatDim<SINGLETS>(),this->responseNRoots_);
      } else if(this->respType_ == RESPONSE_TYPE::PPTDA) {
        T_p1 = resp_p1.template transDen<A_PPTDA_SINGLETS>().block(0,0,
            resp_p1.template nMatDim<A_PPTDA_SINGLETS>(),this->responseNRoots_);
        T_m1 = resp_m1.template transDen<A_PPTDA_SINGLETS>().block(0,0,
            resp_m1.template nMatDim<A_PPTDA_SINGLETS>(),this->responseNRoots_);
      }

      this->singleSlater_undisplaced_->fileio()->out << 
        "  Performing Phase Check on Displaced Transition Vectors" << endl;

//    Phase Check broken as they don't keep phase checked T in this function...
//    this->checkPhase((*this->response_undisplaced_),resp_p1);
//    this->checkPhase((*this->response_undisplaced_),resp_m1);

      cout << endl;
      cout << "  **Checking | T - T' | Before Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " 
           << diffNorm(T_0,T_p1) 
           << endl;  
      cout << "  | T(X,Y) - T(X-DX,Y) | = " 
           << diffNorm(T_0,T_m1) 
           << endl;  
      prettyPrint(cout,T_0.transpose() * T_p1,"+B");
      prettyPrint(cout,T_0.transpose() * T_m1,"-B");

//    this->checkPhase(T_0,T_p1);
//    this->checkPhase(T_0,T_m1);
      TMatrix Inner_0_p1 = T_0.transpose() * T_p1; 
      TMatrix Inner_0_m1 = T_0.transpose() * T_m1; 

      this->checkPhase(T_0,T_p1,Inner_0_p1);
      this->checkPhase(T_0,T_m1,Inner_0_m1);

      cout << "  **Checking | T - T' | After Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " 
           << diffNorm(T_0,T_p1) 
           << endl;  
      cout << "  | T(X,Y) - T(X-DX,Y) | = " 
           << diffNorm(T_0,T_m1) 
           << endl;  
      cout << endl;
 
      prettyPrint(cout,T_0.transpose() * T_p1,"+A");
      prettyPrint(cout,T_0.transpose() * T_m1,"-A");
    }

    Derivatives derv;
    if(this->computeGSGradient) 
      derv.GS_GRAD = this->GSGradient(ss_p1,ss_m1);
    if(this->computeESGradient) 
      derv.ES_GRAD = this->ESGradient(resp_p1,resp_m1);
    
    outGSGrad << std::setprecision(10) << derv.GS_GRAD << endl;
    if(this->computeESGradient){
      outESGrad << std::setprecision(10) 
//              << derv.ES_GRAD[this->responseDiffRoot_] + derv.GS_GRAD 
                << derv.ES_GRAD[1] - derv.ES_GRAD[0]
                << endl;
    }


    this->singleSlater_undisplaced_->fileio()->out
      << "  Summary of Results for IX = " << IX << ":" << endl; 

    this->singleSlater_undisplaced_->fileio()->out
      << "    SCF Energies:" << endl;

    this->singleSlater_undisplaced_->fileio()->out
      << "    E(+) = "
      << std::setprecision(10) <<  ss_p1.totalEnergy << endl
      << "    E(-) = "
      << std::setprecision(10) <<  ss_m1.totalEnergy << endl;
    if(this->computeESGradient){
      this->singleSlater_undisplaced_->fileio()->out 
        << "    ES Energies:" << endl;
      for(auto iRt = 0; iRt < this->responseNRoots_; iRt++){
        if(this->respType_ == RESPONSE_TYPE::CIS){
          this->singleSlater_undisplaced_->fileio()->out 
            << "      W(+," << iRt << ") = "
            << std::setprecision(10) << resp_p1.
              template frequencies<SINGLETS>()(iRt) 
            << endl << "      W(-," << iRt << ") = "
            << std::setprecision(10) << resp_m1.
              template frequencies<SINGLETS>()(iRt) 
            << endl;
        } else if(this->respType_ == RESPONSE_TYPE::PPTDA) {
          this->singleSlater_undisplaced_->fileio()->out 
            << "      W(+," << iRt << ") = "
            << std::setprecision(10) 
            << resp_p1.template frequencies<A_PPTDA_SINGLETS>()(iRt) << endl
            << "      W(-," << iRt << ") = "
            << std::setprecision(10) 
            << resp_m1.template frequencies<A_PPTDA_SINGLETS>()(iRt) << endl;
        }
      }
    }
    if(this->computeES2GSNACME){
//    this->ES2GSNACME(ss_p1,ss_m1,resp_p1,resp_m1,SAO_0_p1,SAO_0_m1,
//      SMO_0_p1,SMO_0_m1);
      derv.ES_GS_NACME = this->ES2GSNACME(ss_p1,ss_m1,T_0,T_p1,T_m1,SAO_0_p1,
	SAO_0_m1,SMO_0_p1,SMO_0_m1);
    }
    if(this->computeES2ESNACME){
      // COPY T's
//    this->ES2ESNACME(ss_p1,ss_m1,resp_p1,resp_m1,SAO_0_p1,SAO_0_m1,
//      SMO_0_p1,SMO_0_m1);
      derv.ES_ES_NACME = this->ES2ESNACME(ss_p1,ss_m1,T_0,T_p1,T_m1,SAO_0_p1,
	SAO_0_m1,SMO_0_p1,SMO_0_m1);
    }
    

    if(this->computeGSGradient)
      this->singleSlater_undisplaced_->fileio()->out 
        << "    GS Gradient = "
        << std::setprecision(10) << derv.GS_GRAD << endl;
    if(this->computeESGradient){
      this->singleSlater_undisplaced_->fileio()->out 
        << "    ES Gradients:" << endl;
      for(auto iRt = 0; iRt < this->responseNRoots_; iRt++)
        this->singleSlater_undisplaced_->fileio()->out 
          << "      W'(" << iRt << ") = "
          << std::setprecision(10) << derv.ES_GRAD(iRt) << endl;
    }

    this->singleSlater_undisplaced_->fileio()->out << endl;
    
   

    this->dervData_.push_back(derv);
    
  }

};

template<typename T>
void NumericalDifferentiation<T>::generateDispGeom(
  Molecule &mol_p1, Molecule &mol_m1, int iXYZ, int iAtm){

  (*mol_p1.cart()) = (*this->molecule_undisplaced_->cart());
  (*mol_m1.cart()) = (*this->molecule_undisplaced_->cart());

  (*mol_p1.cart())(iXYZ,iAtm) += this->step;
  (*mol_m1.cart())(iXYZ,iAtm) -= this->step;

  mol_p1.computeI();
  mol_p1.computeRij();
  mol_p1.computeNucRep();

  mol_m1.computeI();
  mol_m1.computeRij();
  mol_m1.computeNucRep();
};

template<typename T>
void NumericalDifferentiation<T>::computeGS(SingleSlater<T> &ss){

  ss.formGuess();
  ss.formFock();
  ss.computeEnergy();
  ss.SCF();
  ss.computeProperties();
  ss.printProperties();

};

template<typename T>
double NumericalDifferentiation<T>::GSGradient(
  SingleSlater<T> &ss_p1, SingleSlater<T> &ss_m1){
    double scf_p1 = ss_p1.totalEnergy;
    double scf_m1 = ss_m1.totalEnergy;
    double gsdx = (scf_p1 - scf_m1) / (2*this->step);
    return gsdx;
};

template<typename T>
void NumericalDifferentiation<T>::computeES(Response<T> &resp){
  if(this->responseNRoots_ == -1) 
    CErr(
      "Must Set NRoots in NumericalDifferentiation if generating Objects",
      this->singleSlater_undisplaced_->fileio()->out
    );

  if(this->respType_ == RESPONSE_TYPE::NOMETHOD)
    CErr(
      "Must Set RespType in NumericalDifferentiation if generating Objects",
      this->singleSlater_undisplaced_->fileio()->out
    );

  resp.setMeth(this->respType_);
  if(this->singleSlater_undisplaced_->isClosedShell) resp.doSA();
  resp.setNSek(this->responseNRoots_);
  resp.doFull();
  resp.doResponse();

};

template <typename T>
Eigen::VectorXd NumericalDifferentiation<T>::ESGradient(
  Response<T> &resp_p1, Response<T> &resp_m1){
  // This assumes strictly Singlets FIXME

  VectorXd freq_p1, freq_m1;
  if(this->respType_ == RESPONSE_TYPE::CIS){
    freq_p1 = resp_p1.
      template frequencies<SINGLETS>().head(this->responseNRoots_);
    freq_m1 = resp_m1.
      template frequencies<SINGLETS>().head(this->responseNRoots_);
  }else if(this->respType_ == RESPONSE_TYPE::PPTDA){
    freq_p1 = resp_p1.template frequencies<A_PPTDA_SINGLETS>().head(
        this->responseNRoots_);
    freq_m1 = resp_m1.template frequencies<A_PPTDA_SINGLETS>().head(
        this->responseNRoots_);
  }

  Eigen::VectorXd freqDX = (freq_p1 - freq_m1)/(2*this->step);
  return freqDX;
};

template <typename T>
void NumericalDifferentiation<T>::checkDegeneracies(SingleSlater<T> &ss) {

  prettyPrint(cout,*ss.epsA(),"EpsA");
  for(auto iMO = 0; iMO < ss.epsA()->rows() - 1 ; iMO++){
    if( std::abs((*ss.epsA())(iMO) - (*ss.epsA())(iMO+1)) < 1e-8 )
      cout << "WARNING: DEGENERACY IN MOs: " << iMO << " and " << iMO + 1 <<
           endl;
  }
};


#include <numdiff_nacme.h>
