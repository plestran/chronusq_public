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
  
  

  cout << "HERE 1" << endl;
  std::ofstream outGSGrad("geom_gdv.gsgrad");
  std::ofstream outESGrad("geom_gdv.esgrad");
  for(auto iAtm = 0, IX = 0; iAtm < nAtoms; iAtm++)
  for(auto iXYZ = 0; iXYZ < 3     ; iXYZ++, IX++) {
    this->singleSlater_undisplaced_->fileio()->out <<
      "Numerically Evaluating First Derivative IX = " << IX + 1 <<
      " (IATM = " << iAtm + 1 << " , IXYZ = " << iXYZ + 1 <<
      ")" << endl;

  cout << "HERE 2" << endl;
    this->generateDispGeom(mol_p1,mol_m1,iXYZ,iAtm);

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
  cout << "HERE 2" << endl;

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
  cout << "HERE 2" << endl;
    this->computeGS(ss_p1);
    this->computeGS(ss_m1);
  cout << "HERE 2" << endl;

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
        << "Performing Response Calculation at + Displaced Geometry" << endl;
      this->computeES(resp_p1);

      this->singleSlater_undisplaced_->fileio()->out 
        << "Performing Response Calculation at - Displaced Geometry" << endl;
      this->computeES(resp_m1);

    }
  cout << "HERE 2" << endl;

    Derivatives derv;
    if(this->computeGSGradient) derv.GS_GRAD=this->GSGradient(ss_p1,ss_m1);
    if(this->computeESGradient) derv.ES_GRAD=this->ESGradient(resp_p1,resp_m1);
    
    outGSGrad << std::setprecision(10) << derv.GS_GRAD << endl;
    if(this->computeESGradient){
    outESGrad << std::setprecision(10) << derv.ES_GRAD[this->responseDiffRoot_] 
             + derv.GS_GRAD << endl;
    }

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
  Eigen::VectorXd freq_p1=resp_p1.frequencies()[0].head(this->responseNRoots_);
  Eigen::VectorXd freq_m1=resp_m1.frequencies()[0].head(this->responseNRoots_);

  Eigen::VectorXd freqDX = (freq_p1 - freq_m1)/(2*this->step);
  return freqDX;
};
