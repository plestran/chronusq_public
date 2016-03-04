template<typename T>
void NumericalDifferentiation<T>::cartesianDiff(){
  FileIO fileioTmp("NumDiffScr.inp","NumDiffScr.out","NumDiffScr.rst");
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
  
  

  std::ofstream outf("geom_gdv.grad");
//std::ofstream outf2("geom_gdv.energy");
//outf2 << this->singleSlater_undisplaced_->totalEnergy << endl;
  for(auto iAtm = 0, IX = 0; iAtm < nAtoms; iAtm++)
  for(auto iXYZ = 0; iXYZ < 3     ; iXYZ++, IX++) {
    
    (*mol_p1.cart()) = (*this->molecule_undisplaced_->cart());
    (*mol_m1.cart()) = (*this->molecule_undisplaced_->cart());

    (*mol_p1.cart())(iXYZ,iAtm) += step;
    (*mol_m1.cart())(iXYZ,iAtm) -= step;

    mol_p1.computeI();
    mol_p1.computeRij();
    mol_p1.computeNucRep();
    mol_m1.computeI();
    mol_m1.computeRij();
    mol_m1.computeNucRep();

    BasisSet basis_p1; BasisSet basis_m1;
    AOIntegrals aoints_p1; AOIntegrals aoints_m1;
    SingleSlater<double> ss_p1; SingleSlater<double> ss_m1;

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

    ss_p1.setRef(SingleSlater<double>::RHF);
    ss_m1.setRef(SingleSlater<double>::RHF);
 
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
    MOIntegrals<double> moints_p;
    MOIntegrals<double> moints_m;
 
    Response<double> resp_p;
    Response<double> resp_m;
*/
    
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
    ss_p1.formGuess();
    ss_p1.formFock();
    ss_p1.computeEnergy();
    ss_p1.SCF();
    ss_p1.computeProperties();
    ss_p1.printProperties();

/*
    ss_m1.moA()->setZero();
    ss_m1.densityA()->setZero();
*/
    ss_m1.formGuess();
    ss_m1.formFock();
    ss_m1.computeEnergy();
    ss_m1.SCF();
    ss_m1.computeProperties();
    ss_m1.printProperties();

    double scf_p1 = ss_p1.totalEnergy;
    double scf_m1 = ss_m1.totalEnergy;
    double gsdx = (scf_p1 - scf_m1) / (2*step);
    outf << std::setprecision(10) << gsdx << endl;
    
  }

};
