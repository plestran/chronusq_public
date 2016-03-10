NumericalDifferentiation(){
  this->computeGSGradient = false;
  this->computeESGradient = false;
  this->computeES2GSNACME = false;
  this->computeES2ESNACME = false;
  this->doAllCartesianDOF = true;
  this->generateESObjs_   = false;

  molecule_undisplaced_     = NULL;
  singleSlater_undisplaced_ = NULL;
  response_undisplaced_     = NULL;

  this->step = 0.00001;
  this->diffType_ = DiffType::TwoPointSymmetric;
  this->respType_ = RESPONSE_TYPE::NOMETHOD;
  this->responseDiffRoot_ = -1;
  this->responseNRoots_   = -1;
};

NumericalDifferentiation(Molecule *mol) :
NumericalDifferentiation() {
  this->setMolecule(*mol);

};

NumericalDifferentiation(SingleSlater<T> *singleSlater) :
NumericalDifferentiation() {
  this->setSingleSlater(*singleSlater);
};

NumericalDifferentiation(SingleSlater<T> &singleSlater) :
NumericalDifferentiation(&singleSlater){;};

NumericalDifferentiation(Response<T> *resp) :
NumericalDifferentiation() {

  this->setSingleSlater(*resp->singleSlater());

  this->computeESGradient = true;
};

