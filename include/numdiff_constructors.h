NumericalDifferentiation(){
  this->computeGSGradient = false;
  this->computeESGradient = false;
  this->computeES2GSNACME = false;
  this->computeES2ESNACME = false;
  this->doAllCartesianDOF = true;

  molecule_undisplaced_     = NULL;
  singleSlater_undisplaced_ = NULL;
  response_undisplaced_     = NULL;

  this->step = 0.0001;
  this->diffType_ = DiffType::TwoPointSymmetric;
};

NumericalDifferentiation(Molecule *mol) :
NumericalDifferentiation() {

  this->molecule_undisplaced_ = mol;

};

NumericalDifferentiation(SingleSlater<T> *singleSlater) :
NumericalDifferentiation(singleSlater->molecule()) {

  this->singleSlater_undisplaced_ = singleSlater;
  this->computeGSGradient = true;

};
NumericalDifferentiation(SingleSlater<T> &singleSlater) :
NumericalDifferentiation(&singleSlater){;};

NumericalDifferentiation(Response<T> *resp) :
NumericalDifferentiation() {
  this->singleSlater_undisplaced_ = resp->singleSlater();
  this->molecule_undisplaced_     = resp->singleSlater()->molecule();

  this->computeGSGradient = true;
  this->computeESGradient = true;
};
