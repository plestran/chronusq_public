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
