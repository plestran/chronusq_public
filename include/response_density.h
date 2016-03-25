template<typename T>
void Response<T>::formDensity(){
  if(this->iDen_ == TRANSITION)      this->formTransitionDensity();
  else if(this->iDen_ == DIFFERENCE) this->formDifferenceDensity();
};

template<typename T>
void Response<T>::formTransitionDensity(){
  if(this->iClass_ == RESPONSE_CLASS::FOPPA) 
    this->formTransitionDensityFOPPA();
  else if(this->iClass_ == RESPONSE_CLASS::PPRPA)
    this->formTransitionDensityPPRPA();
};

template<typename T>
void Response<T>::formDifferenceDensity(){
  if(this->iClass_ == RESPONSE_CLASS::FOPPA) 
    this->formDifferenceDensityFOPPA();
  else if(this->iClass_ == RESPONSE_CLASS::PPRPA)
    this->formDifferenceDensityPPRPA();
};
