template <typename T>
void RealTime<T>::setIntScheme(const std::string &x){
  if(!x.compare("MMUT") or !x.compare("MODIFIED MIDPOINT") or
     !x.compare("MODIFIEDMIDPOINT"))
    iScheme_ = MMUT;
  else if(!x.compare("EM2") or !x.compare("EXPLICIT MAGNUS 2") or
          !x.compare("EXPLICITMAGNUS2"))
    iScheme_ = ExpMagnus2;
  else
    CErr(x + std::string(" is not a recognized integration scheme."),
     fileio_->out);
};

template <typename T>
void RealTime<T>::setMMUTRstScheme(const std::string &x){
  if(!x.compare("FE") or !x.compare("FORWARD EULER") or
     !x.compare("FORWARDEULER"))
    iRstScheme_ = ForwardEuler;
  else if(!x.compare("EM2") or !x.compare("EXPLICIT MAGNUS 2") or
          !x.compare("EXPLICITMAGNUS2"))
    iRstScheme_ = ExplicitMagnus2;
  else
    CErr(x + std::string(" is not a recognized MMUT restart scheme."),
     fileio_->out);
};

template <typename T>
void RealTime<T>::setEnvlp(const std::string &x) {
  if(!x.compare("CONSTANT") or !x.compare("PLANE WAVE"))
    this->iEnvlp_ = Constant;
  else if(!x.compare("LINRAMP"))
    this->iEnvlp_ = LinRamp;
  else if(!x.compare("GAUSSIAN"))
    this->iEnvlp_ = Gaussian;
  else if(!x.compare("STEP"))
    this->iEnvlp_ = Step;
  else if(!x.compare("DELTA")) {
    this->iEnvlp_ = Step;
    this->tOn_ = 0.0; this->tOff_ = 1e-6;
  } 
  else if(!x.compare("SINSQ"))
    this->iEnvlp_ = SinSq;
  else
    CErr(std::string("RT Envelope ") + x + std::string(" not Recognized"),
      this->fileio_->out);

}
