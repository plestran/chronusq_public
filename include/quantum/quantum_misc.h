template <typename T>
void Quantum<T>::rotateDensities(const std::array<double,3> &n, double theta){
  if(this->nTCS_ != 2) return;

  double sint = std::sin(theta);
  double cost = std::cos(theta);
  double omcost = 1 - cost;
  for(auto i = 0; i < this->onePDMScalar_->rows(); i++)
  for(auto j = 0; j < this->onePDMScalar_->rows(); j++) {
    T x = (*this->onePDMMx_)(i,j);
    T y = (*this->onePDMMy_)(i,j);
    T z = (*this->onePDMMz_)(i,j);

    (*this->onePDMMx_)(i,j)  = (cost + n[0]*n[0] * omcost) * x;
    (*this->onePDMMx_)(i,j) += (n[0]*n[1] * omcost - n[2] * sint) * y;
    (*this->onePDMMx_)(i,j) += (n[0]*n[2] * omcost + n[1] * sint) * z;
    
    (*this->onePDMMy_)(i,j)  = (n[0]*n[1] * omcost + n[2] * sint) * x;
    (*this->onePDMMy_)(i,j) += (cost + n[1]*n[1] * omcost) * y;
    (*this->onePDMMy_)(i,j) += (n[1]*n[2] * omcost - n[0] * sint) * z;

    (*this->onePDMMz_)(i,j)  = (n[0]*n[2] * omcost - n[1] * sint) * x;
    (*this->onePDMMz_)(i,j) += (n[1]*n[2] * omcost + n[0] * sint) * y;
    (*this->onePDMMz_)(i,j) += (cost + n[2]*n[2] * omcost) * z;
  }

};
