template<typename T>
void Response<T>::printInfo(){
  this->fileio_->out << bannerTop << endl;
  
  this->fileio_->out << "Response Settings:" << endl;
  this->fileio_->out << endl;

  this->fileio_->out << std::setw(40) << std::left <<"  Response Function:" << this->methMap_[this->iMeth_] << endl;
  this->fileio_->out << std::setw(40) << std::left <<"    Tamm-Dancoff Approximation:";
  if(this->doTDA_) this->fileio_->out << "Yes";
  else             this->fileio_->out << "No";
  this->fileio_->out << endl;
  this->fileio_->out << std::setw(40) << std::left <<"  Number of Unique Response Matricies:" << iMatIter_.size() << endl;

} // printInfo

