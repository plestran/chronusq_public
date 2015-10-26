template<typename T>
void SingleSlater<T>::formGuess(){
  if(this->guess_ == SAD) this->SADGuess();
  else CErr("Guess NYI",this->fileio_->out);
}
