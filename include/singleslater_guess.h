template<typename T>
void SingleSlater<T>::formGuess(){
  cout << this->guess_ << endl << READ << endl;
  if(this->guess_ == SAD) this->SADGuess();
  else if(this->guess_ == READ) this->READGuess();
  else CErr("Guess NYI",this->fileio_->out);
}
