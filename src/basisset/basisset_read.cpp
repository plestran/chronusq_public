#include <basisset.h>
namespace ChronusQ{
void BasisSet::basisSetRead(FileIO * fileio, Molecule * mol){

  std::string readString;
  
  this->fileio_ = fileio;
 
  this->fileio_->in >> readString; // read the name of the basis set file
  this->findBasisFile(readString); // Try to find the basis set file
  this->parseGlobal();
  this->constructLocal(mol);
  this->makeMapSh2Bf();
  this->makeMapSh2Cen(mol);
  this->makeMapCen2Bf(mol);
  this->printInfo();
  this->renormShells();

}; // basisSetRead

void BasisSet::findBasisFile(std::string fName){
  std::string tmpStr;

  tmpStr = "/" + fName;
  tmpStr.insert(0,BASIS_PATH);
  this->setBasisPath(tmpStr);

  this->basisFile_ = std::unique_ptr<ifstream>(new ifstream(this->basisPath_));

  
  if(!this->basisFile_->fail()){ // Check if file is in BASIS_PATH
    this->fileio_->out << "Reading Basis Set from: " << this->basisPath_ << endl;
  } else {
    this->basisFile_.reset();
    this->setBasisPath(fName);
    this->basisFile_ = std::unique_ptr<ifstream>(new ifstream(this->basisPath_));
    if(!this->basisFile_->fail()){ // Check if file is in PWD
      this->fileio_->out << "Reading Basis Set from: ./" << this->basisPath_ << endl;
    } else CErr("Could not find basis set file \"" + fName + "\"");
  }
};

void BasisSet::parseGlobal(){

  std::string readString;
  std::string nameOfAtom;
  std::string shSymb;
  double      scaleFactor;
  double      dummyDouble;
  int         dummyInt;
  int         contDepth;
  int         nShell_tmp = 0;
  int atomicNumber;
  int indx;
  std::vector<libint2::Shell> tmpShell;

  bool readRec = false;
  bool newRec  = false;
  bool firstRec = true;
  int nEmpty = 0;
  int nComm  = 0;
  int nRec   = 0;

  while(!this->basisFile_->eof()){
    std::getline(*this->basisFile_,readString);
    if(readString.size() == 0)    nEmpty++;
    else if(readString[0] == '!') nComm++;
    else if(!readString.compare("****")){
      std::getline(*this->basisFile_,readString);
      if(readString.size() == 0) { nEmpty++; readRec = false; continue;}
      nRec++;
      readRec = true;
      newRec  = true;
    }

    if(readRec){
      std::istringstream iss(readString);
      std::vector<std::string> tokens(std::istream_iterator<std::string>{iss},
        std::istream_iterator<std::string>{});
      if(newRec){
        if(!firstRec) {
//        cout << " " << nShell_tmp << endl;
          this->refShells_.push_back(ReferenceShell{atomicNumber,indx,tmpShell});
        }
//      for(auto it = tokens.begin(); it != tokens.end(); it++) cout << *it << " ";
//      cout << tokens[0] << " ";
        indx = HashAtom(tokens[0],0);
//      cout << " " << elements[indx].atomicNumber;
        atomicNumber = elements[indx].atomicNumber;
        newRec = false;
        firstRec = false;
        nShell_tmp = 0;
        tmpShell.clear();
      } else {
        contDepth = std::stoi(tokens[1]);
        shSymb    = tokens[0];
        std::vector<double> exp;
        std::vector<double> contPrimary;
        std::vector<double> contSecondary;

        for(auto i = 0; i < contDepth; i++) {
          std::getline(*this->basisFile_,readString);
          std::istringstream iss2(readString);
          std::vector<std::string> tokens2(std::istream_iterator<std::string>{iss2},
            std::istream_iterator<std::string>{});

          exp.push_back(std::stod(tokens2[0]));
          contPrimary.push_back(std::stod(tokens2[1]));
          if(!shSymb.compare("SP"))
            contSecondary.push_back(std::stod(tokens2[2]));
        }

        if(!shSymb.compare("SP")) {
          tmpShell.push_back(
            libint2::Shell{ exp, {{0,this->doSph_,contPrimary}}, {{0,0,0}} }
          );
          tmpShell.push_back(
            libint2::Shell{ exp, {{1,this->doSph_,contSecondary}}, {{0,0,0}} }
          );
        } else {
          tmpShell.push_back(
            libint2::Shell{ exp, {{HashL(shSymb),this->doSph_,contPrimary}}, {{0,0,0}} }
          );
        }

        if(!shSymb.compare("SP")) nShell_tmp += 2;
        else                      nShell_tmp++;
      }
    }
  }
  this->refShells_.push_back(ReferenceShell{atomicNumber,indx,tmpShell}); // Append the last Rec
//cout << endl;
//cout << nEmpty << endl;
//cout << nComm << endl;
//cout << nRec << endl;
 
  for(auto it = this->refShells_.begin(); it != this->refShells_.end(); ++it){
    cout << (*it).atomicNumber << '\t' << (*it).index << endl;
    for(auto jt = (*it).shells.begin(); jt != (*it).shells.end(); ++jt)
      cout << *jt << endl;
    cout << endl;
  }

};

void BasisSet::constructLocal(Molecule * mol){
  cout << mol->nAtoms() << endl << this->refShells_.size() << endl;
  for(auto iAtom = 0; iAtom < mol->nAtoms(); iAtom++){
    bool found = false;
    for(auto iRef = this->refShells_.begin(); iRef != this->refShells_.end(); ++iRef){
      if(mol->index(iAtom) == (*iRef).index){
        for(auto iShell = (*iRef).shells.begin(); iShell != (*iRef).shells.end(); ++iShell)
          this->shells_.push_back(
            libint2::Shell{ 
              iShell->alpha, 
              {iShell->contr[0]}, 
              {{ (*mol->cart())(0,iAtom),
                 (*mol->cart())(1,iAtom),
                 (*mol->cart())(2,iAtom)}}
            }
          );
        found = true;
      }
    }
    if(!found)  
      CErr("Atomic Number " + 
             std::to_string(elements[mol->index(iAtom)].atomicNumber) +
             " not found in current Basis Set",
           this->fileio_->out);
  }
  cout << "Local Shells" << endl;
  for(auto iShell = this->shells_.begin(); iShell != this->shells_.end(); ++iShell)
    cout << *iShell << endl;

  this->computeMeta();
}

void BasisSet::computeMeta(){
  this->nShell_     = this->shells_.size();
  this->nShellPair_ = this->nShell_ * (this->nShell_ + 1) / 2;
  
  for(auto iShell = this->shells_.begin(); iShell != this->shells_.end(); ++iShell){
    this->nBasis_ += (*iShell).size();
    cout << this->nBasis_ << endl;
    auto L = (*iShell).contr[0].l;
    auto shPrim = (*iShell).alpha.size();  
    if( L      > this->maxL_   ) this->maxL_    = L     ;
    if( shPrim > this->maxPrim_) this->maxPrim_ = shPrim;
    this->nPrimitive_ += shPrim * (*iShell).size();
  }

  this->nLShell_ = std::vector<int>(this->maxL_+1,0);
  for(auto shell : this->shells_){
    this->nLShell_[shell.contr[0].l]++;
  }

  cout << "nBasis       " <<  this->nBasis_     << endl; 
  cout << "nPrimitive   " <<  this->nPrimitive_ << endl; 
  cout << "maxPrim      " <<  this->maxPrim_    << endl; 
  cout << "maxL         " <<  this->maxL_       << endl; 
  cout << "nShell       " <<  this->nShell_     << endl; 
  cout << "nShellPair   " <<  this->nShellPair_ << endl; 
  for(auto sh : this->nLShell_) cout << sh << endl;
}

void BasisSet::makeMapSh2Bf(){
  auto n = 0;
  for(auto shell : this->shells_){
     this->mapSh2Bf_.push_back(n);
     n += shell.size();
  }
  this->haveMapSh2Bf = true;
}

void BasisSet::makeMapSh2Cen(Molecule *mol){
  for(auto shell : this->shells_){
    for(auto iAtom = 0; iAtom < mol->nAtoms(); iAtom++){
      std::array<double,3> center = {{ (*mol->cart())(0,iAtom),
                                       (*mol->cart())(1,iAtom),
                                       (*mol->cart())(2,iAtom) }};
      if(shell.O == center){
        this->mapSh2Cen_.push_back(iAtom+1);
        break;
      }
    } 
  }
  this->haveMapSh2Cen = true;
}

void BasisSet::makeMapCen2Bf(Molecule *mol){
  if(!this->haveMapSh2Bf ) this->makeMapSh2Bf();
  if(!this->haveMapSh2Cen) this->makeMapSh2Cen(mol);

/*
  for(auto iAtm = 0; iAtm < mol->nAtoms(); iAtm++){
    for(auto iShell = 0; iShell < this->nShell_; iShell++){
      if(iAtm == this->mapSh2Cen_[iShell]){
        this->mapCen2Bf_.push_back({{ this->mapSh2Bf_[iShell], this->shells_[iShell].size() }});
      }
    }
  }
*/
  cout << "HERE" << endl;
  for(auto iAtm = 0; iAtm < mol->nAtoms(); iAtm++){
    auto nSize = 0;
    for(auto iShell = 0; iShell < this->nShell_; iShell++){
      if((iAtm+1) == this->mapSh2Cen_[iShell]) nSize += this->shells_[iShell].size();
    }
    auto iSt = -1;
    for(auto iShell = 0; iShell < this->nShell_; iShell++){
      cout << this->mapSh2Cen_[iShell] << endl;
      if((iAtm+1) == this->mapSh2Cen_[iShell]){
       iSt = this->mapSh2Bf_[iShell];
       break;
      }
    }
    if(iSt == -1) CErr("Could not find Center in Basis definition",this->fileio_->out);
    this->mapCen2Bf_.push_back({{iSt,nSize}});
  }

  cout << "MC2B" << endl;
  for(auto i : this->mapCen2Bf_)
    cout << i[0] << '\t' << i[1] << endl;
  cout << endl << endl;

  this->haveMapCen2Bf = true;
  
}

template<>
void BasisSet::computeShBlkNorm(bool doBeta, const RealMatrix *DAlpha, 
                                   const RealMatrix *DBeta){
  // If map doesnt exist, make it
  if(!this->haveMapSh2Bf) this->makeMapSh2Bf();

  // Allocate Matricies
  this->shBlkNormAlpha = 
    std::unique_ptr<RealMatrix>(new RealMatrix(this->nShell_,this->nShell_));
  if(doBeta)
    this->shBlkNormBeta = 
      std::unique_ptr<RealMatrix>(new RealMatrix(this->nShell_,this->nShell_));

  for(int s1 = 0; s1 < this->nShell_; s1++) {
    int bf1 = this->mapSh2Bf_[s1];
    int n1  = this->shells_[s1].size();
    for(int s2 = 0; s2 < this->nShell_; s2++) {
      int bf2 = this->mapSh2Bf_[s2];
      int n2  = this->shells_[s2].size();
     
      (*this->shBlkNormAlpha)(s1,s2) = DAlpha->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
      if(doBeta)
        (*this->shBlkNormBeta)(s1,s2) = DBeta->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
    }
  }
} // computeShBlkNorm (TMat = RealMatrix)

template<>
void BasisSet::computeShBlkNorm(bool doBeta, const ComplexMatrix *DAlpha, 
                                   const ComplexMatrix *DBeta){
  // If map doesnt exist, make it
  if(!this->haveMapSh2Bf) this->makeMapSh2Bf();

  // Allocate Matricies
  this->shBlkNormAlpha = 
    std::unique_ptr<RealMatrix>(new RealMatrix(this->nShell_,this->nShell_));
  if(doBeta)
    this->shBlkNormBeta = 
      std::unique_ptr<RealMatrix>(new RealMatrix(this->nShell_,this->nShell_));

  for(int s1 = 0; s1 < this->nShell_; s1++) {
    int bf1 = this->mapSh2Bf_[s1];
    int n1  = this->shells_[s1].size();
    for(int s2 = 0; s2 < this->nShell_; s2++) {
      int bf2 = this->mapSh2Bf_[s2];
      int n2  = this->shells_[s2].size();
     
      (*this->shBlkNormAlpha)(s1,s2) = DAlpha->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
      if(doBeta)
        (*this->shBlkNormBeta)(s1,s2) = DBeta->block(bf1,bf2,n1,n2).lpNorm<Infinity>();
    }
  }
} // computeShBlkNorm (TMat = ComplexMatrix)

void BasisSet::constructExtrn(Molecule * mol, BasisSet *genBasis){
  genBasis->fileio_ = this->fileio_;
  for(auto iAtom = 0; iAtom < mol->nAtoms(); iAtom++){
    bool found = false;
    for(auto iRef = this->refShells_.begin(); iRef != this->refShells_.end(); ++iRef){
      if(mol->index(iAtom) == (*iRef).index){
        for(auto iShell = (*iRef).shells.begin(); iShell != (*iRef).shells.end(); ++iShell)
          genBasis->shells_.push_back(
            libint2::Shell{ 
              iShell->alpha, 
              {iShell->contr[0]}, 
              {{ (*mol->cart())(0,iAtom),
                 (*mol->cart())(1,iAtom),
                 (*mol->cart())(2,iAtom)}}
            }
          );
        found = true;
      }
    }
    if(!found)  
      CErr("Atomic Number " + 
             std::to_string(elements[mol->index(iAtom)].atomicNumber) +
             " not found in current Basis Set",
           this->fileio_->out);
  }
  genBasis->computeMeta();

}

void BasisSet::renormShells(){
  for(auto iShell = this->shells_.begin(); iShell != this->shells_.end(); ++iShell)
    iShell->renorm();
}

} // namespace ChronusQ

