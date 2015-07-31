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
  int         contDepth;
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
          this->refShells_.push_back(ReferenceShell{atomicNumber,indx,tmpShell});
        }
        indx = HashAtom(tokens[0],0);
        atomicNumber = elements[indx].atomicNumber;
        newRec = false;
        firstRec = false;
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
      }
    }
  }
  // Append the last Rec
  this->refShells_.push_back(ReferenceShell{atomicNumber,indx,tmpShell}); 
 
};
} // namespace ChronusQ

