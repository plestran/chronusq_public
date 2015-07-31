#include <basisset.h>
namespace ChronusQ{
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

void BasisSet::renormShells(){
  for(auto iShell = this->shells_.begin(); iShell != this->shells_.end(); ++iShell)
    iShell->renorm();
}


}; //namespace ChronusQ
