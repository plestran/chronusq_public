#include "aointegrals.h"
using ChronusQ::AOIntegrals;

void AOIntegrals::computeOverlapS(){
  double S;
  int i,j,k,ijShell,lA[3],lB[3];
  int bf1,bf2;
//  RealMatrix Overlap(this->nCartBasis_,this->nCartBasis_);
  std::vector<double> tmpS;
  ChronusQ::ShellPair *ijS;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(this->printLevel_>=1) start = std::chrono::high_resolution_clock::now();
  this->overlap_->setZero();
/*
int ll[3]={1,0,0};
int L=1;
int m = 1;
std::complex<double> aaa;
aaa = this->car2sphcoeff(L,m,ll);
cout<<"complex number"<<aaa<<"\t"<<aaa.real()<<"\t"<<aaa.imag()<<endl;
int ll2[3]={0,1,0};
 L=1;
 m = 1;
aaa = this->car2sphcoeff(L,m,ll2);
cout<<"complex number"<<aaa<<endl;
int ll3[3] ={0,0,1};
 L=1;
 m = 1;
aaa = this->car2sphcoeff(L,m,ll3);
cout<<"complex number"<<aaa<<endl;
int ll4[3]={1,0,0};
 L=1;
 m = 0;
aaa = this->car2sphcoeff(L,m,ll4);
cout<<"complex number"<<aaa<<endl;
int ll5[3]={0,1,0};
 L=1;
 m = 0;
aaa = this->car2sphcoeff(L,m,ll5);
cout<<"complex number"<<aaa<<endl;
int ll6[3]={0,0,1};
 L=1;
 m = 0;
aaa = this->car2sphcoeff(L,m,ll6);
cout<<"complex number"<<aaa<<endl;
int ll7[3]={1,0,0};
 L=1;
 m = -1;
aaa = this->car2sphcoeff(L,m,ll7);
cout<<"complex number"<<aaa<<endl;
int ll8[3]={0,1,0};
 L=1;
 m = -1;
aaa = this->car2sphcoeff(L,m,ll8);
cout<<"complex number"<<aaa<<endl;
int ll9[3]={0,0,1};
 L=1;
 m = -1;
aaa = this->car2sphcoeff(L,m,ll9);
cout<<"complex number"<<aaa<<endl;
*/









  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);
    if(ijS->iShell==ijS->jShell) {
      tmpS.resize(ijS->icarsize*ijS->jcarsize);
      for(i=0,bf1=ijS->icarbf_s  ; i<ijS->icarsize; i++,bf1++)  
      for(j=i,bf2=ijS->icarbf_s+i; j<ijS->jcarsize; j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };


        S = this->hRRSab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB);
/*        (*this->overlap_)(bf1,bf2) = S;
        (*this->overlap_)(bf2,bf1) = S;
*/

        if (this->basisSet_->getforceCart()==false) {
//          Overlap(bf1,bf2) =S;
//          Overlap(bf2,bf1) =S;
          tmpS[i*ijS->jcarsize+j] = S;
          tmpS[j*ijS->jcarsize+i] = S;

        }
        else if (this->basisSet_->getforceCart()==true) {
          (*this->overlap_)(bf1,bf2) = S;
          (*this->overlap_)(bf2,bf1) = S;
        }
      };
    } else {
      for(i=0,bf1=ijS->icarbf_s; i<ijS->icarsize; i++,bf1++) 
      for(j=0,bf2=ijS->jcarbf_s; j<ijS->jcarsize; j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        S = this->hRRSab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB);
/*        (*this->overlap_)(bf1,bf2) = S;
        (*this->overlap_)(bf2,bf1) = S;*/
        if (this->basisSet_->getforceCart()==false) {
/*          Overlap(bf1,bf2) =S;
          Overlap(bf2,bf1) =S;
*/
        tmpS.push_back(S);
        }
        else if (this->basisSet_->getforceCart()==true) {
          (*this->overlap_)(bf1,bf2) = S;
          (*this->overlap_)(bf2,bf1) = S;
        }
 
      };
    };
    if (this->basisSet_->getforceCart()==false) {
      std::vector<double> SPH = this->cart2SphTrans(ijS,tmpS.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->overlap_)(bf1,bf2) = SPH[i*ijS->jsphsize+j];
          (*this->overlap_)(bf2,bf1) = SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpS.clear(); 
    }
/*
 //beginning of cartesian to spehrical transform
      if (this->nBasis_==this->nSphBasis_) {
        transdimi = ijS->isphsize;
        transdimj = ijS->jsphsize;
      
        RealMatrix Csp2cai(transdimi ,ijS->iShell->cartesian_l.size());
        RealMatrix Csp2caj(transdimj ,ijS->jShell->cartesian_l.size());
        for( p = 0 ; p<ijS->isphsize; p++ ) {
          for( q = 0 ; q<ijS->icarsize; q++ ) {
            Csp2cai(p,q) = 0;
          }
        }
        for ( p = 0 ; p <ijS->jsphsize; p++ ) {
          for ( q = 0 ; q <ijS->jcarsize; q++ ) {
            Csp2caj(p,q) = 0;
          }
        }
        if (ijS->iShell->l==0) {
          Csp2cai(0,0) = 1;
        }
        else if (ijS->iShell->l==1) {
          Csp2cai(0,0) =1;
          Csp2cai(1,1) =1;
          Csp2cai(2,2) =1;

        }
        else if (ijS->iShell->l==2) {
          Csp2cai(4,0) = sqrt(3)/2;
          Csp2cai(4,3) = -sqrt(3)/2;
      
          Csp2cai(3,2) = sqrt(3);

          Csp2cai(2,0) = -0.5;
          Csp2cai(2,3) = -0.5;
          Csp2cai(2,5) = 1;
  
          Csp2cai(1,4) = sqrt(3);

          Csp2cai(0,1) = sqrt(3);
       
        }
        if (ijS->jShell->l==0) {
          Csp2caj(0,0) = 1;
        }
        else if (ijS->jShell->l==1) {
          Csp2caj(0,0) =1;
          Csp2caj(1,1) =1;
          Csp2caj(2,2) =1;
        }
        else if (ijS->jShell->l==2) {
          Csp2caj(4,0) = sqrt(3)/2;
          Csp2caj(4,3) = -sqrt(3)/2;
      
          Csp2caj(3,2) = sqrt(3);

          Csp2caj(2,0) = -0.5;
          Csp2caj(2,3) = -0.5;
          Csp2caj(2,5) = 1;

          Csp2caj(1,4) = sqrt(3);

          Csp2caj(0,1) = sqrt(3);
        }


        for( i = 0, bf1 = ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++  ) {
          for( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++,bf2++ ) {
            S = 0;
            for( p = 0, cartbf1 = ijS->icarbf_s; p<ijS->icarsize; p++, cartbf1++) {
              for( q = 0, cartbf2 = ijS->jcarbf_s; q<ijS->jcarsize; q++, cartbf2++) {
                S += Csp2cai(i,p)*Csp2caj(j,q)*Overlap(cartbf1,cartbf2);
              }
            }
            (*this->overlap_)(bf1,bf2) = S;
            (*this->overlap_)(bf2,bf1) = S;
          }
        }
      }    // end of cartesian to spherical transform
*/
  };
  prettyPrint(this->fileio_->out,(*this->overlap_),"XSLI Overlap");
//cout<<"overlap finished"<<endl;
};




void AOIntegrals::computeKineticT(){
  double T;
  int i,j,k,ijShell,lA[3],lB[3];
  int bf1,bf2;
//  RealMatrix Kinetic(this->nCartBasis_,this->nCartBasis_);
  std::vector<double> tmpT;
  ChronusQ::ShellPair *ijS;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(this->printLevel_>=1) start = std::chrono::high_resolution_clock::now();
  this->potential_->setZero();
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);

    if(ijS->iShell==ijS->jShell) {
      tmpT.resize( ijS->icarsize*ijS->jcarsize );
      for(i=0,bf1=ijS->icarbf_s  ; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=i,bf2=ijS->icarbf_s+i; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        T = this->RRTab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB);
//        cout<<"VV: "<<bf1<<" "<<bf2<<" "<<V<<endl;
/*        (*this->kinetic_)(bf1,bf2) = T;
        (*this->kinetic_)(bf2,bf1) = T;
*/
 
        if (this->basisSet_->getforceCart()==false) {
          tmpT[i*ijS->jcarsize+j] = T;
          tmpT[j*ijS->jcarsize+i] = T;
        }
        else if (this->basisSet_->getforceCart()==true) {
          (*this->kinetic_)(bf1,bf2) = T;
          (*this->kinetic_)(bf2,bf1) = T;
        }
      }
    } else {
      for(i=0,bf1=ijS->icarbf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jcarbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        T = this->RRTab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB);
//        cout<<"T: "<<bf1<<" "<<bf2<<" "<<T<<endl;
//        (*this->kinetic_)(bf1,bf2) = T;
//        (*this->kinetic_)(bf2,bf1) = T;

        if (this->basisSet_->getforceCart()==false) {
//          Kinetic(bf1,bf2) = T;
//          Kinetic(bf2,bf1) = T;
          tmpT.push_back(T);
        }
        else if (this->basisSet_->getforceCart()==true) {
          (*this->kinetic_)(bf1,bf2) = T;
          (*this->kinetic_)(bf2,bf1) = T;
        }
 
      };
    };
  //beginning of cartesian to spehrical transform
      if (this->basisSet_->getforceCart()==false) {
        std::vector<double> SPH = this->cart2SphTrans(ijS,tmpT.data() );
        for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
          for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
            (*this->kinetic_)(bf1,bf2) = SPH[i*ijS->jsphsize+j];
            (*this->kinetic_)(bf2,bf1) = SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpT.clear(); 
    }

//        transdimi = ijS->isphsize;
//        transdimj = ijS->jsphsize;
      
/*        RealMatrix Csp2cai(ijS->isphsize ,ijS->iShell->cartesian_l.size());
        RealMatrix Csp2caj(ijS->jsphsize ,ijS->jShell->cartesian_l.size());
        for( p = 0 ; p<ijS->isphsize; p++ ) {
          for( q = 0 ; q<ijS->icarsize; q++ ) {
            Csp2cai(p,q) = 0;
          }
        }
        for ( p = 0 ; p <ijS->jsphsize; p++ ) {
          for ( q = 0 ; q <ijS->jcarsize; q++ ) {
            Csp2caj(p,q) = 0;
          }
        }
 
        if (ijS->iShell->l==0) {
          Csp2cai(0,0) = 1;
        }
        else if (ijS->iShell->l==1) {
          Csp2cai(0,0) =1;
          Csp2cai(1,1) =1;
          Csp2cai(2,2) =1;
        }
        else if (ijS->iShell->l==2) {
          Csp2cai(4,0) = sqrt(3)/2;
          Csp2cai(4,3) = -sqrt(3)/2;
      
          Csp2cai(3,2) = sqrt(3);

          Csp2cai(2,0) = -0.5;
          Csp2cai(2,3) = -0.5;
          Csp2cai(2,5) = 1;

          Csp2cai(1,4) = sqrt(3);

          Csp2cai(0,1) = sqrt(3);
        }
        if (ijS->jShell->l==0) {
          Csp2caj(0,0) = 1;
        }
        else if (ijS->jShell->l==1) {
          Csp2caj(0,0) =1;
          Csp2caj(1,1) =1;
          Csp2caj(2,2) =1;
        }
        else if (ijS->jShell->l==2) {
          Csp2caj(4,0) = sqrt(3)/2;
          Csp2caj(4,3) = -sqrt(3)/2;
      
          Csp2caj(3,2) = sqrt(3);

          Csp2caj(2,0) = -0.5;
          Csp2caj(2,3) = -0.5;
          Csp2caj(2,5) = 1;

          Csp2caj(1,4) = sqrt(3);

          Csp2caj(0,1) = sqrt(3);
        }


        for( i = 0, bf1 = ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++  ) {
          for( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++,bf2++ ) {
            T = 0;
            for( p = 0, cartbf1 = ijS->icarbf_s; p<ijS->iShell->cartesian_l.size(); p++, cartbf1++) {
              for( q = 0, cartbf2 = ijS->jcarbf_s; q<ijS->jShell->cartesian_l.size(); q++, cartbf2++) {
                T += Csp2cai(i,p)*Csp2caj(j,q)*Kinetic(cartbf1,cartbf2);
              }
            }
            (*this->kinetic_)(bf1,bf2) = T;
            (*this->kinetic_)(bf2,bf1) = T;
          }
        }
      }    // end of cartesian to spherical transform
*/
  };
  prettyPrint(this->fileio_->out,(*this->kinetic_),"XSLI Kinetic");
};

void AOIntegrals::computePotentialV(){
/*
  double pVp, pVpC, C[3];     //C is the coordinate of atoms
int i,j,k,ijShell,lA[3],lB[3],iPP,iAtom,nu;
  int bf1,bf2;
  ChronusQ::ShellPair *ijS;
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);

      for(i=0,bf1=ijS->ibf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };

        pVp = 0;
        for ( iPP=0; iPP<ijS->nPGTOPair; iPP++) {
          pVpC = 0;
          for (iAtom = 0; iAtom <this->molecularConstants_->nAtom; iAtom++){
             for ( nu = 0; nu <3 ; nu++) {
               C[nu] = this->molecularConstants_->cart[nu][iAtom];
             }
             pVpC += this->molecularConstants_->atomZ[iAtom]*this->hRRiPPVab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB,C,0,iPP);

           }
         pVp += pVpC;
        }  
        (*this->potential_)(bf1,bf2) = -pVp; 
        (*this->potential_)(bf2,bf1) = -pVp;
      };
    };
 
  prettyPrint(this->fileio_->out,(*this->potential_),"potential");
  
  double fobi;
  fobi = this->potential_->frobInner(*this->potential_);
  this->fileio_->out <<"fobinius inner product" << fobi << endl;

};

*/


  double V;
  int i,j,k,ijShell,lA[3],lB[3];
  int bf1,bf2;
//  RealMatrix Potential(this->nCartBasis_,this->nCartBasis_);
  std::vector<double> tmpV;
  ChronusQ::ShellPair *ijS;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(this->printLevel_>=1) start = std::chrono::high_resolution_clock::now();
  this->potential_->setZero();
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
//    cout<<"ShellPairs"<<endl;
    ijS = &(this->shellPairs_[ijShell]);

    if(ijS->iShell==ijS->jShell) {
      tmpV.resize(ijS->icarsize*ijS->jcarsize);
      for(i=0,bf1=ijS->icarbf_s  ; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=i,bf2=ijS->icarbf_s+i; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        V = this->hRRVab(ijS,this->molecularConstants_.get(),ijS->iShell->l,lA,ijS->jShell->l,lB);
//        (*this->potential_)(bf1,bf2) = -V;
//        (*this->potential_)(bf2,bf1) = -V;
//        Potential(bf1,bf2) = -V;
//        Potential(bf2,bf1) = -V;
        if (this->basisSet_->getforceCart() == false ) {
          tmpV[i*ijS->jcarsize+j] = -V;
          tmpV[j*ijS->jcarsize+i] = -V;
        }
        else if (this->basisSet_->getforceCart() ==true) {
          (*this->potential_)(bf1,bf2) = -V; 
          (*this->potential_)(bf2,bf1) = -V;
        }
      }
    } else {
      for(i=0,bf1=ijS->icarbf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jcarbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        V = this->hRRVab(ijS,this->molecularConstants_.get(),ijS->iShell->l,lA,ijS->jShell->l,lB);
        if (this->basisSet_->getforceCart()==false) {
//        (*this->potential_)(bf1,bf2) = -V; 
//        (*this->potential_)(bf2,bf1) = -V;
          tmpV.push_back(-V);
        }
        else if (this->basisSet_->getforceCart() ==true ) {
          (*this->potential_)(bf1,bf2) = -V; 
          (*this->potential_)(bf2,bf1) = -V;
        }
 
      };
    };
    if (this->basisSet_->getforceCart()==false) {
      std::vector<double> SPH = this->cart2SphTrans(ijS,tmpV.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->potential_)(bf1,bf2) = SPH[i*ijS->jsphsize+j];
          (*this->potential_)(bf2,bf1) = SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpV.clear(); 
    }
/*
  //beginning of cartesian to spehrical transform
      if (this->nBasis_==this->nSphBasis_) {
        RealMatrix Csp2cai( ijS->isphsize ,ijS->iShell->cartesian_l.size());
        RealMatrix Csp2caj( ijS->jsphsize ,ijS->jShell->cartesian_l.size());
        for( p = 0 ; p<ijS->isphsize; p++ ) {
          for( q = 0 ; q<ijS->icarsize; q++ ) {
            Csp2cai(p,q) = 0;
          }
        }
        for ( p = 0 ; p <ijS->jsphsize; p++ ) {
          for ( q = 0 ; q <ijS->jcarsize; q++ ) {
            Csp2caj(p,q) = 0;
          }
        }
 
        if (ijS->iShell->l==0) {
          Csp2cai(0,0) = 1;
        }
        else if (ijS->iShell->l==1) {
          Csp2cai(0,0) =1;
          Csp2cai(1,1) =1;
          Csp2cai(2,2) =1;
        }
        else if (ijS->iShell->l==2) {
          Csp2cai(4,0) = sqrt(3)/2;
          Csp2cai(4,3) = -sqrt(3)/2;
      
          Csp2cai(3,2) = sqrt(3);

          Csp2cai(2,0) = -0.5;
          Csp2cai(2,3) = -0.5;
          Csp2cai(2,5) = 1;

          Csp2cai(1,4) = sqrt(3);

          Csp2cai(0,1) = sqrt(3);
        }
        if (ijS->jShell->l==0) {
          Csp2caj(0,0) = 1;
        }
        else if (ijS->jShell->l==1) {
          Csp2caj(0,0) =1;
          Csp2caj(1,1) =1;
          Csp2caj(2,2) =1;
        }
        else if (ijS->jShell->l==2) {
          Csp2caj(4,0) = sqrt(3)/2;
          Csp2caj(4,3) = -sqrt(3)/2;
      
          Csp2caj(3,2) = sqrt(3);

          Csp2caj(2,0) = -0.5;
          Csp2caj(2,3) = -0.5;
          Csp2caj(2,5) = 1;

          Csp2caj(1,4) = sqrt(3);

          Csp2caj(0,1) = sqrt(3);
        }


        for( i = 0, bf1 = ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++  ) {
          for( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++,bf2++ ) {
            V = 0;
            for( p = 0, cartbf1 = ijS->icarbf_s; p<ijS->iShell->cartesian_l.size(); p++, cartbf1++) {
              for( q = 0, cartbf2 = ijS->jcarbf_s; q<ijS->jShell->cartesian_l.size(); q++, cartbf2++) {
                V += Csp2cai(i,p)*Csp2caj(j,q)*Potential(cartbf1,cartbf2);
              }
            }
            (*this->potential_)(bf1,bf2) = V;
            (*this->potential_)(bf2,bf1) = V;
          }
        }
     }
*/
  };
  prettyPrint(this->fileio_->out,(*this->potential_),"XSLI Potential");
};



//SS
void AOIntegrals::computeAngularL(){
  double L[3];
  int i,j,k,ijShell,mu,bf1,bf2,lA[3],lB[3],iPP;
//  RealTensor3d Angular(this->nCartBasis_,this->nCartBasis_,3);
  std::vector<double> tmpLx;
  std::vector<double> tmpLy;
  std::vector<double> tmpLz;
  ChronusQ::ShellPair *ijS;
//  this->angular_->setZero();  //not sure if setZero work
  RealTensor3d OneiOnek(3,3,3); //(1i,1k,mu)
  OneiOnek(0,0,0) = math.zero; 
  OneiOnek(0,0,1) = math.zero; 
  OneiOnek(0,0,2) = math.zero; 
 
  OneiOnek(1,1,0) = math.zero; 
  OneiOnek(1,1,1) = math.zero; 
  OneiOnek(1,1,2) = math.zero; 

  OneiOnek(2,2,0) = math.zero; 
  OneiOnek(2,2,1) = math.zero; 
  OneiOnek(2,2,2) = math.zero; 

  OneiOnek(0,1,0) = math.zero; 
  OneiOnek(0,1,1) = math.zero; 
  OneiOnek(0,1,2) = math.one;

  OneiOnek(0,2,0) = math.zero; 
  OneiOnek(0,2,1) = -math.one;
  OneiOnek(0,2,2) = math.zero; 
 
  OneiOnek(1,0,0) = math.zero; 
  OneiOnek(1,0,1) = math.zero; 
  OneiOnek(1,0,2) = -math.one;

  OneiOnek(1,2,0) = math.one; 
  OneiOnek(1,2,1) = math.zero; 
  OneiOnek(1,2,2) = math.zero;
 
  OneiOnek(2,0,0) = math.zero; 
  OneiOnek(2,0,1) = math.one; 
  OneiOnek(2,0,2) = math.zero; 

 
  OneiOnek(2,1,0) = -math.one; 
  OneiOnek(2,1,1) = math.zero; 
  OneiOnek(2,1,2) = math.zero; 
/*
  RealMatrix OneixBC(3,3);//OneixBC(i,mu)
  OneixBC(0,0) = math.zero;
  OneixBC(0,1) = -ijS->B[2];
  cout<<ijS->B[2]<<"\t";
  OneixBC(0,2) = ijS->B[1];
  cout<<ijS->B[1]<<"\t";
  OneixBC(1,0) = ijS->B[2];
  OneixBC(1,1) = math.zero;
  OneixBC(1,2) = -(ijS->B)[0];
  cout<<ijS->B[0]<<"\t"<<endl;
  OneixBC(2,0) = -ijS->B[1];
  OneixBC(2,1) = ijS->B[0];
  OneixBC(2,2) = math.zero;
*/ 
/*  for(q = 0;q<3;q++) {
    for(r = 0;r<3;r++) {
        cout<<OneixBC(q,r)<<"\t";
    }
    cout<<endl;
  }
  cout<<endl;   
*/    
/*                                    
  RealMatrix OneixAC(3,3);//OneixAC(i,mu)
  OneixAC(0,0) = math.zero;
  OneixAC(0,1) = -ijS->A[2];
  OneixAC(0,2) = ijS->A[1];
  OneixAC(1,0) = ijS->A[2];      
  OneixAC(1,1) = math.zero;         
  OneixAC(1,2) = -ijS->A[0];         
  OneixAC(2,0) = -(ijS->A)[1];      
  OneixAC(2,1) = ijS->A[0];      
  OneixAC(2,2) = math.zero;         
*/
/*  for(q = 0;q<3;q++) {
    for(r = 0;r<3;r++) {
        cout<<OneixAC(q,r)<<"\t";
    }
    cout<<endl;
  }
  cout<<endl;   
*/ 
                                 
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);
/*    if(ijS->iShell==ijS->jShell) {
      for(i=0, bf1 = ijS->ibf_s  ; i < ijS->iShell->cartesian_l.size(); i++,bf1++)
      for(j=i, bf2 = ijS->ibf_s+i; j < ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k = 0; k < 3; k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];          
        }
        for(mu = 0; mu < 3; mu++) {
          L[mu]=0;
          LL[mu]=0;
          for(iPP=0; iPP<ijS->nPGTOPair; iPP++){ 
            L[mu] += this->Labmu(ijS,&OneixAC,&OneixBC,&OneiOnek,ijS->iShell->l,lA,ijS->jShell->l,lB,mu,iPP);
            LL[mu] += this->Labmu(ijS,&OneixBC,&OneixAC,&OneiOnek,ijS->jShell->l,lB,ijS->iShell->l,lA,mu,iPP);
          };
        (*this->angular_)(bf1,bf2,mu) =  L[mu];
        (*this->angular_)(bf2,bf1,mu) =  LL[mu];
        } 
      }
    } else { */
/*
  tmpLx.resize(ijS->icarsize*ijS->jcarsize);
  tmpLy.resize(ijS->icarsize*ijS->jcarsize);
  tmpLz.resize(ijS->icarsize*ijS->jcarsize);
*/
  RealMatrix OneixBC(3,3);//OneixBC(i,mu)
  OneixBC(0,0) = math.zero;
  OneixBC(0,1) = -ijS->B[2];
//  cout<<ijS->B[2]<<"\t";
  OneixBC(0,2) = ijS->B[1];
//  cout<<ijS->B[1]<<"\t";
  OneixBC(1,0) = ijS->B[2];
  OneixBC(1,1) = math.zero;
  OneixBC(1,2) = -(ijS->B)[0];
//  cout<<ijS->B[0]<<"\t"<<endl;
  OneixBC(2,0) = -ijS->B[1];
  OneixBC(2,1) = ijS->B[0];
  OneixBC(2,2) = math.zero;
 
/*  for(q = 0;q<3;q++) {
    for(r = 0;r<3;r++) {
        cout<<OneixBC(q,r)<<"\t";
    }
    cout<<endl;
  }
  cout<<endl;   
*/    
                                    
  RealMatrix OneixAC(3,3);//OneixAC(i,mu)
  OneixAC(0,0) = math.zero;
  OneixAC(0,1) = -ijS->A[2];
  OneixAC(0,2) = ijS->A[1];
  OneixAC(1,0) = ijS->A[2];      
  OneixAC(1,1) = math.zero;         
  OneixAC(1,2) = -ijS->A[0];         
  OneixAC(2,0) = -(ijS->A)[1];      
  OneixAC(2,1) = ijS->A[0];      
  OneixAC(2,2) = math.zero;         




 
       for(i=0,bf1=ijS->icarbf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++)
       for(j=0,bf2=ijS->jcarbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
         for(k=0;k<3;k++){
           lA[k]=ijS->iShell->cartesian_l[i][k];
           lB[k]=ijS->jShell->cartesian_l[j][k];
         }
         for(mu=0;mu<3;mu++) {
           L[mu] = 0;
//           LL[mu] = 0;
           for(iPP=0; iPP<ijS->nPGTOPair; iPP++){
             L[mu] += this->Labmu(ijS,&OneixAC,&OneixBC,&OneiOnek,ijS->iShell->l,lA,ijS->jShell->l,lB,mu,iPP);
//             LL[mu] += this->Labmu(ijS,&OneixBC,&OneixAC,&OneiOnek,ijS->jShell->l,lB,ijS->iShell->l,lA,mu,iPP);
           }
             
//           Angular(bf1,bf2,mu) = L[mu];
//           Angular(bf2,bf1,mu) = -L[mu];

           if (this->basisSet_->getforceCart() == true) {
             (*this->angular_)(bf1,bf2,mu)=L[mu];
             (*this->angular_)(bf2,bf1,mu)=-L[mu];
           }
         }

         if (this->basisSet_->getforceCart()==false) {
           tmpLx.push_back(L[0]);
           tmpLy.push_back(L[1]);
           tmpLz.push_back(L[2]);
       }
     }
// cartesian to spherical transform
    if (this->basisSet_->getforceCart()==false) {
      std::vector<double> SPH = this->cart2SphTrans(ijS,tmpLx.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->angular_)(bf1,bf2,0) = SPH[i*ijS->jsphsize+j];
//cout<<(*this->angular_)(bf1,bf2,0)<<endl;
          (*this->angular_)(bf2,bf1,0) = -SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpLx.clear(); 
      SPH = this->cart2SphTrans(ijS,tmpLy.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->angular_)(bf1,bf2,1) = SPH[i*ijS->jsphsize+j];
          (*this->angular_)(bf2,bf1,1) = -SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpLy.clear(); 
      SPH = this->cart2SphTrans(ijS,tmpLz.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->angular_)(bf1,bf2,2) = SPH[i*ijS->jsphsize+j];
          (*this->angular_)(bf2,bf1,2) = -SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpLz.clear(); 
 
    }

  //beginning of cartesian to spehrical transform
/*      if (this->nBasis_==this->nSphBasis_) {
        RealMatrix Csp2cai(ijS->isphsize ,ijS->iShell->cartesian_l.size());
        RealMatrix Csp2caj(ijS->jsphsize ,ijS->jShell->cartesian_l.size());

        for( mm = 0 ; mm<ijS->isphsize; mm++ ) {
          for( nn = 0 ; nn<ijS->icarsize; nn++ ) {
            Csp2cai(mm,nn) = 0;
          }
        }
        for ( mm = 0 ; mm <ijS->jsphsize; mm++ ) {
          for ( nn = 0 ; nn <ijS->jcarsize; nn++ ) {
            Csp2caj(mm,nn) = 0;
          }
        }
 
        if (ijS->iShell->l==0) {
          Csp2cai(0,0) = 1;
        }
        else if (ijS->iShell->l==1) {
          Csp2cai(0,0) =1;
          Csp2cai(1,1) =1;
          Csp2cai(2,2) =1;
        }
        else if (ijS->iShell->l==2) {
          Csp2cai(4,0) = sqrt(3)/2;
          Csp2cai(4,3) = -sqrt(3)/2;
      
          Csp2cai(3,2) = sqrt(3);

          Csp2cai(2,0) = -0.5;
          Csp2cai(2,3) = -0.5;
          Csp2cai(2,5) = 1;

          Csp2cai(1,4) = sqrt(3);

          Csp2cai(0,1) = sqrt(3);
        }
        if (ijS->jShell->l==0) {
          Csp2caj(0,0) = 1;
        }
        else if (ijS->jShell->l==1) {
          Csp2caj(0,0) =1;
          Csp2caj(1,1) =1;
          Csp2caj(2,2) =1;
        }
        else if (ijS->jShell->l==2) {
          Csp2caj(4,0) = sqrt(3)/2;
          Csp2caj(4,3) = -sqrt(3)/2;
      
          Csp2caj(3,2) = sqrt(3);

          Csp2caj(2,0) = -0.5;
          Csp2caj(2,3) = -0.5;
          Csp2caj(2,5) = 1;

          Csp2caj(1,4) = sqrt(3);

          Csp2caj(0,1) = sqrt(3);
        }

        for( mu = 0 ; mu < 3 ; mu++ ) {
          for( i = 0, bf1 = ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++  ) {
            for( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++,bf2++ ) {
              L[mu] = 0;
              for( mm = 0, cartbf1 = ijS->icarbf_s; mm<ijS->iShell->cartesian_l.size(); mm++, cartbf1++) {
                for( nn = 0, cartbf2 = ijS->jcarbf_s; nn<ijS->jShell->cartesian_l.size(); nn++, cartbf2++) {
                  L[mu] += Csp2cai(i,mm)*Csp2caj(j,nn)*Angular(cartbf1,cartbf2,mu);
                 }
               }
             (*this->angular_)(bf1,bf2,mu) = L[mu];
             (*this->angular_)(bf2,bf1,mu) = -L[mu];
             }
          }
       }  
     }  */   // end of cartesian to spehrical transform
   }
  for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
    RealMap TEMP(&this->angular_->storage()[iXYZ*this->nBasis_*this->nBasis_],this->nBasis_,this->nBasis_);
    prettyPrint(this->fileio_->out,TEMP,"Angular XYZ="+std::to_string(iXYZ));
  }
//  prettyPrint(this->fileio_->out,(*this->kinetic_),"SS Angular");
};


void AOIntegrals::computeSL(){
  double Sl[3],SlC,C[3];
  int i,j,k,ijShell,mu,nu,lA[3],lB[3],iPP,iAtom;
  int bf1,bf2;
//  RealTensor3d SpinOrbit(this->nCartBasis_,this->nCartBasis_,3);
  std::vector<double> tmpSLx;
  std::vector<double> tmpSLy;
  std::vector<double> tmpSLz;
  ChronusQ::ShellPair *ijS;
//cout<<"Here we stop";
// precalculate the tensor
  RealTensor3d OneiOnek(3,3,3); //(1i,1k,mu)
  OneiOnek(0,0,0) = math.zero; 
  OneiOnek(0,0,1) = math.zero; 
  OneiOnek(0,0,2) = math.zero; 
 
  OneiOnek(1,1,0) = math.zero; 
  OneiOnek(1,1,1) = math.zero; 
  OneiOnek(1,1,2) = math.zero; 

  OneiOnek(2,2,0) = math.zero; 
  OneiOnek(2,2,1) = math.zero; 
  OneiOnek(2,2,2) = math.zero; 

  OneiOnek(0,1,0) = math.zero; 
  OneiOnek(0,1,1) = math.zero; 
  OneiOnek(0,1,2) = math.one;

  OneiOnek(0,2,0) = math.zero; 
  OneiOnek(0,2,1) = -math.one;
  OneiOnek(0,2,2) = math.zero; 
 
  OneiOnek(1,0,0) = math.zero; 
  OneiOnek(1,0,1) = math.zero; 
  OneiOnek(1,0,2) = -math.one;

  OneiOnek(1,2,0) = math.one; 
  OneiOnek(1,2,1) = math.zero; 
  OneiOnek(1,2,2) = math.zero;
 
  OneiOnek(2,0,0) = math.zero; 
  OneiOnek(2,0,1) = math.one; 
  OneiOnek(2,0,2) = math.zero; 

 
  OneiOnek(2,1,0) = -math.one; 
  OneiOnek(2,1,1) = math.zero; 
  OneiOnek(2,1,2) = math.zero; 
  
//cout<<"OneiOnek"<<endl;
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {

    ijS = &(this->shellPairs_[ijShell]);
/*    if(ijS->iShell==ijS->jShell) {
      for(i=0,bf1=ijS->ibf_s  ; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=i,bf2=ijS->ibf_s+i; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        V = this->hRRVab(ijS,this->molecularConstants_.get(),ijS->iShell->l,lA,ijS->jShell->l,lB);
        (*this->potential_)(bf1,bf2) = -V;
        (*this->potential_)(bf2,bf1) = -V;
      };
    } else {  */
      for(i=0,bf1=ijS->icarbf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jcarbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        for (mu = 0; mu < 3; mu++) {
          Sl[mu] = 0;
          for (iPP = 0; iPP<ijS->nPGTOPair; iPP++){
            SlC = 0;
            for (iAtom = 0; iAtom <this->molecularConstants_->nAtom; iAtom++){
              // get the coordinate of C
              for (nu = 0; nu<3 ; nu++){
                C[nu] = this->molecularConstants_->cart[nu][iAtom];
              }
              RealMatrix OneixAC(3,3);//OneixAC(i,mu)
              OneixAC(0,0) = math.zero;
              OneixAC(0,1) = -(ijS->A[2] - C[2]);
              OneixAC(0,2) = ijS->A[1] - C[1];
              OneixAC(1,0) = ijS->A[2] - C[2];      
              OneixAC(1,1) = math.zero;         
              OneixAC(1,2) = -(ijS->A[0] - C[0]);         
              OneixAC(2,0) = -((ijS->A)[1] - C[1]);      
              OneixAC(2,1) = ijS->A[0] - C[0];      
              OneixAC(2,2) = math.zero;         


              RealMatrix OneixBC(3,3);//OneixBC(i,mu)
              OneixBC(0,0) = math.zero;
              OneixBC(0,1) = -(ijS->B[2] - C[2]);
              OneixBC(0,2) = ijS->B[1] - C[1];
              OneixBC(1,0) = ijS->B[2] - C[2];
              OneixBC(1,1) = math.zero;
              OneixBC(1,2) = -((ijS->B)[0] - C[0]);
              OneixBC(2,0) = -(ijS->B[1] - C[1]);
              OneixBC(2,1) = ijS->B[0] - C[0];
              OneixBC(2,2) = math.zero;

              SlC += this->molecularConstants_->atomZ[iAtom]*this->Slabmu(ijS,&OneixAC,&OneixBC,&OneiOnek,C,ijS->iShell->l,lA,ijS->jShell->l,lB,mu,0,iPP,iAtom);
              };

            Sl[mu] +=SlC;
            };
//          (*this->SOCoupling_)(bf1,bf2,mu) = Slmu;
//          (*this->SOCoupling_)(bf2,bf1,mu) = -Slmu;
//cout<<"SOCoupling\tbf1\t"<<bf1<<"\tbf2\t"<<bf2<<"\tSOCoupling\t"<<Slmu<<endl;
//            SpinOrbit(bf1,bf2,mu) = Slmu;
//            SpinOrbit(bf2,bf1,mu) = -Slmu;
            if (this->basisSet_->getforceCart() == true) {
              (*this->SOCoupling_)(bf1,bf2,mu) = Sl[mu]; 
              (*this->SOCoupling_)(bf2,bf1,mu) = -Sl[mu];
            }
        }
        if (this->basisSet_->getforceCart() == false) {
          tmpSLx.push_back(Sl[0]);
          tmpSLy.push_back(Sl[1]);
          tmpSLz.push_back(Sl[2]);
        }
      }
    if (this->basisSet_->getforceCart()==false) {
      std::vector<double> SPH = this->cart2SphTrans(ijS,tmpSLx.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->SOCoupling_)(bf1,bf2,0) = SPH[i*ijS->jsphsize+j];
//cout<<(*this->angular_)(bf1,bf2,0)<<endl;
          (*this->SOCoupling_)(bf2,bf1,0) = -SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpSLx.clear(); 
      SPH = this->cart2SphTrans(ijS,tmpSLy.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->SOCoupling_)(bf1,bf2,1) = SPH[i*ijS->jsphsize+j];
          (*this->SOCoupling_)(bf2,bf1,1) = -SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpSLy.clear(); 
      SPH = this->cart2SphTrans(ijS,tmpSLz.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->SOCoupling_)(bf1,bf2,2) = SPH[i*ijS->jsphsize+j];
          (*this->SOCoupling_)(bf2,bf1,2) = -SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmpSLz.clear(); 
 
    }


  //beginning of cartesian to spehrical transform
/*      if (this->nBasis_==this->nSphBasis_) {
        RealMatrix Csp2cai(ijS->isphsize ,ijS->iShell->cartesian_l.size());
        RealMatrix Csp2caj(ijS->jsphsize ,ijS->jShell->cartesian_l.size());

        for( mm = 0 ; mm<ijS->isphsize; mm++ ) {
          for( nn = 0 ; nn<ijS->icarsize; nn++ ) {
            Csp2cai(mm,nn) = 0;
          }
        }
        for ( mm = 0 ; mm <ijS->jsphsize; mm++ ) {
          for ( nn = 0 ; nn <ijS->jcarsize; nn++ ) {
            Csp2caj(mm,nn) = 0;
          }
        }
 
        if (ijS->iShell->l==0) {
          Csp2cai(0,0) = 1;
        }
        else if (ijS->iShell->l==1) {
          Csp2cai(0,0) =1;
          Csp2cai(1,1) =1;
          Csp2cai(2,2) =1;
        }
        else if (ijS->iShell->l==2) {
//gaussian          Csp2cai(3,0) = sqrt(3)/2;
          Csp2cai(3,3) = -sqrt(3)/2;
      
          Csp2cai(1,2) = sqrt(3);

          Csp2cai(0,0) = -0.5;
          Csp2cai(0,3) = -0.5;
          Csp2cai(0,5) = 1;

          Csp2cai(2,4) = sqrt(3);

          Csp2cai(4,1) = sqrt(3);
//gaussian
          Csp2cai(4,0) = sqrt(3)/2;
          Csp2cai(4,3) = -sqrt(3)/2;
      
          Csp2cai(3,2) = sqrt(3);

          Csp2cai(2,0) = -0.5;
          Csp2cai(2,3) = -0.5;
          Csp2cai(2,5) = 1;

          Csp2cai(1,4) = sqrt(3);

          Csp2cai(0,1) = sqrt(3);


        }
        if (ijS->jShell->l==0) {
          Csp2caj(0,0) = 1;
        }
        else if (ijS->jShell->l==1) {
          Csp2caj(0,0) =1;
          Csp2caj(1,1) =1;
          Csp2caj(2,2) =1;
        }
        else if (ijS->jShell->l==2) {
//gaussian          Csp2caj(3,0) = sqrt(3)/2;
          Csp2caj(3,3) = -sqrt(3)/2;
      
          Csp2caj(1,2) = sqrt(3);

          Csp2caj(0,0) = -0.5;
          Csp2caj(0,3) = -0.5;
          Csp2caj(0,5) = 1;

          Csp2caj(2,4) = sqrt(3);

          Csp2caj(4,1) = sqrt(3);
//gaussian
          Csp2caj(4,0) = sqrt(3)/2;
          Csp2caj(4,3) = -sqrt(3)/2;
       
          Csp2caj(3,2) = sqrt(3);

          Csp2caj(2,0) = -0.5;
          Csp2caj(2,3) = -0.5;
          Csp2caj(2,5) = 1;

          Csp2caj(1,4) = sqrt(3);

          Csp2caj(0,1) = sqrt(3);

      }

      for( mu = 0 ; mu < 3 ; mu++ ) {
        for( i = 0, bf1 = ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++  ) {
          for( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++,bf2++ ) {
            Slmu = 0;
            for( mm = 0, cartbf1 = ijS->icarbf_s; mm<ijS->iShell->cartesian_l.size(); mm++, cartbf1++) {
              for( nn = 0, cartbf2 = ijS->jcarbf_s; nn<ijS->jShell->cartesian_l.size(); nn++, cartbf2++) {
                Slmu += Csp2cai(i,mm)*Csp2caj(j,nn)*SpinOrbit(cartbf1,cartbf2,mu);
              }
            }
            (*this->SOCoupling_)(bf1,bf2,mu) = Slmu;
            (*this->SOCoupling_)(bf2,bf1,mu) = -Slmu;
            }
          }
        }
      }     //end of cartesian to spherical transform
*/
    }
  
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
      RealMap TEMP(&this->SOCoupling_->storage()[iXYZ*this->nBasis_*this->nBasis_],this->nBasis_,this->nBasis_);
      double fobi;
      fobi = TEMP.frobInner(TEMP);
       this->fileio_->out <<"fobinius inner product XYZ=\t"<< fobi << endl;
      prettyPrint(this->fileio_->out,TEMP,"SOCoupling XYZ="+std::to_string(iXYZ));
    }
  }


void AOIntegrals::computepVdotp(){
  double pVp, pVpC, C[3];     //C is the coordinate of atoms
  int i,j,k,ijShell,lA[3],lB[3],iPP,iAtom,nu;
//  RealMatrix pVdotp(this->nCartBasis_,this->nCartBasis_);
  std::vector<double> tmppVp;
  int bf1,bf2;
  ChronusQ::ShellPair *ijS;
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);

    if (ijS->iShell==ijS->jShell) { 
      tmppVp.resize(ijS->icarsize*ijS->jcarsize);
      for(i=0,bf1=ijS->icarbf_s  ; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=i,bf2=ijS->icarbf_s+i; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };

        pVp = 0;
        for ( iPP=0; iPP<ijS->nPGTOPair; iPP++) {
          pVpC = 0;
          for (iAtom = 0; iAtom <this->molecularConstants_->nAtom; iAtom++){
             for ( nu = 0; nu <3 ; nu++) {
               C[nu] = this->molecularConstants_->cart[nu][iAtom];
             }
             pVpC += this->molecularConstants_->atomZ[iAtom]*this->pVpab(ijS,C,ijS->iShell->l,lA,ijS->jShell->l,lB,0,iPP,iAtom);

           }
         pVp += pVpC;
        }  
//        (*this->pVp_)(bf1,bf2) = pVp; 
//        (*this->pVp_)(bf2,bf1) = pVp;
//          pVdotp(bf1,bf2) = pVp;
//          pVdotp(bf2,bf1) = pVp;
          if (this->basisSet_->getforceCart() == false ) {
            tmppVp[i*ijS->jcarsize+j] = pVp;
            tmppVp[j*ijS->jcarsize+i] = pVp;
          }
          else if (this->basisSet_->getforceCart()== true) {
            (*this->pVp_)(bf1,bf2) = pVp;  
            (*this->pVp_)(bf2,bf1) = pVp; 
          }
        }
   } else {
      for(i=0,bf1=ijS->icarbf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jcarbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };

        pVp = 0;
        for ( iPP=0; iPP<ijS->nPGTOPair; iPP++) {
          pVpC = 0;
          for (iAtom = 0; iAtom <this->molecularConstants_->nAtom; iAtom++){
             for ( nu = 0; nu <3 ; nu++) {
               C[nu] = this->molecularConstants_->cart[nu][iAtom];
             }
             pVpC += this->molecularConstants_->atomZ[iAtom]*this->pVpab(ijS,C,ijS->iShell->l,lA,ijS->jShell->l,lB,0,iPP,iAtom);

           }
         pVp += pVpC;
        }  
//        (*this->pVp_)(bf1,bf2) = pVp; 
//        (*this->pVp_)(bf2,bf1) = pVp;
//          pVdotp(bf1,bf2) = pVp;
//          pVdotp(bf2,bf1) = pVp;
          if (this->basisSet_->getforceCart()==false ) {
            tmppVp.push_back(pVp);
          }
          else if (this->basisSet_->getforceCart()== true ) {
            (*this->pVp_)(bf1,bf2) = pVp;  
            (*this->pVp_)(bf2,bf1) = pVp; 
          }
      }
      }
    if (this->basisSet_->getforceCart()==false) {
      std::vector<double> SPH = this->cart2SphTrans(ijS,tmppVp.data() );
      for ( i = 0, bf1= ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++) {
        for ( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++, bf2++  ) {
          (*this->pVp_)(bf1,bf2) = SPH[i*ijS->jsphsize+j];
          (*this->pVp_)(bf2,bf1) = SPH[i*ijS->jsphsize+j];
        }
      }
      SPH.clear();
      tmppVp.clear(); 
    }

/*      if (this->nBasis_==this->nSphBasis_) {
        RealMatrix Csp2cai(ijS->isphsize ,ijS->iShell->cartesian_l.size());
        RealMatrix Csp2caj(ijS->jsphsize ,ijS->jShell->cartesian_l.size());

        for( p = 0 ; p<ijS->isphsize; p++ ) {
          for( q = 0 ; q<ijS->icarsize; q++ ) {
            Csp2cai(p,q) = 0;
          }
        }
        for ( p = 0 ; p <ijS->jsphsize; p++ ) {
          for ( q = 0 ; q <ijS->jcarsize; q++ ) {
            Csp2caj(p,q) = 0;
          }
        }
 
        if (ijS->iShell->l==0) {
          Csp2cai(0,0) = 1;
        }
        else if (ijS->iShell->l==1) {
          Csp2cai(0,0) =1;
          Csp2cai(1,1) =1;
          Csp2cai(2,2) =1;
        }
        else if (ijS->iShell->l==2) {
// this transformation is consistent with libint

          Csp2cai(4,0) = sqrt(3)/2;
          Csp2cai(4,3) = -sqrt(3)/2;
      
          Csp2cai(3,2) = sqrt(3);

          Csp2cai(2,0) = -0.5;
          Csp2cai(2,3) = -0.5;
          Csp2cai(2,5) = 1;

          Csp2cai(1,4) = sqrt(3);

          Csp2cai(0,1) = sqrt(3);

      // this transformation is consistent with Gaussian. 
 //for the untransformed cartesian integral, 
 //(1,1,0) (1,0,1) (0,1,1) times sqrt(3) will be equal to the result to gaussian.
 //cartesian gaussians are not normalized properly in the norm constant in libint

//gaussian          Csp2cai(3,0) = sqrt(3)/2;
          Csp2cai(3,3) = -sqrt(3)/2;
      
          Csp2cai(1,2) = sqrt(3);

          Csp2cai(0,0) = -0.5;
          Csp2cai(0,3) = -0.5;
          Csp2cai(0,5) = 1;

          Csp2cai(2,4) = sqrt(3);

          Csp2cai(4,1) = sqrt(3);
 //gaussian
        }
        if (ijS->jShell->l==0) {
          Csp2caj(0,0) = 1;
      }
        else if (ijS->jShell->l==1) {
          Csp2caj(0,0) =1;
          Csp2caj(1,1) =1;
          Csp2caj(2,2) =1;
        }
        else if (ijS->jShell->l==2) {
          Csp2caj(4,0) = sqrt(3)/2;
          Csp2caj(4,3) = -sqrt(3)/2;
       
          Csp2caj(3,2) = sqrt(3);

          Csp2caj(2,0) = -0.5;
          Csp2caj(2,3) = -0.5;
          Csp2caj(2,5) = 1;

          Csp2caj(1,4) = sqrt(3);

          Csp2caj(0,1) = sqrt(3);
//gaussian
          Csp2caj(3,0) = sqrt(3)/2;
          Csp2caj(3,3) = -sqrt(3)/2;
       
          Csp2caj(1,2) = sqrt(3);

          Csp2caj(0,0) = -0.5;
          Csp2caj(0,3) = -0.5;
          Csp2caj(0,5) = 1;

          Csp2caj(2,4) = sqrt(3);

          Csp2caj(4,1) = sqrt(3);
//gaussian
        }


        for( i = 0, bf1 = ijS->isphbf_s; i<ijS->isphsize ; i++,bf1++  ) {
          for( j = 0, bf2 = ijS->jsphbf_s; j<ijS->jsphsize ; j++,bf2++ ) {
            pVp = 0;
            for( p = 0, cartbf1 = ijS->icarbf_s; p<ijS->icarsize; p++, cartbf1++) {
              for( q = 0, cartbf2 = ijS->jcarbf_s; q<ijS->jcarsize; q++, cartbf2++) {
                pVp  += Csp2cai(i,p)*Csp2caj(j,q)*pVdotp(cartbf1,cartbf2);
              }
            }
            (*this->pVp_)(bf1,bf2) = pVp;
            (*this->pVp_)(bf2,bf1) = pVp;
          }
        }
      }    // end of cartesian to spherical transform
*/
    };
 
  prettyPrint(this->fileio_->out,(*this->pVp_),"pVp");
  
  double fobi;
  fobi = this->pVp_->frobInner(*this->pVp_);
  this->fileio_->out <<"fobinius inner product" << std::setprecision(12) << std::scientific << fobi << endl;

};





//SS

//------------------------------------//
// overlap horizontal recursion       //
// (a|b) = (A-B)(a|b-1) + (a+1|b-1)   //
//   LA >= LB                         //
//------------------------------------//
double AOIntegrals::hRRSab(ChronusQ::ShellPair *ijSP,int LA,int *lA,int LB,int *lB) {
  int i,j,iPP;
  double tmpVal = 0.0;
  if(LA==0) {
    // (s|s)
    for(iPP=0;iPP<ijSP->nPGTOPair;iPP++) 
      tmpVal+= ijSP->ss[iPP];
    return tmpVal;
  } else if(LB==0) {
    // (|s)
    for(iPP=0;iPP<ijSP->nPGTOPair;iPP++) 
      tmpVal+= ijSP->ss[iPP]*this->vRRSa0(ijSP,LA,lA,iPP);
    return tmpVal;
  };

  int iWork,lAp1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAp1[iWork]=lA[iWork];     
    lBm1[iWork]=lB[iWork];     
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lAp1[iWork]++;
  lBm1[iWork]--;

  tmpVal = this->hRRSab(ijSP,LA+1,lAp1,LB-1,lBm1);
  tmpVal+= ijSP->AB[iWork]*this->hRRSab(ijSP,LA,lA,LB-1,lBm1);
  return tmpVal;
};

//----------------------------------------------//
//kinetic integrals                             //
//
//[a|T|b]=[a|Tx|b]+[a|Ty|b]+[a|Tz|b]            //
//[a|Ti|b]=beta*(2bi+1)[a|b]-2beta^2 *[a|b+2i]-1/2*bi*(bi-1)[a|b-2i]//
// LA>=LB?                                      //
//S.S.                                          //
//----------------------------------------------//
double AOIntegrals::RRTab(ChronusQ::ShellPair *ijSP,int LA,int *lA,int LB,int *lB) {
  int i,j,iPP;
  double tmpVal = 0.0;
  if(LA==0) {
    // (s|T|s)
    for(iPP=0;iPP<ijSP->nPGTOPair;iPP++) tmpVal+= ijSP->ssT[iPP];
    return tmpVal;
  };

  int iWork,lBm[3];
 
  for(iWork=0;iWork<3;iWork++) { 
    for(i=0;i<3;i++) lBm[i]=lB[i];     
    lBm[iWork] += 2; 

    for(iPP=0;iPP<ijSP->nPGTOPair;iPP++) {
      tmpVal -= math.two*ijSP->beta[iPP]*ijSP->beta[iPP]*this->hRRiPPSab(ijSP,LA,lA,LB+2,lBm,iPP);
      tmpVal += ijSP->beta[iPP]*(2*lB[iWork]+1)*this->hRRiPPSab(ijSP,LA,lA,LB,lB,iPP);
    };
    if(lB[iWork]>=2) {
      lBm[iWork] -= 4; 
      for(iPP=0;iPP<ijSP->nPGTOPair;iPP++) {
        tmpVal -= math.half*lB[iWork]*(lB[iWork]-1)*this->hRRiPPSab(ijSP,LA,lA,LB-2,lBm,iPP);
      };
    };
  };
  return tmpVal;
};
//----------------------------------------------------------//
// overlap vertical recursion                               //
// (a|0) = (P-A)(a-1|0) + halfInvZeta*N_(a-1)*(a-2|0)//
//----------------------------------------------------------//
double AOIntegrals::vRRSa0(ShellPair *ijSP,int LA,int *lA,int iPP) {
  double tmpVal=0.0;
  int iWork;
  int lAm1[3];
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0) iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;

  if(LA==1) return ijSP->PA[iPP][iWork];

  lAm1[iWork]--;
  tmpVal = ijSP->PA[iPP][iWork]*this->vRRSa0(ijSP,LA-1,lAm1,iPP);

  if(LA==2&&lA[iWork]==2) tmpVal += ijSP->halfInvZeta[iPP];
  else if (lA[iWork]>=2) {
    lAm1[iWork]--;
    tmpVal += (lAm1[iWork]+1)*ijSP->halfInvZeta[iPP]*this->vRRSa0(ijSP,LA-2,lAm1,iPP);
  };

  return tmpVal;
};

//---------------------------------------------------------//
// kinetic vertical recursion                              //
// (a|T|b) = (P-B)(a|T|b-1) + halfInvZeta*N_(b-1)(a|T|b-2) //
//         + halfInvZeta*N_(a)(a-1|T|b-1)                  //
//         + 2*Xi*[(a|b) - halfInvZeta*N_(b-1)(a|b-2)]     //
//                                                         //
// (a|T|0) = (P-A)(a-1|T|0) + halfInvZeta*N_(a-1)(a-2|T|0) //
//         + 2*Xi*[(a|0) - halfInvZeta*N_(a-1)(a-2|0)]     //
//---------------------------------------------------------//
double AOIntegrals::vRRTab(ChronusQ::ShellPair *ijSP,
                           int LA,int *lA,int LB,int *lB) {
  double tmpVal = 0.0;
  int i,j,iPP;
  if(LA==0) {
    for(iPP=0; iPP<ijSP->nPGTOPair; iPP++) tmpVal += ijSP->ssT[iPP];
    return tmpVal;
  } else if(LB==0) {
//    cout<<"LA LB: "<<LA<<" "<<LB<<endl;
    for(iPP=0; iPP<ijSP->nPGTOPair; iPP++) tmpVal += this->vRRiPPTab(ijSP,LA,lA,LB,lB,iPP);
    return tmpVal;
  };

  int iWork,lAm1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAm1[iWork]=lA[iWork];
    lBm1[iWork]=lB[iWork];
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lBm1[iWork]--;
  for(iPP=0; iPP<ijSP->nPGTOPair; iPP++)
    tmpVal += ijSP->PB[iPP][iWork]*this->vRRiPPTab(ijSP,LA,lA,LB-1,lBm1,iPP);

  if(lA[iWork]>0) {
    lAm1[iWork]--;
    for(iPP=0; iPP<ijSP->nPGTOPair; iPP++)
      tmpVal += (lAm1[iWork]+1)*ijSP->halfInvZeta[iPP]*this->vRRiPPTab(ijSP,LA-1,lAm1,LB-1,lBm1,iPP);
  };

  if(lB[iWork]>=2) {
    lBm1[iWork]--;
    for(iPP=0; iPP<ijSP->nPGTOPair; iPP++) {
      tmpVal += (lBm1[iWork]+1)*ijSP->halfInvZeta[iPP]*this->vRRiPPTab(ijSP,LA,lA,LB-2,lBm1,iPP);
      tmpVal -= (lBm1[iWork]+1)*math.two*ijSP->Xi[iPP]*ijSP->halfInvZeta[iPP]*this->hRRiPPSab(ijSP,LA,lA,LB-2,lBm1,iPP);
    };
  };

  for(iPP=0; iPP<ijSP->nPGTOPair; iPP++)
    tmpVal += math.two*ijSP->Xi[iPP]*this->hRRiPPSab(ijSP,LA,lA,LB,lB,iPP);

  return tmpVal;
};
//---------------------------------------------------------//
// kinetic vertical recursion (iPP specific)               //
// (a|T|b) = (P-B)(a|T|b-1) + halfInvZeta*N_(b-1)(a|T|b-2) //
//         + halfInvZeta*N_(a)(a-1|T|b-1)                  //
//         + 2*Xi*[(a|b) - halfInvZeta*N_(b-1)(a|b-2)]     //
//                                                         //
// (a|T|0) = (P-A)(a-1|T|0) + halfInvZeta*N_(a-1)(a-2|T|0) //
//         + 2*Xi*[(a|0) - halfInvZeta*N_(a-1)(a-2|0)]     //
//---------------------------------------------------------//
double AOIntegrals::vRRiPPTab(ChronusQ::ShellPair *ijSP,
                           int LA,int *lA,int LB,int *lB,int iPP) {
  double tmpVal = 0.0;
  if(LA==0) {
    tmpVal += ijSP->ssT[iPP];
    return tmpVal;
  } else if(LB==0) {
//    cout<<"LA: "<<LA<<endl;
    int iWork,lAm1[3];
    for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
    if (lA[0]>0)      iWork=0;
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
    lAm1[iWork]--;
//    cout<<"xsli test 1 :"<<iWork<<" "<<ijSP->PA[iPP][iWork]<<endl;
    tmpVal += ijSP->PA[iPP][iWork]*this->vRRiPPTab(ijSP,LA-1,lAm1,LB,lB,iPP);
//    cout<<"xsli test 2 :"<<tmpVal<<endl;
    if(lA[iWork]==2) {
      tmpVal += ijSP->halfInvZeta[iPP]*ijSP->ssT[iPP];
//    cout<<"xsli test 3 :"<<tmpVal<<endl;
      tmpVal -= math.two*ijSP->Xi[iPP]*ijSP->halfInvZeta[iPP]*ijSP->ss[iPP];
//    cout<<"xsli test 4 :"<<tmpVal<<endl;
    } else if(lA[iWork]>2) {
      lAm1[iWork]--;
      tmpVal += (lAm1[iWork]+1)*ijSP->halfInvZeta[iPP]*this->vRRiPPTab(ijSP,LA-2,lAm1,LB,lB,iPP);
      tmpVal -= (lAm1[iWork]+1)*math.two*ijSP->Xi[iPP]*ijSP->halfInvZeta[iPP]*ijSP->ss[iPP]*this->vRRSa0(ijSP,LA-2,lAm1,iPP);
    };
    tmpVal += math.two*ijSP->Xi[iPP]*ijSP->ss[iPP]*this->vRRSa0(ijSP,LA,lA,iPP);
//    cout<<"xsli test 5 :"<<tmpVal<<endl;
    return tmpVal;
  };

  int iWork,lAm1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAm1[iWork]=lA[iWork];
    lBm1[iWork]=lB[iWork];
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lBm1[iWork]--;

  tmpVal += ijSP->PB[iPP][iWork]*this->vRRiPPTab(ijSP,LA,lA,LB-1,lBm1,iPP);

  if(lA[iWork]>0) {
    lAm1[iWork]--;
    tmpVal += (lAm1[iWork]+1)*ijSP->halfInvZeta[iPP]*this->vRRiPPTab(ijSP,LA-1,lAm1,LB-1,lBm1,iPP);
  };

  if(lB[iWork]>=2) {
    lBm1[iWork]--;
    tmpVal += (lBm1[iWork]+1)*ijSP->halfInvZeta[iPP]*this->vRRiPPTab(ijSP,LA,lA,LB-2,lBm1,iPP);
    tmpVal -= (lBm1[iWork]+1)*math.two*ijSP->Xi[iPP]*ijSP->halfInvZeta[iPP]*this->hRRiPPSab(ijSP,LA,lA,LB-2,lBm1,iPP);
  };

  tmpVal += math.two*ijSP->Xi[iPP]*this->hRRiPPSab(ijSP,LA,lA,LB,lB,iPP);

  return tmpVal;
};
//----------------------------------//
// overlap recursion (iPP specific) //
// (a|b) = (A-B)(a|b-1) + (a+1|b-1) //
//   LA >= LB                       //
//----------------------------------//
double AOIntegrals::hRRiPPSab(ChronusQ::ShellPair *ijSP,int LA,int *lA,int LB,int *lB,int iPP) {
  double tmpVal = 0.0;
  if((LA+LB)==0) {
    tmpVal = ijSP->ss[iPP];
    return tmpVal;
  } else if(LB==0) {
    tmpVal = ijSP->ss[iPP]*this->vRRSa0(ijSP,LA,lA,iPP);
    return tmpVal;
  };

  int iWork,lAp1[3],lBm1[3];
  for(iWork=0;iWork<3;iWork++){
    lAp1[iWork]=lA[iWork];     
    lBm1[iWork]=lB[iWork];     
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lAp1[iWork]++;
  lBm1[iWork]--;

  tmpVal = this->hRRiPPSab(ijSP,LA+1,lAp1,LB-1,lBm1,iPP);
  tmpVal+= ijSP->AB[iWork]*this->hRRiPPSab(ijSP,LA,lA,LB-1,lBm1,iPP);
  return tmpVal;
};
//---------------------------------------------------------//
// potential integral horizontal recursion                 //
//  (a|0_c|b) = (a+1|0_c|b-1) + (A-B)(a|0_c|b-1)           //
//   LA >= LB                                              //
//   horizontal recursion doesn't increase (m). and it's only used once, so m=0 here. //
//---------------------------------------------------------//
double AOIntegrals::hRRVab(ChronusQ::ShellPair *ijSP,MolecularConstants *mc,int LA,int *lA,int LB,int *lB){
  int iPP,iAtom,m;
  double tmpVal=0.0;
  double PC[3];
  double rho;

  if(LB==0) {
    // (LA|s)
    for(iPP=0;iPP<ijSP->nPGTOPair;iPP++) 
    for(iAtom=0;iAtom<mc->nAtom;iAtom++) {
//cout<<"number of atom in potential integral\t"<<mc->nAtom<<endl;
      double squarePC=0.0;
      for(m=0;m<3;m++) {
//        PC[m] = ijSP->P[iPP][m] - mc->cart[iAtom][m];
        PC[m] = ijSP->P[iPP][m] - mc->cart[m][iAtom];
        squarePC += PC[m]*PC[m];
      };
//      cout<<LA<<" "<<LB<<" "<<squarePC<<endl;

      double *tmpFmT = new double[ijSP->lTotal+1];
      if (this->useFiniteWidthNuclei == true) {
        rho = ijSP->Zeta[iPP]*this->molecule_->nucShell(iAtom).alpha[0]/(ijSP->Zeta[iPP]+this->molecule_->nucShell(iAtom).alpha[0]);
        this->computeFmTTaylor(tmpFmT,rho*squarePC,ijSP->lTotal,0);
      }
      else if (this->useFiniteWidthNuclei == false ) {
        this->computeFmTTaylor(tmpFmT,ijSP->Zeta[iPP]*squarePC,ijSP->lTotal,0);
//cout<<"not dead yet"<<endl;
      }
      if (LA==0) {
        if ( this->useFiniteWidthNuclei == false ) {
          tmpVal += (static_cast<double>(mc->atomZ[iAtom]))*ijSP->ssV[iPP]*tmpFmT[0];
//cout<<"still good"<<endl;
        }
        else if ( this->useFiniteWidthNuclei == true ) {
          tmpVal += (static_cast<double>(mc->atomZ[iAtom]))*math.two*sqrt(rho/math.pi)*ijSP->ss[iPP]*tmpFmT[0];
        }
//        cout<<this->molecule_->nucShell(iAtom).alpha[0]<<"\t"<<this->molecule_->nucShell(iAtom).contr[0].coeff[0]<<endl;
      }
      else {
        if ( this->useFiniteWidthNuclei == false ) {
//cout<<"before LA= "<<LA<<endl;
//cout<<lA[0]<<";"<<lA[1]<<";"<<lA[2]<<endl;
          tmpVal += mc->atomZ[iAtom]*ijSP->ssV[iPP]*this->vRRVa0(ijSP,tmpFmT,PC,0,LA,lA,iPP,iAtom);
//cout<<"thanks"<<endl;
        } else if ( this->useFiniteWidthNuclei == true ) {
          tmpVal += mc->atomZ[iAtom]*math.two*sqrt(rho/math.pi)*ijSP->ss[iPP]*this->vRRVa0(ijSP,tmpFmT,PC,0,LA,lA,iPP,iAtom);
        }
      }
      delete[] tmpFmT;
    };
  } else {
    int iWork,lAp1[3],lBm1[3];
    for(iWork=0;iWork<3;iWork++) {
      lAp1[iWork]=lA[iWork];
      lBm1[iWork]=lB[iWork];
    };
    if (lB[0]>0) iWork=0;
    else if (lB[1]>0) iWork=1;
    else if (lB[2]>0) iWork=2;

    lAp1[iWork]++;
    lBm1[iWork]--;
    tmpVal = this->hRRVab(ijSP,mc,LA+1,lAp1,LB-1,lBm1);
    tmpVal+= ijSP->AB[iWork]*this->hRRVab(ijSP,mc,LA,lA,LB-1,lBm1);
  };
  return tmpVal;
};
//----------------------------------------------------------------------------//
// potential integral vertical recursion                                      //
//  (a|0_c|0)^(m) = (P-A)(a-1|0_c|0)^(m) - (P-C)(a-1|0_c|0)^(m+1)             //
//                + halfInvZeta*N_(a-1)*[(a-2|0_c|0)^(m)-(a-2|0_c|0)^(m+1)]   //
//  since LB==0, we only decrease a                                           //
//----------------------------------------------------------------------------//
double AOIntegrals::vRRVa0(ChronusQ::ShellPair *ijSP,double *FmT,double *PC,int m,int LA,int *lA,int iPP,int iAtom){

  double rhoovzeta;
  if(LA==0) {
    if (this->useFiniteWidthNuclei == true) {
//cout<<"out here"<<endl;
      rhoovzeta =  this->molecule_->nucShell(iAtom).alpha[0]/(this->molecule_->nucShell(iAtom).alpha[0]+ijSP->Zeta[iPP]);
      return (pow(rhoovzeta,m)*FmT[m]);
    } else if ( this->useFiniteWidthNuclei == false ) {       
//cout<<"no finite"<<endl;
      return FmT[m];  //Z*2*sqrt[zeta/pi]*[s|s]is given in hRRVab, in ssV.
    }
  }

  double tmpVal=0.0;
  int lAm1[3]; //means lA minus 1_i
  int iWork;
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0) iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;
  lAm1[iWork]--;

//cout<<"already in the recursion"<<endl;
if (LA>=-1) {
//  cout<<"LA"<<LA<<endl;
//  cout<<lAm1[0]<<";"<<lAm1[1]<<";"<<lAm1[2]<<"\t"<<lA[0]<<";"<<lA[1]<<";"<<lA[2]<<endl;
}
  tmpVal  = ijSP->PA[iPP][iWork]*this->vRRVa0(ijSP,FmT,PC,m,LA-1,lAm1,iPP,iAtom);
//cout<<"2"<<endl;
  tmpVal -= PC[iWork]*this->vRRVa0(ijSP,FmT,PC,m+1,LA-1,lAm1,iPP,iAtom);
//cout<<"1"<<endl;
  if(lAm1[iWork]>=1){
    lAm1[iWork]--;
    tmpVal += (lAm1[iWork]+1)*ijSP->halfInvZeta[iPP]*(this->vRRVa0(ijSP,FmT,PC,m,LA-2,lAm1,iPP,iAtom)
            - this->vRRVa0(ijSP,FmT,PC,m+1,LA-2,lAm1,iPP,iAtom));
  };
//cout<<"finish the recursion and want to quit here oh yeah"<<endl;
  return tmpVal; 
};

// SS
//----------------------------------------------------------------------------//
// potential integral horizontal recursion iPP specific                       //
//  [a|A(0)|b]^(m) = [a+1i|A(0)|b-1i]^(m) + (Ai-Bi)*[a|A(0)|b-1i]^(m)         //
//  LA >=LB                                                                   //
//  m can be nonzero                                                          //
//----------------------------------------------------------------------------//

double AOIntegrals::hRRiPPVab(ChronusQ::ShellPair *ijSP,int LA,int *lA,int LB,int *lB,double *C,int m,int iPP,int iAtom){
  double tmpVal=0.0;
  double PC[3];
  double rho;
  int k;
  if(LB==0) {
    // (LA|s)
      double squarePC=0.0;
      for(k=0;k<3;k++) {
        PC[k] = ijSP->P[iPP][k] - C[k];
        squarePC += PC[k]*PC[k];
      };
      double *tmpFmT = new double[ijSP->lTotal+m+1];
      if (this->useFiniteWidthNuclei == true) {
        rho = ijSP->Zeta[iPP]*this->molecule_->nucShell(iAtom).alpha[0]/(ijSP->Zeta[iPP]+this->molecule_->nucShell(iAtom).alpha[0]);
        this->computeFmTTaylor(tmpFmT,rho*squarePC,ijSP->lTotal+m,0);
      }
      else if (this->useFiniteWidthNuclei == false ) {
        this->computeFmTTaylor(tmpFmT,ijSP->Zeta[iPP]*squarePC,ijSP->lTotal+m,0);
      }
      if(LA==0) {
        if ( this->useFiniteWidthNuclei == false ) {
          tmpVal = ijSP->ssV[iPP]*tmpFmT[m];
        }
        else if ( this->useFiniteWidthNuclei == true ) {
          
          double  rhoovzeta =  this->molecule_->nucShell(iAtom).alpha[0]/(this->molecule_->nucShell(iAtom).alpha[0]+ijSP->Zeta[iPP]);
          tmpVal = math.two*sqrt(rho/math.pi)*pow(rhoovzeta,m)*ijSP->ss[iPP]*tmpFmT[m];
        }
      }

      else if (LA>0) {
        if ( this->useFiniteWidthNuclei == false ) {
          tmpVal = ijSP->ssV[iPP]*this->vRRVa0(ijSP,tmpFmT,PC,m,LA,lA,iPP,iAtom);
        }
        else if ( this->useFiniteWidthNuclei == true ) {
          tmpVal = math.two*sqrt(rho/math.pi)*ijSP->ss[iPP]*this->vRRVa0(ijSP,tmpFmT,PC,m,LA,lA,iPP,iAtom);
        } 
      }
//cout<<"before delete the tmpFmT"<<endl;
      delete[] tmpFmT;
//cout<<"after delete the tmpFmT"<<endl;
    }
  else if ((LA==0)&&(LB>0)){
//    cout<<"LA==0\tLB>0"<<endl;
      double squarePC=0.0;
      for(k=0;k<3;k++) {
        PC[k] = ijSP->P[iPP][k] - C[k];
        squarePC += PC[k]*PC[k];
      };
      double *tmpFmT = new double[ijSP->lTotal+m+1];
      if (this->useFiniteWidthNuclei == true) {
        rho = ijSP->Zeta[iPP]*this->molecule_->nucShell(iAtom).alpha[0]/(ijSP->Zeta[iPP]+this->molecule_->nucShell(iAtom).alpha[0]);
        this->computeFmTTaylor(tmpFmT,rho*squarePC,ijSP->lTotal+m,0);
        tmpVal =  math.two*sqrt(rho/math.pi)*ijSP->ss[iPP]*this->vRRV0b(ijSP,tmpFmT,PC,m,LB,lB,iPP,iAtom);
      }
      else if (this->useFiniteWidthNuclei == false ) {
        this->computeFmTTaylor(tmpFmT,ijSP->Zeta[iPP]*squarePC,ijSP->lTotal+m,0);

        tmpVal = ijSP->ssV[iPP]*this->vRRV0b(ijSP,tmpFmT,PC,m,LB,lB,iPP,iAtom);
      }

//cout<<"before delete the tmpFmT"<<endl;
      delete[] tmpFmT;
//cout<<"after delete the tmpFmT"<<endl;
    }

  
  else if ((LA>0)&&(LB>0)) {
    int iWork,lAp1[3],lBm1[3];
    for(iWork=0;iWork<3;iWork++) {
      lAp1[iWork]=lA[iWork];
      lBm1[iWork]=lB[iWork];
    };
    if (lB[0]>0) iWork=0;
    else if (lB[1]>0) iWork=1;
    else if (lB[2]>0) iWork=2;

    lAp1[iWork]++;
    lBm1[iWork]--;
    tmpVal = this->hRRiPPVab(ijSP,LA+1,lAp1,LB-1,lBm1,C,m,iPP,iAtom);
    tmpVal+= ijSP->AB[iWork]*this->hRRiPPVab(ijSP,LA,lA,LB-1,lBm1,C,m,iPP,iAtom);
  };
  return tmpVal;
}

//----------------------------------------------------------------------------//
// potential integral vertical recursion                                      //
//  (0|0_c|b)^(m) = (P-B)(0|0_c|b-1i)^(m) - (P-C)(0|0_c|b-1i)^(m+1)           //
//                + halfInvZeta*N_(b-1)*[(0|0_c|b-2i)^(m)-(0|0_c|b-2i)^(m+1)] //
//  since LA==0, we only decrease b                                           //
//----------------------------------------------------------------------------//
double AOIntegrals::vRRV0b(ChronusQ::ShellPair *ijSP,double *FmT,double *PC,int m,int LB,int *lB,int iPP,int iAtom){

  double rhoovzeta;
  if(LB==0) {
     if (this->useFiniteWidthNuclei == true) {
      rhoovzeta =  this->molecule_->nucShell(iAtom).alpha[0]/(this->molecule_->nucShell(iAtom).alpha[0]+ijSP->Zeta[iPP]);
//cout<<"power"<<pow(rhoovzeta,m)<<endl; 
     return (pow(rhoovzeta,m)*FmT[m]);
    } else if ( this->useFiniteWidthNuclei == false ) {       
    
      return FmT[m];  //Z*2*sqrt[zeta/pi]*[s|s]is given in hRRVab, in ssV.
    }
  }

  double tmpVal=0.0;
  int lBm1[3]; //means lA minus 1_i
  int iWork;
  for(iWork=0;iWork<3;iWork++) lBm1[iWork]=lB[iWork];
  if (lB[0]>0) iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lBm1[iWork]--;

  tmpVal  = ijSP->PB[iPP][iWork]*this->vRRV0b(ijSP,FmT,PC,m,LB-1,lBm1,iPP,iAtom);
  tmpVal -= PC[iWork]*this->vRRV0b(ijSP,FmT,PC,m+1,LB-1,lBm1,iPP,iAtom);

  if(lBm1[iWork]>=1){
    lBm1[iWork]--;
    tmpVal += (lBm1[iWork]+1)*ijSP->halfInvZeta[iPP]*(this->vRRV0b(ijSP,FmT,PC,m,LB-2,lBm1,iPP,iAtom) - this->vRRV0b(ijSP,FmT,PC,m+1,LB-2,lBm1,iPP,iAtom));
  };
  return tmpVal; 
};





//SS


//---------------------------------------------------------------------------------//
//angular momentum vertical recursion                                              //
// if LB==0,then [a|L|0]=(Pi-Ai)[a-1i|L|0]]+halfInvZeta*Ni(a-1i)[a-2i|L|0]         //
//                  +zeta_b/zeta*{1i cross(B-C)}_mu*[a-1i|b]                       //
// if LB>0,then  [a|L|b]=(Pi-Bi)[a|L|b-1i]                                         //
//                  +halfInvZeta*Ni(b-1i)[a|L|b-2i]+halfInvZeta*Ni(a)[a-1i|L|b-1i] //
//                  -zeta_a/zeta*{1i cross(A-C)}_mu*[a|b-1i]                       //
//                  -halfInvZeta*Sum_k=x,y,z N_k(a){1i cross 1k}_mu*[a-1k|b-1i]    //
//---------------------------------------------------------------------------------//
double AOIntegrals::Labmu(ChronusQ::ShellPair *ijSP,RealMatrix *OneixAC,RealMatrix *OneixBC,RealTensor3d *OneiOnek,int LA,int *lA,int LB,int *lB,int mu,int iPP){
  double tmpVal=0.0;
  if ((LA+LB)==0){
    tmpVal = ijSP->ssL[iPP][mu];
    return tmpVal;
  };

  int lAm1[3];
  int iWork;
 
  if ((LB==0)&(LA>0)){
   for (iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
   if (lA[0]>0) iWork=0;
   else if (lA[1]>0) iWork=1;
   else if (lA[2]>0) iWork=2;
   lAm1[iWork]--;

   tmpVal  = ijSP->PA[iPP][iWork]*this->Labmu(ijSP,OneixAC,OneixBC,OneiOnek,LA-1,lAm1,LB,lB,mu,iPP);
   tmpVal += ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(*OneixBC)(iWork,mu)*this->hRRiPPSab(ijSP,LA-1,lAm1,LB,lB,iPP); 
//   cout<<(*OneixBC)(iWork,mu)<<endl;//*this->hRRiPPSab(ijSP,LA-1,lAm1,LB,lB,iPP)<<endl; //test 


   if(lAm1[iWork]>=1){
     lAm1[iWork]--;
     tmpVal += (lA[iWork]-1)*ijSP->halfInvZeta[iPP]*(this->Labmu(ijSP,OneixAC,OneixBC,OneiOnek,LA-2,lAm1,LB,lB,mu,iPP));
   };
   return tmpVal;

  }else if ((LA>0)&(LB==0)) cout<<"here"; 
   else if ((LB>0)&(LA>0)){
//   cout<<"here we are"<<endl;
   int lBm1[3],lAm1k[3];
   for(iWork=0;iWork<3;iWork++) {
     lAm1[iWork]=lA[iWork];
     lBm1[iWork]=lB[iWork];
     lAm1k[iWork]=lA[iWork];
   };
   if (lB[0]>0) iWork=0;
   else if (lB[1]>0) iWork=1;
   else if (lB[2]>0) iWork=2;

   lAm1[iWork]--;
   lBm1[iWork]--;
   tmpVal  = ijSP->PB[iPP][iWork]*this->Labmu(ijSP,OneixAC,OneixBC,OneiOnek,LA,lA,LB-1,lBm1,mu,iPP);
   if (lAm1[iWork]>=0){
    tmpVal += ijSP->halfInvZeta[iPP]*lA[iWork]*this->Labmu(ijSP,OneixAC,OneixBC,OneiOnek,LA-1,lAm1,LB-1,lBm1,mu,iPP);
   };
//   double u;
   tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*(*OneixAC)(iWork,mu)*this->hRRiPPSab(ijSP,LA,lA,LB-1,lBm1,iPP);
//   cout<<u<<endl;
//   tmpVal -= u;
   int k;
   for(k=0;k<3;k++) {
     if(lA[k]>0){
      lAm1k[k]--;
      tmpVal -= ijSP->halfInvZeta[iPP]*lA[k]*(*OneiOnek)(iWork,k,mu)*this->hRRiPPSab(ijSP,LA-1,lAm1k,LB-1,lBm1,iPP);
      lAm1k[k]++;
     };
   };
   
   if (lBm1[iWork]>0){
    lBm1[iWork]--;
    tmpVal += (lBm1[iWork]+1)*ijSP->halfInvZeta[iPP]*Labmu(ijSP,OneixAC,OneixBC,OneiOnek,LA,lA,LB-2,lBm1,mu,iPP);
   };
  };
   return tmpVal;
};
     
        
 
//---------------------------------------------------------------------------------//
//spin orbital vertical recursion                                                  //
// if LA==0,LB==0, then [0|S|0]^(m)=4*zetaa*zetab*((A-C)X(B-C))_mu*[s|A(0)|s]^(m+1)//
// if LB==0,then [a|S|0]^(m)=(Pi-Ai)[a-1i|S|0]^(m)-(Pi-Ci)[a-1i|S|0]^(m+1)         //
//                           +halfInvZeta*Ni(a-1i)([a-2i|S|0]^(m)-[a-2i|S|0]^(m+1))//
//                           +2*zeta_b*{1i cross(B-C)}_mu*[a-1i|A(0)|0]^(m+1)      //
// if LB>0,then  [a|S|b]^(m)=(Pi-Bi)[a|S|b-1i]^(m)-(Pi-Ci)[a|S|b-1i]^(m+1)         //
//                  +halfInvZeta*Ni(b-1i)([a|S|b-2i]^(m)-[a|S|b-2i]^(m+1))         //
//                  +halfInvZeta*Ni(a)([a-1i|S|b-1i]^(m)-[a-1i|S|b-1i]^(m+1))      //
//                  -2*zeta_a*{1i cross(A-C)}_mu*[a|A(0)|b-1i]^(m+1)               //
//                  -Sum_k=x,y,z N_k(a)*{1i cross 1k}_mu*[a-1k|A(0)|b-1i]^(m+1)    //
//---------------------------------------------------------------------------------//
double AOIntegrals::Slabmu(ChronusQ::ShellPair *ijSP,RealMatrix *OneixAC,RealMatrix *OneixBC,RealTensor3d *OneiOnek,double *C,int LA,int *lA,int LB,int *lB,int mu,int m,int iPP,int iAtom){
  double tmpVal=0.0;
  double AC[3],BC[3],ACxBC[3];
  int k;

  if ((LA+LB)==0){ 
    for (k = 0; k<3; k++) {
      AC[k] = ijSP->A[k] - C[k];
      BC[k] = ijSP->B[k] - C[k];
    }
  
    ACxBC[0] = AC[1]*BC[2] - AC[2]*BC[1]; 
    ACxBC[1] = AC[2]*BC[0] - AC[0]*BC[2];
    ACxBC[2] = AC[0]*BC[1] - AC[1]*BC[0];
double u;
u = hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+1,iPP,iAtom);
tmpVal = 4*ijSP->zetaa[iPP]*ijSP->zetab[iPP]*ACxBC[mu]*u;        
if (tmpVal>0.001||tmpVal<-0.001) {
 } 
    return tmpVal;
  };

  int lAm1[3];
  int iWork;
 
  double PC[3];
  for (k = 0; k <3; k++) {
    PC[k] = ijSP->P[iPP][k] - C[k];
  } 

  if ((LB==0)&(LA>0)){
   for (iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
   if (lA[0]>0) iWork=0;
   else if (lA[1]>0) iWork=1;
   else if (lA[2]>0) iWork=2;
   lAm1[iWork]--;
   tmpVal  = ijSP->PA[iPP][iWork]*this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-1,lAm1,LB,lB,mu,m,iPP,iAtom);
   tmpVal -= PC[iWork]*this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-1,lAm1,LB,lB,mu,m+1,iPP,iAtom);
   tmpVal += 2*ijSP->zetab[iPP]*(*OneixBC)(iWork,mu)*this->hRRiPPVab(ijSP,LA-1,lAm1,LB,lB,C,m+1,iPP,iAtom); 
   if(lAm1[iWork]>=1){
     lAm1[iWork]--;
     tmpVal += (lA[iWork]-1)*ijSP->halfInvZeta[iPP]*(this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-2,lAm1,LB,lB,mu,m,iPP,iAtom) - this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-2,lAm1,LB,lB,mu,m+1,iPP,iAtom));
   };

   return tmpVal;

  }else if ((LA>0)&(LB==0)) cout<<"here"; 
   else if ((LB>0)&(LA>0)){

//cout<<endl<<"Here we have LA,LB>0"<<endl;

//cout<<lB[0]<<"\t"<<lB[1]<<"\t"<<lB[2]<<endl;
   int lBm1[3],lAm1k[3];
   for(iWork=0;iWork<3;iWork++) {
     lAm1[iWork]=lA[iWork];
     lBm1[iWork]=lB[iWork];
     lAm1k[iWork]=lA[iWork];
   };
   if (lB[0]>0) iWork=0;
   else if (lB[1]>0) iWork=1;
   else if (lB[2]>0) iWork=2;

//cout<<"iWork"<<iWork<<endl;
   lAm1[iWork]--;
   lBm1[iWork]--;
//cout<<"Still working"<<endl;

//cout<<"iWork"<<iWork<<endl;
   tmpVal  = ijSP->PB[iPP][iWork]*this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA,lA,LB-1,lBm1,mu,m,iPP,iAtom);

//cout<<"still working 2"<<endl;
//cout<<"iWork"<<iWork<<endl;
//cout<<"LB-1="<<LB-1<<endl;
//cout<<lBm1[0]<<"\t"<<lBm1[1]<<"\t"<<lBm1[2]<<endl;
   tmpVal -= PC[iWork]*this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA,lA,LB-1,lBm1,mu,m+1,iPP,iAtom);
   if (lAm1[iWork]>=0){
    tmpVal += ijSP->halfInvZeta[iPP]*lA[iWork]*(this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-1,lAm1,LB-1,lBm1,mu,m,iPP,iAtom) - this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-1,lAm1,LB-1,lBm1,mu,m+1,iPP,iAtom));
   };

//cout<<"iWork"<<iWork<<endl;
//cout<<"still working 3"<<endl;
//cout<<lA[0]<<"\t"<<lA[1]<<"\t"<<lA[2]<<endl;
//cout<<lBm1[0]<<"\t"<<lBm1[1]<<"\t"<<lBm1[2]<<endl;
//cout<<"m+1="<<m<<endl;
//cout<<"iWork="<<iWork<<"\t"<<"mu="<<mu<<endl;
//cout<<"OneixAC"<<(*OneixAC)(iWork,mu)<<endl;
   tmpVal -= 2*ijSP->zetaa[iPP]*(*OneixAC)(iWork,mu)*this->hRRiPPVab(ijSP,LA,lA,LB-1,lBm1,C,m+1,iPP,iAtom);
//cout<<"still working 4"<<endl;
   for(k=0;k<3;k++) {
     if(lA[k]>0){
      lAm1k[k]--;
      tmpVal -= lA[k]*(*OneiOnek)(iWork,k,mu)*this->hRRiPPVab(ijSP,LA-1,lAm1k,LB-1,lBm1,C,m+1,iPP,iAtom);
      lAm1k[k]++;
     };
   };
   
   if (lBm1[iWork]>0){
    lBm1[iWork]--;
    tmpVal += (lBm1[iWork]+1)*ijSP->halfInvZeta[iPP]*(this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA,lA,LB-2,lBm1,mu,m,iPP,iAtom)-this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA,lA,LB-2,lBm1,mu,m+1,iPP,iAtom));
   };
  };

//cout<<"abigissuethereisnothingyoucanafford/t/t/tomg"<<endl;
   return tmpVal;
};
   
  
//----------------------------------------------------------------------------------//
//pVp vertical recursion                                                            //
//                                                                                  //
//if LA+LB==0,then                                                                  //
//[s|pVp|s]=(6*zeta_a*zeta_b/zeta+4*zeta_a*zeta_b*(P-A)dot(P-B)[0_A|A(0)|0_B]^(m)   //
//     -[4*zeta_a*zeta_b*(2P-A-B)dot(P-C)+6*zeta_a*zeta_b/zeta]*[0_A|A(0)|0_B]^(m+1)//
//     +4*zeta_a*zeta_b*(P-C)^2*[0_A|A(0)|0_B]^(m+2)                                //
//                                                                                  //
//if LB==0, LA<>0, then                                                             //
//[a|pVp|0]^(m)=(P-A)[a-1i|pVp|0]^(m)-(P-C)[a-1i|pVp|0]^(m+1)                       //
//     +invZeta*Ni(a)*([a-2i|pVp|0]^(m)-[a-2i|pVp|0]^(m+1))                         //
//     -zeta_b/zeta*2*zeta_b*[a|A|1i]^(m)                                           //
//     +zeta_b/zeta*{2*zeta_a*[a+1i|A|0]^(m)-Ni(a)[a-1i|A|0]^(m)}                   //
//     -zeta_a/zeta*2*zeta_b*[a|A|1i]^(m+1)                                         //
//     -zeta_b/zeta*{2*zeta_a*[a+1i|A|0]^(m+1)-Ni(a)*[a-1i|A|0]^(m+1)}              //
//                                                                                  //    
//if LA==0, LB<>0, then                                                             //
//[0|pVp|b]^(m)=(P-B)[0|pVp|b-1i]^(m)-(P-C)[0|pVp|b-1i]^(m+1)                       //
//     +halfinvzeta*Ni(b)*{[0|pVp|b-2i]^(m)-[0|pVp|b-2i]^(m+1)}                     //
//     -zeta_a/zeta*2*zeta_a*[1i|A|b-1i]^(m)                                        //
//     +zeta_a/zeta*{2*zeta_b*[0|A|b]^(m)-Ni(b)[0|A|b-2i]^(m)                       //
//     -zeta_b/zeta*2*zeta_a*[1i|A|b-1i]^(m+1)                                      //
//     -zeta_a/zeta*{2*zeta_b*[a|A|b]^(m+1)-Ni(b)[a|A|b-2i]^(m+1)}                  //
//                                                                                  //
//if LA>0, LB>0, then                                                               //
//[a|pVp|b]^(m)=(P-B)[a|pVp|b-1i]^(m)-(P-C)[a|pVp|b-1i]^(m+1)                       //
//     +halfinvzeta*Ni(b)*{[a|pVp|b-2i]^(m)-[a|pVp|b-2i]^(m+1)}                     //
//     +halfinvzeta*Ni(a)*{[a-1i|pVp|b-1i]^(m)-[a-1i|pVp|b-1i]^(m+1)}               //
//     -zeta_a/zeta*{2*zeta_a*[a+1i|A|b-1i]^(m)-Ni(a)*[a-1i|A|b-1i]^(m)}            //
//     +zeta_a/zeta*{2*zeta_b*[a|A|b]^(m)-Ni(b)[a|A|b-2i]^(m)                       //
//     -zeta_b/zeta*{2*zeta_a*[a+1i|A|b-1i]^(m+1)-Ni(a)[a-1i|A|b-1i]^(m+1)}         //
//     -zeta_a/zeta*{2*zeta_b*[a|A|b]^(m+1)-Ni(b)[a|A|b-2i]^(m+1)}                  //
//----------------------------------------------------------------------------------//
double AOIntegrals::pVpab(ChronusQ::ShellPair *ijSP,double *C,int LA,int *lA,int LB,int *lB,int m,int iPP,int iAtom){
  double tmpVal=0.0;
  double PAPB,PABPC,PCsquare;
  int k;

  if ((LA+LB)==0){
    PAPB = 0.0;
    PABPC = 0.0;
    PCsquare = 0.0;
    for (k = 0; k<3 ; k++) {
      PAPB += (ijSP->P[iPP][k]-ijSP->A[k])*(ijSP->P[iPP][k]-ijSP->B[k]);
      PABPC += (2*ijSP->P[iPP][k]-ijSP->A[k]-ijSP->B[k])*(ijSP->P[iPP][k]-C[k]);
      PCsquare += (ijSP->P[iPP][k]-C[k])*(ijSP->P[iPP][k]-C[k]);
    }
  double u,uone,utwo;
  u = hRRiPPVab(ijSP,LA,lA,LB,lB,C,m,iPP,iAtom);
  uone = hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+1,iPP,iAtom);
  utwo = hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+2,iPP,iAtom);
  tmpVal = (6*ijSP->zetaa[iPP]*ijSP->zetab[iPP]*ijSP->invZeta[iPP]+4*ijSP->zetaa[iPP]*ijSP->zetab[iPP]*PAPB)*u;
  tmpVal -= (PABPC*4*ijSP->zetaa[iPP]*ijSP->zetab[iPP]+6*ijSP->zetaa[iPP]*ijSP->zetab[iPP]*ijSP->invZeta[iPP])*uone;
  tmpVal += 4*ijSP->zetaa[iPP]*ijSP->zetab[iPP]*PCsquare*utwo;

  return tmpVal;
  };
  
  int lAm1[3],lAp1[3],onei[3];
  int iWork;

  double PC[3];
  for (k = 0; k<3 ; k++) {
    PC[k] = ijSP->P[iPP][k] -C[k];
  }

  if ((LB==0)&(LA>0)){
    for (iWork = 0 ; iWork<3 ; iWork++ ) {
      lAm1[iWork] = lA[iWork];
      onei[iWork] = 0;
    }

    if (lA[0]>0) iWork = 0;
    else if (lA[1]>0)  iWork = 1;
    else if (lA[2]>0)  iWork = 2;
    lAm1[iWork]--;
    onei[iWork]=1;

    tmpVal = ijSP->PA[iPP][iWork]*this->pVpab(ijSP,C,LA-1,lAm1,LB,lB,m,iPP,iAtom);
    tmpVal -= PC[iWork]*this->pVpab(ijSP,C,LA-1,lAm1,LB,lB,m+1,iPP,iAtom);
    tmpVal -= ijSP->zetab[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetab[iPP]*this->hRRiPPVab(ijSP,LA-1,lAm1,LB+1,onei,C,m,iPP,iAtom);
    tmpVal += ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(2*ijSP->zetaa[iPP]*this->hRRiPPVab(ijSP,LA,lA,LB,lB,C,m,iPP,iAtom));
    tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetab[iPP]*this->hRRiPPVab(ijSP,LA-1,lAm1,LB+1,onei,C,m+1,iPP,iAtom);
    tmpVal -= ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(2*ijSP->zetaa[iPP]*this->hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+1,iPP,iAtom));

    if (lAm1[iWork]>0) {
      lAm1[iWork]--;
      tmpVal += ijSP->halfInvZeta[iPP]*(lA[iWork]-1)*(this->pVpab(ijSP,C,LA-2,lAm1,LB,lB,m,iPP,iAtom)-this->pVpab(ijSP,C,LA-2,lAm1,LB,lB,m+1,iPP,iAtom));
      tmpVal -= ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(lAm1[iWork]+1)*this->hRRiPPVab(ijSP,LA-2,lAm1,LB,lB,C,m,iPP,iAtom);
      tmpVal += ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(lAm1[iWork]+1)*this->hRRiPPVab(ijSP,LA-2,lAm1,LB,lB,C,m+1,iPP,iAtom);
    }
  }
  else if ((LA==0)&(LB>0)) {
cout<<"here we are wrong"<<endl;
  }

  else if ((LA>0)&(LB)>0) {
    int lBm1[3],lAp1[3];
    for ( iWork = 0 ; iWork<3 ; iWork++) {
      lAm1[iWork] = lA[iWork];
      lBm1[iWork] = lB[iWork];
      lAp1[iWork] = lA[iWork];
    }
    if (lB[0]>0) iWork = 0;
    else if (lB[1]>0) iWork = 1;
    else if (lB[2]>0) iWork = 2;
    
    lAm1[iWork]--;
    lBm1[iWork]--;
    lAp1[iWork]++;

    tmpVal = ijSP->PB[iPP][iWork]*this->pVpab(ijSP,C,LA,lA,LB-1,lBm1,m,iPP,iAtom);
    tmpVal -= PC[iWork]*this->pVpab(ijSP,C,LA,lA,LB-1,lBm1,m+1,iPP,iAtom);
    tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetaa[iPP]*this->hRRiPPVab(ijSP,LA+1,lAp1,LB-1,lBm1,C,m,iPP,iAtom);
    tmpVal += ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetab[iPP]*this->hRRiPPVab(ijSP,LA,lA,LB,lB,C,m,iPP,iAtom);
    tmpVal -= ijSP->zetab[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetaa[iPP]*this->hRRiPPVab(ijSP,LA+1,lAp1,LB-1,lBm1,C,m+1,iPP,iAtom);
    tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetab[iPP]*this->hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+1,iPP,iAtom);

    if (lA[iWork]>0) {
      tmpVal += ijSP->halfInvZeta[iPP]*lA[iWork]* (this->pVpab(ijSP,C,LA-1,lAm1,LB-1,lBm1,m,iPP,iAtom)-this->pVpab(ijSP,C,LA-1,lAm1,LB-1,lBm1,m+1,iPP,iAtom));
      tmpVal += ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*lA[iWork]*this->hRRiPPVab(ijSP,LA-1,lAm1,LB-1,lBm1,C,m,iPP,iAtom);
      tmpVal += ijSP->zetab[iPP]*ijSP->invZeta[iPP]*lA[iWork]*this->hRRiPPVab(ijSP,LA-1,lAm1,LB-1,lBm1,C,m+1,iPP,iAtom);
    }

    if (lBm1[iWork]>0) {
      lBm1[iWork]--;
      tmpVal += ijSP->halfInvZeta[iPP]*(lB[iWork]-1)*(this->pVpab(ijSP,C,LA,lA,LB-2,lBm1,m,iPP,iAtom)-this->pVpab(ijSP,C,LA,lA,LB-2,lBm1,m+1,iPP,iAtom));
      tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*(lB[iWork]-1)*this->hRRiPPVab(ijSP,LA,lA,LB-2,lBm1,C,m,iPP,iAtom);
      tmpVal += ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*(lB[iWork]-1)*this->hRRiPPVab(ijSP,LA,lA,LB-2,lBm1,C,m+1,iPP,iAtom);
    }

  }
  return tmpVal;
} 

/*  
void AOIntegrals::cart2sphtrans(ChronusQ::ShellPair *ijSP, RealMatrix *BeforeTransi, RealMatrix &Aftertransi, double sym) {
int p,q,i,j,bf1,bf2,cartbf1,cartbf2;
double tempVal;
RealMatrix Csp2cai(2*(ijSP->iShell->l)+1,(ijSP->iShell->l+1)*(ijSP->iShell->l+2)/2);

RealMatrix Csp2caj(2*(ijSP->jShell->l)+1,(ijSP->jShell->l+1)*(ijSP->jShell->l+2)/2);
for ( p = 0 ; p <ijSP->isphsize; p++ ) {
  for ( q = 0 ; q<ijSP->icarsize; q++ ) {
    Csp2cai(p,q) = 0;
  }
}
for ( p = 0 ; p <ijSP->jsphsize; p++ ) {
  for ( q = 0 ; q<ijSP->jcarsize; q++ ) {
    Csp2caj(p,q) = 0;
  }
}
 
if (ijSP->iShell->l==0) {
 Csp2cai(0,0) = 1;
}
else if (ijSP->iShell->l==1) {
  Csp2cai(0,0) = 1;
  Csp2cai(1,1) = 1;
  Csp2cai(2,2) = 1;
}
else if (ijSP->iShell->l==2) {
  Csp2cai(4,0) = sqrt(3)/2;
  Csp2cai(4,3) = -sqrt(3)/2;

  Csp2cai(3,2) = sqrt(3);

  Csp2cai(2,0) = -0.5;
  Csp2cai(2,3) = -0.5;
  Csp2cai(2,5) = 1;

  Csp2cai(1,4) = sqrt(3);

  Csp2cai(0,1) = sqrt(3);       
 }

if (ijSP->jShell->l==0) {
  Csp2caj(0,0) = 1;
}
else if (ijSP->jShell->l==1) {
  Csp2caj(0,0) =1;
  Csp2caj(1,1) =1;
  Csp2caj(2,2) =1;
}
else if (ijSP->jShell->l==2) {
  Csp2caj(4,0) = sqrt(3)/2;
  Csp2caj(4,3) = -sqrt(3)/2;

  Csp2caj(3,2) = sqrt(3);

  Csp2caj(2,0) = -0.5;
  Csp2caj(2,3) = -0.5;
  Csp2caj(2,5) = 1;

  Csp2caj(1,4) = sqrt(3);

  Csp2caj(0,1) = sqrt(3);
}

for( i = 0, bf1 = ijSP->isphbf_s; i<ijSP->isphsize ; i++,bf1++  ) {
  for( j = 0, bf2 = ijSP->jsphbf_s; j<ijSP->jsphsize ; j++,bf2++ ) {
    tempVal = 0;
    for( p = 0, cartbf1 = ijSP->icarbf_s; p<ijSP->icarsize; p++, cartbf1++) {
      for( q = 0, cartbf2 = ijSP->jcarbf_s; q<ijSP->jcarsize; q++, cartbf2++) {
        tempVal += Csp2cai(i,p)*Csp2caj(j,q)*(*BeforeTransi)(cartbf1,cartbf2);
      }
    }
    Aftertransi(bf1,bf2) = tempVal;
    Aftertransi(bf2,bf1) = sym*tempVal;
  }
}
   // end of cartesian to spherical transform
} */
   //
std::vector<double> AOIntegrals::cart2SphTrans(ChronusQ::ShellPair *ijSP, double *CART) {
  int L=ijSP->iShell->l;
  int M=ijSP->jShell->l; 
  std::vector<double> SPH((2*L+1)*(2*M+1));
  int CARTL=((L+1)*(L+2)/2);
  int CARTM=((M+1)*(M+2)/2);
  RealMatrix Csp2cai(2*L+1,CARTL);
  RealMatrix Csp2caj(2*M+1,CARTM);
  int i,j,p,q;
  double tempVal;
  int k,m,l[3]; // l is the angular momentum
  std::complex<double> cplxcoeff;  
  double scalecoeff;


  if( L < 2 && M < 2) {
    std::copy(CART,CART+CARTL*CARTM,SPH.begin());
  } else {
   for ( p = 0 ; p <(2*L+1); p++ ) {
     for ( q = 0 ; q<CARTL; q++ ) {
      Csp2cai(p,q) = 0;
     }
   }
   for ( p = 0 ; p <(2*M+1); p++ ) {
     for ( q = 0 ; q<CARTM; q++ ) {
       Csp2caj(p,q) = 0;
     }
   }
   if (L==0) {
     Csp2cai(0,0) = 1;
   }
   else if (L==1) {
     Csp2cai(0,0) = 1;
     Csp2cai(1,1) = 1;
     Csp2cai(2,2) = 1;
   }
   else {
   for ( p = 0 ; p <L+1; p++ ) {
    for ( q = 0 ; q<CARTL; q++ ) {
     for ( k = 0 ; k < 3 ; k++ ){
       l[k] = ijSP->iShell->cartesian_l[q][k];
     }
     m = -L+p;
     if ( m<0 ) {
       cplxcoeff = this->car2sphcoeff(L,m,l);
       Csp2cai(p,q) = sqrt(2)*(-cplxcoeff.imag());

       Csp2cai(2*L-p,q) = sqrt(2)*cplxcoeff.real();
       scalecoeff = sqrt(ChronusQ::doubleFact(L)/
         (ChronusQ::doubleFact(l[0])*ChronusQ::doubleFact(l[1])*ChronusQ::doubleFact(l[2])));
       Csp2cai(p,q) = Csp2cai(p,q)*scalecoeff;
       Csp2cai(2*L-p,q) = Csp2cai(2*L-p,q)*scalecoeff; 
     } else if (m ==0) {
       Csp2cai(p,q) = this->car2sphcoeff(L,m,l).real();
       scalecoeff = sqrt(ChronusQ::doubleFact(L)/
         (ChronusQ::doubleFact(l[0])*ChronusQ::doubleFact(l[1])*ChronusQ::doubleFact(l[2])));
       Csp2cai(p,q) = Csp2cai(p,q)*scalecoeff;
      
     }
   }
  }
  }

/*   } else if (L==2) {
     Csp2cai(4,0) = sqrt(3)/2;
     Csp2cai(4,3) = -sqrt(3)/2;

     Csp2cai(3,2) = sqrt(3);

     Csp2cai(2,0) = -0.5;
     Csp2cai(2,3) = -0.5;
     Csp2cai(2,5) = 1;

     Csp2cai(1,4) = sqrt(3);

     Csp2cai(0,1) = sqrt(3);       
   }
*/ 
   if (M==0) {
     Csp2caj(0,0) = 1;
   }
   else if (M==1) {
     Csp2caj(0,0) =1;
     Csp2caj(1,1) =1;
     Csp2caj(2,2) =1;
   }
   else {
  for ( p = 0 ; p <M+1; p++ ) {
    for ( q = 0 ; q<CARTM; q++ ) {
      for ( k = 0 ; k < 3 ; k++ ){
        l[k] = ijSP->jShell->cartesian_l[q][k];
      }
      m = -M+p;
      if ( m<0 ) {
       cplxcoeff = this->car2sphcoeff(M,m,l);
       Csp2caj(p,q) = sqrt(2)*(-cplxcoeff.imag());
       Csp2caj(2*M-p,q) = sqrt(2)*(cplxcoeff.real());
       scalecoeff = sqrt(ChronusQ::doubleFact(M)/
         (ChronusQ::doubleFact(l[0])*ChronusQ::doubleFact(l[1])*ChronusQ::doubleFact(l[2])));
       Csp2caj(p,q) = Csp2caj(p,q)*scalecoeff;
       Csp2caj(2*M-p,q) = Csp2caj(2*M-p,q)*scalecoeff;
 
     } else if (m ==0) {
       Csp2caj(p,q) = this->car2sphcoeff(M,m,l).real();
       scalecoeff = sqrt(ChronusQ::doubleFact(M)/
         (ChronusQ::doubleFact(l[0])*ChronusQ::doubleFact(l[1])*ChronusQ::doubleFact(l[2])));
       Csp2caj(p,q) = Csp2caj(p,q)*scalecoeff;
 
     }
   }
  }
  }

/*   else if (M==2) {
     Csp2caj(4,0) = sqrt(3)/2;
     Csp2caj(4,3) = -sqrt(3)/2;

     Csp2caj(3,2) = sqrt(3);

     Csp2caj(2,0) = -0.5;
     Csp2caj(2,3) = -0.5;
     Csp2caj(2,5) = 1;

     Csp2caj(1,4) = sqrt(3);

     Csp2caj(0,1) = sqrt(3);
   }
*/
     
   for( i = 0 ; i<2*L+1 ; i++  ) {
     for( j = 0 ; j<2*M+1 ; j++ ) {
       tempVal = 0;
       for( p = 0 ; p<CARTL ; p++ ) {
         for( q = 0 ; q<CARTM ; q++ ) {
           tempVal += Csp2cai(i,p)*Csp2caj(j,q)*CART[p*CARTM+q];
         }
       }
     SPH[i*(2*M+1)+j] = tempVal;
     }
    }   
  } 
// general case
/*  int m;  //m is the label of the component of spherical gaussian
  int k,l[3]; // l is the angular momentum
  std::complex<double> cplxcoeff;  
  double scalecoeff;

  for ( p = 0 ; p <L+1; p++ ) {
    for ( q = 0 ; q<CARTL; q++ ) {
     for ( k = 0 ; k < 3 ; k++ ){
       l[k] = ijSP->iShell->cartesian_l[q][k];
     }
     m = -L+p;
     if ( m<0 ) {
       cplxcoeff = this->car2sphcoeff(L,m,l);
       Csp2cai(p,q) = sqrt(2)*(-cplxcoeff.imag());

       Csp2cai(2*L-p,q) = sqrt(2)*cplxcoeff.real();
       scalecoeff = sqrt(ChronusQ::doubleFact(L)/
         (ChronusQ::doubleFact(l[0])*ChronusQ::doubleFact(l[1])*ChronusQ::doubleFact(l[2])));
       Csp2cai(p,q) = Csp2cai(p,q)*scalecoeff;
       Csp2cai(2*L-p,q) = Csp2cai(2*L-p,q)*scalecoeff; 
     } else if (m ==0) {
       Csp2cai(p,q) = this->car2sphcoeff(L,m,l).real();
       scalecoeff = sqrt(ChronusQ::doubleFact(L)/
         (ChronusQ::doubleFact(l[0])*ChronusQ::doubleFact(l[1])*ChronusQ::doubleFact(l[2])));
       Csp2cai(p,q) = Csp2cai(p,q)*scalecoeff;
      
     }
   }
  }
if (L==1) { 
   for ( p = 0 ; p <(2*L+1); p++ ) {
     for ( q = 0 ; q<CARTL; q++ ) {
      Csp2cai(p,q) = 0;
     }
   }
      Csp2cai(0,0) = 1;
     Csp2cai(1,1) = 1;
     Csp2cai(2,2) = 1;
 
}
  for ( p = 0 ; p <M+1; p++ ) {
    for ( q = 0 ; q<CARTM; q++ ) {
      for ( k = 0 ; k < 3 ; k++ ){
        l[k] = ijSP->jShell->cartesian_l[q][k];
      }
      m = -M+p;
      if ( m<0 ) {
       cplxcoeff = this->car2sphcoeff(M,m,l);
       Csp2caj(p,q) = sqrt(2)*(-cplxcoeff.imag());
       Csp2caj(2*M-p,q) = sqrt(2)*(cplxcoeff.real());
       scalecoeff = sqrt(ChronusQ::doubleFact(M)/
         (ChronusQ::doubleFact(l[0])*ChronusQ::doubleFact(l[1])*ChronusQ::doubleFact(l[2])));
       Csp2caj(p,q) = Csp2caj(p,q)*scalecoeff;
       Csp2caj(2*M-p,q) = Csp2caj(2*M-p,q)*scalecoeff;
 
     } else if (m ==0) {
       Csp2caj(p,q) = this->car2sphcoeff(M,m,l).real();
       scalecoeff = sqrt(ChronusQ::doubleFact(M)/
         (ChronusQ::doubleFact(l[0])*ChronusQ::doubleFact(l[1])*ChronusQ::doubleFact(l[2])));
       Csp2caj(p,q) = Csp2caj(p,q)*scalecoeff;
 
     }
   }
  }
if (M==1) {
   for ( p = 0 ; p <(2*M+1); p++ ) {
     for ( q = 0 ; q<CARTM; q++ ) {
       Csp2caj(p,q) = 0;
     }
   }
      Csp2caj(0,0) =1;
     Csp2caj(1,1) =1;
     Csp2caj(2,2) =1;
} 


  for( i = 0 ; i<2*L+1 ; i++  ) {
   for( j = 0 ; j<2*M+1 ; j++ ) {
     tempVal = 0;
     for( p = 0 ; p<CARTL ; p++ ) {
       for( q = 0 ; q<CARTM ; q++ ) {
         tempVal += Csp2cai(i,p)*Csp2caj(j,q)*CART[p*CARTM+q];
       }
     }
   SPH[i*(2*M+1)+j] = tempVal;
   }
  }  */   
  return SPH;
};

std::complex <double> AOIntegrals::car2sphcoeff(int L,int m,int *l ) {

  int Ltotal;
  std::complex <double> coeff(0.0);
  Ltotal = l[0]+l[1]+l[2];
  double tmp = 0.0;
  if (L!=Ltotal) {
//    coeff = 0.0;
    return  coeff;
  }
  double j;
  j = (double(l[0]+l[1])-std::abs(double(m)))/2;
  if (fmod(j,1)>0) {
//    coeff = 0.0;
    return coeff;
  }
  std::complex<double> sumval(0.0);
  std::complex<double> ttmmpp,sumsumval;
  std::complex<double> pref,absmchooselxm2k,ichoosej;
  int i,k;
  if (Ltotal == L) {
  pref = sqrt(ChronusQ::factorial(l[0]*2)*ChronusQ::factorial(2*l[1])*ChronusQ::factorial(2*l[2])*ChronusQ::factorial(L)*ChronusQ::factorial(L-std::abs(m))
     /(ChronusQ::factorial(2*L)*ChronusQ::factorial(l[0])*ChronusQ::factorial(l[1])*ChronusQ::factorial(l[2])*ChronusQ::factorial(L+std::abs(m))))
     /(ChronusQ::factorial(L)*pow(2,L));
  
  i = 0;
  
  while (i<=double((L-std::abs(m))/2) ) {
    sumsumval = 0.0;
    for ( k = 0 ; k <= j ; k++ ) {
      if (m>=0) {
        ttmmpp = double(std::abs(m)-l[0]+2*k)/2;
      }
      else {
        ttmmpp = -double(std::abs(m)-l[0]+2*k)/2;
      }
      
      if ((std::abs(m)>=(l[0]-2*k))&&((l[0]-2*k)>=0)) {
        absmchooselxm2k =ChronusQ::polyCoeff(std::abs(m),l[0]-2*k);
      }
      else {
        absmchooselxm2k = 0.0;
      }
      sumsumval = sumsumval + ChronusQ::polyCoeff(j,k)*absmchooselxm2k*pow(-1.0,ttmmpp);
    }
    if (i<j||(j<0)) {
       ichoosej = 0.0;
    }
    else {
      ichoosej = ChronusQ::polyCoeff(i,j);
    }
    sumval = sumval + ChronusQ::polyCoeff(L,i)*ichoosej*pow(-1,i)*ChronusQ::factorial(2*L-2*i)/(ChronusQ::factorial(L-std::abs(m)-2*i))*sumsumval;
    i = i + 1;
  }
  coeff = pref * sumval;
  return coeff;
  }
}; 




