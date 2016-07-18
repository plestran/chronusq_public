#include "aointegrals.h"
using ChronusQ::AOIntegrals;

void AOIntegrals::computeOverlapS(){
  double S;
  int i,j,k,ijShell,lA[3],lB[3];
  int bf1,bf2;
  ChronusQ::ShellPair *ijS;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(controls_->printLevel>=1) start = std::chrono::high_resolution_clock::now();
  this->overlap_->setZero();
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);
    if(ijS->iShell==ijS->jShell) {
      for(i=0,bf1=ijS->ibf_s  ; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=i,bf2=ijS->ibf_s+i; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        S = this->hRRSab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB);
        (*this->overlap_)(bf1,bf2) = S;
        (*this->overlap_)(bf2,bf1) = S;
      };
    } else {
      for(i=0,bf1=ijS->ibf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        S = this->hRRSab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB);
        (*this->overlap_)(bf1,bf2) = S;
        (*this->overlap_)(bf2,bf1) = S;
      };
    };
  };
  prettyPrint(this->fileio_->out,(*this->overlap_),"XSLI Overlap");
};

void AOIntegrals::computeKineticT(){
  double T;
  int i,j,k,ijShell,lA[3],lB[3];
  int bf1,bf2;
  ChronusQ::ShellPair *ijS;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(controls_->printLevel>=1) start = std::chrono::high_resolution_clock::now();
  this->potential_->setZero();
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);
    if(ijS->iShell==ijS->jShell) {
      for(i=0,bf1=ijS->ibf_s  ; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=i,bf2=ijS->ibf_s+i; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        T = this->RRTab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB);
//        cout<<"VV: "<<bf1<<" "<<bf2<<" "<<V<<endl;
        (*this->kinetic_)(bf1,bf2) = T;
        (*this->kinetic_)(bf2,bf1) = T;
      };
    } else {
      for(i=0,bf1=ijS->ibf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        T = this->RRTab(ijS,ijS->iShell->l,lA,ijS->jShell->l,lB);
//        cout<<"T: "<<bf1<<" "<<bf2<<" "<<T<<endl;
        (*this->kinetic_)(bf1,bf2) = T;
        (*this->kinetic_)(bf2,bf1) = T;
      };
    };
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
  ChronusQ::ShellPair *ijS;
  std::chrono::high_resolution_clock::time_point start,finish;
  if(controls_->printLevel>=1) start = std::chrono::high_resolution_clock::now();
  this->kinetic_->setZero();
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
//    cout<<"ShellPairs"<<endl;
    ijS = &(this->shellPairs_[ijShell]);
    if(ijS->iShell==ijS->jShell) {
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
    } else {
      for(i=0,bf1=ijS->ibf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        V = this->hRRVab(ijS,this->molecularConstants_.get(),ijS->iShell->l,lA,ijS->jShell->l,lB);
        (*this->potential_)(bf1,bf2) = -V; 
        (*this->potential_)(bf2,bf1) = -V;
      };
    };
  };
  prettyPrint(this->fileio_->out,(*this->potential_),"XSLI Potential");
};



//SS
void AOIntegrals::computeAngularL(){
  double L[3],LL[3];
  int i,j,k,ijShell,mu,bf1,bf2,lA[3],lB[3],iPP,q,r;
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




 
       for(i=0,bf1=ijS->ibf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++)
       for(j=0,bf2=ijS->jbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
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
         (*this->angular_)(bf1,bf2,mu) = L[mu];
         (*this->angular_)(bf2,bf1,mu) = -L[mu];
         }
       }
//     }
   }
  for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
    RealMap TEMP(&this->angular_->storage()[iXYZ*this->nBasis_*this->nBasis_],this->nBasis_,this->nBasis_);
    prettyPrint(this->fileio_->out,TEMP,"Angular XYZ="+std::to_string(iXYZ));
  }
//  prettyPrint(this->fileio_->out,(*this->kinetic_),"SS Angular");
};


void AOIntegrals::computeSL(){
  double Slmu,SlC,C[3];
  int i,j,k,ijShell,mu,nu,lA[3],lB[3],iPP,iAtom;
  int bf1,bf2;
  ChronusQ::ShellPair *ijS;

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
      for(i=0,bf1=ijS->ibf_s; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=0,bf2=ijS->jbf_s; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        for (mu = 0; mu < 3; mu++) {
          Slmu = 0;
//cout<<"Sl[mu]\tbefore\t"<<Slmu<<endl;
          for (iPP = 0; iPP<ijS->nPGTOPair; iPP++){
            SlC = 0;
            for (iAtom = 0; iAtom <this->molecularConstants_->nAtom; iAtom++){
//cout<<"number of atoms"<<this->molecularConstants_->nAtom<<endl;
              // get the coordinate of C
              for (nu = 0; nu<3 ; nu++){
                C[nu] = this->molecularConstants_->cart[nu][iAtom];
//cout<<"C[nu]"<<C[nu]<<endl;
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

//cout<<"OneixBC"<<endl;
//cout<<"atomic number\t"<<this->molecularConstants_->atomZ[iAtom]<<endl;
 
              SlC += this->molecularConstants_->atomZ[iAtom]*this->Slabmu(ijS,&OneixAC,&OneixBC,&OneiOnek,C,ijS->iShell->l,lA,ijS->jShell->l,lB,mu,0,iPP);
//cout<<"SlC:\t"<<SlC<<endl;
              };
//cout<<"SlC\tfinished\tSl:\t"<<Slmu<<"\tmu\t"<<mu<<endl;

//cout<<"SlC\t"<<SlC<<endl;
            Slmu +=SlC;
//cout<<"Sl[mu]:\t"<<Slmu<<endl;
            };
          (*this->SOCoupling_)(bf1,bf2,mu) = Slmu;
          (*this->SOCoupling_)(bf2,bf1,mu) = -Slmu;
//cout<<"SOCoupling\tbf1\t"<<bf1<<"\tbf2\t"<<bf2<<"\tSOCoupling\t"<<Slmu<<endl;
        }
      }
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
  int bf1,bf2;
  ChronusQ::ShellPair *ijS;
  for(ijShell=0;ijShell<this->nShellPair_;ijShell++) {
    ijS = &(this->shellPairs_[ijShell]);
/*    if(ijS->iShell==ijS->jShell) {
      for(i=0,bf1=ijS->ibf_s  ; i<ijS->iShell->cartesian_l.size(); i++,bf1++) 
      for(j=i,bf2=ijS->ibf_s+i; j<ijS->jShell->cartesian_l.size(); j++,bf2++){
        for(k=0;k<3;k++){
          lA[k]=ijS->iShell->cartesian_l[i][k];
          lB[k]=ijS->jShell->cartesian_l[j][k];
        };
        pVp = this->hRRVab(ijS,this->molecularConstants_.get(),ijS->iShell->l,lA,ijS->jShell->l,lB);
        (*this->pVp_)(bf1,bf2) = pVp;
        (*this->pVp_)(bf2,bf1) = pVp;
      };
    } else {    */
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
             pVpC += this->molecularConstants_->atomZ[iAtom]*this->pVpab(ijS,C,ijS->iShell->l,lA,ijS->jShell->l,lB,0,iPP);

           }
         pVp += pVpC;
        }  
        (*this->pVp_)(bf1,bf2) = pVp; 
        (*this->pVp_)(bf2,bf1) = pVp;
      };
    };
 
  prettyPrint(this->fileio_->out,(*this->pVp_),"pVp");
  
  double fobi;
  fobi = this->pVp_->frobInner(*this->pVp_);
  this->fileio_->out <<"fobinius inner product" << fobi << endl;

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
      this->computeFmTTaylor(tmpFmT,ijSP->Zeta[iPP]*squarePC,ijSP->lTotal,0);
      if(LA==0) tmpVal += (static_cast<double>(mc->atomZ[iAtom]))*ijSP->ssV[iPP]*tmpFmT[0];
      else tmpVal += mc->atomZ[iAtom]*ijSP->ssV[iPP]*this->vRRVa0(ijSP,tmpFmT,PC,0,LA,lA,iPP);
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
double AOIntegrals::vRRVa0(ChronusQ::ShellPair *ijSP,double *FmT,double *PC,int m,int LA,int *lA,int iPP){
  if(LA==0) return FmT[m];  //Z*2*sqrt[zeta/pi]*[s|s]is given in hRRVab, in ssV.

  double tmpVal=0.0;
  int lAm1[3]; //means lA minus 1_i
  int iWork;
  for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
  if (lA[0]>0) iWork=0;
  else if (lA[1]>0) iWork=1;
  else if (lA[2]>0) iWork=2;
  lAm1[iWork]--;

  tmpVal  = ijSP->PA[iPP][iWork]*this->vRRVa0(ijSP,FmT,PC,m,LA-1,lAm1,iPP);
  tmpVal -= PC[iWork]*this->vRRVa0(ijSP,FmT,PC,m+1,LA-1,lAm1,iPP);

  if(lAm1[iWork]>=1){
    lAm1[iWork]--;
    tmpVal += (lAm1[iWork]+1)*ijSP->halfInvZeta[iPP]*(this->vRRVa0(ijSP,FmT,PC,m,LA-2,lAm1,iPP)
            - this->vRRVa0(ijSP,FmT,PC,m+1,LA-2,lAm1,iPP));
  };
  return tmpVal; 
};

// SS
//----------------------------------------------------------------------------//
// potential integral horizontal recursion iPP specific                       //
//  [a|A(0)|b]^(m) = [a+1i|A(0)|b-1i]^(m) + (Ai-Bi)*[a|A(0)|b-1i]^(m)         //
//  LA >=LB                                                                   //
//  m can be nonzero                                                          //
//----------------------------------------------------------------------------//

double AOIntegrals::hRRiPPVab(ChronusQ::ShellPair *ijSP,int LA,int *lA,int LB,int *lB,double *C,int m,int iPP){
  double tmpVal=0.0;
  double PC[3];
  int k;
  if(LB==0) {
    // (LA|s)
      double squarePC=0.0;
      for(k=0;k<3;k++) {
        PC[k] = ijSP->P[iPP][k] - C[k];
        squarePC += PC[k]*PC[k];
      };
      double *tmpFmT = new double[ijSP->lTotal+m+1];
      this->computeFmTTaylor(tmpFmT,ijSP->Zeta[iPP]*squarePC,ijSP->lTotal+m,0);
      if(LA==0) {
        tmpVal = ijSP->ssV[iPP]*tmpFmT[m];
//cout<<"ssV\t"<<tmpVal<<endl;
//cout<<"TOTAL ORDER OF BOYS FUNCTION\t"<<ijSP->lTotal+m<<endl;
//cout<<"m equals\t"<<m<<endl;
      }

      else {

//        cout<<"(a|V|0)^(m)"<<endl;
        tmpVal = ijSP->ssV[iPP]*this->vRRVa0(ijSP,tmpFmT,PC,m,LA,lA,iPP);
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
      this->computeFmTTaylor(tmpFmT,ijSP->Zeta[iPP]*squarePC,ijSP->lTotal+m,0);

//        cout<<"(a|V|0)^(m)"<<endl;
        tmpVal = ijSP->ssV[iPP]*this->vRRV0b(ijSP,tmpFmT,PC,m,LB,lB,iPP);

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
    tmpVal = this->hRRiPPVab(ijSP,LA+1,lAp1,LB-1,lBm1,C,m,iPP);
    tmpVal+= ijSP->AB[iWork]*this->hRRiPPVab(ijSP,LA,lA,LB-1,lBm1,C,m,iPP);
  };
  return tmpVal;
}

//----------------------------------------------------------------------------//
// potential integral vertical recursion                                      //
//  (0|0_c|b)^(m) = (P-B)(0|0_c|b-1i)^(m) - (P-C)(0|0_c|b-1i)^(m+1)           //
//                + halfInvZeta*N_(b-1)*[(0|0_c|b-2i)^(m)-(0|0_c|b-2i)^(m+1)] //
//  since LA==0, we only decrease b                                           //
//----------------------------------------------------------------------------//
double AOIntegrals::vRRV0b(ChronusQ::ShellPair *ijSP,double *FmT,double *PC,int m,int LB,int *lB,int iPP){
  if(LB==0) return FmT[m];  //Z*2*sqrt[zeta/pi]*[s|s]is given in hRRVab, in ssV.

  double tmpVal=0.0;
  int lBm1[3]; //means lA minus 1_i
  int iWork;
  for(iWork=0;iWork<3;iWork++) lBm1[iWork]=lB[iWork];
  if (lB[0]>0) iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lBm1[iWork]--;

  tmpVal  = ijSP->PB[iPP][iWork]*this->vRRV0b(ijSP,FmT,PC,m,LB-1,lBm1,iPP);
  tmpVal -= PC[iWork]*this->vRRV0b(ijSP,FmT,PC,m+1,LB-1,lBm1,iPP);

  if(lBm1[iWork]>=1){
    lBm1[iWork]--;
    tmpVal += (lBm1[iWork]+1)*ijSP->halfInvZeta[iPP]*(this->vRRV0b(ijSP,FmT,PC,m,LB-2,lBm1,iPP)
            - this->vRRV0b(ijSP,FmT,PC,m+1,LB-2,lBm1,iPP));
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
double AOIntegrals::Slabmu(ChronusQ::ShellPair *ijSP,RealMatrix *OneixAC,RealMatrix *OneixBC,RealTensor3d *OneiOnek,double *C,int LA,int *lA,int LB,int *lB,int mu,int m,int iPP){
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
u = hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+1,iPP);
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
   tmpVal  = ijSP->PA[iPP][iWork]*this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-1,lAm1,LB,lB,mu,m,iPP);
   tmpVal -= PC[iWork]*this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-1,lAm1,LB,lB,mu,m+1,iPP);
   tmpVal += 2*ijSP->zetab[iPP]*(*OneixBC)(iWork,mu)*this->hRRiPPVab(ijSP,LA-1,lAm1,LB,lB,C,m+1,iPP); 
   if(lAm1[iWork]>=1){
     lAm1[iWork]--;
     tmpVal += (lA[iWork]-1)*ijSP->halfInvZeta[iPP]*(this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-2,lAm1,LB,lB,mu,m,iPP) - this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-2,lAm1,LB,lB,mu,m+1,iPP));
   };

   return tmpVal;

  }else if ((LA>0)&(LB==0)) cout<<"here"; 
   else if ((LB>0)&(LA>0)){
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
   tmpVal  = ijSP->PB[iPP][iWork]*this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA,lA,LB-1,lBm1,mu,m,iPP);

   tmpVal -= PC[iWork]*this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA,lA,LB-1,lBm1,mu,m+1,iPP);
   if (lAm1[iWork]>=0){
    tmpVal += ijSP->halfInvZeta[iPP]*lA[iWork]*(this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-1,lAm1,LB-1,lBm1,mu,m,iPP) - this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA-1,lAm1,LB,lBm1,mu,m+1,iPP));
   };
   tmpVal -= 2*ijSP->zetaa[iPP]*(*OneixAC)(iWork,mu)*this->hRRiPPVab(ijSP,LA,lA,LB-1,lBm1,C,m+1,iPP);
   for(k=0;k<3;k++) {
     if(lA[k]>0){
      lAm1k[k]--;
      tmpVal -= lA[k]*(*OneiOnek)(iWork,k,mu)*this->hRRiPPVab(ijSP,LA-1,lAm1k,LB-1,lBm1,C,m+1,iPP);
      lAm1k[k]++;
     };
   };
   
   if (lBm1[iWork]>0){
    lBm1[iWork]--;
    tmpVal += (lBm1[iWork]+1)*ijSP->halfInvZeta[iPP]*(this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA,lA,LB-2,lBm1,mu,m,iPP)-this->Slabmu(ijSP,OneixAC,OneixBC,OneiOnek,C,LA,lA,LB-2,lBm1,mu,m+1,iPP));
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
double AOIntegrals::pVpab(ChronusQ::ShellPair *ijSP,double *C,int LA,int *lA,int LB,int *lB,int m,int iPP){
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
  u = hRRiPPVab(ijSP,LA,lA,LB,lB,C,m,iPP);
  uone = hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+1,iPP);
  utwo = hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+2,iPP);
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

    tmpVal = ijSP->PA[iPP][iWork]*this->pVpab(ijSP,C,LA-1,lAm1,LB,lB,m,iPP);
    tmpVal -= PC[iWork]*this->pVpab(ijSP,C,LA-1,lAm1,LB,lB,m+1,iPP);
    tmpVal -= ijSP->zetab[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetab[iPP]*this->hRRiPPVab(ijSP,LA-1,lAm1,LB+1,onei,C,m,iPP);
    tmpVal += ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(2*ijSP->zetaa[iPP]*this->hRRiPPVab(ijSP,LA,lA,LB,lB,C,m,iPP));
    tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetab[iPP]*this->hRRiPPVab(ijSP,LA-1,lAm1,LB+1,onei,C,m+1,iPP);
    tmpVal -= ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(2*ijSP->zetaa[iPP]*this->hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+1,iPP));

    if (lAm1[iWork]>=1) {
      lAm1[iWork]--;
      tmpVal += ijSP->halfInvZeta[iPP]*(lA[iWork]-1)*(this->pVpab(ijSP,C,LA-2,lAm1,LB,lB,m,iPP)-this->pVpab(ijSP,C,LA-2,lAm1,LB,lB,m+1,iPP));
      tmpVal -= ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(lAm1[iWork]+1)*this->hRRiPPVab(ijSP,LA-2,lAm1,LB,lB,C,m,iPP);
      tmpVal += ijSP->zetab[iPP]*ijSP->invZeta[iPP]*(lAm1[iWork]+1)*this->hRRiPPVab(ijSP,LA-2,lAm1,LB,lB,C,m+1,iPP);
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

    tmpVal = ijSP->PB[iPP][iWork]*this->pVpab(ijSP,C,LA,lA,LB-1,lBm1,m,iPP);
    tmpVal -= PC[iWork]*this->pVpab(ijSP,C,LA,lA,LB-1,lB,m+1,iPP);
    tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetaa[iPP]*this->hRRiPPVab(ijSP,LA+1,lAp1,LB-1,lBm1,C,m,iPP);
    tmpVal += ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetab[iPP]*this->hRRiPPVab(ijSP,LA,lA,LB,lB,C,m,iPP);
    tmpVal -= ijSP->zetab[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetaa[iPP]*this->hRRiPPVab(ijSP,LA+1,lAp1,LB-1,lBm1,C,m+1,iPP);
    tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*2*ijSP->zetab[iPP]*this->hRRiPPVab(ijSP,LA,lA,LB,lB,C,m+1,iPP);

    if (lAm1[iWork]>=0) {
      tmpVal += ijSP->halfInvZeta[iPP]*lA[iWork]* (this->pVpab(ijSP,C,LA-1,lAm1,LB-1,lBm1,m,iPP)-this->pVpab(ijSP,C,LA-1,lAm1,LB-1,lBm1,m+1,iPP));
      tmpVal += ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*lA[iWork]*this->hRRiPPVab(ijSP,LA-1,lAm1,LB-1,lBm1,C,m,iPP);
      tmpVal += ijSP->zetab[iPP]*ijSP->invZeta[iPP]*lA[iWork]*this->hRRiPPVab(ijSP,LA-1,lAm1,LB-1,lBm1,C,m+1,iPP);
    }

    if (lBm1[iWork]>0) {
      lBm1[iWork]--;
      tmpVal += ijSP->halfInvZeta[iPP]*(lB[iWork]-1)*(this->pVpab(ijSP,C,LA,lA,LB-2,lBm1,m,iPP)-this->pVpab(ijSP,C,LA,lA,LB-2,lBm1,m+1,iPP));
      tmpVal -= ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*(lB[iWork]-1)*this->hRRiPPVab(ijSP,LA,lA,LB-2,lBm1,C,m,iPP);
      tmpVal += ijSP->zetaa[iPP]*ijSP->invZeta[iPP]*(lB[iWork]-1)*this->hRRiPPVab(ijSP,LA,lA,LB-2,lBm1,C,m+1,iPP);
    }

  }
  return tmpVal;
} 
  



