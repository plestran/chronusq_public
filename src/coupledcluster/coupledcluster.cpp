/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#include <coupledcluster.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::MOIntegrals;
using ChronusQ::SingleSlater;
using ChronusQ::CoupledCluster;

void CoupledCluster::iniCoupledCluster( Molecule * molecule, BasisSet * basisSet, MOIntegrals<double> * mointegrals, 
                                FileIO * fileio, Controls * controls, SingleSlater<double> * singleSlater) {
  this->nBasis_         = basisSet->nBasis();
  this->nTCS_           = singleSlater->nTCS();
  this->molecule_       = molecule;
  this->basisSet_       = basisSet;
  this->fileio_         = fileio;
  this->controls_       = controls;
  this->mointegrals_    = mointegrals;
  this->singleSlater_   = singleSlater;
  this->Ref_            = singleSlater->Ref();
  this->nOA_            = this->singleSlater_->nOccA();
  this->nOB_            = this->singleSlater_->nOccB();
  this->nVA_            = this->singleSlater_->nVirA();
  this->nVB_            = this->singleSlater_->nVirB();
  this->nO_             = this->nOA_ + this->nOB_;
  this->nV_             = this->nVA_ + this->nVB_;
  // forms all necessary double bar integrals, replaces them in MOint objects in Dirac notation
  this->mointegrals_->formDBar();
}

double CoupledCluster::CCSD(){
  double ECorr = 0.0;
  double ECCSD = 0.0;
  double EInit = 0.0;
  double Denom = 1.0;
  double E1    = 0.0;
  double E2    = 0.0;
  double E3    = 0.0;
  double E4    = 0.0;
  double E5    = 0.0;

  // One-time intermediates
  RealTensor2d Zbc, Tjk;
  RealTensor4d Sijkl, Yjkbc;
  // Intermediates
  RealTensor2d Hba, Hbj, Hij, Gca, Gik; 
  RealTensor4d Aijkl, Babcd, Hicak;
  // Amplitudes
  RealTensor2d Tia, Wia;
  RealTensor4d Tijab, Tau, Wijab;

  if(this->Ref_ == SingleSlater<double>::TCS){
    // Static (one-time) intermediates
    Sijkl   = RealTensor4d(this->nO_,this->nO_,this->nO_,this->nO_);
    Yjkbc   = RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_);
    Zbc     = RealTensor2d(this->nV_,this->nV_);
    Tjk     = RealTensor2d(this->nO_,this->nO_);
    // Intermediates
    Hba     = RealTensor2d(this->nV_,this->nV_);
    Hbj     = RealTensor2d(this->nV_,this->nO_);
    Hij     = RealTensor2d(this->nO_,this->nO_);
    Gik     = RealTensor2d(this->nO_,this->nO_);
    Gca     = RealTensor2d(this->nV_,this->nV_);
    Aijkl   = RealTensor4d(this->nO_,this->nO_,this->nO_,this->nO_);
    Babcd   = RealTensor4d(this->nV_,this->nV_,this->nV_,this->nV_);
    Hicak   = RealTensor4d(this->nO_,this->nV_,this->nV_,this->nO_);
    // Amplitudes 
    Tia     = RealTensor2d(this->nO_,this->nV_);
    Wia     = RealTensor2d(this->nO_,this->nV_);
    Tijab   = RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_);
    Tau     = RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_);
    Wijab   = RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_);


    // Initialize T1 and T2 amplitudes
    for(auto i = 0; i < this->nO_; i++)
    for(auto a = 0; a < this->nV_; a++) {
      Tia(i,a) = 0.0;
      }
    for(auto i = 0; i < this->nO_; i++) 
    for(auto j = 0; j < this->nO_; j++) 
    for(auto a = 0; a < this->nV_; a++) 
    for(auto b = 0; b < this->nV_; b++) { 
      Denom = (*this->singleSlater_->epsA())(i)
               + (*this->singleSlater_->epsA())(j)
               - (*this->singleSlater_->epsA())(a + this->nO_)
               - (*this->singleSlater_->epsA())(b + this->nO_);
      Tijab(i,j,a,b) = (this->mointegrals_->IJAB(i,j,a,b))/Denom; 
    }

    // Form our first Tau
    for(auto i = 0; i < this->nO_; i++) 
    for(auto j = 0; j < this->nO_; j++) 
    for(auto a = 0; a < this->nV_; a++) 
    for(auto b = 0; b < this->nV_; b++) { 
      Tau(i,j,a,b) = Tijab(i,j,a,b) + (Tia(i,a)*Tia(j,b) - Tia(i,b)*Tia(j,a));
      EInit += (0.25)*(this->mointegrals_->IJAB(i,j,a,b))*Tijab(i,j,a,b);
    }
    // Check... if we made Tau correctly we should get MP2 energy expression back
    cout << "MP2 ENERGY = " << endl; 
    cout << EInit << endl;

    // Sijkl intermediate
    for(auto i = 0; i < this->nO_; i++) 
    for(auto j = 0; j < this->nO_; j++) 
    for(auto k = 0; k < this->nO_; k++) 
    for(auto l = 0; l < this->nO_; l++) { 
      Sijkl(i,j,k,l) += (this->mointegrals_->IJKL(i,j,k,l))*(0.5);
      for(auto c = 0; c < this->nV_; c++) {
        Sijkl(i,j,k,l) -= (this->mointegrals_->IJKA(i,j,l,c))*Tia(k,c);
        for(auto d = 0; d < this->nV_; d++) {
          Sijkl(i,j,k,l) += (this->mointegrals_->IJAB(i,j,c,d))*Tau(k,l,c,d)*(0.25);
        }
      }
    }

   // Zbc intermediate
   for(auto b = 0; b < this->nV_; b++)
   for(auto c = 0; c < this->nV_; c++)
   for(auto k = 0; k < this->nO_; k++)
   for(auto d = 0; d < this->nV_; d++) {
     Zbc(b,c) -= (this->mointegrals_->IABC(k,c,d,b))*Tia(k,d);
     for(auto l = 0; l < this->nO_; l++) {
       Zbc(b,c) -= (this->mointegrals_->IJAB(k,l,d,b))*Tau(k,l,c,d)*(0.5);
     }
   }
   // Tjk intermediate
   for(auto j = 0; j < this->nO_; j++)
   for(auto k = 0; k < this->nO_; k++)
   for(auto l = 0; l < this->nO_; l++)
   for(auto c = 0; c < this->nV_; c++) {
     Tjk(j,k) -= (this->mointegrals_->IJKA)(l,j,k,c)*Tia(l,c);
     for(auto d = 0; d < this->nV_; d++) {
       Tjk(j,k) -= (this->mointegrals_->IJAB(l,j,c,d))*Tau(k,l,c,d)*(0.5);
     }
   }

   // Yjkbc intermediate
   for(auto j = 0; j < this->nO_; j++)
   for(auto k = 0; k < this->nO_; k++)
   for(auto b = 0; b < this->nV_; b++)
   for(auto c = 0; c < this->nV_; c++) {
      Yjkbc(j,k,b,c) -= this->mointegrals_->IAJB(j,c,k,b);
     for(auto d = 0; d < this->nV_; d++) {
       Yjkbc(j,k,b,c) -= (this->mointegrals_->IABC(j,c,d,b))*Tia(k,d);
     }
     for(auto l = 0; l < this->nO_; l++) {
       Yjkbc(j,k,b,c) -= (this->mointegrals_->IJKA(l,j,k,b))*Tia(l,c);
       for(auto d = 0; d < this->nV_; d++) {
         Yjkbc(j,k,b,c) += (this->mointegrals_->IJAB(l,j,d,b))*
                          (Tijab(k,l,c,d) - Tia(l,c)*Tia(k,d));
       }
     }
   }

   // Done with one time intermediates, begin iterations
   for(auto niter = 0; niter < 40; niter++) {
     // Hab
     for(auto b = 0; b < this->nV_; b++)
     for(auto a = 0; a < this->nV_; a++)
     for(auto j = 0; j < this->nO_; j++)
     for(auto k = 0; k < this->nO_; k++)
     for(auto c = 0; c < this->nV_; c++) {
       Hba(b,a) -= (this->mointegrals_->IJAB(j,k,b,c))*Tau(j,k,a,c)*(0.5);
     }

    // Hij
     for(auto i = 0; i < this->nO_; i++)
     for(auto j = 0; j < this->nO_; j++)
     for(auto k = 0; k < this->nO_; k++)
     for(auto b = 0; b < this->nV_; b++)
     for(auto c = 0; c < this->nV_; c++) {
       Hij(i,j) += (this->mointegrals_->IJAB(j,k,b,c))*Tau(i,k,b,c)*(0.5);
     }

   // Hbj
     for(auto b = 0; b < this->nV_; b++)
     for(auto j = 0; j < this->nO_; j++)
     for(auto k = 0; k < this->nO_; k++)
     for(auto c = 0; c < this->nV_; c++) {
       Hbj(b,j) += (this->mointegrals_->IJAB(j,k,b,c))*Tia(k,c);
     }

   // Gca
     for(auto c = 0; c < this->nV_; c++) 
     for(auto a = 0; a < this->nV_; a++) {
       Gca(c,a) += Hba(c,a);
       for(auto k = 0; k < this->nO_; k++) 
       for(auto d = 0; d < this->nV_; d++) {
         Gca(c,a) -= (this->mointegrals_->IABC(k,a,c,d))*Tia(k,d);
       } 
     }


   // Gik
     for(auto i = 0; i < this->nO_; i++)
     for(auto k = 0; k < this->nO_; k++) {
       Gik(i,k) += Hij(i,k);
       for(auto l = 0; l < this->nO_; l++)
       for(auto c = 0; c < this->nV_; c++) {
         Gik(i,k) += (this->mointegrals_->IJKA(k,l,i,c))*Tia(l,c);
       }
     }

   // Aijkl
     for(auto i = 0; i < this->nO_; i++)
     for(auto j = 0; j < this->nO_; j++)
     for(auto k = 0; k < this->nO_; k++)
     for(auto l = 0; l < this->nO_; l++) {
       Aijkl(i,j,k,l) += (this->mointegrals_->IJKL(i,j,k,l));
       for(auto c = 0; c < this->nV_; c++) {
         Aijkl(i,j,k,l) += (this->mointegrals_->IJKA(k,l,i,c))*Tia(j,c) 
                          -(this->mointegrals_->IJKA(k,l,j,c))*Tia(i,c);
         for(auto d = 0; d < this->nV_; d++) {
           Aijkl(i,j,k,l) += (this->mointegrals_->IJAB(k,l,c,d))*Tau(i,j,c,d)*(0.5);
         }
       }
     }

   // Babcd
     
     for(auto c = 0; c < this->nV_; c++)
     for(auto d = 0; d < this->nV_; d++)
     for(auto a = 0; a < this->nV_; a++)
     for(auto b = 0; b < this->nV_; b++) {
       Babcd(a,b,c,d) += (this->mointegrals_->ABCD(a,b,c,d));
       for(auto k = 0; k < this->nO_; k++) {
         Babcd(a,b,c,d) += (this->mointegrals_->IABC(k,a,c,d))*Tia(k,b) 
                          -(this->mointegrals_->IABC(k,b,c,d))*Tia(k,a);
       }
     }

  // Hicak 
     for(auto i = 0; i < this->nO_; i++)
     for(auto c = 0; c < this->nV_; c++)
     for(auto a = 0; a < this->nV_; a++)
     for(auto k = 0; k < this->nO_; k++) {
       Hicak(i,c,a,k) -= (this->mointegrals_->IAJB(i,c,k,a));
       for(auto d = 0; d < this->nV_; d++) {
         Hicak(i,c,a,k) -= (this->mointegrals_->IABC(k,a,d,c))*Tia(i,d);  
       }
       for(auto l = 0; l < this->nO_; l++)  
       for(auto d = 0; d < this->nV_; d++) { 
         Hicak(i,c,a,k) -= (this->mointegrals_->IJAB(k,l,c,d))*
                           (Tijab(i,l,d,a)*(0.5) + Tia(i,d)*Tia(l,a));
       }
       for(auto l = 0; l < this->nO_; l++) { 
         Hicak(i,c,a,k) -= (this->mointegrals_->IJKA(l,k,i,c))*Tia(l,a);
       }
     }

  // Wia time!
     for(auto i = 0; i < this->nO_; i++) 
     for(auto a = 0; a < this->nV_; a++) 
     for(auto j = 0; j < this->nO_; j++) {
       Wia(i,a) -= Hij(i,j)*Tia(j,a);
       for(auto b = 0; b < this->nV_; b++)  
       for(auto c = 0; c < this->nV_; c++) { 
         Wia(i,a) -= (this->mointegrals_->IABC(j,a,b,c))*Tau(i,j,b,c)*(0.5);
       }
     }
     for(auto i = 0; i < this->nO_; i++) 
     for(auto a = 0; a < this->nV_; a++) { 
       for(auto b = 0; b < this->nV_; b++) {
         Wia(i,a) += Hba(b,a)*Tia(i,b);
       }
       for(auto j = 0; j < this->nO_; j++)
       for(auto b = 0; b < this->nV_; b++) {
         Wia(i,a) += Hbj(b,j)*(Tijab(i,j,a,b) + Tia(i,b)*Tia(j,a))
                    -(this->mointegrals_->IAJB(i,b,j,a)*Tia(j,b));
         for(auto k = 0; k < this->nO_; k++) {
           Wia(i,a) -= (this->mointegrals_->IJKA(j,k,i,b))*Tau(j,k,a,b)*(0.5);
         }
       }
     }
     
   // Wijab
     for(auto i = 0; i < this->nO_; i++) 
     for(auto j = 0; j < this->nO_; j++) 
     for(auto a = 0; a < this->nV_; a++) 
     for(auto b = 0; b < this->nV_; b++) { 
       for(auto k = 0; k < this->nO_; k++)  
       for(auto l = 0; l < this->nO_; l++) {
       Wijab(i,j,a,b) += Aijkl(i,j,k,l)*Tau(k,l,a,b)*(0.5);
       }
       E1 += Wijab(i,j,a,b)*Tijab(i,j,a,b);
       for(auto c = 0; c < this->nV_; c++)  
       for(auto d = 0; d < this->nV_; d++) {
         Wijab(i,j,a,b) += Babcd(a,b,c,d)*Tau(i,j,c,d)*(0.5);
       }
       E2 += Wijab(i,j,a,b)*Tijab(i,j,a,b);
       for(auto c = 0; c < this->nV_; c++) {
         Wijab(i,j,a,b) += Gca(c,a)*Tijab(i,j,c,b) -
                           Gca(c,b)*Tijab(i,j,c,b) -
                           (this->mointegrals_->IABC(j,c,a,b))*Tia(i,c) +
                           (this->mointegrals_->IABC(i,c,a,b))*Tia(j,c); 
       }
       E3 += Wijab(i,j,a,b)*Tijab(i,j,a,b);
       for(auto k = 0; k < this->nO_; k++) {
         Wijab(i,j,a,b) += (this->mointegrals_->IJKA(i,j,k,a))*Tia(k,b) -
                           (this->mointegrals_->IJKA(i,j,k,b))*Tia(k,a) -
                           Gik(i,k)*Tijab(k,j,a,b) +
                           Gik(j,k)*Tijab(k,i,a,b); 
       }
       E4 += Wijab(i,j,a,b)*Tijab(i,j,a,b);
       for(auto k = 0; k < this->nO_; k++) 
       for(auto c = 0; c < this->nV_; c++) {
         Wijab(i,j,a,b) += Hicak(i,c,a,k)*Tijab(j,k,b,c) -
                           Hicak(j,c,a,k)*Tijab(i,k,b,c) -
                           Hicak(i,c,b,k)*Tijab(j,k,a,c) +
                           Hicak(j,c,b,k)*Tijab(i,k,a,c) +
                           (this->mointegrals_->IAJB(i,c,k,a))*Tia(j,c)*Tia(k,b) -
                           (this->mointegrals_->IAJB(j,c,k,a))*Tia(i,c)*Tia(k,b) -
                           (this->mointegrals_->IAJB(i,c,k,b))*Tia(j,c)*Tia(k,a) +
                           (this->mointegrals_->IAJB(j,c,k,b))*Tia(i,c)*Tia(k,a); 
       }
       E5 += Wijab(i,j,a,b)*Tijab(i,j,a,b);
     } // end ijab 
 
    // Make new T amplitudes
     for(auto i = 0; i < this->nO_; i++)
     for(auto a = 0; i < this->nO_; i++) {
       Denom = (*this->singleSlater_->epsA())(i)
               - (*this->singleSlater_->epsA())(a + this->nO_);
       Tia(i,a) = Wia(i,a)/Denom;
     }
     for(auto i = 0; i < this->nO_; i++) 
     for(auto j = 0; j < this->nO_; j++) 
     for(auto a = 0; a < this->nV_; a++) 
     for(auto b = 0; b < this->nV_; b++) { 
       Denom = (*this->singleSlater_->epsA())(i)
               + (*this->singleSlater_->epsA())(j)
               - (*this->singleSlater_->epsA())(a + this->nO_)
               - (*this->singleSlater_->epsA())(b + this->nO_);
       Tijab(i,j,a,b) = (Wijab(i,j,a,b) + (this->mointegrals_->IJAB(i,j,a,b)))/Denom; 
     }

     for(auto i = 0; i < this->nO_; i++) 
     for(auto j = 0; j < this->nO_; j++) 
     for(auto a = 0; a < this->nV_; a++) 
     for(auto b = 0; b < this->nV_; b++) { 
       Tau(i,j,a,b) = Tijab(i,j,a,b) + (Tia(i,a)*Tia(j,b) -
                                        Tia(i,b)*Tia(j,a));
       ECCSD += (this->mointegrals_->IJAB(i,j,a,b))*Tau(i,j,a,b);
     }
    
     ECCSD = ECCSD*(0.25);
     cout << "iter # " << niter << "   " << ECCSD << endl;
    



    } // end CCSD iteration
  } // end if TCS
  ECorr = ECCSD;
  return ECorr;
};
