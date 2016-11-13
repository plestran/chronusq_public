/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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

void CoupledCluster::MollerPlesset(){
 /*
   Compute MP2, MP3 and singles contribution to MP4
   Only works for real TCS/GHF at the moment.

   Evaluating MP2, MP3 and MP4 takes us through all
   the possible double bar integrals, so this routine
   is a useful check for AO to MO integral transformation
   correctness.
 */

  if(this->Ref_ == SingleSlater<double>::TCS){
  
  this->fileio_->out << " Moller Plesset Perturbation Theory:" << endl; 
  
    double E2    = 0.0;
    double Denom = 1.0;
  
  // First, calculate MP2 energy and MP1 coefficient
  
    RealTensor4d Aijab1;
    Aijab1 = RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_);
  
    for(auto i = 0; i < this->nO_; i++) 
    for(auto j = 0; j < this->nO_; j++) 
    for(auto a = 0; a < this->nV_; a++) 
    for(auto b = 0; b < this->nV_; b++) { 
      Denom = (*this->singleSlater_->epsA())(a + this->nO_)
               + (*this->singleSlater_->epsA())(b + this->nO_)
               - (*this->singleSlater_->epsA())(i)
               - (*this->singleSlater_->epsA())(j);
      Aijab1(i,j,a,b) = -(this->mointegrals_->IJAB(i,j,a,b))/Denom;
      E2 += (0.25)*(this->mointegrals_->IJAB(i,j,a,b))*Aijab1(i,j,a,b);
      }
  
    this->fileio_->out << "EMP2(corr) = " << std::setprecision(14) << E2 << endl;
  
  // Now try MP3
    double EMP3   = 0.0;
    double EMP3A  = 0.0;
    double EMP3B  = 0.0;
    double EMP3C  = 0.0;
    double Denom1 = 1.0;
    double Denom2 = 1.0;
  
    for(auto i = 0; i < this->nO_; i++) 
    for(auto j = 0; j < this->nO_; j++) 
    for(auto a = 0; a < this->nV_; a++) 
    for(auto b = 0; b < this->nV_; b++)  
    for(auto c = 0; c < this->nV_; c++)  
    for(auto d = 0; d < this->nV_; d++) {
      Denom1 =   (*this->singleSlater_->epsA())(a + this->nO_)
               + (*this->singleSlater_->epsA())(b + this->nO_)
               - (*this->singleSlater_->epsA())(i)
               - (*this->singleSlater_->epsA())(j);
      Denom2 =   (*this->singleSlater_->epsA())(c + this->nO_)
               + (*this->singleSlater_->epsA())(d + this->nO_)
               - (*this->singleSlater_->epsA())(i)
               - (*this->singleSlater_->epsA())(j);
      EMP3A += (this->mointegrals_->IJAB(i,j,a,b))*
               (this->mointegrals_->IJAB(i,j,c,d))*
               (this->mointegrals_->ABCD(a,b,c,d))/
               (Denom1*Denom2);
    }  
  
    for(auto i = 0; i < this->nO_; i++) 
    for(auto j = 0; j < this->nO_; j++) 
    for(auto a = 0; a < this->nV_; a++) 
    for(auto b = 0; b < this->nV_; b++)  
    for(auto k = 0; k < this->nO_; k++)  
    for(auto l = 0; l < this->nO_; l++) {
      Denom1 =   (*this->singleSlater_->epsA())(a + this->nO_)
               + (*this->singleSlater_->epsA())(b + this->nO_)
               - (*this->singleSlater_->epsA())(i)
               - (*this->singleSlater_->epsA())(j);
      Denom2 =   (*this->singleSlater_->epsA())(a + this->nO_)
               + (*this->singleSlater_->epsA())(b + this->nO_)
               - (*this->singleSlater_->epsA())(k)
               - (*this->singleSlater_->epsA())(l);
      EMP3B += (this->mointegrals_->IJAB(i,j,a,b))*
               (this->mointegrals_->IJAB(k,l,a,b))*
               (this->mointegrals_->IJKL(k,l,i,j))/
               (Denom1*Denom2);
    }  
  
    for(auto i = 0; i < this->nO_; i++) 
    for(auto j = 0; j < this->nO_; j++) 
    for(auto a = 0; a < this->nV_; a++) 
    for(auto b = 0; b < this->nV_; b++)  
    for(auto k = 0; k < this->nO_; k++)  
    for(auto c = 0; c < this->nV_; c++) {
      Denom1 =   (*this->singleSlater_->epsA())(a + this->nO_)
               + (*this->singleSlater_->epsA())(b + this->nO_)
               - (*this->singleSlater_->epsA())(i)
               - (*this->singleSlater_->epsA())(j);
      Denom2 =   (*this->singleSlater_->epsA())(c + this->nO_)
               + (*this->singleSlater_->epsA())(b + this->nO_)
               - (*this->singleSlater_->epsA())(k)
               - (*this->singleSlater_->epsA())(j);
      EMP3C -= (this->mointegrals_->IJAB(i,j,a,b))*
               (this->mointegrals_->IJAB(k,j,c,b))*
               (this->mointegrals_->IAJB(k,a,i,c))/
               (Denom1*Denom2);
    }  
    EMP3 = EMP3A/8.0 + EMP3B/8.0 + EMP3C; 
  
    this->fileio_->out << "EMP3(corr) = " << EMP3 << endl;
  
  // Now onto MP4, first compute EMP4(S)
  
    double EUMP4S = 0.0;
  
    RealTensor2d Ubbia, Aia;
    RealTensor4d Ubb2, Aijab2;
    
    Ubbia  = RealTensor2d(this->nO_,this->nV_);
    Aia    = RealTensor2d(this->nO_,this->nV_);
    Ubb2   = RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_);
    Aijab2 = RealTensor4d(this->nO_,this->nO_,this->nV_,this->nV_);
  
    for(auto i = 0; i < this->nO_; i++) 
    for(auto a = 0; a < this->nV_; a++) { 
      for(auto j = 0; j < this->nO_; j++) 
      for(auto b = 0; b < this->nV_; b++)  
      for(auto c = 0; c < this->nV_; c++) {
        Ubbia(i,a) -= (this->mointegrals_->IABC(j,a,b,c))*Aijab1(i,j,b,c)*(0.5);
      }
      for(auto j = 0; j < this->nO_; j++) 
      for(auto k = 0; k < this->nO_; k++)  
      for(auto b = 0; b < this->nV_; b++) {
        Ubbia(i,a) -= (this->mointegrals_->IJKA(j,k,i,b))*Aijab1(j,k,a,b)*(0.5);
      }
      Denom =    (*this->singleSlater_->epsA())(i)
               - (*this->singleSlater_->epsA())(a + this->nO_);
      Aia(i,a) = Ubbia(i,a)/Denom;
    }
    
    for(auto i = 0; i < this->nO_; i++) 
    for(auto a = 0; a < this->nV_; a++) { 
      EUMP4S += Aia(i,a)*Ubbia(i,a);
    }
  
    this->fileio_->out << "EMP4(S)(corr) = " << EUMP4S << endl;
  
  // Now EMP4(D)... NYI 
    double EUMP4D = 0.0;
  } // end TCS logic
}


double CoupledCluster::CCSD(){
  /*
     Formulation based off of Scuseria & Schaefer, J. Chem. Phys. 90 (7),1989 3700

     This simple routine performs just a basic iteration to get CCSD energy and
     amplitudes. The notation should be close to the original paper, but the 
     relevant equation number have been given.
  */
  int    maxIter   = 50;     // Maximum iterations in CCSD

  double ECCSD     = 0.0;    // CCSD correlation energy 
  double ECorr;              // CCSD correlation energy from old iteration.
  double EMP2      = 0.0;    // MP2 correlation energy. We get MP2 for free with CCSD.
  double deltaE;             // Change in energy from
  double threshold = 1e-10;  // Energy convergence threshold
  double Denom;              // Allocate float for denominator in expressions.

  // Intermediates
  RealTensor2d Hba, Hbj, Hij, Gca, Gik; // Eq (5,6,7,9,10) 
  RealTensor4d Aijkl, Babcd, Hicak;     // Eq (11,12,13)
  // Amplitudes 
  RealTensor2d Tia, Wia;                // T1 amplitudes. Wia is 'working' T1.
  RealTensor4d Tijab, Tau, Wijab;       // Tau = Eq(3), Rest are T2 amplitudes

  // Currently only real GHF is supported.
  if(this->Ref_ == SingleSlater<double>::TCS){

    this->fileio_->out << bannerTop << endl;
    this->fileio_->out << "Coupled Cluster Singles and Doubles (CCSD)" << endl;
    this->fileio_->out << bannerMid << endl;

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
    Tia.fill(0.0); // T1 is zero with canonical HF ref
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

    // Initialize Tau
    for(auto i = 0; i < this->nO_; i++) 
    for(auto j = 0; j < this->nO_; j++) 
    for(auto a = 0; a < this->nV_; a++) 
    for(auto b = 0; b < this->nV_; b++) { 
      Tau(i,j,a,b) = Tijab(i,j,a,b) + (Tia(i,a)*Tia(j,b) - Tia(i,b)*Tia(j,a));
      EMP2 += (0.25)*(this->mointegrals_->IJAB(i,j,a,b))*Tijab(i,j,a,b);
    }
    // Check... if we made Tau correctly we should get MP2 energy expression back
    this->fileio_->out << "E(MP2,corr) = " << EMP2 << endl; 

    // Top of CCSD amplitude iteration
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(12) << "Iteration" 
                       << std::setw(20) << "Energy, Corr (au)" 
                       << std::setw(20) << "DeltaE, Corr (au)" << endl; 
    this->fileio_->out << bannerMid << endl;

    for(auto niter = 0; niter < maxIter; niter++) {
      ECorr = ECCSD;
      // Hba, Eq(5)
      Hba.fill(0.0);
      for(auto b = 0; b < this->nV_; b++)
      for(auto a = 0; a < this->nV_; a++)
      for(auto j = 0; j < this->nO_; j++)
      for(auto k = 0; k < this->nO_; k++)
      for(auto c = 0; c < this->nV_; c++) {
        Hba(b,a) -= (this->mointegrals_->IJAB(j,k,b,c))*Tau(j,k,a,c)*(0.5);
      }

     // Hij, Eq(6)
      Hij.fill(0.0);
      for(auto i = 0; i < this->nO_; i++)
      for(auto j = 0; j < this->nO_; j++)
      for(auto k = 0; k < this->nO_; k++)
      for(auto b = 0; b < this->nV_; b++)
      for(auto c = 0; c < this->nV_; c++) {
        Hij(i,j) += (this->mointegrals_->IJAB(j,k,b,c))*Tau(i,k,b,c)*(0.5);
      }

    // Hbj, Eq(7)
      Hbj.fill(0.0);
      for(auto b = 0; b < this->nV_; b++)
      for(auto j = 0; j < this->nO_; j++)
      for(auto k = 0; k < this->nO_; k++)
      for(auto c = 0; c < this->nV_; c++) {
        Hbj(b,j) += (this->mointegrals_->IJAB(j,k,b,c))*Tia(k,c);
      }

    // Gca, Eq(9)
      Gca.fill(0.0);
      for(auto c = 0; c < this->nV_; c++) 
      for(auto a = 0; a < this->nV_; a++) {
        Gca(c,a) += Hba(c,a);
        for(auto k = 0; k < this->nO_; k++) 
        for(auto d = 0; d < this->nV_; d++) {
          Gca(c,a) -= (this->mointegrals_->IABC(k,a,c,d))*Tia(k,d);
        } 
      }

    // Gik, Eq(10)
      Gik.fill(0.0);
      for(auto i = 0; i < this->nO_; i++)
      for(auto k = 0; k < this->nO_; k++) {
        Gik(i,k) += Hij(i,k);
        for(auto l = 0; l < this->nO_; l++)
        for(auto c = 0; c < this->nV_; c++) {
          Gik(i,k) += (this->mointegrals_->IJKA(k,l,i,c))*Tia(l,c);
        }
      }

    // Aijkl, Eq(11)
      Aijkl.fill(0.0);
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

    // Babcd, Eq(12)
      Babcd.fill(0.0); 
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

   // Hicak, Eq(13) 
      Hicak.fill(0.0);
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

   // Wia, rearranged version of Eq(4)
      Wia.fill(0.0);
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
                     -(this->mointegrals_->IAJB(i,b,j,a))*Tia(j,b);
          for(auto k = 0; k < this->nO_; k++) {
            Wia(i,a) -= (this->mointegrals_->IJKA(j,k,i,b))*Tau(j,k,a,b)*(0.5);
          }
        }
      }
      
    // Wijab, rearranged version of Eq(8)
      Wijab.fill(0.0);
      for(auto i = 0; i < this->nO_; i++) 
      for(auto j = 0; j < this->nO_; j++) 
      for(auto a = 0; a < this->nV_; a++) 
      for(auto b = 0; b < this->nV_; b++) { 
        for(auto k = 0; k < this->nO_; k++)  
        for(auto l = 0; l < this->nO_; l++) {
        Wijab(i,j,a,b) += Aijkl(i,j,k,l)*Tau(k,l,a,b)*(0.5);
        }
        for(auto c = 0; c < this->nV_; c++)  
        for(auto d = 0; d < this->nV_; d++) {
          Wijab(i,j,a,b) += Babcd(a,b,c,d)*Tau(i,j,c,d)*(0.5);
        }
        for(auto c = 0; c < this->nV_; c++) {
          Wijab(i,j,a,b) += Gca(c,a)*Tijab(i,j,c,b) -
                            Gca(c,b)*Tijab(i,j,c,a) -
                            (this->mointegrals_->IABC(j,c,a,b))*Tia(i,c) +
                            (this->mointegrals_->IABC(i,c,a,b))*Tia(j,c); 
        }
        for(auto k = 0; k < this->nO_; k++) {
          Wijab(i,j,a,b) += (this->mointegrals_->IJKA(i,j,k,a))*Tia(k,b) -
                            (this->mointegrals_->IJKA(i,j,k,b))*Tia(k,a) -
                            Gik(i,k)*Tijab(k,j,a,b) +
                            Gik(j,k)*Tijab(k,i,a,b); 
        }
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
      }  

      // Make new T amplitudes
      // T1
      Tia.fill(0.0);
      for(auto i = 0; i < this->nO_; i++)
      for(auto a = 0; a < this->nV_; a++) {
        Denom = (*this->singleSlater_->epsA())(i)
                - (*this->singleSlater_->epsA())(a + this->nO_);
        Tia(i,a) = Wia(i,a)/Denom;
      }

      // T2
      Tijab.fill(0.0);
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


      // Update Tau and evaluate E(corr) for CCSD
      ECCSD = 0.0;
      for(auto i = 0; i < this->nO_; i++) 
      for(auto j = 0; j < this->nO_; j++) 
      for(auto a = 0; a < this->nV_; a++) 
      for(auto b = 0; b < this->nV_; b++) { 
        Tau(i,j,a,b) = Tijab(i,j,a,b) + (Tia(i,a)*Tia(j,b) -
                                         Tia(i,b)*Tia(j,a));
        ECCSD += (this->mointegrals_->IJAB(i,j,a,b))*Tau(i,j,a,b);
      }

      // GHF reference, scale by factor of 1/4, evaluate change in energy
      ECCSD = ECCSD*(0.25);
      deltaE = ECCSD - ECorr; 

      // Print to stdout
      this->fileio_->out << " Iter #:" <<  std::setw(4) << niter + 1 
                         << std::setw(20) << ECCSD 
                         << std::setw(20) << deltaE << endl;

      // threshold, right now we do energy
      if (std::abs(deltaE) < threshold) break;

      } // end CCSD iteration
    this->fileio_->out << bannerEnd << endl;
   } // end if TCS
   ECorr = ECCSD;
   return ECorr;
};
