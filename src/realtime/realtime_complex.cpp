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
#include <realtime.h>
#include <aointegrals.h>
using ChronusQ::AOIntegrals;
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::FileIO;
using ChronusQ::SingleSlater;
using ChronusQ::RealTime;

namespace ChronusQ {

template<>
void RealTime<dcomplex>::iniDensity() {
  bool inOrthoBas;
  bool idempotent;

// Form the orthonormal transformation matrices
  if (this->typeOrtho_ == 1) {  
   // Lowdin transformation 
   // V1 = S^(-1/2)
   // V2 = S^(1/2)
    (*this->oTrans1_).real() = (*this->aointegrals_->overlap_).pow(-0.5);
    (*this->oTrans2_).real() = (*this->aointegrals_->overlap_).pow(0.5);
    if(this->controls_->printLevel>3) {
      prettyPrint(this->fileio_->out,(*this->oTrans1_),"S^(-1/2)");
      prettyPrint(this->fileio_->out,(*this->oTrans2_),"S^(1/2)");
    }
  }
  else if (this->typeOrtho_ == 2) {  
  // Cholesky transformation
    CErr("Cholesky orthogonalization NYI",this->fileio_->out);
  }
  else if (this->typeOrtho_ == 3) {  	
    CErr("Canonical orthogonalization NYI",this->fileio_->out);
  // Canonical orthogonalization
  // V1 = U*s^(-1/2)
  // V2 = S*V1
  }

// Form the initial density
  if (this->initDensity_ == 0) { 
// Use converged ground-state density
    inOrthoBas = false;
    idempotent = true;
  }
  else if (this->initDensity_ == 1) { 
// Form the initial density by swaping MOs
    inOrthoBas = false;
    idempotent = true;
    if (this->swapMOA_ != 0) {
      int iA = ((this->swapMOA_)/1000), jA = ((this->swapMOA_)%1000); // MOs to swap
      this->fileio_->out<<"\nAlpha MOs swapped: "<<iA<<" <-> "<<jA<<endl;
      if(this->controls_->printLevel>3) {
        prettyPrint(this->fileio_->out,(*this->ssPropagator_->moA()),"Initial Alpha MO");
      }
      this->ssPropagator_->moA()->col(jA-1).swap(this->ssPropagator_->moA()->col(iA-1));
    }
    if (this->swapMOB_ != 0 && !this->isClosedShell_) {
      int iB = (this->swapMOB_/1000), jB = (this->swapMOB_%1000); // MOs to swap
      this->fileio_->out<<"\nBeta MOs swapped: "<<iB<<" <-> "<<jB<<endl;
      if(this->controls_->printLevel>3) {
        prettyPrint(this->fileio_->out,(*this->ssPropagator_->moB()),"Initial Beta MO");
      }
      this->ssPropagator_->moB()->col(jB-1).swap(this->ssPropagator_->moB()->col(iB-1));
    }
    this->ssPropagator_->formDensity();

    // Transform the ground state MO to orthonormal basis
    this->initMOA_->setZero();
    *this->initMOA_= *this->groundState_->moA();
    *this->initMOA_ = (*this->oTrans2_)*(*this->initMOA_);
    if(!this->isClosedShell_) {
      this->initMOB_->setZero();
      *this->initMOB_ = *this->groundState_->moB();
      *this->initMOB_ = (*this->oTrans2_)*(*this->initMOB_);
    }
  }
  else if (this->initDensity_ == 2) { 
// Read in the AO density from checkpoint file
    CErr("Read in the AO density from checkpint file NYI",this->fileio_->out);
  }
  else if (this->initDensity_ == 3) { 
// Read in the orthonormal density from checkpoint file
    CErr("Read in the orthonormal density from checkpoint file NYI",this->fileio_->out);
  }

  if (!inOrthoBas) { 
// Transform density from AO to orthonormal basis
    *this->POA_    = (*this->oTrans2_)*(*this->ssPropagator_->densityA())*(*this->oTrans2_);
    *this->POAsav_ = *this->POA_;
    if(!this->isClosedShell_) {
      *this->POB_    = (*this->oTrans2_)*(*this->ssPropagator_->densityB())*(*this->oTrans2_);
      *this->POBsav_ = *this->POB_;
    }
  }
  else { 
// Transform density from orthonormal to AO basis
    *this->ssPropagator_->densityA() = (*this->oTrans1_)*(*this->POAsav_)*(*this->oTrans1_);
    if(!this->isClosedShell_) *this->ssPropagator_->densityB() = (*this->oTrans1_)*(*this->POB_)*(*this->oTrans1_);
  }
};

template<>
void RealTime<dcomplex>::formEDField() {  
  int IEnvlp   = this->controls_->rtEnvelope_;
  double Ex    = this->controls_->rtField_[0]; 
  double Ey    = this->controls_->rtField_[1]; 
  double Ez    = this->controls_->rtField_[2]; 
  double TOn   = (this->controls_->rtTOn_)/phys.AuToFs;
  double TOff  = (this->controls_->rtTOff_)/phys.AuToFs;
  double Omega = this->controls_->rtFreq_/phys.eVPerHartree;
  double Phase = this->controls_->rtPhase_;
  double Time  = currentTime_;
  double OmegT;
  double TMax;

  if (IEnvlp == 1) { 
  //   Constant envelope (plane wave)
    if(Time >= TOn && Time <= TOff) {
      OmegT = Omega*(Time-TOn) + Phase;
      (*this->EDField_)[0] = Ex*std::cos(OmegT);
      (*this->EDField_)[1] = Ey*std::cos(OmegT);
      (*this->EDField_)[2] = Ez*std::cos(OmegT);
    } else {
      (*this->EDField_)[0] = 0.0;
      (*this->EDField_)[1] = 0.0;
      (*this->EDField_)[2] = 0.0;
    }
  } 
  else if (IEnvlp == 2) { 
  /*   
       Linearly ramping up to the maximum in the first cycle
       Constant envelope afterwards
       Linearly ramping off to zero
  */
    if(Time >= TOn && Time <= TOff) {
      OmegT = Omega*(Time-TOn) + Phase;
    } else {
      (*this->EDField_)[0] = 0.0;
      (*this->EDField_)[1] = 0.0;
      (*this->EDField_)[2] = 0.0;
    }

  }
  else if (IEnvlp == 3) { 
  /*
       Gaussian envelope: E(t) = E * exp(-(a*(t-t0))^2) * Sin(wt)
           a  = the range of frequency (FWHM) 
           t0 = the time when the amplitude reaches maximum 
                  the default for t0 = sqrt(LN(1000))/a 
                  (at t=0, the amplitude is 1/1000 times maximum. 
                  This ensures a smooth turning-on of the field)             
        Dipole length approximation only
  */
    if(Time >= TOn && Time <= TOff) {
      OmegT = Omega*(Time-TOn) + Phase;
    } else {
      (*this->EDField_)[0] = 0.0;
      (*this->EDField_)[1] = 0.0;
      (*this->EDField_)[2] = 0.0;
    }

  }
  else if (IEnvlp == 4) { 
  /*
        Step function
  */
    if(Time >= TOn && Time <= TOff) {
      (*this->EDField_)[0] = Ex;
      (*this->EDField_)[1] = Ey;
      (*this->EDField_)[2] = Ez;
    } else {
      (*this->EDField_)[0] = 0.0;
      (*this->EDField_)[1] = 0.0;
      (*this->EDField_)[2] = 0.0;
    }

  }
  else if (IEnvlp == 5) { 
  /*
        Sine-square envelope: 
          A(r,t) = (E0/w)*(sin(Pi*t/T))^2*sin(k*r-w*t) 
          E(r,t) = E0*(sin(Pi*t/T))^2*cos(k*r-w*t)
                    -(E0*Pi/w*T)*sin(2*Pi*t/T)*sin(k*r-w*t)
        only dipole length/velocity and quardrupole velocity
  */
    if(Time >= TOn && Time <= TOff) {
      OmegT = Omega*(Time-TOn) + Phase;
    } else {
      (*this->EDField_)[0] = 0.0;
      (*this->EDField_)[1] = 0.0;
      (*this->EDField_)[2] = 0.0;
    }
  }
};

template<>
void RealTime<dcomplex>::printRT() {
  double EDx = (*this->ssPropagator_->dipole())(0)/phys.debye;
  double EDy = (*this->ssPropagator_->dipole())(1)/phys.debye;
  double EDz = (*this->ssPropagator_->dipole())(2)/phys.debye;

  if (currentTime_ <= 0.0) {
    this->fileio_->out << bannerMid << endl;
    this->fileio_->out << std::setw(11) << std::right << "Time" 
                       << std::setw(1)  << " "
                       << std::setw(16) << std::right << "Energy"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << "Ex Dipole"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " Ey Dipole" 
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " Ez Dipole"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " Etot Dipole"
                       << std::endl;
    this->fileio_->out << std::setw(11) << std::right << "    (au)" 
                       << std::setw(1)  << " "
                       << std::setw(16) << std::right << "    (Eh)"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " (debye)"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " (debye)" 
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " (debye)"
                       << std::setw(1)  << " "
                       << std::setw(12) << std::right << " (debye)"
                       << std::endl;
    this->fileio_->out << bannerMid << endl;
  }
  this->fileio_->out << std::setw(11) << std::right << std::setprecision(4) << currentTime_
                     << std::setw(1)  << " "
                     << std::setw(16) << std::right << std::setprecision(10) << this->ssPropagator_->totalEnergy
                     << std::setw(1)  << " "
                     << std::setw(12) << std::right << std::setprecision(6) << EDx 
                     << std::setw(1)  << " "
                     << std::setw(12) << std::right << std::setprecision(6) << EDy
                     << std::setw(1)  << " "
                     << std::setw(12) << std::right << std::setprecision(6) << EDz 
                     << std::setw(1)  << " "
                     << std::setw(12) << std::right << std::setprecision(6) <<
                        std::sqrt( std::pow(EDx,2.0) +
                                   std::pow(EDy,2.0) +
                                   std::pow(EDz,2.0))
                     << std::endl;
};


template<>
void RealTime<dcomplex>::formUTrans() {  
//
// Form the unitary transformation matrix: 
// U = exp(-i*dT*F)
//
  if (this->methFormU_ == 1) { 
   //  Eigen-decomposition
    ComplexMatrix EVec(this->nBasis_,this->nBasis_);
    RealMatrix 	  EVal(this->nBasis_,1);

    Eigen::SelfAdjointEigenSolver<ComplexMatrix> sys(*this->ssPropagator_->fockA());
    EVec = sys.eigenvectors();
    EVal = sys.eigenvalues();
    this->uTransA_->setZero();
    for (int i=0;i<this->nBasis_;i++) {
      (*this->uTransA_)(i,i) = dcomplex(cos(deltaT_*EVal(i,0)),-sin(deltaT_*EVal(i,0)));
    }
    *this->uTransA_ = EVec*(*this->uTransA_)*(EVec.adjoint());
    if(!this->isClosedShell_) {
      Eigen::SelfAdjointEigenSolver<ComplexMatrix> sys(*this->ssPropagator_->fockB());
      EVec = sys.eigenvectors();
      EVal = sys.eigenvalues();
      this->uTransB_->setZero();
      for (int i=0;i<this->nBasis_;i++) {
        (*this->uTransB_)(i,i) = dcomplex(cos(deltaT_*EVal(i,0)),-sin(deltaT_*EVal(i,0)));
      }
    *this->uTransB_ = EVec*(*this->uTransB_)*(EVec.adjoint());
    }
  }
  else if (this->methFormU_ == 2) { 
  // Taylor expansion
    *this->scratch_ = -math.ii*deltaT_*(*this->ssPropagator_->fockA());
    *this->uTransA_ = (*this->scratch_).exp();
    if(!this->isClosedShell_) {
      *this->scratch_ = -math.ii*deltaT_*(*this->ssPropagator_->fockB());
      *this->uTransB_ = (*this->scratch_).exp();
    }
  }
//    prettyPrint(this->fileio_->out,(*this->uTransA_),"uTransA");
//    if(!this->isClosedShell_) prettyPrint(this->fileio_->out,(*this->uTransB_),"uTransB");
};
  
template<>
void RealTime<dcomplex>::doPropagation() {
  long int iStep;
  bool checkFP = false;

  currentTime_ = 0.0;
  for (iStep=0; iStep<=this->maxSteps_; iStep++) {
    //this->fileio_->out<<"\nStep "<<iStep<<":\n"<<endl;
    //this->fileio_->out<<std::right<<std::setw(20)<<"Time = "<<std::setw(15)<<currentTime_<<std::setw(5)<<" a.u. "<<endl;

    if (iStep==0) deltaT_ = this->stepSize_;
    else deltaT_ = 2.0*(this->stepSize_);

    *this->scratch_ = *this->POA_;
    *this->POA_	    = *this->POAsav_;
    *this->POAsav_  = *this->scratch_;
    if (!this->isClosedShell_) {
      *this->scratch_ = *this->POB_;
      *this->POB_     = *this->POBsav_;
      *this->POBsav_  = *this->scratch_;
    }
//  Print 
    if(this->controls_->printLevel>=1) {
      //prettyPrintComplex(this->fileio_->out,(*this->ssPropagator_->densityA()),"Alpha AO Density");
      //if(!this->isClosedShell_) prettyPrintComplex(this->fileio_->out,(*this->ssPropagator_->densityB()),"Beta AO Density");

//  Form AO Fock matrix
    this->formEDField();
    this->ssPropagator_->setField(*this->EDField_);
    this->ssPropagator_->formFock();
    this->ssPropagator_->computeEnergy();
    this->ssPropagator_->computeMultipole();
    this->printRT();


//  Transform Fock from AO to orthonormal basis
    *this->ssPropagator_->fockA() = (*this->oTrans1_)*(*this->ssPropagator_->fockA())*(*this->oTrans1_);
    if (!this->isClosedShell_) {
      *this->ssPropagator_->fockB() = (*this->oTrans1_)*(*this->ssPropagator_->fockB())*(*this->oTrans1_);
    }

//  Form the unitary propagation matrix
    this->formUTrans();

    if ((this->initDensity_==0) && checkFP) {
//    Check [F,P] for converged density, should equal to zero
      *this->scratch_ = (*this->ssPropagator_->fockA())*(*this->POAsav_) - (*this->POAsav_)*(*this->ssPropagator_->fockA());
      prettyPrint(this->fileio_->out,(*this->scratch_),"[FOA,POA]");
      if (!this->isClosedShell_) {
        *this->scratch_ = (*this->ssPropagator_->fockB())*(*this->POBsav_) - (*this->POBsav_)*(*this->ssPropagator_->fockB());
        prettyPrint(this->fileio_->out,(*this->scratch_),"[FOB,POB]");
      }
    }

//  Propagate the density matrix
    *this->POA_ = (*this->uTransA_)*(*this->POA_)*((*this->uTransA_).adjoint());
    if (!this->isClosedShell_) {
      *this->POB_ = (*this->uTransB_)*(*this->POB_)*((*this->uTransB_).adjoint());
    }

//  Transform density matrix from orthonormal to AO basis
    *this->ssPropagator_->densityA() = (*this->oTrans1_)*(*this->POA_)*(*this->oTrans1_);
    if (!this->isClosedShell_) {
      *this->ssPropagator_->densityB() = (*this->oTrans1_)*(*this->POB_)*(*this->oTrans1_);
    }

//  Advance step
    currentTime_ += this->stepSize_;
    };
  }
};

} // namespace ChronusQ
