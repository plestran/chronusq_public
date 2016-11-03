/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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

//#include <singleslater/singleslater_oldscf.h>
#include <singleslater/singleslater_levelshift.h>
#include <singleslater/singleslater_scfutils.h>



template<typename T>
void SingleSlater<T>::SCF3(){
  // Compute the 1-Body Integrals if they haven't already been 
  // computed
  if(!this->aointegrals_->haveAOOneE) 
    this->aointegrals_->computeAOOneE();

  this->initSCFMem3();

  // Print the SCF Header
  if(this->printLevel_ > 0)
    this->printSCFHeader(this->fileio_->out);
  
  this->computeEnergy();
  this->isConverged = false;
  std::size_t iter;
  bool doLevelShift;
  if(this->printLevel_ > 0){
    this->fileio_->out << "    *** INITIAL GUESS ENERGY = " 
      << this->totalEnergy_ << " Eh ***" << endl;
    
    double HOMO_LUMO_GAP = 
      (*this->epsA_)(this->nO_) - (*this->epsA_)(this->nO_-1);

    if(HOMO_LUMO_GAP < 1e-10)
      this->fileio_->out << "    *** WARNING: DEGENERATE HOMO-LUMO GAP ***"
                         << endl;
  }

  this->doIncFock_ = this->isPrimary and !this->isDFT;
  this->doIncFock_ = false;

  for(iter = 0; iter < this->maxSCFIter_; iter++){
    /*
    this->computeSExpect(*this->aointegrals_->overlap_);
    this->printSExpect();
    */

    // Write the current matricies to disk: 
    //   F(k-1), D(k), C(k),
    //   FO(k-1), DO(k)
    this->writeSCFFiles();

//  if(iter != 0 and this->doIncFock_)
//    this->copyDeltaDtoD();

    // Form the AO Fock Matrix: F(k)[D(k)]
    this->formFock(this->doIncFock_ && iter != 0);

   
//  if(iter != 0 and this->doIncFock_){
//    this->incPT();
//    this->copyDOldtoD();
//  } 
   

    // Orthonormalize the Fock: F(k) -> FO(k)
    this->orthoFock3();
    // Form the orthonormal product of fock and density:
    //   FPO = FO(k) * DO(k), FPO -> FP
    this->formFP();

    // Copy over the necessary moieties for DIIS Extrapolation:
    //   FO(k), DO(k), FP - PF
    auto IDIISIter = iter - this->iDIISStart_;
    if(this->doDIIS){
      this->cpyFockDIIS(IDIISIter);
      this->cpyDenDIIS(IDIISIter);
      this->genDIISCom(IDIISIter); 
    }

    // DIIS Extrapolation of the Fock
    //   Extrapolates the AO Fock and then transforms into orthonormal
    //   basis, i.e. both F(k) and FO(k) are now overwritten with the
    //   DIIS extrapolated quantities.
    if(this->doDIIS and IDIISIter > 0 and this->diisAlg_ != NO_DIIS) { 
      if(this->diisAlg_ == CDIIS)
        this->CDIIS4(std::min(IDIISIter+1,std::size_t(this->nDIISExtrap_)));
    }

    // Damping
    if(this->doDamp){
      this->fockDamping();
    }
/*
    // Level shifting
    if(this->doLevelShift){
      this->levelShift2();
    }
*/

/*
    doLevelShift = this->nLevelShift_ != 0;
    doLevelShift = doLevelShift and iter >= this->iStartLevelShift_;
    doLevelShift = 
      doLevelShift and iter < this->iStartLevelShift_ + 15;

    if(doLevelShift) this->levelShift2();
*/

    if(this->doDMS){
      this->formDMSErr(IDIISIter);
    }

    if(this->doDMS and IDIISIter > 0){
      // DIIS (DMS) Extrapolation of the Density
      this->DMSExtrap(std::min(IDIISIter+1,std::size_t(this->nDIISExtrap_)));
    } else {

      // Copies over the ortho focks to the MO storage
      //   C(k) = FO(k)
      this->populateMO4Diag();

      
      // Diagonalizes the orthonormal fock matricies
      // FO(k) {stored in C} -> C(k+1) {in orthonormal basis}
      this->diagFock2();

      // This stupidly computes the orthonormal density and stores it in the
      // AO storage
      //   DO(k+1) = C(k+1) * C(k+1)**H {stored in D storage}
      this->formDensity();

    }

    // This copy operation negates the mispopulation of AO storage
    //   DO(k+1) = D(k+1)
    this->cpyAOtoOrthoDen();

    // Transform D into the AO basis
    //   DO(k+1) -> D(k+1)
    this->unOrthoDen3();
    if(this->printLevel_ > 3) this->printDensity();

    // The generation of the commutator for DIIS also checks for convergece
    // i.e. if [F,P] is converged, this is a necessary and sufficient 
    // condition for SCF convergence and will cause problems in the DIIS
    // extrapolation if not accounted for (ill posed DIIS problem)
    if(this->isConverged)
      this->fileio_->out << "   *** SCF Converged by [F,P] ***" << endl;

    // Evaluate standard convergence critera (RMS density diffence,
    //   energy change, etc)
    SCFConvergence CONVER = this->evalConver3();

    // Print line to output file regarding the SCF iteration
    if(this->printLevel_ > 0)
      this->printSCFIter(iter,CONVER.EDelta,CONVER.PSRMS,
        CONVER.PMRMS);

    this->nSCFIter++;

    // If convergence is reached, break out of SCF loop
    if(this->isConverged) break;
  }

  // Deallocate all of the SCF temporaries
  this->cleanupSCFMem3();

  // Transform the MO coefficients into the AO basis as all of the
  //   arithmetic was done in the orthonormal basis
  //   C = X*C
  this->backTransformMOs();

  // Force a certain phase critera on the MO coefficients (which
  //   have gauge freedom), s.t. the "largest" element of each MO
  //   coefficient will be real and positive. 
  this->fixPhase();

  // Print out information regarding SCF not converging
  if(!this->isConverged)
    CErr("SCF Failed to converge within MAXITER iterations",
        this->fileio_->out);


  // Print out SCF energy results...
  if(this->printLevel_ > 0 ){
    this->fileio_->out 
      << endl << "SCF Completed: E(" 
      << this->SCFTypeShort_ << ") = ";
    this->fileio_->out 
      << std::fixed << std::setprecision(10) 
      << this->totalEnergy_ << "  Eh after  " << iter + 1 
      << "  SCF Iterations" << endl;
    this->fileio_->out << bannerEnd <<endl;
  }
};

