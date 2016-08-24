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
      << this->totalEnergy << " Eh ***" << endl;
  }

  this->doIncFock_ = this->isPrimary and !this->isDFT;
  this->doIncFock_ = false;

  for(iter = 0; iter < this->maxSCFIter_; iter++){
/*
    DeltaDScalar_->read(this->NBSqScratch_->data(),H5PredType<T>());
    prettyPrint(cout,*this->onePDMA_,"Den at top " +std::to_string(iter));
    prettyPrint(cout,*this->fockA_,"Fock at top " +std::to_string(iter));
    prettyPrint(cout,*this->PTA_,"PT at top " +std::to_string(iter));
    prettyPrint(cout,*this->NBSqScratch_,"Delta at top "+std::to_string(iter));
*/
    this->copyDen();
    this->copyPT();
    if(iter != 0 and this->doIncFock_)
      this->copyDeltaDtoD();

 // prettyPrint(cout,*this->onePDMA_,"Den before Fock " + std::to_string(iter));
    this->formFock(this->doIncFock_ && iter != 0);

    if(iter != 0 and this->doIncFock_){
 //   prettyPrint(cout,*this->fockA_,"Fock Before at " + std::to_string(iter));
      this->incPT();
 //   prettyPrint(cout,*this->fockA_,"Fock After at " + std::to_string(iter));
      this->copyDOldtoD();
    } else {
 //   prettyPrint(cout,*this->fockA_,"Fock at " + std::to_string(iter));
    }
 // prettyPrint(cout,*this->onePDMA_,"Den at " + std::to_string(iter));
   

/*
    doLevelShift = this->nLevelShift_ != 0;
    doLevelShift = doLevelShift and iter >= this->iStartLevelShift_;
    doLevelShift = 
      doLevelShift and iter < this->iStartLevelShift_ + 15;

    if(doLevelShift) this->levelShift2();
*/
    this->orthoFock3();
    this->formFP();

    auto IDIISIter = iter - this->iDIISStart_;
    if(this->doDIIS){
      this->cpyFockDIIS(IDIISIter);
      this->cpyDenDIIS(IDIISIter);
      this->genDIISCom(IDIISIter); 
    }

    // DIIS Extrapolation of the Fock
    if(this->doDIIS and IDIISIter > 0) 
      this->CDIIS4(std::min(IDIISIter+1,std::size_t(this->nDIISExtrap_)));
//  prettyPrint(cout,*this->fockA_,"Fock After DIIS at " + std::to_string(iter));

    if(this->doDMS){
      this->formDMSErr(IDIISIter);
    }

    if(this->doDMS and IDIISIter > 0){
      // DIIS (DMS) Extrapolation of the Density
      this->DMSExtrap(std::min(IDIISIter+1,std::size_t(this->nDIISExtrap_)));
    } else {
      // Copies over the ortho focks to the MO storage
      this->populateMO4Diag();

      // Diagonalizes the orthonormal fock matricies
      this->diagFock2();

      // This stupidly computes the orthonormal density and stores it in the
      // AO storage
      this->formDensity();

    }
    // This copy operation negates the problem above
    this->cpyAOtoOrthoDen();

    // Transform D into the AO basis
    this->unOrthoDen3();

    // Calls formFock
    SCFConvergence CONVER = this->evalConver3();

    if(this->printLevel_ > 0)
      this->printSCFIter(iter,CONVER.EDelta,CONVER.PSRMS,
        CONVER.PMRMS);

    this->nSCFIter++;

    if(this->isConverged) break;
  }

  this->cleanupSCFMem3();
  this->backTransformMOs();
  this->fixPhase();

  if(!this->isConverged)
    CErr("SCF Failed to converge within MAXITER iterations",
        this->fileio_->out);
  if(this->printLevel_ > 0 ){
    this->fileio_->out 
      << endl << "SCF Completed: E(" 
      << this->SCFTypeShort_ << ") = ";
    this->fileio_->out 
      << std::fixed << std::setprecision(10) 
      << this->totalEnergy << "  Eh after  " << iter + 1 
      << "  SCF Iterations" << endl;
    this->fileio_->out << bannerEnd <<endl;
  }
};

