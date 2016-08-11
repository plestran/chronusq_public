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
#include <singleslater.h>
#include<dft.h>
namespace ChronusQ {

template<>
void SingleSlater<dcomplex>::formVXC_new(){;};

template<>
void SingleSlater<double>::formVXC_new(){
//Timing
//  this->screenVxc = false;

  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point finish;
  std::chrono::duration<double> duration_formVxc;
  std::chrono::duration<double> T1(0.0);
  std::chrono::duration<double> T2(0.0);
  std::chrono::duration<double> T3(0.0);
  std::chrono::duration<double> T4(0.0);
  std::chrono::duration<double> T5(0.0);
  std::chrono::duration<double> T6(0.0);
  
  std::vector<std::chrono::duration<double>> 
    TF(this->dftFunctionals_.size(),std::chrono::duration<double>(0.0));
/*
  if(this->printLevel_ >= 3) {
    start = std::chrono::high_resolution_clock::now();
  }
*/
//  bool isGGA   = true;
//  bool isGGA   = false;
  RealMatrix SCRATCH2(this->nBasis_,this->nBasis_);
  VectorXd   SCRATCH1(this->nBasis_);
  RealMatrix SCRATCH2X(this->nBasis_,this->nBasis_);
  RealMatrix SCRATCH2Y(this->nBasis_,this->nBasis_);
  RealMatrix SCRATCH2Z(this->nBasis_,this->nBasis_);
  VectorXd   SCRATCH1X(this->nBasis_);
  VectorXd   SCRATCH1Y(this->nBasis_);
  VectorXd   SCRATCH1Z(this->nBasis_);
  int NDer = 0;
  if(this->isGGA) NDer = 1; 

  std::array<double,3>  drhoT = {0.0,0.0,0.0}; ///< array TOTAL density gradient components
  std::array<double,3>  drhoS = {0.0,0.0,0.0}; ///< array SPIN  density gradient components
  std::array<double,3>  drhoA = {0.0,0.0,0.0}; ///< array ALPHA  density gradient components
  std::array<double,3>  drhoB = {0.0,0.0,0.0}; ///< array BETA  density gradient components
  RealVecMap GradRhoT(&drhoT[0],3);
  RealVecMap GradRhoS(&drhoS[0],3);
  RealVecMap GradRhoA(&drhoA[0],3);
  RealVecMap GradRhoB(&drhoB[0],3);

  double rhoA;
  double rhoB;
  double gammaAA = 0.0;
  double gammaBB = 0.0;
  double gammaAB = 0.0;
  std::vector<bool> shMap(this->basisset_->nShell()+1);
  int NSkip2(0);
  int NSkip3(0);
  int NSkip4(0);
  int NSkip5(0);
  bool doTimings(false);

  auto Newstart = std::chrono::high_resolution_clock::now();
  auto Newend = std::chrono::high_resolution_clock::now();

  std::vector<std::size_t> closeShells;
  VectorXd OmegaA(this->nBasis_), OmegaB(this->nBasis_);
  VectorXd DENCOL(this->nBasis_);
  double fact = 2.0;
  double *SCRATCH1DATA = SCRATCH1.data();
  double *SCRATCH1XDATA = SCRATCH1X.data();
  double *SCRATCH1YDATA = SCRATCH1Y.data();
  double *SCRATCH1ZDATA = SCRATCH1Z.data();
  double *OmegaADATA = OmegaA.data();
  double *OmegaBDATA = OmegaB.data();

  std::vector<std::size_t> shSizes;
  for(auto iSh = 0; iSh < this->basisset_->nShell(); iSh++)
    shSizes.push_back(this->basisset_->shells(iSh).size());

// If screen in not on loop over all shell for each atoms, no need to
// figure which are the close shell to each point
  if (!this->screenVxc){
    closeShells.resize(this->basisset_->nShell());
    std::iota(closeShells.begin(),closeShells.end(),0);
    cout << "MAch. Prec " << std::numeric_limits<double>::epsilon() << endl;
    for(auto iShell : closeShells) {
       cout << "iShell " << iShell << " closeS[] " << closeShells[iShell] << endl; 
    }
  }
  auto valVxc = [&](std::size_t iAtm, ChronusQ::IntegrationPoint &pt, 
  KernelIntegrand<double> &result) -> void {


    cartGP &GP = pt.pt;


    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();
    // For each new sphere, determine list of close basis shells (if screen ON)
    if(pt.I == 0 && this->screenVxc) {
      closeShells.clear();
      double DX = (*this->molecule_->cart())(0,iAtm) - bg::get<0>(GP);
      double DY = (*this->molecule_->cart())(1,iAtm) - bg::get<1>(GP);
      double DZ = (*this->molecule_->cart())(2,iAtm) - bg::get<2>(GP);

      double R = DX*DX + DY*DY + DZ*DZ;
      R = std::sqrt(R);

      for(auto jAtm = 0; jAtm < this->molecule_->nAtoms(); jAtm++){
        double RAB = (*this->molecule_->rIJ())(iAtm,jAtm);
        double DIST = std::abs(RAB - R);
        for(auto iSh = 0; iSh < this->basisset_->nShell(); iSh++){
          if(this->basisset_->mapSh2Cen(iSh)-1 != jAtm) continue;
          if(this->basisset_->radCutSh()[iSh] > DIST)
            closeShells.push_back(iSh);
        }
      }
    }
    if(doTimings){ 
      Newend = std::chrono::high_resolution_clock::now();
      T1 += Newend - Newstart;
    }


    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();

    // Compute value of "close" basis functions at the given point
    for(auto iShell : closeShells) {
      int shSize= shSizes[iShell];
      int iSt = this->basisset_->mapSh2Bf(iShell);
      double * buff = this->basisset_->basisDEval(NDer,this->basisset_->shells(iShell),&pt.pt);

      std::memcpy(SCRATCH1DATA + iSt,buff,shSize*sizeof(double));
      if(NDer>0){
        double * ds1EvalX = buff + shSize;
        double * ds1EvalY = ds1EvalX + shSize;
        double * ds1EvalZ = ds1EvalY + shSize;
        std::memcpy(SCRATCH1XDATA + iSt,ds1EvalX,shSize*sizeof(double));
        std::memcpy(SCRATCH1YDATA + iSt,ds1EvalY,shSize*sizeof(double));
        std::memcpy(SCRATCH1ZDATA + iSt,ds1EvalZ,shSize*sizeof(double));
      }
    }
    if(doTimings){
      Newend = std::chrono::high_resolution_clock::now();
      T2 += Newend - Newstart;
    }

    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();

    // Determine if we computed Zeros
    double S1Norm = SCRATCH1.norm();
    double S1XNorm = SCRATCH1X.norm();
    double S1YNorm = SCRATCH1Y.norm();
    double S1ZNorm = SCRATCH1Z.norm();

    if(doTimings){
      Newend = std::chrono::high_resolution_clock::now();
      T3 += Newend - Newstart;
    }

    if(this->screenVxc && S1Norm < this->epsScreen) {NSkip4++; return;}
//A    else if(!this->screenVxc && S1Norm < this->epsScreen) {cout << "screenOFF 0" <<endl;}


    double rhoT(0.0);
    double rhoS(0.0);
    double Tt(0.0), Ts(0.0);
    double Pt(0.0), Ps(0.0);
    GradRhoT.setZero();
    GradRhoS.setZero();

    // Evaluate density and optionally density gradient at the given point
      
    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();
    double * DENT, *DENS;
    // Loop over close shells "I"
    for(auto iShell : closeShells) {
      int iSz= shSizes[iShell];
      int iSt = this->basisset_->mapSh2Bf(iShell);
      for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
        Tt = 0.0; Ts = 0.0;
        if(this->nTCS_ == 1 && this->isClosedShell){
          DENT = this->onePDMA_->data() + iBf*this->nBasis_;
        } else {
          DENT = this->onePDMScalar_->data() + iBf*this->nBasis_;
          DENS = this->onePDMMz_->data() + iBf*this->nBasis_;
        }

        // Loop over close shells "J"
        for(auto jShell : closeShells) {
          int jSz= shSizes[jShell];
          int jSt = this->basisset_->mapSh2Bf(jShell);
          for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
            Pt = fact * DENT[jBf];
            if (this->screenVxc) {
              if(std::abs(Pt) > this->epsScreen){
                Tt += Pt * SCRATCH1DATA[jBf];
              } 
            } else {
//               cout << "SCREEN OFF 2" <<endl;
               Tt += Pt * SCRATCH1DATA[jBf];
            } //Screening
          } // jBf
        } // jShell      

        if(this->nTCS_ == 1 && !this->isClosedShell){
          for(auto jShell : closeShells) {
            int jSz= shSizes[jShell];
            int jSt = this->basisset_->mapSh2Bf(jShell);
            for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
              Ps = fact * DENS[jBf];
              if (this->screenVxc) {
                if( std::abs(Ps) > this->epsScreen){
                  Ts += Ps * SCRATCH1DATA[jBf];
                }
              } else {
//                cout << "SCREEN OFF 3" <<endl;
                Ts += Ps * SCRATCH1DATA[jBf];
              }//Screening
            } // jBf
          } // jShell      
        } //UKS
        rhoT += Tt * SCRATCH1DATA[iBf];
        if(this->nTCS_ == 1 && !this->isClosedShell){
          rhoS += Ts * SCRATCH1DATA[iBf];
        }  //UKS
        if(NDer > 0) {
          drhoT[0] += Tt * SCRATCH1XDATA[iBf];
          drhoT[1] += Tt * SCRATCH1YDATA[iBf];
          drhoT[2] += Tt * SCRATCH1ZDATA[iBf];
          if(this->nTCS_ == 1 && !this->isClosedShell){
            drhoS[0] += Ts * SCRATCH1XDATA[iBf];
            drhoS[1] += Ts * SCRATCH1YDATA[iBf];
            drhoS[2] += Ts * SCRATCH1ZDATA[iBf];
          } //UKS
        } //GGA
      } // iBf
    } // iShell
    if(doTimings){
      Newend = std::chrono::high_resolution_clock::now();
      T4 += Newend - Newstart;
    }

    rhoT *= 0.5;
    rhoS *= 0.5;
    rhoA = 0.5 * (rhoT + rhoS);
    rhoB = 0.5 * (rhoT - rhoS);
    if( NDer > 0 ){
      GradRhoA.noalias() = 0.5 * (GradRhoT + GradRhoS);
      GradRhoB.noalias() = 0.5 * (GradRhoT - GradRhoS);
      gammaAA = GradRhoA.dot(GradRhoA);           
      gammaBB = GradRhoB.dot(GradRhoB);           
      gammaAB = GradRhoA.dot(GradRhoB);           
    }
    
//these if statements prevent numerical instability with zero guesses
    if ((std::abs(rhoT)) <= std::numeric_limits<double>::epsilon() ){return;}
    if (rhoT    <= 0.0 ) {CErr("Numerical noise in the density");}
    if(this->screenVxc && rhoT < this->epsScreen){return;} 

    // Evaluate density functional
    DFTFunctional::DFTInfo kernelXC;
    for(auto i = 0; i < this->dftFunctionals_.size(); i++){
      if(doTimings) 
        Newstart = std::chrono::high_resolution_clock::now();
//      if (NDer > 0) {
      kernelXC += this->dftFunctionals_[i]->eval(
          rhoA,rhoB,gammaAA,gammaAB,gammaBB);
//      } else {
//        kernelXC += this->dftFunctionals_[i]->eval(rhoA, rhoB);
//      }
      if(doTimings){
        Newend = std::chrono::high_resolution_clock::now();
        TF[i] += Newend - Newstart;
      }
    } // loop over kernels


    // Evaluate Functional derivatives

    if(doTimings) Newstart = std::chrono::high_resolution_clock::now();
      if(NDer > 0 ) {
        for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
          double GA(GradRhoA(iXYZ)), GB(GradRhoB(iXYZ));
            GradRhoA(iXYZ) = pt.weight * 
              ( 2.0 * GA * kernelXC.ddgammaAA + GB * kernelXC.ddgammaAB);
            GradRhoB(iXYZ) = pt.weight * 
              ( 2.0 * GB * kernelXC.ddgammaBB + GA * kernelXC.ddgammaAB);
        } // Grad
        OmegaA.setZero(); 
        for(auto iShell : closeShells) {
          int iSz= shSizes[iShell];
          int iSt = this->basisset_->mapSh2Bf(iShell);
       
          for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
            OmegaADATA[iBf] += GradRhoA(0) * SCRATCH1XDATA[iBf] +
                               GradRhoA(1) * SCRATCH1YDATA[iBf] +
                               GradRhoA(2) * SCRATCH1ZDATA[iBf];
          } // iBf
        } // iShell
        if(this->nTCS_ == 1 && !this->isClosedShell){
          OmegaB.setZero(); 
          for(auto iShell : closeShells) {
             int iSz= shSizes[iShell];
             int iSt = this->basisset_->mapSh2Bf(iShell);
             for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
                OmegaBDATA[iBf] += GradRhoB(0) * SCRATCH1XDATA[iBf] + 
                                   GradRhoB(1) * SCRATCH1YDATA[iBf] + 
                                   GradRhoB(2) * SCRATCH1ZDATA[iBf];  
             } // iBf
          } // iShell
        } //UKS
      } //GGA
    result.Energy += pt.weight * kernelXC.eps;
    kernelXC.ddrhoA *= pt.weight;
    kernelXC.ddrhoB *= pt.weight;

    for(auto iShell : closeShells) {
      int iSz= shSizes[iShell];
      int iSt = this->basisset_->mapSh2Bf(iShell);

      for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
        double Ta = kernelXC.ddrhoA*SCRATCH1DATA[iBf] + OmegaADATA[iBf];
        for(auto jShell : closeShells) {
          int jSz= shSizes[jShell];
          int jSt = this->basisset_->mapSh2Bf(jShell);
          double *VXCDATA = result.VXCA.data() + iBf*this->nBasis_;
          for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
            if(jBf < iBf) continue;
            VXCDATA[jBf] += Ta*SCRATCH1DATA[jBf] + SCRATCH1DATA[iBf]*OmegaADATA[jBf]; 
          } // jBf
        } // jShell
      } // iBf
    } // iShell

    if(this->nTCS_ == 1 && !this->isClosedShell){
      for(auto iShell : closeShells) {
        int iSz= shSizes[iShell];
        int iSt = this->basisset_->mapSh2Bf(iShell);
  
        for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
          double Tb = kernelXC.ddrhoB*SCRATCH1DATA[iBf] + OmegaBDATA[iBf];
          for(auto jShell : closeShells) {
            int jSz= shSizes[jShell];
            int jSt = this->basisset_->mapSh2Bf(jShell);
            double *VXCDATB = result.VXCB.data() + iBf*this->nBasis_;
            for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
              if(jBf < iBf) continue;
              VXCDATB[jBf] += Tb*SCRATCH1DATA[jBf] + SCRATCH1DATA[iBf]*OmegaBDATA[jBf]; 
            } // jBf
          } // jShell
        } // iBf
      } // iShell
    } // UKS

    if(doTimings){
      Newend = std::chrono::high_resolution_clock::now();
      T5 += Newend - Newstart;
    }

  };

//  ChronusQ::AtomicGrid AGrid(100,302,ChronusQ::GRID_TYPE::EULERMAC,
//      ChronusQ::GRID_TYPE::LEBEDEV,ChronusQ::ATOMIC_PARTITION::BECKE,
//      this->molecule_->cartArray(),this->molecule_->rIJ(),0,1.0,false);


//if(this->isGGA) cout << "GGA ON " << this->isGGA <<endl ; 
//if(!this->isGGA) cout << "GGA OFF " << this->isGGA <<endl ; 

/*
  ChronusQ::AtomicGrid AGrid(this->nRadDFTGridPts_,this->nAngDFTGridPts_,
      ChronusQ::GRID_TYPE::GAUSSCHEBFST,ChronusQ::GRID_TYPE::LEBEDEV,
      ChronusQ::ATOMIC_PARTITION::BECKE,this->molecule_->cartArray(),
      this->molecule_->rIJ(),0,this->epsScreen,1e6,1.0,false);

*/   
  ChronusQ::AtomicGrid AGrid(this->nRadDFTGridPts_,this->nAngDFTGridPts_,
      static_cast<ChronusQ::GRID_TYPE>(this->dftGrid_),ChronusQ::GRID_TYPE::LEBEDEV,
      static_cast<ChronusQ::ATOMIC_PARTITION>(this->weightScheme_),this->molecule_->cartArray(),
      this->molecule_->rIJ(),0,this->epsScreen,1e6,1.0,false);
  KernelIntegrand<double> res(this->vXA_->cols());
  //Screaning based on cutoff
//  the radCutoof vector is populated in singleSlater_real_guess.cpp

  this->totalEx    = 0.0;
  this->vXA()->setZero();   // Set to zero every occurence of the SCF

  std::vector<double> atomRadCutoff(this->molecule_->nAtoms(),0.0);

  if (this->screenVxc) {
    for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
      for(auto iSh = 0; iSh < this->basisset_->nShell(); iSh++){
        if(this->basisset_->mapSh2Cen(iSh)-1 != iAtm) continue;
        if(this->basisset_->radCutSh()[iSh] > atomRadCutoff[iAtm])
          atomRadCutoff[iAtm] = this->basisset_->radCutSh()[iSh];
//A        atomRadCutoff[iAtm] = 1.e100;
      }
    } 
  }
  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    AGrid.center() = iAtm;
//    cout << "Cutoff Rad = " << atomRadCutoff[iAtm] << "for iAt= " << iAtm << endl;
    if (this->screenVxc) {
      AGrid.setRadCutOff(atomRadCutoff[iAtm]);
      AGrid.findNearestNeighbor();
    } else {
      AGrid.SwitchScreen();
    }
    AGrid.scalingFactor()=0.5 *
      elements[this->molecule_->index(iAtm)].sradius/phys.bohr;

    std::function<void(IntegrationPoint&,KernelIntegrand<double>&)>
      wrapper = std::bind(valVxc,iAtm,std::placeholders::_1,
                  std::placeholders::_2);

    AGrid.integrate<KernelIntegrand<double>>(wrapper,res);
  };
  (*this->vXA_) = 4*math.pi*res.VXCA;
  (*this->vXA_) = this->vXA_->selfadjointView<Lower>();
  this->totalEx = 4*math.pi*res.Energy;
  if(!this->isClosedShell && this->nTCS_ != 2){
    (*this->vXB_) = 4*math.pi*res.VXCB;
    (*this->vXB_) = this->vXB_->selfadjointView<Lower>();
  }

  if(doTimings) {
    cout << "T1 = " << T1.count() << endl;
    cout << "T2 = " << T2.count() << endl;
    cout << "T3 = " << T3.count() << endl;
    cout << "T4 = " << T4.count() << endl;
    cout << "T5 = " << T5.count() << endl;
    cout << "T6 = " << T6.count() << endl;
  }
//cout << "NSkip2 = " << NSkip2 << endl;
//cout << "NSkip3 = " << NSkip3 << endl;
//cout << "NSkip4 = " << NSkip4 << endl;
//cout << "NSkip5 = " << NSkip5 << endl;
//if(doTimings)
//  for(auto i : TF) cout << "TF " << i.count() << endl;

  if(this->printLevel_ >= 3) {
    finish = std::chrono::high_resolution_clock::now();
    duration_formVxc = finish - start;
    prettyPrint(this->fileio_->out,(*this->vXA()),"LDA Vxc alpha");
    if(!this->isClosedShell && this->nTCS_ != 2){
      prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vxc beta");
    }
    this->fileio_->out << "VXC Energy= " <<  this->totalEx << endl, 
    this->fileio_->out << endl << "CPU time for VXC integral:  "
                       << duration_formVxc.count() << " seconds." 
                       << endl;
  }
}

} // Namespace ChronusQ
