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
#include <aointegrals.h>
using ChronusQ::AOIntegrals;

void AOIntegrals::OneEDriver(libint2::Operator iType) {
  std::vector<RealMap> mat;
  int NB = this->nBasis_;
  int NBSq = NB*NB;

  if(iType == libint2::Operator::overlap){
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
  } else if(iType == libint2::Operator::kinetic) {
    mat.push_back(RealMap(this->kinetic_->data(),NB,NB));
  } else if(iType == libint2::Operator::nuclear) {
    mat.push_back(RealMap(this->potential_->data(),NB,NB));
  } else if(iType == libint2::Operator::emultipole1) {
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
    for(auto i = 0, IOff=0; i < 3; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecDipole_->storage()[IOff],NB,NB));
  } else if(iType == libint2::Operator::emultipole2) {
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
    for(auto i = 0, IOff=0; i < 3; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecDipole_->storage()[IOff],NB,NB));
    for(auto i = 0, IOff=0; i < 6; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecQuadpole_->storage()[IOff],NB,NB));
  } else if(iType == libint2::Operator::emultipole3) {
    mat.push_back(RealMap(this->overlap_->data(),NB,NB));
    for(auto i = 0, IOff=0; i < 3; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecDipole_->storage()[IOff],NB,NB));
    for(auto i = 0, IOff=0; i < 6; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecQuadpole_->storage()[IOff],NB,NB));
    for(auto i = 0, IOff=0; i < 10; i++,IOff+=NBSq)
      mat.push_back(RealMap(&this->elecOctpole_->storage()[IOff],NB,NB));
  } else {
    CErr("libint2::Operator type not recognized",this->fileio_->out);
  }
 
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  // Define integral Engine
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(iType,this->basisSet_->maxPrim(),this->basisSet_->maxL(),0);
  engines[0].set_precision(0.0);

  // If engine is V, define nuclear charges
  if(iType == libint2::Operator::nuclear){
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(int i = 0; i < this->molecule_->nAtoms(); i++) {
      q.push_back(
        {
          static_cast<double>((*this->molecularConstants_).atomZ[i]), 
          {
            {
	      (*this->molecularConstants_).cart[0][i],
	      (*this->molecularConstants_).cart[1][i],
	      (*this->molecularConstants_).cart[2][i]
	    }
	  }
	}
      );
    }
    engines[0].set_params(q);
  }
  for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];

  if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(); 
#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
    for(auto s1=0l, s12=0l; s1 < this->basisSet_->nShell(); s1++){
      int bf1_s = this->basisSet_->mapSh2Bf(s1);
      int n1  = this->basisSet_->shells(s1).size();
      for(int s2=0; s2 <= s1; s2++, s12++){
        if(s12 % nthreads != thread_id) continue;
        int bf2_s = this->basisSet_->mapSh2Bf(s2);
        int n2  = this->basisSet_->shells(s2).size();
  
        const double* buff = engines[thread_id].compute(
          this->basisSet_->shells(s1),
          this->basisSet_->shells(s2)
        );

        int IOff = 0;
        for(auto nMat = 0; nMat < mat.size(); nMat++) {
          Eigen::Map<
            const Eigen::Matrix<double,Dynamic,Dynamic,Eigen::RowMajor>>
            bufMat(&buff[IOff],n1,n2);

          mat[nMat].block(bf1_s,bf2_s,n1,n2) = bufMat;
//        for(auto i = 0, bf1 = bf1_s; i < n1; i++, bf1++) 
//        for(auto j = 0, bf2 = bf2_s; j < n2; j++, bf2++){            
//          mat[nMat](bf1,bf2) = bufMat(i,j);
//        }
          IOff += n1*n2;
        }
      }
    }
  } // end openmp parallel
  for(auto nMat = 0; nMat < mat.size(); nMat++) 
    mat[nMat] = mat[nMat].selfadjointView<Lower>();
}

double AOIntegrals::formBeckeW(cartGP gridPt, int iAtm){
//     Generate Frisch (not-normalized yet) Weights (if frischW) according to 
//     the partition schems in (Chem. Phys. Let., 257, 213-223 (1996)) 
//     using Eq. 11 and 14
//
//     Note these Weights have to be normailzed 
//     Generate Becke not-normalized Weights (if becke) according to the 
//     partition schems in (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii 
//     Fuzzi Cells
//
//     Note these Weights have to be normailzed (see normBeckeW) 
       int nAtom = this->molecule_->nAtoms();
       double WW = 1.0;
       double tmp;
       double muij;   ///< elliptical coordinate (ri -rj / Rij)
       cartGP rj;     ///< Cartisian position of Atom j
       cartGP ri;     ///< Cartisian position of Atom i
       ri.set<0>((*this->molecule_->cart())(0,iAtm) );
       ri.set<1>((*this->molecule_->cart())(1,iAtm) );
       ri.set<2>((*this->molecule_->cart())(2,iAtm) );
       for(auto jAtm = 0; jAtm < nAtom; jAtm++){
         if(jAtm != iAtm){
           muij = 0.0;
//       Vector rj (Atoms (j.ne.i) position)
           rj.set<0>((*this->molecule_->cart())(0,jAtm));
           rj.set<1>((*this->molecule_->cart())(1,jAtm));
           rj.set<2>((*this->molecule_->cart())(2,jAtm));
//       Coordinate of the Grid point in elleptical (Eq. 11) 
           muij = (boost::geometry::distance(gridPt,ri) - 
                   boost::geometry::distance(gridPt,rj))/
                   (*this->molecule_->rIJ())(iAtm,jAtm) ;
//       Do the product over all atoms i .ne. j (Eq. 13 using definition Eq. 21 with k=3)
////           if (weightScheme_ == FRISCH) 
////             WW *= 0.5*(1.0-this->twodgrid_->frischpol(muij,0.64));
////           else if (weightScheme_ == BECKE)  
             WW *= 0.5 * 
                   (1.0-this->twodgrid_->voronoii(
                          this->twodgrid_->voronoii(
                            this->twodgrid_->voronoii(muij))));
           }
         }
       return WW;
}; //End formBeckeW


double AOIntegrals::normBeckeW(cartGP gridPt){
//     Normalization of Becke/Frisch Weights
//     (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii Fuzzi Cells
//     Eq. 22
       int   nAtom = this->molecule_->nAtoms();
       double norm = 0.0;
       for(auto iAtm = 0; iAtm < nAtom; iAtm++){
         norm += this->formBeckeW(gridPt,iAtm);
         }
       return norm ;
}; //End normBeckeW



void AOIntegrals::computeAORcrossDel(){
// build R cross Del matrices (Numerical)
// It also creates the grid and can perform also other
// integrals (i.e. <mu/Del/nu> if required.

//  Hard coded for now (we need to set the flag in aointegrals)
//  I also hard coded the type of grid (GaussChebyshev1stGridInf).
//  I hard coded the weight scheme in formBeckeW.
  bool screenVxc = false;   //Screening logical
  double epsScreen     = 1.0e-10; // screening constant
  double epsConv       = 1.0e-7;  // average values used for screening a entire shell
  int maxiter       = 50;         // Max number for figuring the cut off radius for a given shell (given an apsConv in radCut)
  int nRadDFTGridPts_ = 100;     // Number of Radial grid points
  int nAngDFTGridPts_ = 302;    // Number of Angular grid points
  int ngpts;                   // Total Number of grid points
  int nDer = 1 ; //Order of differenzation (In this case it will tell to compute also Del/mu
  double val;  // tmp for loading mu(r) or Del mu(r)
  double x;   // tmp r_x
  double y;   // tmp r_y
  double z;   // tmp r_z
// End Parameters
  this->basisSet_->radcut(epsScreen, maxiter, epsConv);
  ngpts = nRadDFTGridPts_*nAngDFTGridPts_;  // Number of grid point for each center
  OneDGrid * Rad ;                              // Pointer for Radial Grid
  LebedevGrid GridLeb(nAngDFTGridPts_);   // Angular Grid
  Rad = new  EulerMaclaurinGrid(nRadDFTGridPts_,0.0,1.0);   
  GridLeb.genGrid();                            
  auto NBSQ = this->nBasis_*this->nBasis_;
  RealMap rdotpX(&this->RcrossDel_->storage()[0],this->nBasis_,this->nBasis_);
  RealMap rdotpY(&this->RcrossDel_->storage()[NBSQ],this->nBasis_,this->nBasis_);
  RealMap rdotpZ(&this->RcrossDel_->storage()[2*NBSQ],this->nBasis_,this->nBasis_);

  std::vector<RealSparseMatrix>  sparsedmudX_; ///< basis derivative (x) (NbasisxNbasis)
  std::vector<RealSparseMatrix>  sparsedmudY_; ///< basis derivative (y) (NbasisxNbasis)
  std::vector<RealSparseMatrix>  sparsedmudZ_; ///< basis derivative (z) (NbasisxNbasis) 
  std::vector<RealSparseMatrix> sparseMap_;     // BasisFunction Map 
  std::vector<RealSparseMatrix> sparseWeights_; // Weights Map
  std::vector<RealSparseMatrix> sparseDoRho_; // Evaluate density Map
  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    sparseMap_.push_back(RealSparseMatrix(this->nBasis_,ngpts));
    sparseWeights_.push_back(RealSparseMatrix(ngpts,1));
    sparseDoRho_.push_back(RealSparseMatrix(ngpts,1));
//  Derivative
    sparsedmudX_.push_back(RealSparseMatrix(this->nBasis_,ngpts));
    sparsedmudY_.push_back(RealSparseMatrix(this->nBasis_,ngpts));
    sparsedmudZ_.push_back(RealSparseMatrix(this->nBasis_,ngpts));
//  Mapping over iAtom
    RealSparseMatrix *Map        = &sparseMap_[iAtm];
    RealSparseMatrix *WeightsMap = &sparseWeights_[iAtm];
    RealSparseMatrix *DoRhoMap   = &sparseDoRho_[iAtm];
    RealSparseMatrix *MapdX      = &sparsedmudX_[iAtm];
    RealSparseMatrix *MapdY      = &sparsedmudY_[iAtm];
    RealSparseMatrix *MapdZ      = &sparsedmudZ_[iAtm];
    Rad->genGrid(); 
    Rad->atomGrid((elements[this->molecule_->index(iAtm)].sradius)) ;  
    TwoDGrid Raw3Dg(ngpts,Rad,&GridLeb);             
    Raw3Dg.centerGrid(
      (*this->molecule_->cart())(0,iAtm),
      (*this->molecule_->cart())(1,iAtm),
      (*this->molecule_->cart())(2,iAtm)
    );
    std::vector<bool> mapRad_(this->basisSet_->nShell()+1);
    for (auto ipts =0; ipts < ngpts; ipts++){
      cartGP pt = Raw3Dg.gridPtCart(ipts); 
//    auto mapRad_ = this->basisSet_->MapGridBasis(pt);
      this->basisSet_->MapGridBasis(mapRad_,pt);
      // Evaluate each Becke fuzzy call weight, normalize it and muliply by 
      //   the Raw grid weight at that point
      auto bweight = (this->formBeckeW(pt,iAtm)) 
                     / (this->normBeckeW(pt)) ;
      x = bg::get<0>(pt) ;
      y = bg::get<1>(pt) ;
      z = bg::get<2>(pt) ;
      if (screenVxc) {
        if (mapRad_[0] || (bweight < epsScreen)) continue;
        }
         WeightsMap->insert(ipts,0) = Raw3Dg.getweightsGrid(ipts) * bweight;
         for (int s1=0; s1 < this->basisSet_->nShell(); s1++){
           if (screenVxc) {
             if (!mapRad_[s1+1]) continue;
             }
           int bf1_s = this->basisSet_->mapSh2Bf(s1);
           auto shSize = this->basisSet_->shells(s1).size(); 
           libint2::Shell shTmp = this->basisSet_->shells(s1);
           double * s1Eval = this->basisSet_->basisDEval(nDer,shTmp,&pt);
           double * ds1EvalX = s1Eval + shSize;
           double * ds1EvalY = ds1EvalX + shSize;
           double * ds1EvalZ = ds1EvalY + shSize;
           for (auto mu =0; mu < shSize; mu++){ 
             val = s1Eval[mu];
             if (std::abs(val) > epsScreen){
               Map->insert(bf1_s+mu,ipts) = val;
               }
             if (nDer ==1) {
               val = ds1EvalX[mu];
               if (std::abs(val) > epsScreen){
                 MapdX->insert(bf1_s+mu,ipts) = val;
                 }
               val = ds1EvalY[mu];
               if (std::abs(val) > epsScreen){
                 MapdY->insert(bf1_s+mu,ipts) = val;
                 }
               val = ds1EvalZ[mu];
               if (std::abs(val) > epsScreen){
                 MapdZ->insert(bf1_s+mu,ipts) = val;
                 }
               }
             } //loop over basis (within a given ishell)
         } //loop over shells
//            (*rdotpX) += (*WeightsMap).coeff(ipts,0)*(*Map).col(ipts)*(*Map).col(ipts).transpose();
//          Figure how populate RcrossDel_ tensor 
            rdotpX +=  (*WeightsMap).coeff(ipts,0)*(
                          ( (*Map).col(ipts)*(*MapdZ).col(ipts).transpose()*y) -
                          ( (*Map).col(ipts)*(*MapdY).col(ipts).transpose()*z) );

            rdotpY +=  (*WeightsMap).coeff(ipts,0)*(
                          ( (*Map).col(ipts)*(*MapdX).col(ipts).transpose()*z) -
                          ( (*Map).col(ipts)*(*MapdZ).col(ipts).transpose()*x) );

            rdotpZ +=  (*WeightsMap).coeff(ipts,0)*(
                          ( (*Map).col(ipts)*(*MapdY).col(ipts).transpose()*x) -
                          ( (*Map).col(ipts)*(*MapdX).col(ipts).transpose()*y) );
         if (screenVxc && Map->col(ipts).norm() > epsScreen){
         DoRhoMap->insert(ipts,0) = 2;}
       }//End ipts
    }//End iAtm
//       Finishing the numerical integration (int * 4*pi)
         rdotpX *= 4.0*math.pi;
         rdotpY *= 4.0*math.pi;
         rdotpZ *= 4.0*math.pi;
//       printing
//         prettyPrint(cout,rdotpX,"Numeric <dipole vel> - x comp");
//         prettyPrint(cout,rdotpY,"Numeric <dipole vel> - y comp");
//         prettyPrint(cout,rdotpZ,"Numeric <dipole vel> - z comp");
//end for now (comment later)
//  cout << "Call HERE " <<endl;
//  CErr();
} 

void AOIntegrals::computeAOOneE(){
  // Collect Relevant data into a struct (odd, but convienient) 
  this->iniMolecularConstants();
/*
  BasisSet newBasis;
  newBasis.communicate(*this->fileio_);
  this->basisSet_->genUCvomLocal(&newBasis);
  newBasis.makeMaps(this->molecule_);
  if(this->isPrimary) {
    cout << "Old Basis" << endl;
    for(auto i = 0 ; i < this->basisSet_->nShell(); i++)
      cout << this->basisSet_->shells(i) << endl;
    cout << "New Basis" << endl;
    for(auto i = 0 ; i < newBasis.nShell(); i++)
      cout << newBasis.shells(i) << endl;
  }
*/

  // Start timer for one-electron integral evaluation
  auto oneEStart = std::chrono::high_resolution_clock::now();

  // Compute and time overlap integrals
  auto OStart = std::chrono::high_resolution_clock::now();
  if(this->maxMultipole_ >= 3)      OneEDriver(libint2::Operator::emultipole3);
  else if(this->maxMultipole_ == 2) OneEDriver(libint2::Operator::emultipole2);
  else if(this->maxMultipole_ == 1) OneEDriver(libint2::Operator::emultipole1);
  else OneEDriver(libint2::Operator::overlap);
  if(this->maxMultipole_ == 4) this->computeAORcrossDel();
  auto OEnd = std::chrono::high_resolution_clock::now();



  // Compute and time kinetic integrals
  auto TStart = std::chrono::high_resolution_clock::now();
  OneEDriver(libint2::Operator::kinetic);
  auto TEnd = std::chrono::high_resolution_clock::now();

  // Compute and time nuclear attraction integrals (negative sign is factored in)
  auto VStart = std::chrono::high_resolution_clock::now();

// this->useFiniteWidthNuclei = true;
  if(this->isPrimary && this->useFiniteWidthNuclei) 
    this->finiteWidthPotential();
  else                
    OneEDriver(libint2::Operator::nuclear);

  auto VEnd = std::chrono::high_resolution_clock::now();

// add DKH correction to kinetic energy
//  if (this->isPrimary) this->DKH0();

// Build Core Hamiltonian
  (*this->oneE_) = (*this->kinetic_) + (*this->potential_);

  // Get end time of one-electron integral evaluation
  auto oneEEnd = std::chrono::high_resolution_clock::now();

  // Compute Orthonormal transformation matricies
  this->computeOrtho();

 // this->formP2Transformation();

  if(this->printLevel_ >= 2) this->printOneE();

  // Compute time differenes
  this->OneED = oneEEnd - oneEStart;
  this->SED = OEnd - OStart;
  this->TED = TEnd - TStart;
  this->VED = VEnd - VStart;
  this->haveAOOneE = true;
  this->breakUpMultipole();
  if(this->isPrimary) this->writeOneE();
  /*
  if(this->isPrimary) {
    for(auto ix = 0, ioff = 0; ix < 3; ix++, ioff += this->nBasis_*this->nBasis_){
      RealMap A(&this->RcrossDel_->storage()[ioff],this->nBasis_,this->nBasis_);
      prettyPrint(cout,A,"R:"+std::to_string(ix));
    }
  }
  */
}

void AOIntegrals::computeSchwartz(){
  if(getRank() == 0) {
    RealMatrix *ShBlk; 
    this->schwartz_->setZero();
 
    // Define Integral Engine
    libint2::Engine engine(
        libint2::Operator::coulomb,this->basisSet_->maxPrim(),
        this->basisSet_->maxL(),0);

    engine.set_precision(0.); // Don't screen primitives during schwartz
 
    auto start =  std::chrono::high_resolution_clock::now();
    for(int s1=0; s1 < this->basisSet_->nShell(); s1++){
      int n1  = this->basisSet_->shells(s1).size();
      for(int s2=0; s2 <= s1; s2++){
        int n2  = this->basisSet_->shells(s2).size();
 
        const auto* buff = engine.compute(
          this->basisSet_->shells(s1),
          this->basisSet_->shells(s2),
          this->basisSet_->shells(s1),
          this->basisSet_->shells(s2)
        );
 
        
        ShBlk = new RealMatrix(n1,n2);
        ShBlk->setZero();
 
        for(auto i = 0, ij = 0; i < n1; i++)
        for(auto j = 0; j < n2; j++, ij++){
          (*ShBlk)(i,j) = buff[ij*n1*n2 + ij];
        }
 
        (*this->schwartz_)(s1,s2) = std::sqrt(ShBlk->lpNorm<Infinity>());
        
        delete ShBlk;
      }
    }
    auto finish =  std::chrono::high_resolution_clock::now();
    this->SchwartzD = finish - start;
    (*this->schwartz_) = this->schwartz_->selfadjointView<Lower>();
  }
#ifdef CQ_ENABLE_MPI
  MPI_Bcast(this->schwartz_->data(),this->schwartz_->size(),MPI_DOUBLE,0,
    MPI_COMM_WORLD);
#endif

  this->haveSchwartz = true;
}
void AOIntegrals::computeAOTwoE(){
  if(!this->haveSchwartz) this->computeSchwartz();
  if(getRank() != 0) return;

#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(libint2::Operator::coulomb,
      this->basisSet_->maxPrim(),this->basisSet_->maxL(),0);
  engines[0].set_precision(std::numeric_limits<double>::epsilon());

  for(int i=1; i<nthreads; i++) engines[i] = engines[0];
  if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(); 

  this->aoERI_->fill(0.0);

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
  for(int s1 = 0, s1234=0; s1 < this->basisSet_->nShell(); s1++) {
    int bf1_s = this->basisSet_->mapSh2Bf(s1);
    int n1    = this->basisSet_->shells(s1).size();
    for(int s2 = 0; s2 <= s1; s2++) {
      int bf2_s = this->basisSet_->mapSh2Bf(s2);
      int n2    = this->basisSet_->shells(s2).size();
      for(int s3 = 0; s3 <= s1; s3++) {
        int bf3_s = this->basisSet_->mapSh2Bf(s3);
        int n3    = this->basisSet_->shells(s3).size();
        int s4_max = (s1 == s3) ? s2 : s3;
        for(int s4 = 0; s4 <= s4_max; s4++, s1234++) {

          if(s1234 % nthreads != thread_id) continue;

          int bf4_s = this->basisSet_->mapSh2Bf(s4);
          int n4    = this->basisSet_->shells(s4).size();
    
          // Schwartz and Density screening
          if((*this->schwartz_)(s1,s2) * (*this->schwartz_)(s3,s4)
              < this->thresholdSchwartz_ ) continue;
 
          const double* buff = engines[thread_id].compute(
            this->basisSet_->shells(s1),
            this->basisSet_->shells(s2),
            this->basisSet_->shells(s3),
            this->basisSet_->shells(s4));


          std::vector<std::array<int,4>> lower;
          std::vector<std::array<int,4>> upper;
          for(int i = 0, bf1 = bf1_s, ijkl = 0 ; i < n1; ++i, bf1++) 
          for(int j = 0, bf2 = bf2_s           ; j < n2; ++j, bf2++) 
          for(int k = 0, bf3 = bf3_s           ; k < n3; ++k, bf3++) 
          for(int l = 0, bf4 = bf4_s           ; l < n4; ++l, bf4++, ++ijkl) {
            (*this->aoERI_)(bf1,bf2,bf3,bf4) = buff[ijkl];
            (*this->aoERI_)(bf1,bf2,bf4,bf3) = buff[ijkl];
            (*this->aoERI_)(bf2,bf1,bf3,bf4) = buff[ijkl];
            (*this->aoERI_)(bf2,bf1,bf4,bf3) = buff[ijkl];
            (*this->aoERI_)(bf3,bf4,bf1,bf2) = buff[ijkl];
            (*this->aoERI_)(bf4,bf3,bf1,bf2) = buff[ijkl];
            (*this->aoERI_)(bf3,bf4,bf2,bf1) = buff[ijkl];
            (*this->aoERI_)(bf4,bf3,bf2,bf1) = buff[ijkl];
            /*
            if(this->nTCS_ == 2){
              (*this->aoERI_)(bf1+1,bf2+1,bf3+1,bf4+1) = buff[ijkl];
              (*this->aoERI_)(bf1+1,bf2+1,bf4+1,bf3+1) = buff[ijkl];
              (*this->aoERI_)(bf2+1,bf1+1,bf3+1,bf4+1) = buff[ijkl];
              (*this->aoERI_)(bf2+1,bf1+1,bf4+1,bf3+1) = buff[ijkl];
              (*this->aoERI_)(bf3+1,bf4+1,bf1+1,bf2+1) = buff[ijkl];
              (*this->aoERI_)(bf4+1,bf3+1,bf1+1,bf2+1) = buff[ijkl];
              (*this->aoERI_)(bf3+1,bf4+1,bf2+1,bf1+1) = buff[ijkl];
              (*this->aoERI_)(bf4+1,bf3+1,bf2+1,bf1+1) = buff[ijkl];

              (*this->aoERI_)(bf1+1,bf2+1,bf3,bf4)     = buff[ijkl];
              (*this->aoERI_)(bf1+1,bf2+1,bf4,bf3)     = buff[ijkl];
              (*this->aoERI_)(bf2+1,bf1+1,bf3,bf4)     = buff[ijkl];
              (*this->aoERI_)(bf2+1,bf1+1,bf4,bf3)     = buff[ijkl];
              (*this->aoERI_)(bf3+1,bf4+1,bf1,bf2)     = buff[ijkl];
              (*this->aoERI_)(bf4+1,bf3+1,bf1,bf2)     = buff[ijkl];
              (*this->aoERI_)(bf3+1,bf4+1,bf2,bf1)     = buff[ijkl];
              (*this->aoERI_)(bf4+1,bf3+1,bf2,bf1)     = buff[ijkl];

              (*this->aoERI_)(bf1,bf2,bf3+1,bf4+1)     = buff[ijkl];
              (*this->aoERI_)(bf1,bf2,bf4+1,bf3+1)     = buff[ijkl];
              (*this->aoERI_)(bf2,bf1,bf3+1,bf4+1)     = buff[ijkl];
              (*this->aoERI_)(bf2,bf1,bf4+1,bf3+1)     = buff[ijkl];
              (*this->aoERI_)(bf3,bf4,bf1+1,bf2+1)     = buff[ijkl];
              (*this->aoERI_)(bf4,bf3,bf1+1,bf2+1)     = buff[ijkl];
              (*this->aoERI_)(bf3,bf4,bf2+1,bf1+1)     = buff[ijkl];
              (*this->aoERI_)(bf4,bf3,bf2+1,bf1+1)     = buff[ijkl];
            }
            */
	  }

	}
      }
    }
  }
  } // OMP Parallel
  this->haveAOTwoE = true;

//this->fileio_->out << "Two-Electron Integrals (ERIs)" << endl;
//for(auto i = 0; i < this->nBasis_; i++)
//for(auto j = 0; j < this->nBasis_; j++)
//for(auto k = 0; k < this->nBasis_; k++)
//for(auto l = 0; l < this->nBasis_; l++){
//  this->fileio_->out << "(" << i << "," << j << "|" << k << "," << l << ")  ";
//  this->fileio_->out << (*this->aoERI_)(i,j,k,l) << endl;
//};
}


void AOIntegrals::breakUpMultipole(){
  int NB = this->nBasis_;
  this->elecDipoleSep_.clear();
  this->elecQuadpoleSep_.clear();
  this->elecOctpoleSep_.clear();
  if(this->maxMultipole_ >= 1)
    for(int ixyz = 0, iBuf = 0; ixyz < 3; ++ixyz, iBuf += NB*NB)
      this->elecDipoleSep_.emplace_back(
        ConstRealMap(&this->elecDipole_->storage()[iBuf],NB,NB)
      );
  if(this->maxMultipole_ >= 2)
    for(int jxyz = 0, iBuf = 0; jxyz < 3; ++jxyz               )
    for(int ixyz = jxyz       ; ixyz < 3; ++ixyz, iBuf += NB*NB)
      this->elecQuadpoleSep_.emplace_back(
        ConstRealMap(&this->elecQuadpole_->storage()[iBuf],NB,NB)
      );
  if(this->maxMultipole_ >= 3)
    for(int kxyz = 0, iBuf = 0; kxyz < 3; ++kxyz               )
    for(int jxyz = kxyz       ; jxyz < 3; ++jxyz               )
    for(int ixyz = jxyz       ; ixyz < 3; ++ixyz, iBuf += NB*NB)
      this->elecOctpoleSep_.emplace_back(
        ConstRealMap(&this->elecOctpole_->storage()[iBuf],NB,NB)
      );
  if(this->maxMultipole_ >= 4)
    for(int ixyz = 0, iBuf = 0; ixyz < 3; ++ixyz, iBuf += NB*NB)
      this->orbitalDipoleSep_.emplace_back(
        ConstRealMap(&this->RcrossDel_->storage()[iBuf],NB,NB)
      );
};

void AOIntegrals::finiteWidthPotential() {
#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
#else
  int nthreads = 1;
#endif

  this->molecule_->generateFiniteWidthNuclei();
  this->potential_->setZero();

  // Define integral Engine
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(libint2::Operator::coulomb,
      this->basisSet_->maxPrim(), this->basisSet_->maxL(),0);
  engines[0].set_precision(0.0);

  for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];

  if(!this->basisSet_->haveMapSh2Bf) this->basisSet_->makeMapSh2Bf(); 
#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif
    for(auto s1=0l, s12=0l; s1 < this->basisSet_->nShell(); s1++){
      int bf1_s = this->basisSet_->mapSh2Bf(s1);
      int n1  = this->basisSet_->shells(s1).size();
      for(int s2=0; s2 <= s1; s2++, s12++){
        if(s12 % nthreads != thread_id) continue;
        int bf2_s = this->basisSet_->mapSh2Bf(s2);
        int n2  = this->basisSet_->shells(s2).size();
  
        for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
          const double* buff = engines[thread_id].compute(
            this->basisSet_->shells(s1),
            this->basisSet_->shells(s2),
            this->molecule_->nucShell(iAtm),
            libint2::Shell::unit()
          );

          Eigen::Map<
            const Eigen::Matrix<double,Dynamic,Dynamic,Eigen::RowMajor>>
            bufMat(buff,n1,n2);

          this->potential_->block(bf1_s,bf2_s,n1,n2) -= bufMat;
          /*
          for(auto i = 0, bf1 = bf1_s; i < n1; i++, bf1 += this->nTCS_) 
          for(auto j = 0, bf2 = bf2_s; j < n2; j++, bf2 += this->nTCS_){            
            (*this->potential_)(bf1,bf2) -= bufMat(i,j);
            if(this->nTCS_ == 2) 
              (*this->potential_)(bf1+1,bf2+1) -= bufMat(i,j);
          }
          */

        }
      }
    }
  } // end openmp parallel

  (*this->potential_) = this->potential_->selfadjointView<Lower>();
};

void AOIntegrals::computeLowdin() {
  char JOBZ = 'V';
  char UPLO = 'L';
  int INFO;

  int LWORK = 3 * this->nBasis_;
  std::vector<double> ovlpEigValues(this->nBasis_);
  std::vector<double> ovlpEigVectrs(this->nBasis_ * this->nBasis_);
  std::vector<double> ovlpCpy(this->nBasis_ * this->nBasis_);
  std::vector<double> WORK(LWORK);

  std::memcpy(&ovlpEigVectrs[0],this->overlap_->data(),
      this->nBasis_*this->nBasis_*sizeof(double));

  dsyev_(&JOBZ,&UPLO,&this->nBasis_,&ovlpEigVectrs[0],&this->nBasis_,
      &ovlpEigValues[0],&WORK[0],&LWORK,&INFO);

  std::memcpy(&ovlpCpy[0],&ovlpEigVectrs[0],
      this->nBasis_*this->nBasis_*sizeof(double));

  RealMap S(&ovlpCpy[0],      this->nBasis_,this->nBasis_);
  RealMap V(&ovlpEigVectrs[0],this->nBasis_,this->nBasis_);

  for(auto i = 0; i < this->nBasis_; i++)
    S.col(i) /= std::sqrt(ovlpEigValues[i]);

  (*this->ortho1_) = S * V.transpose();

  for(auto i = 0; i < this->nBasis_; i++)
    S.col(i) *= ovlpEigValues[i];

  (*this->ortho2_) = S * V.transpose();
//prettyPrint(cout,(*this->ortho1_),"XN");
};
