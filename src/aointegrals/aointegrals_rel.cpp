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
#include <quantum.h>
#include <aointegrals.h>
#include <grid2.h>
using ChronusQ::AOIntegrals;

void AOIntegrals::formP2Transformation(){
  if(!this->isPrimary) return;
  this->basisSet_->makeMapPrim2Bf();
  auto unContractedShells = this->basisSet_->uncontractBasis();
  int nUncontracted = 0;
  for(auto i : unContractedShells) nUncontracted += i.size();

  int nUnSq = nUncontracted * nUncontracted;

  RealMap SUncontracted(
    this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap TUncontracted(
    this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap VUncontracted(
    this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  
  //General Scratch
  RealMap nUnSqScratch(
    this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  ComplexMap C4nUnSqScratch(
    this->memManager_->malloc<dcomplex>(4*nUnSq),2*nUncontracted,2*nUncontracted);

  libint2::Engine engineS(
      libint2::Operator::overlap,1,this->basisSet_->maxL(),0);
  libint2::Engine engineT(
      libint2::Operator::kinetic,1,this->basisSet_->maxL(),0);
  libint2::Engine engineV(
      libint2::Operator::nuclear,1,this->basisSet_->maxL(),0);
  libint2::Engine engineC(
      libint2::Operator::coulomb,1,this->basisSet_->maxL(),0);

  engineS.set_precision(0.0);
  engineT.set_precision(0.0);
  engineV.set_precision(0.0);
  engineC.set_precision(0.0);

  // Loop through and uncontract S and T
  //
  for(auto s1 = 0l, s12 = 0l, bf1_s = 0l; s1 < unContractedShells.size();
      bf1_s += unContractedShells[s1].size(), ++s1){
    int n1 = unContractedShells[s1].size();

  for(auto s2 = 0l, bf2_s = 0l; s2 < unContractedShells.size();
      bf2_s += unContractedShells[s2].size(), ++s2, ++s12){
    int n2 = unContractedShells[s2].size();

    const double * buffS = engineS.compute(
        unContractedShells[s1],
        unContractedShells[s2]
    );
    const double * buffT = engineT.compute(
        unContractedShells[s1],
        unContractedShells[s2]
    );
    Eigen::Map<
      const Eigen::Matrix<double,Dynamic,Dynamic,Eigen::RowMajor>>
      bufMatS(buffS,n1,n2);
    Eigen::Map<
      const Eigen::Matrix<double,Dynamic,Dynamic,Eigen::RowMajor>>
      bufMatT(buffT,n1,n2);

    SUncontracted.block(bf1_s,bf2_s,n1,n2) = bufMatS;
    TUncontracted.block(bf1_s,bf2_s,n1,n2) = bufMatT;
    
    //------------------------------------------------
    //Finite Width Nuclei
    
    for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
      const double* buff = engineC.compute(
        unContractedShells[s1],
        unContractedShells[s2],
        this->molecule_->nucShell(iAtm),
        libint2::Shell::unit()
      );

      Eigen::Map<
        const Eigen::Matrix<double,Dynamic,Dynamic,Eigen::RowMajor>>
        bufMatV(buff,n1,n2);

      VUncontracted.block(bf1_s,bf2_s,n1,n2) -= bufMatV;
    }
    //------------------------------------------------

    //Point Nuclei
    /* 
    const double * buffV = engineV.compute(
        unContractedShells[s1],
        unContractedShells[s2]
    );

    Eigen::Map<
      const Eigen::Matrix<double,Dynamic,Dynamic,Eigen::RowMajor>>
      bufMatV(buffV,n1,n2);

    VUncontracted.block(bf1_s,bf2_s,n1,n2) = bufMatV;
    */
    //------------------------------------------------

  } // s2
  } // s1

  SUncontracted = SUncontracted.selfadjointView<Lower>();
  TUncontracted = TUncontracted.selfadjointView<Lower>();
  
  if(this->printLevel_ >= 3){
    prettyPrintSmart(this->fileio_->out,SUncontracted,"S uncontracted");
    prettyPrintSmart(this->fileio_->out,TUncontracted,"T uncontracted");
  }

  // ------------------------------------------------------------------------
  // Probably should move this to be in the function that makes the mapPrim2Bf. 
  // This is necessary to get the correct mapping.
  RealMap SnonRel(this->memManager_->malloc<double>(this->nBasis_*this->nBasis_),
    this->nBasis_,this->nBasis_);
  SnonRel = (*this->basisSet_->mapPrim2Bf()) * SUncontracted
	* (*this->basisSet_->mapPrim2Bf()).transpose();
  for (auto row = 0; row < this->basisSet_->nBasis(); row++){
      (*this->basisSet_->mapPrim2Bf()).block(row,0,1,nUncontracted) /=
        std::sqrt(SnonRel(row,row)); //scale by appropriate factor
    }
  this->memManager_->free(SnonRel.data(),this->nBasis_*this->nBasis_);
  // -----------------------------------------------------------------------

  RealMap SUn(this->memManager_->malloc<double>(nUnSq),
    nUncontracted,nUncontracted);
  SUn = SUncontracted.real(); // Save S for later

  char JOBU = 'O', JOBVT = 'N';
  int LWORK = 6*nUncontracted;
  int INFO;
  std::vector<double> ovlpEigValues(nUncontracted);
  std::vector<double> WORK(LWORK);

  // Get eigenvalues for overlap matrix S (via SVD) 
  // Store the orthogonal transformation U back in S

  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,SUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],SUncontracted.data(),&nUncontracted,
      SUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);

  // Count linear dependencies here
  int nZero = 0;
  for(auto iS = 0; iS < nUncontracted; iS++)
    if(std::abs(ovlpEigValues[iS]) < 1e-12) nZero++;

  cout << "NZERO " << nZero << endl;
  
  // Normalize
  // What happens when we divide by zero?
  for(auto iS = 0; iS < nUncontracted; iS++){
    SUncontracted.col(iS) /= std::sqrt(ovlpEigValues[iS]);
  }

  // Put the kinetic energy in the orthonormal basis 

  nUnSqScratch = SUncontracted.transpose() * TUncontracted;
  TUncontracted = nUnSqScratch * SUncontracted;

  // Get rid of the linear dependencies?
  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,TUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],TUncontracted.data(),&nUncontracted,
      TUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);

  // Form the K transformation matrix from both pieces
  // (diagonalizes S) and (diagonalizes T). See eq. 12 in Rieher's paper from 2013

  if(this->printLevel_ >= 3){
    prettyPrintSmart(this->fileio_->out,SUncontracted,"S transformation");
    prettyPrintSmart(this->fileio_->out,TUncontracted,"T transformation");
    prettyPrintSmart(this->fileio_->out,TUncontracted.transpose()*TUncontracted, "T^{dagger} T");
  }

  RealMap UK(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  UK = SUncontracted * TUncontracted;

  if(this->printLevel_ >= 3){ 
    prettyPrintSmart(this->fileio_->out,UK,"Uk transformation");
  }

  // Now we transform V to V' 
  RealMap P2_Potential(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);

  nUnSqScratch = UK.transpose() * VUncontracted;
  P2_Potential = nUnSqScratch * UK;

  this->memManager_->free(SUncontracted.data(),nUnSq);
  this->memManager_->free(TUncontracted.data(),nUnSq);
  this->memManager_->free(VUncontracted.data(),nUnSq);
  // Next we need to get W' (from pVp integrals)

  RealVecMap SCRATCH1UnContracted(this->memManager_->malloc<double>(nUncontracted),nUncontracted);
  RealMap SCRATCH2UnContracted(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealVecMap SCRATCHDXUnContracted(this->memManager_->malloc<double>(nUncontracted),nUncontracted);
  RealVecMap SCRATCHDYUnContracted(this->memManager_->malloc<double>(nUncontracted),nUncontracted);
  RealVecMap SCRATCHDZUnContracted(this->memManager_->malloc<double>(nUncontracted),nUncontracted);

  // compute pVp numerically
  auto PVP = [&](IntegrationPoint pt, std::vector<RealMatrix> &result) {
    for(auto iShell = 0, b_s = 0; iShell < unContractedShells.size();

         b_s += unContractedShells[iShell].size(),++iShell) {
      int size= unContractedShells[iShell].size();

      double * 
        buff=this->basisSet_->basisDEval(1,unContractedShells[iShell], &pt.pt);

      RealMap bMap( buff         ,size,1);
      RealMap dxMap(buff + size  ,size,1);
      RealMap dyMap(buff + 2*size,size,1);
      RealMap dzMap(buff + 3*size,size,1);

      SCRATCH1UnContracted.block( b_s,0,size,1) = bMap;
      SCRATCHDXUnContracted.block(b_s,0,size,1) = dxMap;
      SCRATCHDYUnContracted.block(b_s,0,size,1) = dyMap;
      SCRATCHDZUnContracted.block(b_s,0,size,1) = dzMap;

     // delete [] buff;
    };

    std::vector<std::pair<double,std::array<double,3>>> q;
    q.push_back(
      {1.0, {{bg::get<0>(pt.pt),bg::get<1>(pt.pt),bg::get<2>(pt.pt)}}});
    engineV.set_params(q);

    for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
      // Gaussian Nuclei
      const double * gamma = engineV.compute(this->molecule_->nucShell(iAtm),
          libint2::Shell::unit());
                                                    
    //SCRATCH2UnContracted.noalias() = 
    //  SCRATCH1UnContracted * SCRATCH1UnContracted.transpose();
    //result[0].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDXUnContracted * SCRATCHDXUnContracted.transpose();
      result[1].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDXUnContracted * SCRATCHDYUnContracted.transpose();
      result[2].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDXUnContracted * SCRATCHDZUnContracted.transpose();
      result[3].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDYUnContracted * SCRATCHDXUnContracted.transpose();
      result[4].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDYUnContracted * SCRATCHDYUnContracted.transpose();
      result[5].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDYUnContracted * SCRATCHDZUnContracted.transpose();
      result[6].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDZUnContracted * SCRATCHDXUnContracted.transpose();
      result[7].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDZUnContracted * SCRATCHDYUnContracted.transpose();
      result[8].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDZUnContracted * SCRATCHDZUnContracted.transpose();
      result[9].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;
    }

  };

  auto atomicCenters = this->molecule_->cartArray();
  this->basisSet_->radcut(1.0e-10, 50, 1.0e-7);
  std::vector<double> atomRadCutoff(this->molecule_->nAtoms(),0.0);
  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    for(auto iSh = 0; iSh < this->basisSet_->nShell(); iSh++){
      if(this->basisSet_->mapSh2Cen(iSh)-1 != iAtm) continue;
      if(this->basisSet_->radCutSh()[iSh] > atomRadCutoff[iAtm])
        atomRadCutoff[iAtm] = this->basisSet_->radCutSh()[iSh];
    }
  } 

  AtomicGrid AGrid(100,590,GAUSSCHEBFST,LEBEDEV,BECKE,atomicCenters,
    this->molecule_->rIJ(),0,-1,1e6,1.0,false);

  std::vector<RealMatrix> numPot(10,RealMatrix::Zero(nUncontracted,
        nUncontracted));

  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    AGrid.center() = iAtm;
    AGrid.setRadCutOff(atomRadCutoff[iAtm]);
    AGrid.findNearestNeighbor();
    AGrid.scalingFactor()=
      0.5*elements[this->molecule_->index(iAtm)].sradius/phys.bohr;
    AGrid.integrate<std::vector<RealMatrix>>(PVP,numPot);
  };

  for(auto i = 0; i < 10; i++) numPot[i] *= 4 * math.pi;

  this->memManager_->free(SCRATCH1UnContracted.data(),nUncontracted);
  this->memManager_->free(SCRATCH2UnContracted.data(),nUnSq);
  this->memManager_->free(SCRATCHDXUnContracted.data(),nUncontracted);
  this->memManager_->free(SCRATCHDYUnContracted.data(),nUncontracted);
  this->memManager_->free(SCRATCHDZUnContracted.data(),nUncontracted);
  
  RealMap PVPS(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap PVPX(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap PVPY(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap PVPZ(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  
  PVPS = numPot[1] + numPot[5] + numPot[9]; 
  PVPX = numPot[6] - numPot[8];
  PVPY = numPot[7] - numPot[3];
  PVPZ = numPot[2] - numPot[4]; 

  // Apply the Uk unitary transformation ( Uk' * PVP * Uk)
  nUnSqScratch = UK.adjoint() * PVPS;
  PVPS = nUnSqScratch * UK;
  nUnSqScratch = UK.adjoint() * PVPX;
  PVPX = nUnSqScratch * UK;
  nUnSqScratch = UK.adjoint() * PVPY;
  PVPY = nUnSqScratch * UK;
  nUnSqScratch = UK.adjoint() * PVPZ;
  PVPZ = nUnSqScratch * UK;

  if(this->printLevel_ >= 3){
    prettyPrintSmart(cout,P2_Potential,"V");
    prettyPrintSmart(cout,PVPS,"dot(P,VP)");
    prettyPrintSmart(cout,PVPX,"cross(P,VP) X");
    prettyPrintSmart(cout,PVPY,"cross(P,VP) Y");
    prettyPrintSmart(cout,PVPZ,"cross(P,VP) Z");
    cout << "|V| = " << P2_Potential.squaredNorm() << endl;
    cout << "|dot(P,VP)| = " << PVPS.squaredNorm() << endl;
    cout << "|cross(P,VP) X| = " << PVPX.squaredNorm() << endl;
    cout << "|cross(P,VP) Y| = " << PVPY.squaredNorm() << endl;
    cout << "|cross(P,VP) Z| = " << PVPZ.squaredNorm() << endl;
  }

  RealVecMap PMap(&ovlpEigValues[0],nUncontracted);

  if(this->printLevel_ >= 3){
    prettyPrintSmart(cout,PMap,"P^2");
  }

  PMap = 2*PMap;
  PMap = PMap.cwiseSqrt();
  PMap = PMap.cwiseInverse();

  nUnSqScratch = PMap.asDiagonal() * PVPS;
  PVPS = nUnSqScratch * PMap.asDiagonal();
  nUnSqScratch = PMap.asDiagonal() * PVPX;
  PVPX = nUnSqScratch * PMap.asDiagonal();
  nUnSqScratch = PMap.asDiagonal() * PVPY;
  PVPY = nUnSqScratch * PMap.asDiagonal();
  nUnSqScratch = PMap.asDiagonal() * PVPZ;
  PVPZ = nUnSqScratch * PMap.asDiagonal();

  if(this->printLevel_ >= 3){
    cout << "|dot(P,VP)| = " << PVPS.squaredNorm() << endl;
    cout << "|cross(P,VP) X| = " << PVPX.squaredNorm() << endl;
    cout << "|cross(P,VP) Y| = " << PVPY.squaredNorm() << endl;
    cout << "|cross(P,VP) Z| = " << PVPZ.squaredNorm() << endl;
    prettyPrintSmart(cout,PVPS,"scaled dot(P,VP)");
    prettyPrintSmart(cout,PVPX,"scaled cross(P,VP) X");
    prettyPrintSmart(cout,PVPY,"scaled cross(P,VP) Y");
    prettyPrintSmart(cout,PVPZ,"scaled cross(P,VP) Z");
  }

  // Apply Boettger Scaling Here?
  if(this->twoEFudge == 2){
  
    if(this->printLevel_ >= 2){
      (this->fileio_->out) << "Applying 2e fudge factor to PVP integrals" << endl;
    }
  // Q(l) as an function:
  //    double QofL(double L) { return L*(L+1)*(2*L+1)/3; };
  
  // Loop over all shells
    int idxRow = 0;
    for (auto iShell = 0; iShell < this->basisSet_->nShell(); ++iShell) {
        int nIPrims; //allocate space
        int nJPrims;
        int idxCol = 0;
      for (auto jShell = 0; jShell < this->basisSet_->nShell(); ++jShell) {
        auto iL = this->basisSet_->shells(iShell).contr[0].l;
        auto jL = this->basisSet_->shells(jShell).contr[0].l;
          
        nIPrims = this->basisSet_->shells(iShell).size() *
            this->basisSet_->shells(iShell).contr[0].coeff.size();
        nJPrims = this->basisSet_->shells(jShell).size() *
            this->basisSet_->shells(jShell).contr[0].coeff.size();
    //    cout << "nIPrims: " << nIPrims  << " , nJPrims: " << nJPrims << endl;
        if (iL > 0 and jL > 0) {
          double fudgeFactor = (iL)*(iL+1)*(2*iL+1)/3;
          fudgeFactor *= (jL)*(jL+1)*(2*jL+1)/3;

          fudgeFactor /= this->molecule_->atomicZ(this->basisSet_->mapSh2Cen(iShell)-1); 
          fudgeFactor /= this->molecule_->atomicZ(this->basisSet_->mapSh2Cen(jShell)-1); 
          fudgeFactor = sqrt(fudgeFactor);
          
          PVPX.block(idxRow,idxCol,nIPrims,nJPrims) -= 
            PVPX.block(idxRow,idxCol,nIPrims,nJPrims) * fudgeFactor;
          
          PVPY.block(idxRow,idxCol,nIPrims,nJPrims) -= 
            PVPY.block(idxRow,idxCol,nIPrims,nJPrims) * fudgeFactor;
          
          PVPZ.block(idxRow,idxCol,nIPrims,nJPrims) -= 
            PVPZ.block(idxRow,idxCol,nIPrims,nJPrims) * fudgeFactor;
        }
        idxCol += nJPrims;
      } //loop jShells
      idxRow += nIPrims;
    }//loop iShells
  }


  ComplexMap W(this->memManager_->malloc<dcomplex>(4*nUnSq),
    2*nUncontracted,2*nUncontracted);
  W.block(0,0,nUncontracted,nUncontracted).real() = PVPS;
  W.block(0,0,nUncontracted,nUncontracted).imag() = PVPZ;
  W.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = 
    PVPS;
  W.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).imag() = 
   - PVPZ;
  W.block(0,nUncontracted,nUncontracted,nUncontracted).real() = PVPY;
  W.block(0,nUncontracted,nUncontracted,nUncontracted).imag() = PVPX;
  W.block(nUncontracted,0,nUncontracted,nUncontracted).real() = -PVPY;
  W.block(nUncontracted,0,nUncontracted,nUncontracted).imag() = PVPX;

  // No longer need PVP integrals, so free memory here.
  this->memManager_->free(PVPS.data(),nUnSq); 
  this->memManager_->free(PVPX.data(),nUnSq); 
  this->memManager_->free(PVPY.data(),nUnSq); 
  this->memManager_->free(PVPZ.data(),nUnSq); 
  
  // -------------------------------------------------
  // Put all pieces of core Hamiltonian in block form:
  // [ V'      cp    ]
  // [    	 	   ]
  // [ cp   W'-2mc^2 ]
  // -------------------------------------------------

  PMap = PMap.cwiseInverse(); //switch from p^-1 back to p
  
  ComplexMap CORE_HAMILTONIAN(this->memManager_->malloc<dcomplex>(16*nUnSq),
    4*nUncontracted,4*nUncontracted);

  CORE_HAMILTONIAN.block(0,0,nUncontracted,nUncontracted).real() = P2_Potential;
  CORE_HAMILTONIAN.block(
      nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = 
    P2_Potential;

  CORE_HAMILTONIAN.block(
      2*nUncontracted,2*nUncontracted,2*nUncontracted,2*nUncontracted)
    = W;
  CORE_HAMILTONIAN.block(
      2*nUncontracted,2*nUncontracted,2*nUncontracted,2*nUncontracted)
    -= 2*phys.SPEED_OF_LIGHT*phys.SPEED_OF_LIGHT * 
    ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted);


  CORE_HAMILTONIAN.block(
      0,2*nUncontracted,nUncontracted,nUncontracted).real()
    = phys.SPEED_OF_LIGHT * PMap.asDiagonal();
  CORE_HAMILTONIAN.block(
      nUncontracted,3*nUncontracted,nUncontracted,nUncontracted).real()
    = phys.SPEED_OF_LIGHT * PMap.asDiagonal();
  CORE_HAMILTONIAN.block(
      2*nUncontracted,0,2*nUncontracted,2*nUncontracted)
    = CORE_HAMILTONIAN.block(
      0,2*nUncontracted,2*nUncontracted,2*nUncontracted);

  // ------------------------------
  // Diagonalize 
  // ------------------------------

  Eigen::SelfAdjointEigenSolver<ComplexMatrix> es;
  es.compute(CORE_HAMILTONIAN);
  
  RealVecMap HEV(this->memManager_->malloc<double>(4*nUncontracted),4*nUncontracted);
  ComplexMap HEVx(this->memManager_->malloc<dcomplex>(16*nUnSq),
    4*nUncontracted,4*nUncontracted);
  
  HEV = es.eigenvalues();
  HEVx = es.eigenvectors();

  // Now free CORE_HAMILTONIAN after diagonalization
  this->memManager_->free(CORE_HAMILTONIAN.data(),16*nUnSq);

  // Print out the energies (eigenvalues) and eigenvectors
  if(this->printLevel_ >= 2){
    prettyPrintSmart(this->fileio_->out,HEV,"HEV");
    prettyPrintSmart(this->fileio_->out,HEVx,"HEVc");
  }
  
  // Grab C_L (+) and C_S (+) - the large and small components
  // of the electronic (positive energy) solutions
  //
  ComplexMap L(this->memManager_->malloc<dcomplex>(4*nUnSq),
    2*nUncontracted,2*nUncontracted);
  ComplexMap S(this->memManager_->malloc<dcomplex>(4*nUnSq),
    2*nUncontracted,2*nUncontracted);
  L = HEVx.block(0,2*nUncontracted,2*nUncontracted,2*nUncontracted);
  S = HEVx.block(2*nUncontracted,2*nUncontracted,2*nUncontracted,2*nUncontracted);

  this->memManager_->free(HEV.data(),4*nUncontracted);
  this->memManager_->free(HEVx.data(),16*nUnSq);

  /* Maybe get inverse from SVD?
  Eigen::JacobiSVD<ComplexMatrix> 
    svd(L,Eigen::ComputeThinU | Eigen::ComputeThinV);

  VectorXd SigmaL = svd.singularValues();
  ComplexMatrix SVL = svd.matrixU();
  */

  ComplexMap X(this->memManager_->malloc<dcomplex>(4*nUnSq),
    2*nUncontracted,2*nUncontracted);
  X = S * L.inverse();

  this->memManager_->free(L.data(),4*nUnSq);
  this->memManager_->free(S.data(),4*nUnSq);

  // Print out X and its squared norm
  if(this->printLevel_ >= 4){
    prettyPrintSmart(this->fileio_->out,X,"X");
    (this->fileio_->out) << X.squaredNorm() << endl;
  }
  
  // Calculate Y = sqrt(1 + X'X)
  // Also known as the 'renormalization matrix' R
  ComplexMap Y(this->memManager_->malloc<dcomplex>(4*nUnSq),
     2*nUncontracted,2*nUncontracted);
  Y = (ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) 
     + X.adjoint() * X).pow(-0.5);

  // Print out Y and its squared norm
  if(this->printLevel_ >= 4){
    prettyPrintSmart(this->fileio_->out,Y,"Y");
    (this->fileio_->out) << Y.squaredNorm() << endl;
  }

  // Get the momentum p = PMapC (as 2n by 2n matrix)
  ComplexMap PMapC(this->memManager_->malloc<dcomplex>(4*nUnSq),2*nUncontracted,2*nUncontracted);
  PMapC.setZero(); //Important to zero out the matrix first!

  PMapC.block(0,0,nUncontracted,nUncontracted).real() = PMap.asDiagonal();
  PMapC.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = PMap.asDiagonal();

  /* ADDITIONAL DEBUG
  // Compute p^2 and then the relativistic kinetic energy
  // Here we have       __________________
  //              T = \/m^2c^4 + c^2 * p^2  - mc^2
  //
  // Where m is the rest mass of the electron (1 in atomic units)             
  //
  ComplexMap KinEn(this->memManager_->malloc<dcomplex>(4*nUnSq),
    2*nUncontracted,2*nUncontracted);

  C4nUnSqScratch = PMapC.cwiseProduct(PMapC);
  KinEn = C4nUnSqScratch * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT;
  KinEn += ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted)* phys.SPEED_OF_LIGHT *
	phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT; 
  KinEn = KinEn.cwiseSqrt();
  KinEn -= ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT;

  if(this->printLevel_ >= 2){
    prettyPrintSmart(this->fileio_->out,KinEn,"Relativistic Kinetic Energy");
  }
  */

  // Get V' (in p-space)
  ComplexMap HCore2C(this->memManager_->malloc<dcomplex>(4*nUnSq),
    2*nUncontracted,2*nUncontracted);
  HCore2C.setZero();
  HCore2C.block(0,0,nUncontracted,nUncontracted).real() = P2_Potential;
  HCore2C.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = P2_Potential;

  if(this->printLevel_ >= 2){
    prettyPrintSmart(this->fileio_->out,HCore2C,"V prime (p space)");
  } 
  // Calculate the 2-component core Hamiltonian in the uncontracted basis
  //  αα | αβ
  //  -------
  //  βα | ββ  

  C4nUnSqScratch = phys.SPEED_OF_LIGHT * PMapC * X;
  HCore2C = HCore2C + C4nUnSqScratch;
  C4nUnSqScratch = phys.SPEED_OF_LIGHT * X.adjoint() * PMapC;
  HCore2C = HCore2C + C4nUnSqScratch;
  C4nUnSqScratch = 2 * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT * 
	ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted);
  C4nUnSqScratch = W - C4nUnSqScratch;
  C4nUnSqScratch = X.adjoint() * C4nUnSqScratch;
  C4nUnSqScratch = C4nUnSqScratch * X;
  HCore2C = HCore2C + C4nUnSqScratch;
  HCore2C = HCore2C * Y;
  HCore2C = Y * HCore2C;
  
  this->memManager_->free(PMapC.data(),4*nUnSq);
  this->memManager_->free(X.data(),4*nUnSq);
  this->memManager_->free(Y.data(),4*nUnSq);

  if(this->printLevel_ >= 3){
    prettyPrintSmart(this->fileio_->out,HCore2C,"Transformed 2c-Hamilotnian (in p space) ");
  }

  /*ADDITIONAL DEBUG
  ComplexMap Veff(this->memManager_->malloc<dcomplex>(4*nUnSq),
  2*nUncontracted,2*nUncontracted);
  Veff = HCore2C - KinEn;

  if(this->printLevel_ >= 3){
    prettyPrintSmart(this->fileio_->out,Veff,"V effective (p space)");
  }
  */

  RealMap Hs(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap Hx(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap Hy(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap Hz(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);

  Hs = 0.5 * (HCore2C.block(0,0,nUncontracted,nUncontracted).real()
	+ HCore2C.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real());
  Hz = 0.5 * (HCore2C.block(0,0,nUncontracted,nUncontracted).imag()
	- HCore2C.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).imag());
  Hx = 0.5 * (HCore2C.block(0,nUncontracted,nUncontracted,nUncontracted).imag()
	+ HCore2C.block(nUncontracted,0,nUncontracted,nUncontracted).imag());
  Hy = 0.5 * (HCore2C.block(0,nUncontracted,nUncontracted,nUncontracted).real()
	- HCore2C.block(nUncontracted,0,nUncontracted,nUncontracted).real());

  this->memManager_->free(HCore2C.data(),4*nUnSq);

  if(this->printLevel_ >= 3){
    prettyPrintSmart(this->fileio_->out,Hs,"Hs (p space)");
    prettyPrintSmart(this->fileio_->out,Hz,"Hz (p space)");
    prettyPrintSmart(this->fileio_->out,Hx,"Hx (p space)");
    prettyPrintSmart(this->fileio_->out,Hy,"Hz (p space)");
  }

  nUnSqScratch = Hs * UK.adjoint() * SUn;
  Hs = SUn * UK * nUnSqScratch; 
  nUnSqScratch = Hz * UK.adjoint() * SUn;
  Hz = SUn * UK * nUnSqScratch; 
  nUnSqScratch = Hx * UK.adjoint() * SUn;
  Hx = SUn * UK * nUnSqScratch; 
  nUnSqScratch = Hy * UK.adjoint() * SUn;
  Hy = SUn * UK * nUnSqScratch; 

  if(this->printLevel_ >= 3){
    prettyPrintSmart(this->fileio_->out,Hs,"Hs (r space)");
    prettyPrintSmart(this->fileio_->out,Hz,"Hz (r space)");
    prettyPrintSmart(this->fileio_->out,Hx,"Hx (r space)");
    prettyPrintSmart(this->fileio_->out,Hy,"Hy (r space)");
  }

  /*ADDITIONAL DEBUG HERE
  if(this->printLevel_ >= 3){
    RealMatrix TCon = 0.5 * (KinEn.block(0,0,nUncontracted,nUncontracted).real() + 
      KinEn.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real());
    RealMatrix VCon = 0.5 * (Veff.block(0,0,nUncontracted,nUncontracted).real() + 
      Veff.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real());

    prettyPrintSmart(this->fileio_->out,TCon,"T relativistic (p space)");
    prettyPrintSmart(this->fileio_->out,VCon,"V relativistic (p space)");
    
    nUnSqScratch = TCon * UK.adjoint() * SUn;
    TCon = SUn * UK * nUnSqScratch;
    nUnSqScratch = VCon * UK.adjoint() * SUn;
    VCon = SUn * UK * nUnSqScratch;

    prettyPrintSmart(this->fileio_->out,TCon,"T relativistic (r space)");
    prettyPrintSmart(this->fileio_->out,VCon,"V relativistic (r space)");

    TCon = (*this->basisSet_->mapPrim2Bf()) * TCon * (*this->basisSet_->mapPrim2Bf()).transpose();;
    VCon = (*this->basisSet_->mapPrim2Bf()) * VCon * (*this->basisSet_->mapPrim2Bf()).transpose();;

    prettyPrintSmart(this->fileio_->out,TCon,"T (contracted basis)");
    prettyPrintSmart(this->fileio_->out,VCon,"V (contracted basis)");
  }
  
  this->memManager_->free(Veff.data(),4*nUnSq);
  */

  // Recontract the basis

  RealMap CoreS(this->memManager_->malloc<double>(this->nBasis_ * this->nBasis_),
    this->nBasis_,this->nBasis_);
  RealMap CoreX(this->memManager_->malloc<double>(this->nBasis_ * this->nBasis_),
    this->nBasis_,this->nBasis_);
  RealMap CoreY(this->memManager_->malloc<double>(this->nBasis_ * this->nBasis_),
    this->nBasis_,this->nBasis_);
  RealMap CoreZ(this->memManager_->malloc<double>(this->nBasis_ * this->nBasis_),
    this->nBasis_,this->nBasis_);

  CoreS = (*this->basisSet_->mapPrim2Bf()) * Hs * (*this->basisSet_->mapPrim2Bf()).transpose();
  CoreZ = (*this->basisSet_->mapPrim2Bf()) * Hz * (*this->basisSet_->mapPrim2Bf()).transpose();
  CoreX = (*this->basisSet_->mapPrim2Bf()) * Hx * (*this->basisSet_->mapPrim2Bf()).transpose();
  CoreY = (*this->basisSet_->mapPrim2Bf()) * Hy * (*this->basisSet_->mapPrim2Bf()).transpose();

  this->memManager_->free(Hs.data(),nUnSq);
  this->memManager_->free(Hx.data(),nUnSq);
  this->memManager_->free(Hy.data(),nUnSq);
  this->memManager_->free(Hz.data(),nUnSq);

  if(this->twoEFudge == 1){
    if (this->printLevel_ >= 2){
        (this->fileio_->out) << "Applying 2e fudge factor to spin-orbit Hamiltonian" << endl;
    }
    // Q(l) as an function:
    //    double QofL(double L) { return L*(L+1)*(2*L+1)/3; };
    
    // Loop over all shells
    for (auto iShell = 0; iShell < this->basisSet_->nShell(); ++iShell) {
    for (auto jShell = 0; jShell < this->basisSet_->nShell(); ++jShell) {
  
      auto iL = this->basisSet_->shells(iShell).contr[0].l;
      auto jL = this->basisSet_->shells(jShell).contr[0].l;
      // Calculate Fudge Factor term by Term
      double fudgeFactor = (iL)*(iL+1)*(2*iL+1)/3;
      fudgeFactor *= (jL)*(jL+1)*(2*jL+1)/3;
      fudgeFactor /= this->molecule_->atomicZ(this->basisSet_->mapSh2Cen(iShell)-1); 
      fudgeFactor /= this->molecule_->atomicZ(this->basisSet_->mapSh2Cen(jShell)-1); 
      fudgeFactor = sqrt(fudgeFactor);
  
      CoreX.block(this->basisSet_->mapSh2Bf(iShell),
      this->basisSet_->mapSh2Bf(jShell),this->basisSet_->shells(iShell).size(),
      this->basisSet_->shells(jShell).size()) -=
       (CoreX.block(this->basisSet_->mapSh2Bf(iShell),
       this->basisSet_->mapSh2Bf(jShell),this->basisSet_->shells(iShell).size(),
       this->basisSet_->shells(jShell).size()) * fudgeFactor);     
     
      CoreY.block(this->basisSet_->mapSh2Bf(iShell),
      this->basisSet_->mapSh2Bf(jShell),this->basisSet_->shells(iShell).size(),
      this->basisSet_->shells(jShell).size()) -=
       (CoreY.block(this->basisSet_->mapSh2Bf(iShell),
       this->basisSet_->mapSh2Bf(jShell),this->basisSet_->shells(iShell).size(),
       this->basisSet_->shells(jShell).size()) * fudgeFactor);
    
      CoreZ.block(this->basisSet_->mapSh2Bf(iShell),
      this->basisSet_->mapSh2Bf(jShell),this->basisSet_->shells(iShell).size(),
      this->basisSet_->shells(jShell).size()) -=
       (CoreZ.block(this->basisSet_->mapSh2Bf(iShell),
       this->basisSet_->mapSh2Bf(jShell),this->basisSet_->shells(iShell).size(),
       this->basisSet_->shells(jShell).size()) * fudgeFactor);
    }
    }//loop shells
  } //twoE Fudge
  
  if(this->printLevel_ >= 2) {
    prettyPrintSmart(this->fileio_->out,CoreS,"X2C core Hamiltonian (scalar)");
    prettyPrintSmart(this->fileio_->out,CoreZ,"X2C core Hamiltonian (mz)");
    prettyPrintSmart(this->fileio_->out,CoreX,"X2C core Hamiltonian (mx)");
    prettyPrintSmart(this->fileio_->out,CoreY,"X2C core Hamiltonian (my)");
  }

  *this->coreH_ = CoreS;
  *this->oneEmx_ = CoreX;
  *this->oneEmy_ = CoreY;
  *this->oneEmz_ = CoreZ;

  // Free all memory before return
  this->memManager_->free(CoreS.data(),this->nBasis_*this->nBasis_);
  this->memManager_->free(CoreX.data(),this->nBasis_*this->nBasis_);
  this->memManager_->free(CoreY.data(),this->nBasis_*this->nBasis_);
  this->memManager_->free(CoreZ.data(),this->nBasis_*this->nBasis_);
  
  this->memManager_->free(nUnSqScratch.data(),nUnSq);
  this->memManager_->free(C4nUnSqScratch.data(),4*nUnSq);

  this->memManager_->free(SUn.data(),nUnSq);
  this->memManager_->free(UK.data(),nUnSq);
  
  // Might be able to free these earlier
  this->memManager_->free(P2_Potential.data(),nUnSq);
  this->memManager_->free(W.data(),4*nUnSq);

  //NOT CURRENTLY USED
  //this->memManager_->free(KinEn.data(),4*nUnSq);


/*
//ADDITIONAL DEBUG HERE
if(this->printLevel_ >= 4){
  ComplexMatrix TCSham(2*nBasis_,2*nBasis_);
  
  TCSham.block(0,0,nBasis_,nBasis_).real() = CoreS;
  TCSham.block(nBasis_,nBasis_,nBasis_,nBasis_).real() = CoreS;
  TCSham.block(0,0,nBasis_,nBasis_).imag() = CoreZ;
  TCSham.block(nBasis_,nBasis_,nBasis_,nBasis_).imag() = -1*CoreZ;
  TCSham.block(0,nBasis_,nBasis_,nBasis_).imag() = CoreX;
  TCSham.block(nBasis_,0,nBasis_,nBasis_).imag() = CoreX;
  TCSham.block(0,nBasis_,nBasis_,nBasis_).real() = CoreY;
  TCSham.block(nBasis_,0,nBasis_,nBasis_).real() = -1*CoreY;
 
  prettyPrintSmart(this->fileio_->out,TCSham.real(),"Two component Hamiltonian (real)");
  prettyPrintSmart(this->fileio_->out,TCSham.imag(),"Two component Hamiltonian (imag)");

  // put spin as fastest running index
  ComplexMatrix SpinHam(2*nBasis_,2*nBasis_);

  for (auto row = 0; row < nBasis_; row++){
    for (auto col = 0; col < nBasis_; col++){
        SpinHam(2*row,2*col) = TCSham(row,col);
        SpinHam(2*row,2*col+1) = TCSham(row,col+nBasis_);
        SpinHam(2*row+1,2*col) = TCSham(row+nBasis_,col);
        SpinHam(2*row+1,2*col+1) = TCSham(row+nBasis_,col+nBasis_);
        }
    }  
  prettyPrintSmart(this->fileio_->out,SpinHam.real(),"Spin-Blocked 2c-Hamiltonian (real)");
  prettyPrintSmart(this->fileio_->out,SpinHam.imag(),"Spin-Blocked 2c-Hamiltonian (imag)");
  }
//END OF DEBUG
*/
}

