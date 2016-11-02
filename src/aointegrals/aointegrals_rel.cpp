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
  int nthreads = omp_get_max_threads();
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
  RealMap nUnSqScratch2(
    this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  ComplexMap C4nUnSqScratch(
    this->memManager_->malloc<dcomplex>(4*nUnSq),2*nUncontracted,2*nUncontracted);
  ComplexMap C4nUnSqScratch2(
    this->memManager_->malloc<dcomplex>(4*nUnSq),2*nUncontracted,2*nUncontracted);

  double * SCRUnCon = this->memManager_->malloc<double>(nUncontracted*this->nBasis_);
  RealMap SCRUnConMap(SCRUnCon,nUncontracted,this->nBasis_);
  RealMap SCRConUnMap(SCRUnCon,this->nBasis_,nUncontracted);
  

  libint2::Engine engineS(
      libint2::Operator::overlap,1,this->basisSet_->maxL(),0);
  libint2::Engine engineT(
      libint2::Operator::kinetic,1,this->basisSet_->maxL(),0);
  libint2::Engine engineC(
      libint2::Operator::coulomb,1,this->basisSet_->maxL(),0);


  engineS.set_precision(0.0);
  engineT.set_precision(0.0);
  engineC.set_precision(0.0);

  std::vector<libint2::Engine> engineV(nthreads);
  engineV[0] = libint2::Engine(libint2::Operator::nuclear,1,this->basisSet_->maxL(),0);
  engineV[0].set_precision(0.0);

  for(auto i = 1; i < nthreads; i++) engineV[i] = engineV[0];

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

//SnonRel = (*this->basisSet_->mapPrim2Bf()) * SUncontracted
//  * (*this->basisSet_->mapPrim2Bf()).transpose();
  SCRConUnMap.noalias() = (*this->basisSet_->mapPrim2Bf()) * SUncontracted;	
  SnonRel.noalias() = SCRConUnMap * (*this->basisSet_->mapPrim2Bf()).transpose();

  for (auto row = 0; row < this->basisSet_->nBasis(); row++){
      (*this->basisSet_->mapPrim2Bf()).block(row,0,1,nUncontracted) /=
        std::sqrt(SnonRel(row,row)); //scale by appropriate factor
    }
  this->memManager_->free(SnonRel.data(),this->nBasis_*this->nBasis_);
  // -----------------------------------------------------------------------

  RealMap SUn(this->memManager_->malloc<double>(nUnSq),
    nUncontracted,nUncontracted);
  SUn = SUncontracted; // Save S for later

  char JOBU = 'O', JOBVT = 'N';
  int LWORK = 10*nUncontracted;
  int INFO;
//std::vector<double> ovlpEigValues(nUncontracted);
//std::vector<double> WORK(LWORK);
  double * ovlpEigValues = 
    this->memManager_->malloc<double>(nUncontracted);
  double * WORK = this->memManager_->malloc<double>(LWORK);

  // Get eigenvalues for overlap matrix S (via SVD) 
  // Store the orthogonal transformation U back in S

  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,SUncontracted.data(),
      &nUncontracted,ovlpEigValues,SUncontracted.data(),&nUncontracted,
      SUncontracted.data(),&nUncontracted,WORK,&LWORK,&INFO);

  // Count linear dependencies here
  int nZero = 0;
  for(auto iS = 0; iS < nUncontracted; iS++)
    if(std::abs(ovlpEigValues[iS]) < 1e-12) nZero++;

//cout << "NZERO " << nZero << endl;
  
  // Normalize
  // What happens when we divide by zero?
  for(auto iS = 0; iS < nUncontracted; iS++){
    SUncontracted.col(iS) /= std::sqrt(ovlpEigValues[iS]);
  }

  // Put the kinetic energy in the orthonormal basis 

  nUnSqScratch.noalias() = SUncontracted.transpose() * TUncontracted;
  TUncontracted.noalias() = nUnSqScratch * SUncontracted;

  // Get rid of the linear dependencies?
  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,TUncontracted.data(),
      &nUncontracted,ovlpEigValues,TUncontracted.data(),&nUncontracted,
      TUncontracted.data(),&nUncontracted,WORK,&LWORK,&INFO);


  // Deallocate workspace
  this->memManager_->free(WORK,LWORK);

  // Form the K transformation matrix from both pieces
  // (diagonalizes S) and (diagonalizes T). See eq. 12 in Rieher's paper from 2013

  if(this->printLevel_ >= 3){
    prettyPrintSmart(this->fileio_->out,SUncontracted,"S transformation");
    prettyPrintSmart(this->fileio_->out,TUncontracted,"T transformation");
    prettyPrintSmart(this->fileio_->out,TUncontracted.transpose()*TUncontracted, "T^{dagger} T");
  }

  RealMap UK(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  UK.noalias() = SUncontracted * TUncontracted;

  if(this->printLevel_ >= 3){ 
    prettyPrintSmart(this->fileio_->out,UK,"Uk transformation");
  }

  // Now we transform V to V' 
  RealMap P2_Potential(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);

  nUnSqScratch.noalias() = UK.transpose() * VUncontracted;
  P2_Potential.noalias() = nUnSqScratch * UK;

  this->memManager_->free(SUncontracted.data(),nUnSq);
  this->memManager_->free(TUncontracted.data(),nUnSq);
  this->memManager_->free(VUncontracted.data(),nUnSq);
  // Next we need to get W' (from pVp integrals)

  std::vector<RealVecMap> SCRATCHDXUnContracted;
  std::vector<RealVecMap> SCRATCHDYUnContracted;
  std::vector<RealVecMap> SCRATCHDZUnContracted;

  for(auto i = 0; i < nthreads; i++) {
    SCRATCHDXUnContracted.emplace_back(this->memManager_->malloc<double>(nUncontracted),nUncontracted);
    SCRATCHDYUnContracted.emplace_back(this->memManager_->malloc<double>(nUncontracted),nUncontracted);
    SCRATCHDZUnContracted.emplace_back(this->memManager_->malloc<double>(nUncontracted),nUncontracted);

    SCRATCHDXUnContracted.back().setZero();
    SCRATCHDYUnContracted.back().setZero();
    SCRATCHDZUnContracted.back().setZero();
  }
  // compute pVp numerically
  auto PVP = [&](IntegrationPoint pt, std::vector<std::vector<RealMap>> &result) {
    int thread_id = omp_get_thread_num();
    for(auto iShell = 0, b_s = 0; iShell < unContractedShells.size();
         b_s += unContractedShells[iShell].size(),++iShell) {
      int size= unContractedShells[iShell].size();

      double * 
        buff=this->basisSet_->basisDEval(1,unContractedShells[iShell], &pt.pt);

      RealMap bMap( buff         ,size,1);
      RealMap dxMap(buff + size  ,size,1);
      RealMap dyMap(buff + 2*size,size,1);
      RealMap dzMap(buff + 3*size,size,1);

//      SCRATCH1UnContracted[thread_id].block( b_s,0,size,1) = bMap;
      SCRATCHDXUnContracted[thread_id].block(b_s,0,size,1) = dxMap;
      SCRATCHDYUnContracted[thread_id].block(b_s,0,size,1) = dyMap;
      SCRATCHDZUnContracted[thread_id].block(b_s,0,size,1) = dzMap;

     // delete [] buff;
    };

    std::vector<std::pair<double,std::array<double,3>>> q;
    q.push_back(
      {1.0, {{bg::get<0>(pt.pt),bg::get<1>(pt.pt),bg::get<2>(pt.pt)}}});
    engineV[thread_id].set_params(q);

    for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
      // Gaussian Nuclei
      //
/*
      const double * gamma = engineV.compute(this->molecule_->nucShell(iAtm),
          libint2::Shell::unit());
    //SCRATCH2UnContracted.noalias() = 
    //  SCRATCH1UnContracted * SCRATCH1UnContracted.transpose();
    //result[0].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

      SCRATCH2UnContracted.noalias() = 
        SCRATCHDXUnContracted * SCRATCHDXUnContracted.transpose();
      result[1].noalias() += (pt.weight * (*gamma)) *  SCRATCH2UnContracted;

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
      result[5].noalias() += (pt.weight * (*gamma))* SCRATCH2UnContracted;

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
*/
//New AP
    const double * gamma = engineV[thread_id].compute(this->molecule_->nucShell(iAtm),
          libint2::Shell::unit());
    double gamw = (*gamma) * pt.weight;
// Screening ????
//    if (std::abs(gamw) < 1.e-13) continue ;
    auto ShUnSize = unContractedShells.size();
    auto iSt = 0;
    for(auto iShell=0; iShell < ShUnSize; iShell++){ 
      int iSz= unContractedShells[iShell].size();
      for(auto iBf = iSt; iBf < (iSt + iSz); iBf++){
    auto jSt = 0;
    for(auto jShell=0; jShell < ShUnSize; jShell++){ 
          int jSz= unContractedShells[jShell].size();
          double *pVpSDATA = result[thread_id][0].data() + iBf*nUncontracted;
          double *pVpXDATA = result[thread_id][1].data() + iBf*nUncontracted;
          double *pVpYDATA = result[thread_id][2].data() + iBf*nUncontracted;
          double *pVpZDATA = result[thread_id][3].data() + iBf*nUncontracted;
          double *SCRATCHDXUnContractedData = SCRATCHDXUnContracted[thread_id].data();
          double *SCRATCHDYUnContractedData = SCRATCHDYUnContracted[thread_id].data();
          double *SCRATCHDZUnContractedData = SCRATCHDZUnContracted[thread_id].data();
          for(auto jBf = jSt; jBf < (jSt + jSz); jBf++){
            if(jBf < iBf) continue;
// PvP Scalar
           pVpSDATA[jBf] += gamw * SCRATCHDXUnContractedData[iBf] * SCRATCHDXUnContractedData[jBf];
           pVpSDATA[jBf] += gamw * SCRATCHDYUnContractedData[iBf] * SCRATCHDYUnContractedData[jBf];
	   pVpSDATA[jBf] += gamw * SCRATCHDZUnContractedData[iBf] * SCRATCHDZUnContractedData[jBf];
// PvpX
	   pVpXDATA[jBf] -= gamw * SCRATCHDYUnContractedData[iBf] * SCRATCHDZUnContractedData[jBf];
           pVpXDATA[jBf] += gamw * SCRATCHDZUnContractedData[iBf] * SCRATCHDYUnContractedData[jBf];
// PvpY
	   pVpYDATA[jBf] -= gamw * SCRATCHDZUnContractedData[iBf] * SCRATCHDXUnContractedData[jBf];
           pVpYDATA[jBf] += gamw * SCRATCHDXUnContractedData[iBf] * SCRATCHDZUnContractedData[jBf];
// PvpZ
	   pVpZDATA[jBf] -= gamw * SCRATCHDXUnContractedData[iBf] * SCRATCHDYUnContractedData[jBf];
           pVpZDATA[jBf] += gamw * SCRATCHDYUnContractedData[iBf] * SCRATCHDXUnContractedData[jBf];
          } // jBf
        jSt += unContractedShells[jShell].size();
        } // jShell
      } // iBf
     iSt += unContractedShells[iShell].size() ;
    } // iShell

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
/*
  std::vector<RealMatrix> numPot(10,RealMatrix::Zero(nUncontracted,
        nUncontracted));
*/
//  std::vector<RealMatrix> numPot(4,RealMatrix::Zero(nUncontracted,
//        nUncontracted));
  std::vector<std::vector<RealMap> > numPot(nthreads);
  for(auto i = 0; i < nthreads; i++) {
    for(auto j = 0; j < 4; j++) {
      numPot[i].emplace_back(this->memManager_->malloc<double>(nUnSq),
        nUncontracted,nUncontracted);
      numPot[i].back().setZero();
    }
  }

  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    cout << iAtm << endl;
    AGrid.center() = iAtm;
    AGrid.setRadCutOff(atomRadCutoff[iAtm]);
    AGrid.findNearestNeighbor();
    AGrid.scalingFactor()=
      0.5*elements[this->molecule_->index(iAtm)].sradius/phys.bohr;
    AGrid.integrate<std::vector<std::vector<RealMap>>>(PVP,numPot);
  };

  for(auto i = 0; i < nthreads; i++) {
    this->memManager_->free(SCRATCHDXUnContracted[i].data(),nUncontracted);
    this->memManager_->free(SCRATCHDYUnContracted[i].data(),nUncontracted);
    this->memManager_->free(SCRATCHDZUnContracted[i].data(),nUncontracted);
  }
  
/*
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
*/

  RealMap PVPS(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap PVPX(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap PVPY(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  RealMap PVPZ(this->memManager_->malloc<double>(nUnSq),nUncontracted,nUncontracted);
  PVPS.setZero();
  PVPX.setZero();
  PVPY.setZero();
  PVPZ.setZero();
   for(auto i = 0; i < nthreads; i++) {
     PVPS += numPot[i][0];
     PVPX += numPot[i][1];
     PVPY += numPot[i][2];
     PVPZ += numPot[i][3];
   }
   double fact4pi = 4 * math.pi; 
   PVPS.triangularView<Upper>() =  PVPS.transpose(); // Sym
   PVPX.triangularView<Upper>() = -PVPX.transpose(); // AntiSymm
   PVPY.triangularView<Upper>() = -PVPY.transpose(); // AntiSymm
   PVPZ.triangularView<Upper>() = -PVPZ.transpose(); // AntiSymm
   PVPS *= fact4pi;
   PVPX *= fact4pi;
   PVPY *= fact4pi;
   PVPZ *= fact4pi;


  // Apply the Uk unitary transformation ( Uk' * PVP * Uk)
  nUnSqScratch.noalias() = UK.adjoint() * PVPS;
  PVPS.noalias() = nUnSqScratch * UK;
  nUnSqScratch.noalias() = UK.adjoint() * PVPX;
  PVPX.noalias() = nUnSqScratch * UK;
  nUnSqScratch.noalias() = UK.adjoint() * PVPY;
  PVPY.noalias() = nUnSqScratch * UK;
  nUnSqScratch.noalias() = UK.adjoint() * PVPZ;
  PVPZ.noalias() = nUnSqScratch * UK;


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

  RealVecMap PMap(ovlpEigValues,nUncontracted);

  if(this->printLevel_ >= 3){
    prettyPrintSmart(cout,PMap,"P^2");
  }

  PMap = 2*PMap;
  PMap = PMap.cwiseSqrt();
  PMap = PMap.cwiseInverse();

  nUnSqScratch.noalias() = PMap.asDiagonal() * PVPS;
  PVPS.noalias() = nUnSqScratch * PMap.asDiagonal();
  nUnSqScratch.noalias() = PMap.asDiagonal() * PVPX;
  PVPX.noalias() = nUnSqScratch * PMap.asDiagonal();
  nUnSqScratch.noalias() = PMap.asDiagonal() * PVPY;
  PVPY.noalias() = nUnSqScratch * PMap.asDiagonal();
  nUnSqScratch.noalias() = PMap.asDiagonal() * PVPZ;
  PVPZ.noalias() = nUnSqScratch * PMap.asDiagonal();

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

  RealVecMap HEV(this->memManager_->malloc<double>(4*nUncontracted),4*nUncontracted);
  ComplexMap HEVx(this->memManager_->malloc<dcomplex>(16*nUnSq),
    4*nUncontracted,4*nUncontracted);

  char JOBZ = 'V';
  char UPLO = 'U';
  int FourUnCon = 4*nUncontracted;
  int LRWORK = 3*FourUnCon - 2;
  dcomplex * CWORK = this->memManager_->malloc<dcomplex>(LWORK);
  double   * RWORK = this->memManager_->malloc<double>(LRWORK);
  zheev_(&JOBZ,&UPLO,&FourUnCon,CORE_HAMILTONIAN.data(),&FourUnCon,HEV.data(),CWORK,&LWORK,RWORK,&INFO);
 
  // FIXME: This copy is not needed and HEVx is not needed as a result
  HEVx = CORE_HAMILTONIAN;

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

  // Replaces L with L^{-1}
  int TwoUnCon = 2 * nUncontracted;
  std::vector<int> iPiv(TwoUnCon);
  zgetrf_(&TwoUnCon,&TwoUnCon,L.data(),&TwoUnCon,&iPiv[0],&INFO);
  zgetri_(&TwoUnCon,L.data(),&TwoUnCon,&iPiv[0],CWORK,&LWORK,&INFO);


  X.noalias() = S * L;

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

  // FIXME: (DBWY) Need to compute the sqrt smartly here. I'm 90% sure that
  // you can use an SVD here as it constitues the "symmetric / Lowdin" orthogonalization
  // step here...
  Y = (ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) 
     + X.adjoint() * X).pow(-0.5);

/*
  C4nUnSqScratch  = (ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) + X.adjoint() * X);

  prettyPrintSmart(cout,C4nUnSqScratch.pow(-0.5),"True");
  double * SingVal = this->memManager_->malloc<double>(TwoUnCon);
  JOBU = 'O';JOBVT = 'A';
  zgesvd_(&JOBU,&JOBVT,&TwoUnCon,&TwoUnCon,C4nUnSqScratch.data(),&TwoUnCon,SingVal,C4nUnSqScratch.data(),&TwoUnCon,C4nUnSqScratch2.data(),&TwoUnCon,CWORK,&LWORK,RWORK,&INFO);

  Y = C4nUnSqScratch * C4nUnSqScratch2;
  prettyPrintSmart(cout,Y,"Attempt");
  CErr();
*/

  // Free up CWORK and RWORK
  this->memManager_->free(CWORK,LWORK);
  this->memManager_->free(RWORK,LRWORK);

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

  // FIXME: (DBWY) Can you put some docs here that explain the GEMM calls?
  C4nUnSqScratch.noalias() = phys.SPEED_OF_LIGHT * PMapC * X;
  HCore2C.noalias() += C4nUnSqScratch;
  C4nUnSqScratch.noalias() = phys.SPEED_OF_LIGHT * X.adjoint() * PMapC;
  HCore2C.noalias() += C4nUnSqScratch;
  C4nUnSqScratch.noalias() = 2 * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT * 
	ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted);
  C4nUnSqScratch2.noalias() = W - C4nUnSqScratch;
  C4nUnSqScratch.noalias() = X.adjoint() * C4nUnSqScratch2;
  C4nUnSqScratch2.noalias() = C4nUnSqScratch * X;
  HCore2C.noalias() += C4nUnSqScratch2;
  C4nUnSqScratch.noalias() = HCore2C * Y;
  HCore2C.noalias() = Y * C4nUnSqScratch;
  
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

  Hs.noalias() = 0.5 * (HCore2C.block(0,0,nUncontracted,nUncontracted).real()
	+ HCore2C.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real());
  Hz.noalias() = 0.5 * (HCore2C.block(0,0,nUncontracted,nUncontracted).imag()
	- HCore2C.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).imag());
  Hx.noalias() = 0.5 * (HCore2C.block(0,nUncontracted,nUncontracted,nUncontracted).imag()
	+ HCore2C.block(nUncontracted,0,nUncontracted,nUncontracted).imag());
  Hy.noalias() = 0.5 * (HCore2C.block(0,nUncontracted,nUncontracted,nUncontracted).real()
	- HCore2C.block(nUncontracted,0,nUncontracted,nUncontracted).real());

  this->memManager_->free(HCore2C.data(),4*nUnSq);

  if(this->printLevel_ >= 3){
    prettyPrintSmart(this->fileio_->out,Hs,"Hs (p space)");
    prettyPrintSmart(this->fileio_->out,Hz,"Hz (p space)");
    prettyPrintSmart(this->fileio_->out,Hx,"Hx (p space)");
    prettyPrintSmart(this->fileio_->out,Hy,"Hz (p space)");
  }

/*
  nUnSqScratch = Hs * UK.adjoint() * SUn;
  Hs = SUn * UK * nUnSqScratch; 
  nUnSqScratch = Hz * UK.adjoint() * SUn;
  Hz = SUn * UK * nUnSqScratch; 
  nUnSqScratch = Hx * UK.adjoint() * SUn;
  Hx = SUn * UK * nUnSqScratch; 
  nUnSqScratch = Hy * UK.adjoint() * SUn;
  Hy = SUn * UK * nUnSqScratch; 
*/

  // FIXME?: (Possibly) we always compute SUn * UK, we could store it once
  nUnSqScratch.noalias() = Hs * UK.adjoint();
  nUnSqScratch2.noalias() = nUnSqScratch * SUn;
  nUnSqScratch.noalias() = SUn * UK;
  Hs.noalias() = nUnSqScratch * nUnSqScratch2; 

  nUnSqScratch.noalias() = Hz * UK.adjoint();
  nUnSqScratch2.noalias() = nUnSqScratch * SUn;
  nUnSqScratch.noalias() = SUn * UK;
  Hz.noalias() = nUnSqScratch * nUnSqScratch2; 

  nUnSqScratch.noalias() = Hx * UK.adjoint();
  nUnSqScratch2.noalias() = nUnSqScratch * SUn;
  nUnSqScratch.noalias() = SUn * UK;
  Hx.noalias() = nUnSqScratch * nUnSqScratch2; 

  nUnSqScratch.noalias() = Hy * UK.adjoint();
  nUnSqScratch2.noalias() = nUnSqScratch * SUn;
  nUnSqScratch.noalias() = SUn * UK;
  Hy.noalias() = nUnSqScratch * nUnSqScratch2; 

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

/*
  CoreS = (*this->basisSet_->mapPrim2Bf()) * Hs * (*this->basisSet_->mapPrim2Bf()).transpose();
  CoreZ = (*this->basisSet_->mapPrim2Bf()) * Hz * (*this->basisSet_->mapPrim2Bf()).transpose();
  CoreX = (*this->basisSet_->mapPrim2Bf()) * Hx * (*this->basisSet_->mapPrim2Bf()).transpose();
  CoreY = (*this->basisSet_->mapPrim2Bf()) * Hy * (*this->basisSet_->mapPrim2Bf()).transpose();
*/
  SCRConUnMap.noalias() = (*this->basisSet_->mapPrim2Bf()) * Hs;
  CoreS.noalias() = SCRConUnMap * (*this->basisSet_->mapPrim2Bf()).transpose();
  SCRConUnMap.noalias() = (*this->basisSet_->mapPrim2Bf()) * Hz;
  CoreZ.noalias() = SCRConUnMap * (*this->basisSet_->mapPrim2Bf()).transpose();
  SCRConUnMap.noalias() = (*this->basisSet_->mapPrim2Bf()) * Hx;
  CoreX.noalias() = SCRConUnMap * (*this->basisSet_->mapPrim2Bf()).transpose();
  SCRConUnMap.noalias() = (*this->basisSet_->mapPrim2Bf()) * Hy;
  CoreY.noalias() = SCRConUnMap * (*this->basisSet_->mapPrim2Bf()).transpose();

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
  this->memManager_->free(nUnSqScratch2.data(),nUnSq);
  this->memManager_->free(C4nUnSqScratch.data(),4*nUnSq);
  this->memManager_->free(C4nUnSqScratch2.data(),4*nUnSq);
  this->memManager_->free(SCRUnCon,nUncontracted*this->nBasis_);

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

