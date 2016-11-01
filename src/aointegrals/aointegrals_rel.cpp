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
  this->basisSet_->makeMapPrim2Bf();
  if(!this->isPrimary) return;
  auto unContractedShells = this->basisSet_->uncontractBasis();
  int nUncontracted = 0;
  for(auto i : unContractedShells) nUncontracted += i.size();

  RealMatrix SUncontracted(nUncontracted,nUncontracted);
  RealMatrix TUncontracted(nUncontracted,nUncontracted);
  RealMatrix VUncontracted(nUncontracted,nUncontracted);

if(this->printLevel_ >= 2){
  RealMatrix TCpy(*this->kinetic_);
  prettyPrintSmart(this->fileio_->out,TCpy,"T (non-rel)");
  prettyPrintSmart(this->fileio_->out,*this->overlap_,"overlap (non-rel)");
}

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

  cout << "Calculating uncontracted Libint Ints...";
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

  } // s2
  } // s1
  cout << "done!" << endl;

  SUncontracted = SUncontracted.selfadjointView<Lower>();
  TUncontracted = TUncontracted.selfadjointView<Lower>();

  RealMatrix SnonRel = (*this->basisSet_->mapPrim2Bf()) * SUncontracted
	* (*this->basisSet_->mapPrim2Bf()).transpose();
  if(this->printLevel_ >= 3){
  prettyPrintSmart(this->fileio_->out,SnonRel,"S nonRel (take 1)");
    }
  for (auto row = 0; row < this->basisSet_->nBasis(); row++){
      (*this->basisSet_->mapPrim2Bf()).block(row,0,1,nUncontracted) /=
        std::sqrt(SnonRel(row,row)); //scale by appropraite factor
    }

if(this->printLevel_ >= 2){
  SnonRel = (*this->basisSet_->mapPrim2Bf()) * SUncontracted
    * (*this->basisSet_->mapPrim2Bf()).transpose();
  prettyPrintSmart(this->fileio_->out,SnonRel,"S nonRel");
}

if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,SUncontracted,"S uncontracted");
  prettyPrintSmart(this->fileio_->out,TUncontracted,"T uncontracted");

  RealMatrix TnonRel = (*this->basisSet_->mapPrim2Bf()) * TUncontracted
	* (*this->basisSet_->mapPrim2Bf()).transpose();  
  prettyPrintSmart(this->fileio_->out,TnonRel,"T nonRel");
}

  RealMatrix SUn(nUncontracted,nUncontracted);
  SUn = SUncontracted.real(); // Save S for later

  char JOBU = 'O', JOBVT = 'N';
  int LWORK = 6*nUncontracted;
  int INFO;
  std::vector<double> ovlpEigValues(nUncontracted);
  std::vector<double> WORK(LWORK);

// Get eigenvalues for overlap matrix S (via SVD) 
// Store the orthogonal transformation U back in S

  cout << "Performing SVD on uncontracted overlap...";
  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,SUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],SUncontracted.data(),&nUncontracted,
      SUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);
  cout << "done!" << endl;

// normalize
// What happens when we divide by zero?
  for(auto iS = 0; iS < nUncontracted; iS++){
    SUncontracted.col(iS) /= std::sqrt(ovlpEigValues[iS]);
  }

// Count linear dependencies here?
  int nZero = 0;
  for(auto iS = 0; iS < nUncontracted; iS++)
    if(std::abs(ovlpEigValues[iS]) < 1e-12) nZero++;

  cout << "NZERO " << nZero << endl;

// Put the kinetic energy in the orthonormal basis 

  RealMatrix TMP = SUncontracted.transpose() * TUncontracted;
  TUncontracted = TMP * SUncontracted;

  cout << "Performing SVD on uncontracted T...";
// Get rid of the linear dependencies?
  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,TUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],TUncontracted.data(),&nUncontracted,
      TUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);

  cout << "done!" << endl;
// Form the K transformation matrix from both pieces
// (diagonalizes S) and (diagonalizes T). See eq. 12 in Rieher's paper from 2013

if(this->printLevel_ >= 3){
  prettyPrintSmart(this->fileio_->out,SUncontracted,"S transformation");
  prettyPrintSmart(this->fileio_->out,TUncontracted,"T transformation");
  prettyPrintSmart(this->fileio_->out,TUncontracted.transpose()*TUncontracted, "T^{dagger} T");
}

  RealMatrix UK = SUncontracted * TUncontracted;

if(this->printLevel_ >= 3){ 
  prettyPrintSmart(this->fileio_->out,UK,"U (K) transformation");
}

// Now we transform V to V' 
  TMP = UK.transpose() * VUncontracted;
  RealMatrix P2_Potential = TMP * UK;


// Next we need to get W' (from pVp integrals)


  VectorXd SCRATCH1UnContracted(nUncontracted);
  RealMatrix SCRATCH2UnContracted(nUncontracted,nUncontracted);
  VectorXd SCRATCHDXUnContracted(nUncontracted);
  VectorXd SCRATCHDYUnContracted(nUncontracted);
  VectorXd SCRATCHDZUnContracted(nUncontracted);

  size_t tracker(0);
// compute pVp numerically
  auto PVP = [&](IntegrationPoint pt, std::vector<RealMatrix> &result) {
    tracker++;
//  cout << tracker << endl;
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
    const double * gamma = engineV.compute(this->molecule_->nucShell(iAtm),
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
          double *pVpSDATA = result[0].data() + iBf*nUncontracted;
          double *pVpXDATA = result[1].data() + iBf*nUncontracted;
          double *pVpYDATA = result[2].data() + iBf*nUncontracted;
          double *pVpZDATA = result[3].data() + iBf*nUncontracted;
          double *SCRATCHDXUnContractedData = SCRATCHDXUnContracted.data();
          double *SCRATCHDYUnContractedData = SCRATCHDYUnContracted.data();
          double *SCRATCHDZUnContractedData = SCRATCHDZUnContracted.data();
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
  std::vector<RealMatrix> numPot(4,RealMatrix::Zero(nUncontracted,
        nUncontracted));

  cout << "Performing numerical PVP integrals..." << endl;
  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    tracker = 0;
    cout << "IAtm = " << iAtm << endl;
    AGrid.center() = iAtm;
    AGrid.setRadCutOff(atomRadCutoff[iAtm]);
    AGrid.findNearestNeighbor();
    AGrid.scalingFactor()=
      0.5*elements[this->molecule_->index(iAtm)].sradius/phys.bohr;
    AGrid.integrate<std::vector<RealMatrix>>(PVP,numPot);
  };
  cout << "done!" << endl;

  
/*
  for(auto i = 0; i < 10; i++) numPot[i] *= 4 * math.pi;
  // Scalar = mu(x)nu(x) + mu(y)nu(y) +  mu(z)nu(z) 
  RealMatrix PVPS = numPot[1] + numPot[5] + numPot[9]; 
  // X = mu(y)nu(z) - mu(z)nu(y)
  RealMatrix PVPX = numPot[6] - numPot[8];
  // Y = mu(z)nu(x) - mu(x)nu(z)
  RealMatrix PVPY = numPot[7] - numPot[3];
  // Z = mu(x)nu(y) - mu(y)nu(x)
  RealMatrix PVPZ = numPot[2] - numPot[4]; 
*/

  double fact4pi = 4 * math.pi; 
  for(auto i = 0; i <=3 ; i++) numPot[i] *= fact4pi;
   numPot[0].triangularView<Upper>() =  numPot[0].transpose(); // Sym
   numPot[1].triangularView<Upper>() = -numPot[1].transpose(); // AntiSymm
   numPot[2].triangularView<Upper>() = -numPot[2].transpose(); // AntiSymm
   numPot[3].triangularView<Upper>() = -numPot[3].transpose(); // AntiSymm
   RealMatrix PVPS  = numPot[0];
   RealMatrix PVPX  = numPot[1];
   RealMatrix PVPY  = numPot[2];
   RealMatrix PVPZ  = numPot[3];
   prettyPrint(this->fileio_->out,PVPS,"pVpSnum");
   prettyPrint(this->fileio_->out,PVPX,"pVpXnum");
   prettyPrint(this->fileio_->out,PVPY,"pVpYnum");
   prettyPrint(this->fileio_->out,PVPZ,"pVpZnum");
// Apply the Uk unitary transformation ( Uk' * PVP * Uk)

  TMP = UK.adjoint() * PVPS;
  PVPS = TMP * UK;
  TMP = UK.adjoint() * PVPX;
  PVPX = TMP * UK;
  TMP = UK.adjoint() * PVPY;
  PVPY = TMP * UK;
  TMP = UK.adjoint() * PVPZ;
  PVPZ = TMP * UK;

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


  TMP = PMap.asDiagonal() * PVPS;
  PVPS = TMP * PMap.asDiagonal();
  TMP = PMap.asDiagonal() * PVPX;
  PVPX = TMP * PMap.asDiagonal();
  TMP = PMap.asDiagonal() * PVPY;
  PVPY = TMP * PMap.asDiagonal();
  TMP = PMap.asDiagonal() * PVPZ;
  PVPZ = TMP * PMap.asDiagonal();

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


  ComplexMatrix W(2*nUncontracted,2*nUncontracted);
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

//W = -W;
//P2_Potential = - P2_Potential;
//W = W.conjugate();

// -------------------------------------------------
// Put all pieces of core Hamiltonian in block form:
// [ V'      cp    ]
// [    	 	   ]
// [ cp   W'-2mc^2 ]
// -------------------------------------------------

  cout << "Forming 4C Core Hamiltonian...";
  PMap = PMap.cwiseInverse(); //switch from p^-1 back to p
  ComplexMatrix CORE_HAMILTONIAN(4*nUncontracted,4*nUncontracted);

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

   cout << "done!" << endl;
// ------------------------------
// Diagonalize 
// ------------------------------

if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,CORE_HAMILTONIAN,"H");
  (this->fileio_->out) << "|H|" << CORE_HAMILTONIAN.squaredNorm() << endl;
}

  cout << "Diagonalizing 4C Core Hamiltonian...";
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> es;
  es.compute(CORE_HAMILTONIAN);
  RealMatrix HEV= es.eigenvalues();
  ComplexMatrix HEVx= es.eigenvectors();
  cout << "done!" << endl;

// Print out the energies (eigenvalues) and eigenvectors
if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,HEV,"HEV");
  prettyPrintSmart(this->fileio_->out,HEVx,"HEVc");
}

// Grab C_L (+) and C_S (+) - the large and small components
// of the electronic (positive energy) solutions
//
// NOT SURE IF THESE ARE CORRECT!!!! (seems to be wrong 'column')
  ComplexMatrix L = 
    HEVx.block(0,2*nUncontracted,2*nUncontracted,2*nUncontracted);
  ComplexMatrix S = 
    HEVx.block(2*nUncontracted,2*nUncontracted,2*nUncontracted,2*nUncontracted);

  cout << "SVD of L...";
// Do we even use this SVD in calculating L inverse?
  Eigen::JacobiSVD<ComplexMatrix> 
    svd(L,Eigen::ComputeThinU | Eigen::ComputeThinV);

  VectorXd SigmaL = svd.singularValues();
  ComplexMatrix SVL = svd.matrixU();

  ComplexMatrix X = S * L.inverse(); //See above!
  cout << "done!" << endl;

// Print out X and its squared norm
if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,X,"X");
  (this->fileio_->out) << X.squaredNorm() << endl;
}
  
// Calculate Y = sqrt(1 + X'X)
// Also known as the 'renormalization matrix' R
  ComplexMatrix Y = 
    (ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) 
     + X.adjoint() * X).pow(-0.5);

// Print out Y and its squared norm
if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,Y,"Y");
  (this->fileio_->out) << Y.squaredNorm() << endl;
}

// Get PMapC = p (as 2n by 2n matrix)
  ComplexMatrix PMapC(2*nUncontracted,2*nUncontracted);
  PMapC.block(0,0,nUncontracted,nUncontracted).real() = PMap.asDiagonal();
  PMapC.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = PMap.asDiagonal();


// Compute p^2 and then the relativistic kinetic energy
// Here we have       __________________
//              T = \/m^2c^4 + c^2 * p^2  - mc^2
//
// Where m is the rest mass of the electron (1 in atomic units)             
//
  ComplexMatrix P2MapC = PMapC.cwiseProduct(PMapC);
  ComplexMatrix KinEn = P2MapC * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT;
  KinEn += ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted)* phys.SPEED_OF_LIGHT *
	phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT; 
  KinEn = KinEn.cwiseSqrt();
  KinEn -= ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT;

if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,KinEn,"Relativistic Kinetic Energy");
}


// Get P2_PotC == V'
  ComplexMatrix P2_PotC(2*nUncontracted,2*nUncontracted);
  P2_PotC.block(0,0,nUncontracted,nUncontracted).real() = P2_Potential;
  P2_PotC.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = P2_Potential;

if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,P2_PotC,"V prime (p space)");
}

// Calculate the 2-component core Hamiltonian in the uncontracted basis
//  αα | αβ
//  -------
//  βα | ββ
  
  ComplexMatrix HCore(2*nUncontracted,2*nUncontracted);
  HCore = P2_PotC;

  ComplexMatrix TEMP(2*nUncontracted,2*nUncontracted);

  TEMP = phys.SPEED_OF_LIGHT * PMapC * X;
  HCore = HCore + TEMP;
  TEMP = phys.SPEED_OF_LIGHT * X.adjoint() * PMapC;
  HCore = HCore + TEMP;
  TEMP = 2 * phys.SPEED_OF_LIGHT * phys.SPEED_OF_LIGHT * 
	ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted);
  TEMP = W - TEMP;
  TEMP = X.adjoint() * TEMP;
  TEMP = TEMP * X;
  HCore = HCore + TEMP;
  HCore = HCore * Y;
  HCore = Y * HCore;

if(this->printLevel_ >= 3){
  prettyPrintSmart(this->fileio_->out,HCore,"Transformed HCore (in p space) ");
}

  ComplexMatrix Veff(2*nUncontracted,2*nUncontracted);
  Veff = HCore - KinEn;

  RealMatrix Hs(nUncontracted,nUncontracted);
  RealMatrix Hz(nUncontracted,nUncontracted);
  RealMatrix Hx(nUncontracted,nUncontracted);
  RealMatrix Hy(nUncontracted,nUncontracted);

  Hs = 0.5 * (HCore.block(0,0,nUncontracted,nUncontracted).real()
	+ HCore.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real());
  Hz = 0.5 * (HCore.block(0,0,nUncontracted,nUncontracted).imag()
	- HCore.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).imag());
  Hx = 0.5 * (HCore.block(0,nUncontracted,nUncontracted,nUncontracted).imag()
	+ HCore.block(nUncontracted,0,nUncontracted,nUncontracted).imag());
  Hy = 0.5 * (HCore.block(0,nUncontracted,nUncontracted,nUncontracted).real()
	- HCore.block(nUncontracted,0,nUncontracted,nUncontracted).real());

if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,Hs,"Hs (p space)");
  prettyPrintSmart(this->fileio_->out,Hz,"Hz (p space)");
  prettyPrintSmart(this->fileio_->out,Hx,"Hx (p space)");
  prettyPrintSmart(this->fileio_->out,Hy,"Hz (p space)");
}

  RealMatrix rTEMP(nUncontracted,nUncontracted);

  rTEMP = Hs * UK.adjoint() * SUn;
  Hs = SUn * UK * rTEMP; 
  rTEMP = Hz * UK.adjoint() * SUn;
  Hz = SUn * UK * rTEMP; 
  rTEMP = Hx * UK.adjoint() * SUn;
  Hx = SUn * UK * rTEMP; 
  rTEMP = Hy * UK.adjoint() * SUn;
  Hy = SUn * UK * rTEMP; 

if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,Hs,"Hs (r space)");
  prettyPrintSmart(this->fileio_->out,Hz,"Hz (r space)");
  prettyPrintSmart(this->fileio_->out,Hx,"Hx (r space)");
  prettyPrintSmart(this->fileio_->out,Hy,"Hy (r space)");
}

     
//  prettyPrintSmart(this->fileio_->out,Veff,"Veff (p space)");
/* 
  ComplexMatrix SUK(2*nUncontracted,2*nUncontracted);
  SUK.block(0,0,nUncontracted,nUncontracted) = SUn * UK;  
  SUK.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted) = SUn * UK;

  TEMP = SUK * Veff;
  Veff = TEMP * SUK.adjoint();  
   
  prettyPrintSmart(this->fileio_->out,Veff,"Veff (r space)");

  TEMP = SUK * KinEn;
  KinEn = TEMP * SUK.adjoint();

  prettyPrintSmart(this->fileio_->out,KinEn,"Trel (r space)");
*/

// --------------------------------------------------------

//  TEMP = SUK * HCore;
//  HCore = TEMP * SUK.adjoint();

//  cout << HCore.squaredNorm() << " HCore norm" << endl;

//  prettyPrintSmart(this->fileio_->out,HCore,"HCore (r space)");




// Recontract the basis
//  prettyPrintSmart(this->fileio_->out,*this->basisSet_->mapPrim2Bf(),"primitive transformation matrix");
  RealMatrix IPrim2Bf = (*this->basisSet_->mapPrim2Bf()).transpose();

  RealMatrix TCon = 0.5 * (KinEn.block(0,0,nUncontracted,nUncontracted).real() + 
	KinEn.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real());
  RealMatrix VCon = 0.5 * (Veff.block(0,0,nUncontracted,nUncontracted).real() + 
	Veff.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real());
if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,TCon,"Trel (p space)");
  prettyPrintSmart(this->fileio_->out,VCon,"Vrel (p space)");
 }
  rTEMP = TCon * UK.adjoint() * SUn;
  TCon = SUn * UK * rTEMP;
  rTEMP = VCon * UK.adjoint() * SUn;
  VCon = SUn * UK * rTEMP;
if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,TCon,"Trel (r space)");
  prettyPrintSmart(this->fileio_->out,VCon,"Vrel (r space)");
 }

  TCon = (*this->basisSet_->mapPrim2Bf()) * TCon * IPrim2Bf;
  VCon = (*this->basisSet_->mapPrim2Bf()) * VCon * IPrim2Bf;

if(this->printLevel_ >= 2){
  prettyPrintSmart(this->fileio_->out,TCon,"TCon");
  prettyPrintSmart(this->fileio_->out,VCon,"VCon");
 }
 
  RealMatrix CoreS = (*this->basisSet_->mapPrim2Bf()) * Hs * IPrim2Bf;
  RealMatrix CoreZ = (*this->basisSet_->mapPrim2Bf()) * Hz * IPrim2Bf;
  RealMatrix CoreX = (*this->basisSet_->mapPrim2Bf()) * Hx * IPrim2Bf;
  RealMatrix CoreY = (*this->basisSet_->mapPrim2Bf()) * Hy * IPrim2Bf;

if(this->printLevel_ >= 2) {
  prettyPrintSmart(this->fileio_->out,CoreS,"Core (scalar)");
  prettyPrintSmart(this->fileio_->out,CoreZ,"Core (mz)");
  prettyPrintSmart(this->fileio_->out,CoreX,"Core (mx)");
  prettyPrintSmart(this->fileio_->out,CoreY,"Core (my)");
}

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

  *this->coreH_ = CoreS;
  *this->oneEmx_ = CoreX;
  *this->oneEmy_ = CoreY;
  *this->oneEmz_ = CoreZ;

  ComplexMatrix TCSham(2*nBasis_,2*nBasis_);

/*
  std::vector<std::reference_wrapper<RealMatrix>> mats;
  
  mats.emplace_back(CoreS);
  mats.emplace_back(CoreZ);
  mats.emplace_back(CoreY);
  mats.emplace_back(CoreX);

  Quantum<double>::spinGather(TCSham,mats);
*/
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

/*
  es.compute(TCSham);
  HEV= es.eigenvalues();
//  HEVx= es.eigenvectors();

  prettyPrintSmart(cout,HEV,"HEV");
//  prettyPrintSmart(cout,HEVx,"HEVc");
*/
}

