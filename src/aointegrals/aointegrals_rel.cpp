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

  } // s2
  } // s1

  SUncontracted = SUncontracted.selfadjointView<Lower>();
  TUncontracted = TUncontracted.selfadjointView<Lower>();

  char JOBU = 'O', JOBVT = 'N';
  int LWORK = 6*nUncontracted;
  int INFO;
  std::vector<double> ovlpEigValues(nUncontracted);
  std::vector<double> WORK(LWORK);

  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,SUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],SUncontracted.data(),&nUncontracted,
      SUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);


  for(auto iS = 0; iS < nUncontracted; iS++){
    SUncontracted.col(iS) /= std::sqrt(ovlpEigValues[iS]);
  }


  int nZero = 0;
  for(auto iS = 0; iS < nUncontracted; iS++)
    if(std::abs(ovlpEigValues[iS]) < 1e-6) nZero++;

  cout << "NZERO " << nZero << endl;

  RealMatrix TMP = SUncontracted.transpose() * TUncontracted;
  TUncontracted = TMP * SUncontracted;

  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,TUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],TUncontracted.data(),&nUncontracted,
      TUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);

  RealMatrix UK = SUncontracted * TUncontracted;

  TMP = UK.transpose() * VUncontracted;;
  RealMatrix P2_Potential = TMP * UK;

  VectorXd SCRATCH1UnContracted(nUncontracted);
  RealMatrix SCRATCH2UnContracted(nUncontracted,nUncontracted);
  VectorXd SCRATCHDXUnContracted(nUncontracted);
  VectorXd SCRATCHDYUnContracted(nUncontracted);
  VectorXd SCRATCHDZUnContracted(nUncontracted);

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

      delete [] buff;
    };

    std::vector<std::pair<double,std::array<double,3>>> q;
    q.push_back(
      {1.0, {{bg::get<0>(pt.pt),bg::get<1>(pt.pt),bg::get<2>(pt.pt)}}});
    engineV.set_params(q);

    for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
      // Gaussian Nuclei
      const double * gamma = engineV.compute(this->molecule_->nucShell(iAtm),
          libint2::Shell::unit());
                                                    
      SCRATCH2UnContracted.noalias() = 
        SCRATCH1UnContracted * SCRATCH1UnContracted.transpose();
      result[0].noalias() += (pt.weight * (*gamma)) * SCRATCH2UnContracted;

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
  AtomicGrid AGrid(100,590,GAUSSCHEBFST,LEBEDEV,BECKE,atomicCenters,0,1.0,
      false);
  std::vector<RealMatrix> numPot(10,RealMatrix::Zero(nUncontracted,
        nUncontracted));

  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    AGrid.center() = iAtm;
    AGrid.scalingFactor()=
      0.5*elements[this->molecule_->index(iAtm)].sradius/phys.bohr;
    AGrid.integrate<std::vector<RealMatrix>>(PVP,numPot);
  };

  for(auto i = 0; i < 10; i++) numPot[i] *= 4 * math.pi;

  RealMatrix PVPS = numPot[1] + numPot[5] + numPot[9]; 
  RealMatrix PVPX = numPot[6] - numPot[8];
  RealMatrix PVPY = numPot[7] - numPot[3];
  RealMatrix PVPZ = numPot[2] - numPot[4]; 


  TMP = UK.transpose() * PVPS;
  PVPS = TMP * UK;
  TMP = UK.transpose() * PVPX;
  PVPX = TMP * UK;
  TMP = UK.transpose() * PVPY;
  PVPY = TMP * UK;
  TMP = UK.transpose() * PVPZ;
  PVPZ = TMP * UK;

//prettyPrint(cout,P2_Potential,"V");
//prettyPrint(cout,PVPS,"dot(P,VP)");
//prettyPrint(cout,PVPX,"cross(P,VP) X");
//prettyPrint(cout,PVPY,"cross(P,VP) Y");
//prettyPrint(cout,PVPZ,"cross(P,VP) Z");
//cout << "|V| = " << P2_Potential.squaredNorm() << endl;
//cout << "|dot(P,VP)| = " << PVPS.squaredNorm() << endl;
//cout << "|cross(P,VP) X| = " << PVPX.squaredNorm() << endl;
//cout << "|cross(P,VP) Y| = " << PVPY.squaredNorm() << endl;
//cout << "|cross(P,VP) Z| = " << PVPZ.squaredNorm() << endl;

  RealVecMap PMap(&ovlpEigValues[0],nUncontracted);
//prettyPrint(cout,PMap,"P^2");

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

//cout << "|dot(P,VP)| = " << PVPS.squaredNorm() << endl;
//cout << "|cross(P,VP) X| = " << PVPX.squaredNorm() << endl;
//cout << "|cross(P,VP) Y| = " << PVPY.squaredNorm() << endl;
//cout << "|cross(P,VP) Z| = " << PVPZ.squaredNorm() << endl;
//prettyPrint(cout,PVPS,"scaled dot(P,VP)");
//prettyPrint(cout,PVPX,"scaled cross(P,VP) X");
//prettyPrint(cout,PVPY,"scaled cross(P,VP) Y");
//prettyPrint(cout,PVPZ,"scaled cross(P,VP) Z");

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

  PMap = PMap.cwiseInverse();
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


//prettyPrintComplex(cout,CORE_HAMILTONIAN,"H");
//cout << "|H|" << CORE_HAMILTONIAN.squaredNorm();

  Eigen::SelfAdjointEigenSolver<ComplexMatrix> es;
  es.compute(CORE_HAMILTONIAN);
  RealMatrix HEV= es.eigenvalues();
  ComplexMatrix HEVx= es.eigenvectors();
//prettyPrint(cout,HEV,"HEV");
//prettyPrintComplex(cout,HEVx,"HEVc");

  ComplexMatrix L = 
    HEVx.block(0,2*nUncontracted,2*nUncontracted,2*nUncontracted);
  ComplexMatrix S = 
    HEVx.block(2*nUncontracted,2*nUncontracted,2*nUncontracted,2*nUncontracted);

  Eigen::JacobiSVD<ComplexMatrix> 
    svd(L,Eigen::ComputeThinU | Eigen::ComputeThinV);

  VectorXd SigmaL = svd.singularValues();
  ComplexMatrix SVL = svd.matrixU();

  ComplexMatrix X = S * L.inverse();
//prettyPrintComplex(cout,X,"X");
//cout << X.squaredNorm() << endl;

  
  ComplexMatrix Y = 
    (ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) 
     + X.adjoint() * X).pow(-0.5);
//prettyPrintComplex(cout,Y,"Y");
//cout << Y.squaredNorm() << endl;


};
