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

  // Loop through and uncontract all basis functions
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

  RealMatrix SUn(nUncontracted,nUncontracted);
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

// normalize
// What happens when we divide by zero?
  for(auto iS = 0; iS < nUncontracted; iS++){
    SUncontracted.col(iS) /= std::sqrt(ovlpEigValues[iS]);
  }

// Count linear dependencies here?
  int nZero = 0;
  for(auto iS = 0; iS < nUncontracted; iS++)
    if(std::abs(ovlpEigValues[iS]) < 1e-6) nZero++;

  cout << "NZERO " << nZero << endl;

  prettyPrint(this->fileio_->out,SUncontracted.transpose()*SUncontracted, "S^{dagger} S");
// Put the kinetic energy in the orthonormal basis 

  RealMatrix TMP = SUncontracted.transpose() * TUncontracted;
  TUncontracted = TMP * SUncontracted;

// Get rid of the linear dependencies?
  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,TUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],TUncontracted.data(),&nUncontracted,
      TUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);

// Form the K transformation matrix from both pieces
// (diagonalizes S) and (diagonalizes T). See eq. 12 in Rieher's paper from 2013
  prettyPrint(this->fileio_->out,SUncontracted,"S transformation");
  prettyPrint(this->fileio_->out,TUncontracted,"T transformation");
  prettyPrint(this->fileio_->out,TUncontracted.transpose()*TUncontracted, "T^{dagger} T");


  RealMatrix UK = SUncontracted * TUncontracted;
 
  prettyPrint(this->fileio_->out,UK,"U (K) transformation");

  prettyPrint(this->fileio_->out,UK.transpose()*UK, "U^{dagger} U");
  
// Now we transform V to V' 
  TMP = UK.transpose() * VUncontracted;
  RealMatrix P2_Potential = TMP * UK;


// Next we need to get W' (from pVp integrals)


  VectorXd SCRATCH1UnContracted(nUncontracted);
  RealMatrix SCRATCH2UnContracted(nUncontracted,nUncontracted);
  VectorXd SCRATCHDXUnContracted(nUncontracted);
  VectorXd SCRATCHDYUnContracted(nUncontracted);
  VectorXd SCRATCHDZUnContracted(nUncontracted);

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
  AtomicGrid AGrid(100,590,GAUSSCHEBFST,LEBEDEV,BECKE,atomicCenters,this->molecule_->rIJ(),1.0,
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

// Apply the Uk unitary transformation ( Uk' * PVP * Uk)
  TMP = UK.adjoint() * PVPS;
  PVPS = TMP * UK;
  TMP = UK.adjoint() * PVPX;
  PVPX = TMP * UK;
  TMP = UK.adjoint() * PVPY;
  PVPY = TMP * UK;
  TMP = UK.adjoint() * PVPZ;
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

// -------------------------------------------------
// Put all pieces of core Hamiltonian in block form:
// [ V'      cp    ]
// [	 	   ]
// [ cp   W'-2mc^2 ]
// -------------------------------------------------

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

// ------------------------------
// Diagonalize 
// ------------------------------

//prettyPrintComplex(cout,CORE_HAMILTONIAN,"H");
//cout << "|H|" << CORE_HAMILTONIAN.squaredNorm();

  Eigen::SelfAdjointEigenSolver<ComplexMatrix> es;
  es.compute(CORE_HAMILTONIAN);
  RealMatrix HEV= es.eigenvalues();
  ComplexMatrix HEVx= es.eigenvectors();

// Print out the energies (eigenvalues) and eigenvectors
//  prettyPrint(cout,HEV,"HEV");
//  prettyPrintComplex(cout,HEVx,"HEVc");

// Grab C_L (+) and C_S (+) - the large and small components
// of the electronic (positive energy) solutions
//
//
// NOT SURE IF THESE ARE CORRECT!!!! (seems to be wrong 'column')
  ComplexMatrix L = 
    HEVx.block(0,2*nUncontracted,2*nUncontracted,2*nUncontracted);
  ComplexMatrix S = 
    HEVx.block(2*nUncontracted,2*nUncontracted,2*nUncontracted,2*nUncontracted);

// Do we even use this SVD in calculating L inverse?
  Eigen::JacobiSVD<ComplexMatrix> 
    svd(L,Eigen::ComputeThinU | Eigen::ComputeThinV);

  VectorXd SigmaL = svd.singularValues();
  ComplexMatrix SVL = svd.matrixU();

  ComplexMatrix X = S * L.inverse(); //See above!

// Print out X and its squared norm
//  prettyPrintComplex(cout,X,"X");
//  cout << X.squaredNorm() << endl;

  
// Calculate Y = sqrt(1 + X'X)
// Also known as the 'renormalization matrix' R
  ComplexMatrix Y = 
    (ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) 
     + X.adjoint() * X).pow(-0.5);

// Print out Y and its squared norm
//  prettyPrintComplex(cout,Y,"Y");
//  cout << Y.squaredNorm() << endl;


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

  prettyPrint(this->fileio_->out,KinEn,"Relativistic Kinetic Energy");
  cout << KinEn.squaredNorm() << " Kintetic energy norm" << endl;

// Get P2_PotC == V'
  ComplexMatrix P2_PotC(2*nUncontracted,2*nUncontracted);
  P2_PotC.block(0,0,nUncontracted,nUncontracted).real() = P2_Potential;
  P2_PotC.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = P2_Potential;

  prettyPrint(this->fileio_->out,P2_PotC,"V prime");

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
  prettyPrint(this->fileio_->out,HCore,"HCore += TEMP");
  HCore = HCore * Y;
  HCore = Y * HCore;
  prettyPrint(this->fileio_->out,HCore,"Y * HCore * Y");

  cout << HCore.squaredNorm() << " HCore norm" << endl;
  ComplexMatrix Veff(2*nUncontracted,2*nUncontracted);
  Veff = HCore - KinEn;

     
  prettyPrint(this->fileio_->out,Veff,"Veff (p space)");
 
  ComplexMatrix SUK(2*nUncontracted,2*nUncontracted);
  SUK.block(0,0,nUncontracted,nUncontracted) = SUn * UK;  
  SUK.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted) = SUn * UK;

  TEMP = SUK * Veff;
  Veff = TEMP * SUK.adjoint();  
   
  prettyPrint(this->fileio_->out,Veff,"Veff (r space)");

  // All correct to here (up to a sign/permutation)
/*
  RealMatrix rTEMP(nUncontracted,nUncontracted);
  rTEMP = HCore.block(0,0,nUncontracted,nUncontracted).real() * UK.adjoint() * SUn;
  HCore.block(0,0,nUncontracted,nUncontracted).real() = SUn * UK * rTEMP; 
  rTEMP = HCore.block(nUncontracted,0,nUncontracted,nUncontracted).real() * UK.adjoint() * SUn;
  HCore.block(nUncontracted,0,nUncontracted,nUncontracted).real() = SUn * UK * rTEMP; 
  rTEMP = HCore.block(0,nUncontracted,nUncontracted,nUncontracted).real() * UK.adjoint() * SUn;
  HCore.block(0,nUncontracted,nUncontracted,nUncontracted).real() = SUn * UK * rTEMP; 
  rTEMP = HCore.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() * UK.adjoint() * SUn;
  HCore.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = SUn * UK * rTEMP; 
*/

  TEMP = SUK * HCore;
  HCore = TEMP * SUK.adjoint();

  cout << HCore.squaredNorm() << " HCore norm" << endl;

  prettyPrint(this->fileio_->out,HCore,"HCore (r space)");


  RealMatrix SpinCore(2*nUncontracted,2*nUncontracted);
  for (auto row = 0; row < nUncontracted; row++) {
	for (auto col = 0; col < nUncontracted; col++) {
	    SpinCore(2*row,2*col) = HCore(row,col).real(); // αα
	    SpinCore(2*row,2*col+1) = HCore(row,col+nUncontracted).real(); // αβ
	    SpinCore(2*row+1,2*col) = HCore(row+nUncontracted,col).real(); // βα
	    SpinCore(2*row+1,2*col+1) = HCore(row+nUncontracted,col+nUncontracted).real(); // ββ
	    }
	}

  cout << SpinCore.squaredNorm() << " Spin-Core norm" << endl;


// Now split these blocks into the scalar and magnetization parts
// Hs =  αα + ββ
// Hz =  αα - ββ
// Hx =  αβ + βα
// Hy = (αβ - βα)i
  RealMatrix Hs(nUncontracted,nUncontracted);
  RealMatrix Hz(nUncontracted,nUncontracted);
  RealMatrix Hx(nUncontracted,nUncontracted);
  RealMatrix Hy(nUncontracted,nUncontracted);
  
  std::vector<std::reference_wrapper<RealMatrix>> mats;

  // v.push_back(A) -> v.end() = A
  // v.emplace_back(A) -> v[end](A)

  mats.emplace_back(Hs);
  mats.emplace_back(Hz);
  mats.emplace_back(Hy);
  mats.emplace_back(Hx);

  Quantum<double>::spinScatter(SpinCore,mats);


  cout << Hs.squaredNorm() << " Hs norm" << endl;
  cout << Hz.squaredNorm() << " Hz norm" << endl;
  cout << Hx.squaredNorm() << " Hx norm" << endl;
  cout << Hy.squaredNorm() << " Hy norm" << endl;

/*
  Hs = 0.5 * (HCore.block(0,0,nUncontracted,nUncontracted).real()
	+ HCore.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real());
  Hz = 0.5 * (HCore.block(0,0,nUncontracted,nUncontracted).imag()
	- HCore.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).imag());
  Hx = 0.5 * (HCore.block(0,nUncontracted,nUncontracted,nUncontracted).imag()
	+ HCore.block(nUncontracted,0,nUncontracted,nUncontracted).imag());
  Hy = 0.5 * (HCore.block(0,nUncontracted,nUncontracted,nUncontracted).real()
	- HCore.block(nUncontracted,0,nUncontracted,nUncontracted).real());
*/
/*
  RealMatrix rTEMP(nUncontracted,nUncontracted);
  rTEMP = UK * Hs;
  Hs = rTEMP * UK.transpose();
  rTEMP = UK * Hz;
  Hz = rTEMP * UK.transpose();
  rTEMP = UK * Hx;
  Hx = rTEMP * UK.transpose();
  rTEMP = UK * Hy;
  Hy = rTEMP * UK.transpose();
*/
  prettyPrint(this->fileio_->out,Hs,"Hcore (scalar)");
  prettyPrint(this->fileio_->out,Hz,"Hcore (mz)");
  prettyPrint(this->fileio_->out,Hx,"Hcore (mx)");
  prettyPrint(this->fileio_->out,Hy,"Hcore (my)");




// Recontract the basis
  RealMatrix IPrim2Bf = (*this->basisSet_->mapPrim2Bf()).transpose();
  
  Hs = Hs * IPrim2Bf;
  Hs = (*this->basisSet_->mapPrim2Bf()) * Hs;
  Hz = Hz * IPrim2Bf;
  Hz = (*this->basisSet_->mapPrim2Bf()) * Hz;
  Hx = Hx * IPrim2Bf;
  Hx = (*this->basisSet_->mapPrim2Bf()) * Hx;
  Hy = Hy * IPrim2Bf;
  Hy = (*this->basisSet_->mapPrim2Bf()) * Hy;

  prettyPrint(this->fileio_->out,Hs,"Hs (scalar)");
  prettyPrint(this->fileio_->out,Hz,"Hz (mz)");
  prettyPrint(this->fileio_->out,Hx,"Hx (mx)");
  prettyPrint(this->fileio_->out,Hy,"Hy (my)");


  CErr();


/* THIS IS OLD AND DOESN'T WORK 
//  This doesn't seem to work...  
//  RealMatrix CUK = UK * (*this->basisSet_->mapPrim2Bf());

//
// What is the point of the matrix YCUK?
  ComplexMatrix YCUK(2*nUncontracted,this->nBasis_);
//  YCUK = Y * CUK;

  cout << "HERE" << endl;
  cout << YCUK.rows() << " " << YCUK.cols() << endl;
  cout << Y.rows() << " " << Y.cols() << endl;
  cout << CUK.rows() << " " << CUK.cols() << endl;
  cout << "HERE" << endl;

  ComplexMatrix PMapC(2*nUncontracted,2*nUncontracted);
  PMapC.block(0,0,nUncontracted,nUncontracted).real() = PMap.asDiagonal();
  PMapC.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = PMap.asDiagonal();
  ComplexMatrix P2MapC = PMapC.cwiseProduct(PMapC);
  ComplexMatrix P2_PotC(2*nUncontracted,2*nUncontracted);
  P2_PotC.block(0,0,nUncontracted,nUncontracted).real() = P2_Potential;
  P2_PotC.block(nUncontracted,nUncontracted,nUncontracted,nUncontracted).real() = P2_Potential;

  cout << "HERE" << endl;
  ComplexMatrix HCore = YCUK.inverse() * (
    P2_PotC+
    phys.SPEED_OF_LIGHT * PMapC * X +
    phys.SPEED_OF_LIGHT * X.adjoint() * PMapC +
    X.adjoint() * (W - 2*phys.SPEED_OF_LIGHT*phys.SPEED_OF_LIGHT * ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) ) * X -
    phys.SPEED_OF_LIGHT * 
      (phys.SPEED_OF_LIGHT*phys.SPEED_OF_LIGHT * ComplexMatrix::Identity(2*nUncontracted,2*nUncontracted) - P2MapC).pow(0.5)
    ) * YCUK;
  cout << "|HC|" << HCore.norm() << endl;

*/

}

