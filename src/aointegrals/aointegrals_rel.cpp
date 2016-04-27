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

void AOIntegrals::formP2Transformation(){

  if(!this->isPrimary) return;
  auto unContractedShells = this->basisSet_->uncontractBasis();
  int nUncontracted = 0;
  for(auto i : unContractedShells) nUncontracted += i.size();

  RealMatrix SUncontracted(nUncontracted,nUncontracted);
  RealMatrix TUncontracted(nUncontracted,nUncontracted);

  libint2::Engine engineS(
      libint2::Operator::overlap,1,this->basisSet_->maxL(),0);
  libint2::Engine engineT(
      libint2::Operator::kinetic,1,this->basisSet_->maxL(),0);

  engineS.set_precision(0.0);
  engineT.set_precision(0.0);
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

  } // s2
  } // s1

  SUncontracted = SUncontracted.selfadjointView<Lower>();
  TUncontracted = TUncontracted.selfadjointView<Lower>();

  char JOBU = 'O', JOBVT = 'N';
  int LWORK = 6*nUncontracted;
  int INFO;
  std::vector<double> ovlpEigValues(nUncontracted);
  std::vector<double> WORK(LWORK);

  prettyPrint(cout,SUncontracted,"S");
  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,SUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],SUncontracted.data(),&nUncontracted,
      SUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);
  prettyPrint(cout,SUncontracted,"SVD");


  for(auto iS = 0; iS < nUncontracted; iS++){
    SUncontracted.col(iS) /= std::sqrt(ovlpEigValues[iS]);
  }

  prettyPrint(cout,SUncontracted,"SVD Scaled");

  int nZero = 0;
  for(auto iS = nUncontracted - 1; iS >= 0; iS--)
    if(std::abs(ovlpEigValues[iS]) < 1e-6) nZero++;

  cout << "S" << endl;
  for(auto iS = 0; iS < nUncontracted; iS++)
    cout << ovlpEigValues[iS] << endl;

  prettyPrint(cout,TUncontracted,"T");
  cout << "NZERO " << nZero << endl;
//RealMatrix TMP = SUncontracted * TUncontracted;
//TUncontracted = TMP * SUncontracted.transpose();

  RealMatrix TMP = SUncontracted.transpose() * TUncontracted;
  TUncontracted = TMP * SUncontracted;

  dgesvd_(&JOBU,&JOBVT,&nUncontracted,&nUncontracted,TUncontracted.data(),
      &nUncontracted,&ovlpEigValues[0],TUncontracted.data(),&nUncontracted,
      TUncontracted.data(),&nUncontracted,&WORK[0],&LWORK,&INFO);
  cout << "T" << endl;
  for(auto iT = 0; iT < nUncontracted; iT++)
    cout << ovlpEigValues[iT] << endl;
};
