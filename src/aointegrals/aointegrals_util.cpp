#include <aointegrals.h>
using ChronusQ::AOIntegrals;
using ChronusQ::BasisSet;

RealMatrix AOIntegrals::genSpx(BasisSet &obs1, BasisSet &obs2){
  if(obs1.nShell() != obs2.nShell()) {
    cout << "genSpx cannot take two basis sets of different size" 
         << endl;
    std::exit(1);
  }
  RealMatrix S(obs1.nBasis(),obs2.nBasis());

  libint2::OneBodyEngine engine(libint2::OneBodyEngine::overlap,
    obs1.maxPrim(),obs1.maxL(),0);

  S.setZero();
  for(auto iSh = 0; iSh < obs1.nShell(); iSh++)
  for(auto jSh = 0; jSh < obs2.nShell(); jSh++){
    auto s1 = obs1.mapSh2Bf(iSh);
    auto s2 = obs1.mapSh2Bf(jSh);
    auto n1 = obs1.shells(iSh).size();
    auto n2 = obs2.shells(jSh).size();
    const auto *buff = engine.compute(
      obs1.shells(iSh),obs2.shells(jSh)
    );

    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,
                                   Eigen::Dynamic,
                                   Eigen::RowMajor>
      > sMap(buff,n1,n2);
    S.block(s1,s2,n1,n2) = sMap;
  }

/*
  // Renormalize
  Eigen::VectorXd norms1(obs1.nbf());
  for(auto iSh = 0; iSh < obs1.size(); iSh++){
    auto s1 = shell2bf[iSh];
    auto n1 = obs1[iSh].size();
    const auto *buff = engine.compute(obs1[iSh],obs1[iSh]);
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,
                                   Eigen::Dynamic,
                                   Eigen::RowMajor>
      > sMap(buff,n1,n1);

    for(auto iBf = s1, iloc = 0ul; iBf < s1 + n1; iBf++, iloc++)
      norms1(iBf) = std::sqrt(sMap(iloc,iloc));
  }

  Eigen::VectorXd norms2(obs2.nbf());
  for(auto iSh = 0; iSh < obs2.size(); iSh++){
    auto s2 = shell2bf[iSh];
    auto n2 = obs2[iSh].size();
    const auto *buff = engine.compute(obs2[iSh],obs2[iSh]);
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,
                                   Eigen::Dynamic,
                                   Eigen::RowMajor>
      > sMap(buff,n2,n2);

    for(auto iBf = s2, iloc = 0ul; iBf < s2 + n2; iBf++, iloc++)
      norms2(iBf) = std::sqrt(sMap(iloc,iloc));
  }
  
  for(auto iBf = 0; iBf < obs1.nbf(); iBf++)
  for(auto jBf = 0; jBf < obs2.nbf(); jBf++) {
    S(iBf,jBf) = S(iBf,jBf) / (norms1(iBf)*norms2(jBf));
  }
*/

  return S;
};
