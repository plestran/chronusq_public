#include <grid2.h>
namespace ChronusQ {
  void Grid2::printGrid(std::ostream &out){
    /*
    out << "ListPointPlot3D[" << endl << "{";
    for(auto iPt = 0; iPt < this->nPts_; ++iPt){
      out << "{" << 
        bg::get<0>(operator[](iPt).pt) << "," <<
        bg::get<1>(operator[](iPt).pt) << "," <<
        bg::get<2>(operator[](iPt).pt) <<
        "}," << endl;
    };
    out << "}]" << endl;
    */
  };


  double AtomicGrid::evalPartitionWeight(cartGP& pt){
    // Evaluate current centers unnormalized partition weight
    // Generate Becke not-normalized Weights (if becke) 
    // according to the partition schems in
    // (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii Fuzzi 
    // Cells
    // Note these Weights have to be normailzed (see normBeckeW) 

    //VectorXd rA(3), rB(3), rAB(3);
    VectorXd rA(3), rB(3);
    auto h = [](double x) -> double{
    // Eq. 19
      return 1.5 * x - 0.5 * x * x * x;
    };

    auto gBecke = [&](double x) -> double {
      return h(h(h(x)));
    };

    // Frisch Weights
    double alpha = 0.64;
    auto z = [&](double x) -> double {
      double tmp,tmp2,tmp3;
      tmp  = x / alpha;
      tmp2 = tmp * tmp;

      tmp3 =  35 * tmp;
      tmp3 -= 35 * tmp * tmp2;
      tmp3 += 21 * tmp * tmp2 * tmp2;
      tmp3 -=  5 * tmp * tmp2 * tmp2 * tmp2;
   
      return tmp3 / 16.0; 
    };

    auto gFrisch = [&](double x) -> double {
      if(       x <= -alpha ) return -1.0;
      else if ( x >= alpha  ) return  1.0;
      else                    return  z(x);
    };


    auto g = [&](double x) -> double {
      if(     this->partitionScheme_ == BECKE)  return gBecke(x);
      else if(this->partitionScheme_ == FRISCH) return gFrisch(x);
    };


    for(auto iCenter = 0; iCenter < this->centers_.size(); iCenter++){
      this->partitionScratch_[iCenter] = 1.0;
      rA(0) = bg::get<0>(pt) - this->centers_[iCenter][0];
      rA(1) = bg::get<1>(pt) - this->centers_[iCenter][1];
      rA(2) = bg::get<2>(pt) - this->centers_[iCenter][2];
      for(auto jCenter = 0; jCenter < this->centers_.size(); jCenter++){
        if(iCenter == jCenter) continue;
        rB(0) = bg::get<0>(pt) - this->centers_[jCenter][0];
        rB(1) = bg::get<1>(pt) - this->centers_[jCenter][1];
        rB(2) = bg::get<2>(pt) - this->centers_[jCenter][2];
      //rAB(0) = this->centers_[iCenter][0] - this->centers_[jCenter][0];
      //rAB(1) = this->centers_[iCenter][1] - this->centers_[jCenter][1];
      //rAB(2) = this->centers_[iCenter][2] - this->centers_[jCenter][2];
      //double mu = (rA.norm() - rB.norm()) / rAB.norm();
        double mu = (rA.norm() - rB.norm())/(*this->rIJ_)(iCenter,jCenter);
        this->partitionScratch_[iCenter] *= 0.5 * (1.0 - g(mu));
      }
    }
    // Normalization
    double sum = 0.0;
    for(auto iCenter = 0; iCenter < this->centers_.size(); iCenter++)
      sum += this->partitionScratch_[iCenter];
    return this->partitionScratch_[this->centerIndx_] / sum;
  };
};
