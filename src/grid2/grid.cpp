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
    double P = 1.0;
    VectorXd rA(3), rB(3), rAB(3);

    auto h = [](double x) -> double{
      return 1.5 * x - 0.5 * x * x * x;
    };

    // Computes Unnormalized P for Grind Center and the contribution of
    // the grid center to all of the other partition weights
    rA(0) = bg::get<0>(pt) - this->center_[0];
    rA(1) = bg::get<1>(pt) - this->center_[1];
    rA(2) = bg::get<2>(pt) - this->center_[2];
    for(auto iCenter = 0; iCenter < this->otherCenters_.size(); iCenter++){
      this->partitionScratch_[iCenter] = 1.0;
      rB(0) = bg::get<0>(pt) - this->otherCenters_[iCenter][0];
      rB(1) = bg::get<1>(pt) - this->otherCenters_[iCenter][1];
      rB(2) = bg::get<2>(pt) - this->otherCenters_[iCenter][2];
      rAB(0) = this->center_[0] - this->otherCenters_[iCenter][0];
      rAB(1) = this->center_[1] - this->otherCenters_[iCenter][1];
      rAB(2) = this->center_[2] - this->otherCenters_[iCenter][2];

      double mu = (rA.norm() - rB.norm()) / rAB.norm();
      //cout << "MU" << mu << endl;
      //cout << "HHH" << h(h(h(mu))) << endl;
      P *= 0.5 * (1.0 - h(h(h(mu))));
      this->partitionScratch_[iCenter] *= 0.5 * (1.0 - h(h(h(-mu))));
    };

    for(auto iCenter = 0; iCenter < this->otherCenters_.size(); iCenter++){
      rA(0) = bg::get<0>(pt) - this->otherCenters_[iCenter][0];
      rA(1) = bg::get<1>(pt) - this->otherCenters_[iCenter][1];
      rA(2) = bg::get<2>(pt) - this->otherCenters_[iCenter][2];
      for(auto jCenter = 0; jCenter < this->otherCenters_.size(); jCenter++){
        if(iCenter == jCenter) continue;

//        cout << "HERE1" << endl;
        rB(0) = bg::get<0>(pt) - this->otherCenters_[jCenter][0];
        rB(1) = bg::get<1>(pt) - this->otherCenters_[jCenter][1];
        rB(2) = bg::get<2>(pt) - this->otherCenters_[jCenter][2];
        rAB(0) = 
          this->otherCenters_[iCenter][0] - this->otherCenters_[jCenter][0];
        rAB(1) = 
          this->otherCenters_[iCenter][1] - this->otherCenters_[jCenter][1];
        rAB(2) = 
          this->otherCenters_[iCenter][2] - this->otherCenters_[jCenter][2];
     
        double mu = (rA.norm() - rB.norm()) / rAB.norm();
        this->partitionScratch_[iCenter] *= 0.5 * (1.0 - h(h(h(mu))));
      }
    }

    double sum = P;
    for(auto iCenter = 0; iCenter < this->otherCenters_.size(); iCenter++)
      sum += this->partitionScratch_[iCenter];

    //cout << "WEIGHT " << P / sum << endl;
    return P / sum;

  };
};
