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
    VectorXd rA(3), rB(3), rAB(3);

    auto h = [](double x) -> double{
      return 1.5 * x - 0.5 * x * x * x;
    };


    for(auto iCenter = 0; iCenter < this->centers_.size(); iCenter++){
      this->partitionScratch_[iCenter] = 1.0;
      rA(0) = bg::get<0>(pt) - this->centers_[iCenter][0];
      rA(1) = bg::get<1>(pt) - this->centers_[iCenter][1];
      rA(2) = bg::get<2>(pt) - this->centers_[iCenter][2];
      for(auto jCenter = 0; jCenter < this->centers_.size(); jCenter++){
        if(iCenter == jCenter) continue;

//        cout << "HERE_PART i=" << iCenter << " j= "<< jCenter <<endl;
        rB(0) = bg::get<0>(pt) - this->centers_[jCenter][0];
        rB(1) = bg::get<1>(pt) - this->centers_[jCenter][1];
        rB(2) = bg::get<2>(pt) - this->centers_[jCenter][2];
        rAB(0) = this->centers_[iCenter][0] - this->centers_[jCenter][0];
        rAB(1) = this->centers_[iCenter][1] - this->centers_[jCenter][1];
        rAB(2) = this->centers_[iCenter][2] - this->centers_[jCenter][2];
     
        double mu = (rA.norm() - rB.norm()) / rAB.norm();
        this->partitionScratch_[iCenter] *= 0.5 * (1.0 - h(h(h(mu))));
      }
    }

    double sum = 0.0;
    for(auto iCenter = 0; iCenter < this->centers_.size(); iCenter++)
      sum += this->partitionScratch_[iCenter];

    //cout << "WEIGHT " << P / sum << endl;
//    cout << "END Point " <<endl;
    return this->partitionScratch_[this->centerIndx_] / sum;

  };
};
