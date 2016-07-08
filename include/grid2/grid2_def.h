#ifndef INCLUDED_GRID2_DEF_H
#define INCLUDED_GRID2_DEF_H

#include <global.h>
#include <cerr.h>

namespace ChronusQ {

enum GRID_TYPE {
  GAUSSCHEBFST,
  GAUSSCHEBSND,
  EULERMAC,
  LEBEDEV
};

struct IntegrationPoint {
  cartGP pt;
  double weight;
  bool evalpt;
  IntegrationPoint(double pt_ = 0, double weight_ = 0) : 
    pt(pt_), weight(weight_), evalpt(true){ };
  IntegrationPoint(cartGP pt_ = cartGP(0.0,0.0,0.0), double weight_ = 0) : 
    pt(pt_), weight(weight_), evalpt(true){ };
  IntegrationPoint(const IntegrationPoint &other):
    pt(other.pt),
    weight(other.weight),
    evalpt(other.evalpt){ };
};

class Grid2 {
protected:
    size_t nPts_;  ///< number of grid points       
    std::vector<cartGP> gPoints_;
    std::vector<double> weights_;
    bool onTheFly_;
    bool haveGPs_;
public:
    Grid2(size_t npts = 0, bool onTheFly = true) : nPts_(npts),
      onTheFly_(onTheFly), haveGPs_(false){ 
      
    };

    virtual ~Grid2(){ };

    virtual IntegrationPoint operator[](size_t) = 0;
    size_t npts(){return this->nPts_;};
    virtual void generateGridPoints() {
      this->gPoints_.reserve(this->nPts_);
      this->weights_.reserve(this->nPts_);

      for(auto iPt = 0; iPt < nPts_; iPt++) {
        IntegrationPoint tmp = operator[](iPt);
        this->gPoints_.push_back(tmp.pt);
        this->weights_.push_back(tmp.weight);
      }
      this->haveGPs_ = true;
    };

    // Integrate Functions
    template<typename T>
    inline T integrate(std::function< T(cartGP) > func) {

      cout << "Integrate 1" << endl;
      T result = (*this)[0].weight * func((*this)[0].pt);

      for(auto iPt = 1; iPt < this->nPts_; iPt++){
        IntegrationPoint intPoint = (*this)[iPt];
        cartGP gp = intPoint.pt;
        result += intPoint.weight * func(gp);
      }
      return result;
    };

    template <typename T>
    inline void integrate(std::function< T(IntegrationPoint) > func,
        T& result) {

      cout << "Integrate 2" << endl;
      result = func((*this)[0]);
      for(auto iPt = 1; iPt < this->nPts_; iPt++)
        result += func((*this)[iPt]);

    };

    template <typename T>
    inline void integrate(std::function< void(IntegrationPoint,T&) > func,
        T& result) {
      std::size_t NSkip(0);
      for(auto iPt = 0; iPt < this->nPts_; iPt++){
        IntegrationPoint IP((*this)[iPt]);
        if(IP.evalpt){
          func(IP,result);
        } else NSkip++;
      }
      cout << "NSkip : " << NSkip << endl;
    };


    // Print function
    void printGrid(std::ostream&);
}; // class Grid2
}
#endif

