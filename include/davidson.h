#include <global.h>

namespace ChronusQ {
  template <typename TMat>
  class Davidson {
    int     n_;          // Dimension of the problem (LDA)
    std::shared_ptr<TMat> mat_;        // The full matrix to be diagonalized (?)
    bool                  hermetian_;  // Whether or not the problem is hemetian

    std::shared_ptr<TMat> guess_;      // Guess vectors
    std::shared_ptr<TMat> eigenvalues_;
    std::shared_ptr<TMat> eigenvector_;


    int     maxSubSpace_; // Maximum iterative subspace
    int     maxIter_;     // Maximum number of iterations (micro)
    int     MaxIter_;     // Maximum number of iterations (macro) 
    int          nSek_;        // Number of desired roots
    int          nGuess_;      // Number of given (or created) guesses

    void runMicro(ostream &output=cout);
    bool converged_;


  public:
    inline void run(ostream &output=cout) {
      time_t currentTime;
      std::chrono::high_resolution_clock::time_point start,finish;
      std::chrono::duration<double> elapsed;
      output << bannerTop << endl << endl;
      time(&currentTime);
      output << "Davidson Diagonalization Started: " << ctime(&currentTime) << endl;
      this->printInfo(output);
      start = std::chrono::high_resolution_clock::now();
      this->runMicro(output);
      finish = std::chrono::high_resolution_clock::now();
      elapsed = finish-start;
      time(&currentTime);
      output << "Davidson Diagonalization Finished: " << ctime(&currentTime);
      output << "Time Elapsed: " << elapsed.count() << " sec" << endl;
      output << bannerEnd << endl << endl;;
    };
    Davidson() {
      this->mat_         = nullptr;
      this->guess_       = nullptr;
      this->eigenvalues_ = nullptr;
      this->eigenvector_ = nullptr;
      this->maxSubSpace_ = 250;
      this->maxIter_     = 20;
      this->MaxIter_     = 20;
      this->nSek_        = 0;
      this->nGuess_      = 0;
      this->n_           = 0;
      this->converged_   = false;
      this->hermetian_   = false;
    }

    Davidson(std::shared_ptr<TMat> A, int nSek) {
      this->maxSubSpace_ = 250;
      this->maxIter_     = 128;
      this->MaxIter_     = 20;
      this->converged_   = false;
      this->mat_    = A;
      this->nSek_   = nSek;
      this->nGuess_ = 2*nSek;
      this->n_      = A->cols();
/*
      this->guess_  = 
        std::unique_ptr<TMat>(new TMat(this->n_,this->nGuess_));
      this->eigenvalues_ = 
        std::unique_ptr<TMat>(new TMat(this->nSek_,1));
      this->eigenvector_ = 
        std::unique_ptr<TMat>(new TMat(this->n_,this->nSek_));
*/
      this->guess_  = 
        std::make_shared<TMat>(this->n_,this->nGuess_);
      this->eigenvalues_ = 
        std::make_shared<TMat>(this->nSek_,1);
      this->eigenvector_ = 
        std::make_shared<TMat>(this->n_,this->nSek_);
      this->hermetian_ = true; // Only supports Hermetian for time being

      (*this->guess_) = TMat::Identity(this->n_,this->nGuess_); // Identity guess (primitive)
    }
    ~Davidson(){;};
    
    inline void printInfo(ostream &output=cout) {
      output << bannerTop << endl;
      output << "Davidson Diagonalization Settings:" << endl << endl;
 
      output << std::setw(50) << std::left << "  Dimension of the Matrix:" << this->n_ << endl;
      output << std::setw(50) << std::left << "  Number of Desired Roots:" << this->nSek_ << endl;
      output << std::setw(50) << std::left << "  Number of Initial Guess Vectors:" << this->nGuess_;
      if(this->nGuess_ == 2*this->nSek_) output << "    (default 2*NSek)";
      output << endl;
      output << std::setw(50) << std::left << "  Maximum Dimension of Iterative Subspace:"
           << this->maxSubSpace_ << endl;
      output << std::setw(50) << std::left << "  Maximum Number of Micro Iterations:"
           << this->maxIter_ << endl;
      output << std::setw(50) << std::left << "  Maximum Number of Macro Iterations:"
           << this->MaxIter_ << endl;
      output << std::setw(50) << std::left << "  Using an Hermetian algorithm?:";
      if(this->hermetian_) output << "Yes";
      else output << "No";
      output << endl;
 
      output << endl << bannerEnd << endl;
    }
    
  }; // class Davidson
} // namespace ChronusQ
