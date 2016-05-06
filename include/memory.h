#include <global.h>
#include <boost/pool/simple_segregated_storage.hpp>
#include <boost/pool/object_pool.hpp>

#ifndef INCLUDED_CQMEM
#define INCLUDED_CQMEM
namespace ChronusQ {
class CQMemManager : public boost::simple_segregated_storage<std::size_t> {
  std::size_t N_;
  std::size_t NAlloc_;
  std::size_t NBlockSize_;
  std::vector<char> V_;

  inline void fixBlockNumber(){ 
    if(N_ % NBlockSize_ != 0) N_ += NBlockSize_ - (N_ % NBlockSize_);
  };

  public:
    bool isAllocated;

    CQMemManager(std::size_t N = 0, std::size_t BlockSize = 256) : 
      boost::simple_segregated_storage<std::size_t>(), 
      N_(N), NBlockSize_(BlockSize), NAlloc_(0), isAllocated(false){ 
        if(N_ != 0 && BlockSize != 0) this->allocMem(); 
      }

    ~CQMemManager(){ };

    void allocMem(){
      if(isAllocated) return;

      this->fixBlockNumber();
      this->V_ = std::vector<char>(N_);
      this->add_block(&this->V_.front(),this->V_.size(),this->NBlockSize_);
      this->isAllocated = true;
//    cout << "Creating CQ Memory Partition of " << N_ << " Bytes " <<
//     " startring at "; 
//    cout << (int*) &this->V_[0] << endl;
    } 

    template<typename T>
    T * malloc(std::size_t n){
//    cout << "Allocating " << n << " works of " << typeid(T).name()
//      << " data, which is ";
      std::size_t nBlocks = ( (n-1) * sizeof(T) ) / this->NBlockSize_ + 1;
//    cout << nBlocks << " blocks of data" << endl;
      this->NAlloc_ += nBlocks * this->NBlockSize_;

      if(this->NAlloc_ > this->N_) {
        std::bad_alloc excp;
        throw excp;
      };

      void * ptr = 
          boost::simple_segregated_storage<std::size_t>::malloc_n(
            nBlocks,this->NBlockSize_);
//    cout << ptr << endl;
//    cout << static_cast<T*>(ptr) << endl;

      return static_cast<T*>(ptr);
//    return static_cast<T*>(
//        boost::simple_segregated_storage<std::size_t>::malloc_n(
//          1,nBlocks * this->NBlockSize_)
//        );
    };

    template<typename T>
    void free( T * ptr, std::size_t n){
//    cout << "Freeing " << n << " works of " << typeid(T).name()
//      << " data, which is ";
      std::size_t nBlocks = ( (n-1) * sizeof(T) ) / this->NBlockSize_ + 1;
//      cout << nBlocks << " blocks of data" << endl;
      this->NAlloc_ -= nBlocks * this->NBlockSize_;
      boost::simple_segregated_storage<std::size_t>::free_n(
          ptr,nBlocks,this->NBlockSize_);
    }


    inline std::size_t NBlocksAllocated(){ 
      return this->NAlloc_ / this->NBlockSize_; 
    };

    inline std::size_t NBlocksFree(){ 
      return (this->N_ - this->NAlloc_) / this->NBlockSize_; 
    };

    void printSummary(std::ostream &out){
      out << "Memory Allocation Summary:" << std::endl;

      out << std::setw(30) << "Total Memory Allocated: ";
      out << std::setw(15) << this->N_ << " B / "; 
      out << std::setw(8) << std::fixed << this->N_ / 1e3 << " kB / "; 
      out << std::setw(8) << std::fixed << this->N_ / 1e6 << " MB / "; 
      out << std::setw(8) << std::fixed << this->N_ / 1e9 << " GB"; 

      out << std::endl;

      out << std::setw(30) << "Block Size: ";
      out << std::setw(15) << this->NBlockSize_ << " B"; 
      out << " (" << this->N_ / this->NBlockSize_ << " Blocks)";

      out << std::endl;

      out << std::setw(30) << "Reserved Memory: ";
      out << std::setw(15) << this->NAlloc_ << " B"; 
      out << " (" << this->NAlloc_ / this->NBlockSize_ << " Blocks)";

      out << std::endl;

      out << std::setw(30) << "Free Memory: ";
      out << std::setw(15) << this->N_ - this->NAlloc_ << " B"; 
      out << " (" << (this->N_ -this->NAlloc_) / this->NBlockSize_ 
        << " Blocks)";

      out << std::endl;
      
    }

    inline void setTotalMem(std::size_t N){ if(!isAllocated) this->N_ = N; };
    inline void setBlockSize(std::size_t N){ if(!isAllocated) this->NBlockSize_ = N; };
};
}; // namespace ChronusQ
#endif
