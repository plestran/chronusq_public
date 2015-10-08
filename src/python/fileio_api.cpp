#include <workers.h>
#include <global.h>

using ChronusQ::SingleSlater;
using ChronusQ::Molecule;
using ChronusQ::Atoms;
using ChronusQ::FileIO;
using ChronusQ::BasisSet;
using ChronusQ::Controls;
using ChronusQ::AOIntegrals;

namespace ChronusQ {
  void FileIO::write(std::string str){
    this->out << str << std::endl; 
  }
};


