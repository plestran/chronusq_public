#include <pythonapi.h>
#include <global.h>
namespace ChronusQ {
  void Wrapper_CErr_Default(FileIO &fileio) {
    CErr(fileio.out);
  }
  void Wrapper_CErr_Message(FileIO &fileio, std::string msg){
    CErr(msg,fileio.out);
  }
}
