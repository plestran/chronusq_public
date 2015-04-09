#include "matrix.h"
using ChronusQ::Matrix;
//-------------//
// constructor //
//-------------//
namespace ChronusQ {
template<>
Matrix<double>::Matrix(int rows, int cols, char *nm, char *format) {
  // Initialize the Internal Pointers to NULL (to avoid problems later)
  fixMem();

  // Name the matrix
  strcpy(name_,nm);

  // Determine Storage Format
  strupr(format);
  if (!strcmp(format,"STD")) format_ = 0;
  else if(!strcmp(format,"LT")) format_ = 1;

  // Get Dimensions
  if (rows==0||cols==0) throw 3000; // No Zero matricies allowed
  rows_ = rows;
  cols_ = cols;
  if(format_==1 && rows_!=cols_) throw 3001; // Packed Maricies must be square
  vectorized_ = false;
  realCmplxEig_ = false;

  // Determine storage length (size of the data_ array)
  if(format_==0) len_ = rows_*cols_;
  else if(format_==1) len_ = rows_*(rows_+1)/2;

  // Allocate internal storage (data_)
  data_ = new (nothrow) double[len_];
  if (data_==NULL) throw 3002; // If there's no space, throw an error to be handled

  // Set up initial job variables
  haveEigen_ = 0;  // No eigensystem upon initialisation (i.e. not diagonalized)
  JOBVL_ = 'V';    // Compute left eigenvectors by default
  JOBVR_ = 'V';    // Compute right eigenvectors by default
  trans_ = 'N';    // Untransposed by default

  // Determine symmetry
  if(format_==0){      symm_ = 'G';}
  else if(format_==1){ symm_ = 'S';};

  // Try to determine full storage size of the matrix (buggy, hasn't been updated)
  char testChar;
  int  testInt;
  double testDouble;
  size_ = MAXNAMELEN*sizeof(testChar) + 6*sizeof(testInt)/sizeof(testChar) + len_*sizeof(testDouble)/sizeof(testChar);
  if(haveEigen_) size_ = size_ + (rows_+1)*rows_*sizeof(testDouble)/sizeof(testChar);
};
}
