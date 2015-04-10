/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explictly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
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
