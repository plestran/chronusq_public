/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
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
#include "cerr.h"

void CErr(int msg) {
  if(msg==3000) {
    cout << "FATAL: Cannot construct Matrix object with any zero dimension (Matrix.Matrix)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3001) {
    cout << "FATAL: Cannot store Matrix object as a packed lower triangle for non-square matricies (Matrix.Matrix)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3002) {
    cout << "FATAL: Error in Memory allocation for Matrix object (Matrix.Matrix)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3003) {
    cout << "FATAL: Error in Memory Allocation for Matrix object (Matrix.ioRead)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3004) {
    cout << "FATAL: Eigenvectors have not been calculated for this Matrix object (Matrix.ioRead)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3005) {
    cout << "FATAL: Eigenvalues have not been calculated for this Matrix object (Matrix.ioRead)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3006) {
    cout << "FATAL: Scalar Product Dimensions must match (Matrix.scalarProd)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3007) {
    cout << "FATAL: Trace only available for square Matrix Objects (Matrix.trace)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3008) {
    cout << "FATAL: Matrix dimensions are not compatible for Matrix Multiplication (Matrix.*)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3009) {
    cout << "FATAL: Left Eigenvectors have not been calculated for this Matrix object (Matrix.ioRead)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3010) {
    cout << "FATAL: Right Eigenvectors have not been calculated for this Matrix object (Matrix.ioRead)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3011) {
    cout << "FATAL: Matrix diagonalization only available for square matricies (Matrix.diag)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3012) {
    cout << "FATAL: Matrix diagonalization only available for type 'G' (Matrix.diag)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3013) {
    cout << "FATAL: Error Allocating Space for Eigenvetors (Matrix.diag)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3014) {
    cout << "FATAL: Error Allocating Space for Eigenalues (Matrix.diag)" << endl;
    exit(EXIT_FAILURE);
  } else if(msg==3016) {
    cout << "FATAL: Error Allocating Space for CPY for DGEEV (Matrix.diag)" << endl;
    exit(EXIT_FAILURE);
  };
}
