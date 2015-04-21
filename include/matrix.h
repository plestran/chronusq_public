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
#ifndef INCLUDED_MATRIX
#define INCLUDED_MATRIX
#include "global.h"
#include "fileio.h"
#include "tools.h"
#include "clapack.h"
#include "cerr.h"

/****************************/
/* Error Messages 3000-3999 */
/****************************/
namespace ChronusQ {

template <class T>
class Matrix {
  /*****************
   *  Member Data  *
   *****************/
  int       rows_;             // Number of rows (leading dimension) to be used by routines
  int       cols_;             // Number of columns (slowest running index) to be used by routies
  int       rows_t_;           // Number of rows (leading dimension) upon allocation (permanent)
  int       cols_t_;           // Number of columns (slowest running index) upon allocation (permanent)
  int       len_;              // Number of elements being stored in the matrix (length of data array)
  int       format_;           // Storage format                                  
                               //   0:  Standard dense matrix (Column-Major)      
                               //   1:  Packed lower triangle (Column-Major)      
                               //   2:  Packed upper triangle (Column-Major) (NYI)
                               //   3:  Banded                (Column-Major) (NYI)
  int       size_;             // Size of the matrix object (all internal storage, buggy)
  int       haveEigen_;        // Level of eigensystem storage (compounds)
                               //   0:  No eigensystem is stored
                               //   1:  Eigenvalues are stored
                               //   2:  Right Eigenvectors are stored
                               //   3:  Left Eigenvectors are stored  (If applicable)
                               //   4:  Eigenvalues are broken into   (If applicable)
  bool      vectorized_;       // Boolean as to whether or not to treat the matrix as a vector
  bool      realCmplxEig_;     // Complex eigenvalues of a real matrix
  char      name_[MAXNAMELEN]; // Name of the Matrix
  char      trans_;            // Transpose flag, used for LAPACK routines
  char      symm_;             // Symmetry of Matrix
                               //  G: Nonsymmetric
                               //  S: Symmetric                        (NYI)
                               //  A: Antisymmetric                    (NYI)
                               //  H: Hermetian     (only for complex) (NYI)
                               //  X: Antihermetian (only for complex) (NYI)
  char      JOBVL_;            // Whether or not to compute Left Eigenvectors  (LAPACK)
  char      JOBVR_;            // Whether or not to compute Right Eigenvectors (LAPACK)
  T         *data_;            // Internal element storage
  T         *eigenvector_;     // Default eigenvector storage (will point to eigenvector_r_)
  T         *eigenvector_r_;   // Default right eigenvector storage
  T         *eigenvector_l_;   // Default left eigenvector storage (will point to eigenvector_r_ for 'S' and 'H')
  T         *eigenvalue_;      // Detault eigenvalue storage (will point to eigenvalue_re_ for double)
  double    *eigenvalue_re_;   // Real part of the eigenvalue vector
  double    *eigenvalue_im_;   // Imag part of the eigenvalue vector
  dcomplex  *eigenvaluez_;     // Complex cast of the eigenvalue vector

  // Private Functions
  inline void fixMem(){
    data_ = NULL;
    eigenvector_ = NULL;
    eigenvalue_ = NULL;
    eigenvector_r_ = NULL;
    eigenvector_l_ = NULL;
    eigenvalue_re_ = NULL;
    eigenvalue_im_ = NULL;
    eigenvaluez_ = NULL;
  };
  friend void doTNT(Matrix<T>*,const Matrix<T>*, const Matrix<T>*);
  friend void doTTN(Matrix<T>*,const Matrix<T>*, const Matrix<T>*);
  void doDiag(Matrix<T>*);
  void allocEigen();
  void cleanEigen(){
    delete[] eigenvector_;
    delete[] eigenvalue_;
    delete[] eigenvector_l_;
    delete[] eigenvector_r_;
    delete[] eigenvalue_re_;
    delete[] eigenvalue_im_;
    delete[] eigenvaluez_;
  }
  T**  allocEigScr(bool,int&,Matrix*);
  int  getLWORK(bool);
public:

  // Constructor
  Matrix(int, int, char *nm="No Name", char *format="STD");
  // Destrictor
  ~Matrix() {
/*  // This shouldn't be a necessary check, but I'll leave it here
    // incase we run into problems later (should be fixed by fixMem)
    if(         data_!=NULL) delete[] data_;
    if(  eigenvector_!=NULL) delete[] eigenvector_;
    if(   eigenvalue_!=NULL) delete[] eigenvalue_;
    if(eigenvector_l_!=NULL) delete[] eigenvector_l_;
    if(eigenvector_r_!=NULL) delete[] eigenvector_r_;
    if(eigenvalue_re_!=NULL) delete[] eigenvalue_re_;
    if(eigenvalue_im_!=NULL) delete[] eigenvalue_im_;
    if(  eigenvaluez_!=NULL) delete[] eigenvaluez_;
*/
    delete[] data_;
    delete[] eigenvector_;
    delete[] eigenvalue_;
    delete[] eigenvector_l_;
    delete[] eigenvector_r_;
    delete[] eigenvalue_re_;
    delete[] eigenvalue_im_;
    delete[] eigenvaluez_;
  };

  /***************
   *  Matrix IO  *
   ***************/
  // **DBWY Keep in mind that we want to print these to file at some point**
  inline void printDim(ostream &output=cout) {output<<this->name_<<" is a "<<this->rows_<<"(row)x"<<this->cols_<<"(column) matrix."<<endl;};
  void printAll(int list=5, ostream &output=cout);
  inline void printMem(){
    cout << endl;
    cout << "Memory Allocation Dump for Matrix Object: " << this->name_ << endl;
    cout << std::setw(15) << std::left << "data_" << this->data_;
    if(data_==NULL) cout << " (Not Allocated)"; cout << endl;
    cout << std::setw(15) << std::left << "eigenvector_" << this->eigenvector_;
    if(eigenvector_==NULL) cout << " (Not Allocated)"; cout << endl;
    cout << std::setw(15) << std::left << "eigenvalue_" << this->eigenvalue_ ;
    if(eigenvalue_==NULL) cout << " (Not Allocated)"; cout << endl;
    cout << std::setw(15) << std::left << "eigenvector_l_" << this->eigenvector_l_ ;
    if(eigenvector_l_==NULL) cout << " (Not Allocated)"; cout << endl;
    cout << std::setw(15) << std::left << "eigenvector_r_" << this->eigenvector_r_ ;
    if(eigenvector_r_==NULL) cout << " (Not Allocated)"; cout << endl;
    cout << std::setw(15) << std::left << "eigenvalue_re_" << this->eigenvalue_re_ ;
    if(eigenvalue_re_==NULL) cout << " (Not Allocated)"; cout << endl;
    cout << std::setw(15) << std::left << "eigenvalue_im_" << this->eigenvalue_im_ ;
    if(eigenvalue_im_==NULL) cout << " (Not Allocated)"; cout << endl;
    cout << std::setw(15) << std::left << "eigenvaluez_" << this->eigenvaluez_  ;
    if(eigenvaluez_==NULL) cout << " (Not Allocated)"; cout << endl;
    
  }
  // read|write scratch|binary files
  void ioRead (ChronusQ::FileIO*,int,char*,int charOffset=-1); // NYI Complex
  void ioWrite(ChronusQ::FileIO*,int,char*,int charOffset=-1); // NYI Complex

  // Zero out the Matrix data_
  inline void clearAll(ostream &output=cout) { memset(data_,0,len_*sizeof(T));};

  /*****************************************
   *  Getters (Access to private storage)  *
   *****************************************/
  inline int    len(){            return this->len_;          }; // Return the number of elements stored for matrix
  inline int   rows(){            return this->rows_;         }; // Return the number of rows being used in current treatment
  inline int   cols(){            return this->cols_;         }; // Return the number of cols being used in current treatment
  inline char  symm(){            return this->symm_;         }; // Return the symmetry currently being used
  inline T    *data(){            return this->data_;         }; // Return the internal storage array        (template)
  inline int   size(){            return this->size_;         }; // Return the size of the matrix object     (buggy)
  inline int format(){            return this->format_;       }; // Return the storage format

/*
  inline T   *eigenvector(){      if(this->eigenvector_  !=NULL) return this->eigenvector_;   else throw 3102;}; // Return the defualt eigenvector array    
  inline T *eigenvector_r(){      if(this->eigenvector_r_!=NULL) return this->eigenvector_r_; else throw 3101;}; // Return the right eigenvector array      
  inline T *eigenvector_l(){      if(this->eigenvector_l_!=NULL) return this->eigenvector_l_; else throw 3100;}; // Return the left eigenvector array       
  inline T    *eigenvalue(){      if(this->eigenvalue_   !=NULL) return this->eigenvalue_;    else throw 3103;}; // Return the default eigenvalue array     
  inline double* eigenvalue_re(){ if(this->eigenvalue_re_!=NULL) return this->eigenvalue_re_; else throw 3104;}; // Return the real part of the eigenvalues 
  inline double* eigenvalue_im(){ if(this->eigenvalue_im_!=NULL) return this->eigenvalue_im_; else throw 3105;}; // Return the imag part of the eigenvalues 
*/
  inline dcomplex* eigenvaluez(){  if(this->eigenvaluez_  !=NULL) return this->eigenvaluez_;   else throw 3106;}; // Return the eigenvalues in complex format

  inline T *eigenvector_r(){                                // Return the right eigenvector array
    if(haveEigen_==0) throw 3000;
    if(this->symm_=='G') return this->eigenvector_r_;
    if(this->symm_=='S'||this->symm_=='H') return this->eigenvector_;
  };
  inline T *eigenvector_l(){                                // Return the left eigenvector array
    if(haveEigen_==0) throw 3000;
    if(this->symm_=='G') return this->eigenvector_l_;
    if(this->symm_=='S'||this->symm_=='H') return this->eigenvector_;
  };
  inline T *eigenvector(){ return this->eigenvector_r();};  // Return the defualt eigenvector array

  T *eigenvalue();                                   // Return the default eigenvalue array
  double* eigenvalue_re();
  double* eigenvalue_im();
    
    

    

  /************************************
   *  Setters (Edit private storage)  *
   ************************************/
  template <class U> void cpyData(U*);                                          // Copy contents of an array to internal storage    (template)
  inline             void setName(char *x){ strcpy(this->name_,x);     };       // Set the name of the matrix object
//inline             void setSymm(char  x){ this->symm_ = x; this->fixSymm();}; // Set the symmetry of the matrix object
  inline             void setSymm(char  x){ this->symm_ = x;};                   // Set the symmetry of the matrix object
  template <class U> void  setDag(U*);                                          // Set the diagonal elements of the matrix to array (template)

  // Overload of (i,j) operator
  inline T& operator()(int row, int col){  
    if (row>this->rows_||col>this->cols_||row<0||col<0) throw 3003;
    if (this->format_==0) return this->data_[col*(this->rows_)+row];
    else if (this->format_==1&&row>=col) return this->data_[col*(this->rows_)-col*(col-1)/2+row-col];
    else if (this->format_==1&&row<col) return this->data_[row*(this->rows_)-row*(row-1)/2+col-row];
  };


  /******************************/
  /* Mathematical Manipulations */
  /******************************/
  // Structs
  struct Product{ const Matrix *a; const Matrix *b;};  // Struct for Multiplication
  struct     Sum{ const Matrix *a; const Matrix *b;};  // Struct for Addition
  struct    Diff{ const Matrix *a; const Matrix *b;};  // Struct for Subtraction
  struct     TTN{ const Matrix *a; const Matrix *x;};  // Struct for (XT)AX transform
  struct     TNT{ const Matrix *a; const Matrix *x;};  // Struct for XA(XT) transform
  struct     PWR{ const Matrix *a; const T *x;};  // Struct for A^x
  struct     EXP{ const Matrix *a;};                   // Struct for exp(A)


  // Operator Overload
  void    operator=(const Matrix*);                    // Assignment operator for Copy (Matrix)
  void    operator=(T *);                              // Assignment operator for Copy (Array)
  void    operator+=(const Matrix*);                   // Increment function
  void    operator-=(const Matrix *);                  // Decrement function
  void    operator*=(const Matrix *);                  //                              (NYI)
  void    operator=(const     Sum&);                   // Assignment operator for Summation
  void    operator=(const    Diff&);                   // Assignment operator for Subtraction
  void    operator=(const Product&);                   // Assignment operator for Multiplication
  void    operator=(const     TTN&);                   // Assignment operator for (XT)AX transform
  void    operator=(const     TNT&);                   // Assignment operator for XA(XT) transform
  void    operator=(const     PWR&);                   // Assignment operator for A^x 
  void    operator=(const     EXP&);                   // Assignment operator for exp(A)
  Sum     operator+(const Matrix&) const;              // Summation
  Diff    operator-(const Matrix&) const;              // Subtraction
  Product operator*(const Matrix&) const;              // Multiplication
  PWR     operator^(const T&     ) const;              // A^x

  // Function overload
  friend EXP exp(const Matrix&); // exp(A)

  // Misc functions
  inline void scale(T x){ for(int i=0;i<this->len_;i++) this->data_[i] = this->data_[i]*x;};
  void add(const Matrix*, const Matrix*);           // Addition function            (Obsolete)
  void sub(const Matrix*, const Matrix*);           // Subtraction function         (Obsolete)
  void scaleDag(T);
  void scaleODag(T);
  void inversion();         // Inversion           (NYI)
//void diag();              // Diagonalization
  void diag(Matrix *B=NULL);      // GEP Diagonalization
  T trace();                // Trace
  T scalarProd(Matrix*);    // Scalar Product
  inline void transpose(){
    if(this->trans_=='N'){ this->trans_='T';}
    else{                  this->trans_='N';};
  };
  inline void adjoint(){
    if(this->trans_=='N'){ this->trans_='C';}
    else{                  this->trans_='N';};
  };
  void transposeHard();
  void adjointHard();
  inline void vectorize(){
    this->rows_t_ = this->rows_;
    this->cols_t_ = this->cols_;
    this->rows_ = this->len_;
    this->cols_ = 1;
  };
  inline void unvectorize(){
    this->rows_ = this->rows_t_;
    this->cols_ = this->cols_t_;
  };
  void pack();
  void unpack();
  void norm_col();
  TNT transNT(const Matrix &);
  TTN transTN(const Matrix &);
  void eSort();
  T infNorm();
//void cleanEigen();
//void allocEigen();
//int  iniScratch();
//int eigLWORK(int, char); 
  //dbwye

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagMatrix);
  void mpiRecv(int,int tag=tagMatrix);
//dbwys
  /************************/
  /*  Gen Matrix Routines */
  /************************/
//dbwye
};

} // namespace ChronusQ
#endif
