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
#ifndef INCLUDED_GRID
#define INCLUDED_GRID
#include <global.h>
namespace ChronusQ {

// Classes
class Grid {
   protected:
      int                            nPts_;  ///< number of grid points       
   public:
      Grid(int npts = 0){
        this->nPts_ = npts;
      };
      virtual void                   genGrid() = 0; ///< virtual function to generate the grid points
      virtual void                   printGrid() = 0; ///< virtual function to print the grid points
      virtual void                   transformPts() = 0; ///< virtual function to transform the integral interval
      virtual double                 integrate()    = 0; ///< virtual function to integrate
      int npts(){return this->nPts_;};
  }; // class Grid

class OneDGrid : public Grid {
   protected:
      double * gridPts_;
      sph2GP * grid2GPts_;
      double * weights_;
      double   f_val(double);
      double   f2_val(double, double);
      double   norm_;
      std::array<double,2> range_;
      bool intas2GPt_ = false;
   public:
      OneDGrid(
        int npts = 0, double beg = 0.0, double end = 0.0):
        Grid(npts){
          this->range_ = {beg, end};
          this->norm_ = (end-beg)/2.0;
      };
// deconstructor
      ~OneDGrid(){
         if(intas2GPt_){
         delete [] this->grid2GPts_;
         }else{
         delete [] this->gridPts_;
         }
         delete [] this->weights_;
//       cout << "Deleting" <<endl;
       };
// access to protected data
       inline double * gridPts(){ return this->gridPts_;};
       inline double gridPts(int i){ return this->gridPts_[i];};
       inline double * weights(){ return this->weights_;};
       inline sph2GP * grid2GPts(){return this->grid2GPts_;};
       inline sph2GP grid2GPts(int i){return this->grid2GPts_[i];};
       inline double norm(){ return this->norm_;};
       double integrate();
       void printGrid();
}; // Class OneGrid (one dimensional grid)

class TwoDGrid : public Grid {
      protected:
            OneDGrid * Gr_;        ///< pointer to Radial OneD grid
            OneDGrid * Gs_;        ///< pointer to Angular OneD grid
            double *  GridCarX_;  ///<  Cartesian X component of Grid points
            double *  GridCarY_;  ///<  Cartesian Y component of Grid points
            double *  GridCarZ_;  ///<  Cartesian Z component of Grid points
            double   *   weightsGrid_; ///< weights
      public:
//      Constructor
        TwoDGrid(
          int npts = 0,OneDGrid *Gr = NULL, OneDGrid *Gs = NULL):
          Grid(npts){
//      Pointers
          this->Gr_ =  Gr;
          this->Gs_ =  Gs;
          this->GridCarX_     = new double [npts];    ///< x component
          this->GridCarY_     = new double [npts]; 
          this->GridCarZ_     = new double [npts]; 
          this->weightsGrid_  = new double [npts];
          };
//    Function Declaration //
//    Grid ones 
      double integrate();
      void   printGrid();
      void   genGrid();
      void   centerGrid (double cx, double cy, double cz);
      void   transformPts();
//    TwdDgrid ones
      double * BuildDensity(double * Sum,double * Buff,int n1, int n2);
      void iniTwoDGrid(int Ngridr, int NLeb);
      double voronoii(double mu);
      double step_fun(double mu);
      inline double * weightsGrid(){ return this->weightsGrid_;};
      inline double getweightsGrid(int i){ return this->weightsGrid_[i];};
      inline sph3GP gridPt(int i, int j){
        sph3GP x(bg::get<0>(Gs_->grid2GPts(j)),bg::get<1>(Gs_->grid2GPts(j)),Gr_->gridPts(i));
        return x;
      };
      inline cartGP gridPtCart(int ipts){
        cartGP pt ( this->GridCarX_[ipts],this->GridCarY_[ipts],this->GridCarZ_[ipts]);
        return pt;
      };
      inline void SetgridPtCart(int ipts, double x, double y, double z){
        this->GridCarX_[ipts] = x;
        this->GridCarY_[ipts] = y;
        this->GridCarZ_[ipts] = z;
      };
//    Deconstructors //
      ~TwoDGrid(){
        delete [] this->weightsGrid_;
        delete [] this->GridCarX_;
        delete [] this->GridCarY_;
        delete [] this->GridCarZ_;
     };
}; //   Class TwoDGrid

class LebedevGrid : public OneDGrid {
    public:
      LebedevGrid(
        int npts = 0):
        OneDGrid(npts,0.0,0.0){
          this->grid2GPts_ = new sph2GP[this->nPts_];  //< Lebedev polar coordinates [
          this->weights_   = new double[this->nPts_];
          this->intas2GPt_ = true;
        };
      void transformPts();
      void genGrid();
      void gen6_A1(int num, double a, double v);
      void gen12_A2(int num, double a, double v);
      void gen8_A3(int num, double a, double v);
      void gen24_Cn(int num, double q, double v);
      void gen24_Bn(int num, double l, double v);
      void gen48_Dn(int num, double u, double r, double v);
  }; // class LebedevGrid

  class GaussChebyshev1stGrid : public OneDGrid {
    public:
      GaussChebyshev1stGrid(
        int npts = 0, double beg = 0.0, double end = 0.0):
        OneDGrid(npts,beg,end){
          this->gridPts_ = new double[this->nPts_];        ///< Grid Points
          this->weights_ = new double[this->nPts_];        ///< Weights
        };
  // Class Functions
    void genGrid();                                      
    void transformPts();
  }; // class GaussChebyshev1stGrid

  class GaussChebyshev1stGridInf : public OneDGrid {
    public:
      GaussChebyshev1stGridInf(
        int npts = 0, double beg = 0.0, double end = 0.0):
        OneDGrid(npts,beg,end){
          this->gridPts_ = new double[this->nPts_];        ///< Grid Points
          this->weights_ = new double[this->nPts_];        ///< Weights
        };
  // Class Functions
    void genGrid();                                      
    void transformPts();
  }; // class GaussChebyshev1stGridInf

}; // namespace ChronusQ
#endif
