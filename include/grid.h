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
       cout << "Deleting" <<endl;
       };
// access to protected data
       inline double * gridPts(){ return this->gridPts_;};
       inline double * weights(){ return this->weights_;};
       inline sph2GP * grid2GPts(){return this->grid2GPts_;};
       inline double norm(){ return this->norm_;};
       double integrate();
       void printGrid();
}; // Class OneGrid (one dimensional grid)

class TwoDGrid : public Grid {
      protected:
            OneDGrid * Gr_;
            OneDGrid * Gs_;
      double   fsphe(double,double,double);
//            int * Gsnpts_;
      public:
        TwoDGrid(OneDGrid *Gr, OneDGrid *Gs){
        this->Gr_ =  Gr;
        this->Gs_ =  Gs;
          };
      double integrate();
      void printGrid();
      void genGrid();
      void transformPts();
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
      void gen24_Cn(int num, double a, double v);
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

}; // namespace ChronusQ
