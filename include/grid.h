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
      virtual void                   transformPts() = 0; ///< virtual function to transform the integral interval
      virtual double                 integrate()    = 0; ///< virtual function to integrate
  }; // class Grid

class OneDGrid : public Grid {
   protected:
      double * gridPts_;
      double * weights_;
      double   f_val(double);
      double   norm_;
      std::array<double,2> range_;
   public:
      OneDGrid(
        int npts = 0, double beg = 0.0, double end = 0.0):
        Grid(npts){
          this->range_ = {beg, end};
          this->norm_ = (end-beg)/2.0;
      };
// deconstructor
      ~OneDGrid(){
         delete [] this->gridPts_;
         delete [] this->weights_;
       cout << "Deliting" <<endl;
       };
// access to protected data
       inline double * gridPts(){ return this->gridPts_;};
       inline double * weights(){ return this->weights_;};
       inline double norm(){ return this->norm_;};
       double integrate();
}; // Class OneGrid (one dimensional grid)

class twoDGrid : public Grid {
   protected:
      sph2GP * grid2GPts_;
      double * weights_;
   public:
      twoDGrid(
       int npts = 0):
        Grid(npts){
        };
// deconstructor
      ~twoDGrid(){
         delete [] this->grid2GPts_;
//         boost::geometry::clear(this->grid2GPts_);
         delete [] this->weights_;
       cout << "Deliting" <<endl;
       };
// access to protected data
//       inline double * get_grid2GPts_elev(){ 
//     elevation from [0,2PI]
//       double elev = bg::get<0>(this->grid2GPts_);
//       return &elev;};
//       inline double * get_grid2GPts_azim(){ 
//     azimuthal from [0,PI]
//       double azim = bg::get<1>(this->grid2GPts_);
//       return &azim;};
// Functions to inizialize data
//       void * set_grid2GPPts_elev(double * val){
//        bg::set<0>((this->grid2GPts_), (* val));}; 
//       void * set_grid2GPPts_azim(double * val){
//        bg::set<1>((this->grid2GPts_), (* val));}; 
//   2D integration
       double integrate();
       void   transformPts();
}; // Class twoD Grid


class LebedevGrid : public twoDGrid {
    public:
      LebedevGrid(
        int npts = 0):
        twoDGrid(npts){
          this->grid2GPts_ = new sph2GP[this->nPts_];  //< Lebedev polar coordinates [
          this->weights_   = new double[this->nPts_];
        };
      void genGrid();
      void gen6_A1(int num, long double a, long double v);
      void gen12_A2(int num, long double a, long double v);
  }; // class LebedevGrid

/*
class threeDGrid : public Grid {
   protected:
      bg::model::point< double, 3,bg::cs::cartesian> * gridPts_;
//      double * Pcart_;
      double * weights_;
      double   f_val(double);
      double   norm_;
      std::array<double,2> range_;
//   public:
      threeDGrid(
       int npts = 0, double beg = 0.0, double end = 0.0):
        Grid(npts){
          this->range_ = {beg, end};
          this->norm_ = (end-beg)/2.0;
      };
// deconstructor
      ~threeDGrid(){
         delete [] this->gridPts_;
         delete [] this->weights_;
       cout << "Deliting" <<endl;
       };
// access to protected data
       inline double * get_gridPts_X(){ 
       double xval = bg::get<0>(this->gridPts_);
       return &xval;};
       inline double * get_gridPts_Y(){ 
       double yval = bg::get<1>(this->gridPts_);
       return &yval;};
       inline double * get_gridPts_Z(){ 
       double zval = bg::get<2>(this->gridPts_);
       return &zval;};
       inline double * weights(){ return this->weights_;};
       inline double norm(){ return this->norm_;};
// Functions to inizialize data
       void * setX_gridPts(double * val){
        bg::set<0>((this->gridPts_), (* val));}; 
       void * setY_gridPts(double * val){
        bg::set<1>((this->gridPts_), (* val));}; 
       void * setZ_gridPts(double * val){
        bg::set<2>((this->gridPts_), (* val));}; 
//   3D integration
       double integrate();
}; // Class 3DGrid (three dimensional grid)
*/

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
//    double integrate();
//    double f_val( double rad);
  }; // class GaussChebyshev1stGrid

}; // namespace ChronusQ
