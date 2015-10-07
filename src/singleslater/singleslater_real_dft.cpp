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
#include <singleslater.h>
namespace ChronusQ {
//----------------------------//
// form the Vxc matrix        //
//----------------------------//

template<>
void SingleSlater<double>::formVXC(){
//  Right now is just a place holder (we print the overlap, the actual
//  integration is perfomed in the grid class. see scr/grid/grid.cpp
//  Variable Definitions:
    int nAtom    = this->molecule_->nAtoms();    // Number of Atoms
    int nRad     = 100;                          // Number of Radial grid points for each center
    int nAng     = 194;                          // Number of Angular grid points for each center (only certain values are allowed - see grid.h)
    int npts     = nRad*nAng;                    // Total Number of grid point for each center
    int ipts     = 0;                            
    double weight= 0.0;                            
    double Cx = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));    //TF LDA Prefactor
    double val = 4.0*math.pi*Cx;
//  Generating grids (Raw grid, it has to be centered and integrated over each center
    GaussChebyshev1stGridInf Rad(nRad,0.0,1.0);
    LebedevGrid GridLeb(nAng);
    Rad.genGrid();
    GridLeb.genGrid();
    TwoDGrid Raw3Dg(npts,&Rad,&GridLeb);
    this->vXCA()->setZero(); 
    cout << "Erased Vxc term " << endl;
//  Loop over each centers (Atoms) (I think can be distribuited over different cpus)
    for(int iAtm = 0; iAtm < nAtom; iAtm++){
      Raw3Dg.centerGrid((*this->molecule_->cart())(0,iAtm),(*this->molecule_->cart())(1,iAtm),(*this->molecule_->cart())(2,iAtm));
     for(int ipts = 0; ipts < npts; ipts++){
        weight = Raw3Dg.getweightsGrid(ipts)  
                 * (this->formBeckeW((Raw3Dg.gridPtCart(ipts)),iAtm))
                 / (this->normBeckeW(Raw3Dg.gridPtCart(ipts))) ;
        this->buildVxc((Raw3Dg.gridPtCart(ipts)),weight);
        }
    } // end loop natoms
    (*this->vXCA()) =  val * (*this->vXCA());
  double Energy;
  Energy = (*this->vXCA_).frobInner(this->densityA_->conjugate());
  std::cout.precision(10);
  cout << " E_XC = " << Energy <<endl;
  cout << "Single Slater Numeric : Print" <<endl;
  cout << (*this->vXCA_)  << endl;
}; //End

template<>
void SingleSlater<double>::EnVXC(){
  double Energy;
  double resLDA = -11.611162519357;
  Energy = (*this->vXCA_).frobInner(this->densityA_->conjugate());
  std::cout.precision(10);
  cout << " E_XC = " << Energy <<endl;
  cout << "LDA Err " << (Energy-resLDA) << endl;
  cout << "Single Slater Numeric : Print" <<endl;
  cout << (*this->vXCA_)  << endl;
};



} // Namespace ChronusQ
