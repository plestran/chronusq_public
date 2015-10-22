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
    int nAtom    = this->molecule_->nAtoms();    // Number of Atoms
    int nRad     = 100;                          // Number of Radial grid points for each center
    int nAng     = 302;                          // Number of Angular grid points for each center (only certain values are allowed - see grid.h)
    int npts     = nRad*nAng;                    // Total Number of grid point for each center
    double weight= 0.0;                            
    double CxVx = -(std::pow((3.0/math.pi),(1.0/3.0)));    //TF LDA Prefactor
    double CxEn =  (3.0/4.0);    //TF LDA Prefactor
    double val = 4.0*math.pi*CxVx;
    this->totalEx = 0.0;
//  Generating grids (Raw grid, it has to be centered and integrated over each center and centered over each atom)
    GaussChebyshev1stGridInf Rad(nRad,0.0,1.0);   // Radial Grid
    LebedevGrid GridLeb(nAng);                    // Angular Grid
    GridLeb.genGrid();
//    Rad.genGrid();                               
//    TwoDGrid Raw3Dg(npts,&Rad,&GridLeb);          // Final Raw (not centered) 3D grid (Radial times Angular grid)
    this->vXCA()->setZero();                      // Set to zero every occurence of the SCF
    this->vCorA()->setZero();                      // Set to zero every occurence of the SCF
    if(!this->isClosedShell && this->Ref_ != TCS) this->vXCB()->setZero();
//  Loop over each centers (Atoms) (I think can be distribuited over different cpus)
    for(int iAtm = 0; iAtm < nAtom; iAtm++){
//    cout << "Atmoic Radius Slater " << elements[this->molecule_->index(iAtm)].sradius << endl; 
//    The Radial grid is generated and scaled for each atom
      Rad.genGrid();
      Rad.scalePts((elements[this->molecule_->index(iAtm)].sradius)) ; 
      TwoDGrid Raw3Dg(npts,&Rad,&GridLeb);          // Final Raw (not centered) 3D grid (Radial times Angular grid)
      //Center the Grid at iAtom
      Raw3Dg.centerGrid((*this->molecule_->cart())(0,iAtm),(*this->molecule_->cart())(1,iAtm),(*this->molecule_->cart())(2,iAtm));
//    Loop over grid points
      for(int ipts = 0; ipts < npts; ipts++){
//    Evaluate each Becke fuzzy call weight, normalize it and muliply by the Rwa grid weight at that point
        weight = Raw3Dg.getweightsGrid(ipts)  
                 * (this->formBeckeW((Raw3Dg.gridPtCart(ipts)),iAtm))
                 / (this->normBeckeW(Raw3Dg.gridPtCart(ipts))) ;
//    Build the Vxc for the ipts grid point (Vxc will be ready at the end of the two loop
        this->buildVxc((Raw3Dg.gridPtCart(ipts)),weight);
        } //end loop over Raw grid points
    } // end loop natoms
//  Finishing the Vxc using the TF factor and the integration prefactor over a solid sphere
    (*this->vXCA()) =  val * (*this->vXCA());
    this->totalEx   =  val * CxEn * (this->totalEx);
    if(!this->isClosedShell && this->Ref_ != TCS) (*this->vXCA()) =  std::pow(2.0,(1.0/3.0)) * (*this->vXCA()) ;
    if(!this->isClosedShell && this->Ref_ != TCS) (*this->vXCB()) =  std::pow(2.0,(1.0/3.0)) * val * (*this->vXCB()) ;
    if(!this->isClosedShell && this->Ref_ != TCS) this->totalEx = std::pow(2.0,(1.0/3.0)) * this->totalEx;
    prettyPrint(this->fileio_->out,(*this->vXCA()),"LDA Vx alpha");
    if(!this->isClosedShell && this->Ref_ != TCS) prettyPrint(this->fileio_->out,(*this->vXCB()),"LDA Vx beta");
    this->fileio_->out << "Total LDA Ex =" << this->totalEx <<endl;
//    cout << "TotalEx "  << this->totalEx  <<endl;
//    cout << "Factor "  << CxVx*CxEn*std::pow(2.0,(1.0/3.0))  <<endl;
//    cout << "Vx_A" <<endl;
//    cout << (*this->vXCA()) <<endl;
//    cout << "Vx_B" <<endl;
//    cout << (*this->vXCB()) <<endl;

//  Comment to avoid the printing
/*
  double Energy;
  Energy = CxEn*((*this->vXCA_).frobInner(this->densityA_->conjugate()));
  std::cout.precision(10);
  cout << " E_XC = " << Energy <<endl;
  if(!this->isClosedShell && this->Ref_ != TCS) {
  (*this->vXCB()) =  (*this->vXCA());
  double EnergyB;
  EnergyB = CxEn*((*this->vXCB_).frobInner(this->densityB_->conjugate()));
  std::cout.precision(10);
  cout << " E_XCiB = " << EnergyB <<endl;
  cout << " E_Xalpha_beta = " << Energy+EnergyB <<endl;
  }
*/
}; //End

} // Namespace ChronusQ
