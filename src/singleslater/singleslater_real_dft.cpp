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
    this->totalEcorr = 0.0;
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
    (*this->vCorA()) =  4.0*math.pi * (*this->vCorA());
    this->totalEcorr   =  4.0*math.pi * (this->totalEcorr);
    if(!this->isClosedShell && this->Ref_ != TCS) (*this->vXCA()) =  std::pow(2.0,(1.0/3.0)) * (*this->vXCA()) ;
    if(!this->isClosedShell && this->Ref_ != TCS) (*this->vXCB()) =  std::pow(2.0,(1.0/3.0)) * val * (*this->vXCB()) ;
    if(!this->isClosedShell && this->Ref_ != TCS) this->totalEx = std::pow(2.0,(1.0/3.0)) * this->totalEx;
    prettyPrint(this->fileio_->out,(*this->vXCA()),"LDA Vx alpha");
    if(!this->isClosedShell && this->Ref_ != TCS) prettyPrint(this->fileio_->out,(*this->vXCB()),"LDA Vx beta");
    this->fileio_->out << "Total LDA Ex =" << this->totalEx <<endl;
    cout << "Total LDA Ex =" << this->totalEx << " Total VWN Corr= " << this->totalEcorr <<endl;
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

template<>
void SingleSlater<double>::formVWNPara(double rho){
// Parameter for the fit (according Eq 4.4 Vosko Can. J. Phys. 1980
   double A  = 0.0621814; // In the text page 1207 (A^P)
   double b  = 3.72744;  // Caption Table 5
   double c  = 12.9352;  // Caption Table 5
   double x0 = -0.10498; // Caption Table 5
   double b1 = (b*x0 - c)/(c*x0); 
   double b2 = (x0 - b)/(c*x0); 
   double b3 = (-1.0)/(c*x0); 
   double Q =std::pow((4.0*c - b*b),(1.0/2.0)); 
   double r_s = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/3.0));
   double r_s_sqrt = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/6.0));
   double r_s_32 = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/2.0));
   double X = r_s + b*(r_s_sqrt) + c; 
   double X_x0 = x0*x0 + b*(x0) + c; 
   this->eps_corr = 0.0;
   this->mu_corr  = 0.0;
   this->eps_corr = A *
          ( std::log(r_s/X) + 
            2.0*b*(std::atan(Q/(2.0*r_s_sqrt + b)))/Q  -
            b*x0*(std::log( (std::pow((r_s_sqrt-x0),2.0))/X ))/X_x0 +
            2.0*(b + 2.0*x0)*(std::atan(Q/(2.0*r_s_sqrt + b)))/ Q           );
   this->mu_corr = -A * ( (1.0 + b1*r_s_sqrt)/(1.0 + b1*r_s_sqrt + b2*r_s + b3*r_s_32)) / 3.0 ;
   this->mu_corr += this->eps_corr ;
//   cout << "rho " << rho << " eps "<< this->eps_corr << " mu " << this->mu_corr <<endl;
};

template<>
void SingleSlater<double>::formVWNFerr(double rho_A, double rho_B){
};

} // Namespace ChronusQ
