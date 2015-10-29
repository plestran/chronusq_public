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
    double CxVx = -(std::pow((3.0/math.pi),(1.0/3.0)));    //TF LDA Prefactor (for Vx)
    double CxEn =  (3.0/4.0);                              //TF LDA Prefactor to finish the X-Energy
    double val = 4.0*math.pi*CxVx;
    this->totalEx = 0.0;                                  // Total Exchange Energy
    this->totalEcorr = 0.0;                               // Total Correlation Energy
    this->cor_type = 1;                                  // Define Correlation Type (1 for VWN5 and 2 for VWN3)
//  Generating grids (Raw grid, it has to be centered and integrated over each center and centered over each atom)
    GaussChebyshev1stGridInf Rad(nRad,0.0,1.0);   // Radial Grid
    LebedevGrid GridLeb(nAng);                    // Angular Grid
    GridLeb.genGrid();                            // Generate Angular Grid (it alwais the same for all atoms)
    this->vXA()->setZero();                       // Set to zero every occurence of the SCF
    this->vCorA()->setZero();                     // Set to zero every occurence of the SCF
    if(!this->isClosedShell && this->Ref_ != TCS) this->vXB()->setZero();
    if(!this->isClosedShell && this->Ref_ != TCS) this->vCorB()->setZero();
//  Loop over each centers (Atoms) (I think can be distribuited over different cpus)
    for(int iAtm = 0; iAtm < nAtom; iAtm++){
//    The Radial grid is generated and scaled for each atom
      Rad.genGrid();
      Rad.scalePts((elements[this->molecule_->index(iAtm)].sradius)) ;  // Scale the grid according the Atomic Bragg-Slater Radius 
      TwoDGrid Raw3Dg(npts,&Rad,&GridLeb);                              // Final Raw (not centered) 3D grid (Radial times Angular grid)
      //Center the Grid at iAtom
      Raw3Dg.centerGrid((*this->molecule_->cart())(0,iAtm),(*this->molecule_->cart())(1,iAtm),(*this->molecule_->cart())(2,iAtm));
//    Loop over grid points
      for(int ipts = 0; ipts < npts; ipts++){
//    Evaluate each Becke fuzzy call weight, normalize it and muliply by the Raw grid weight at that point
        weight = Raw3Dg.getweightsGrid(ipts)  
                 * (this->formBeckeW((Raw3Dg.gridPtCart(ipts)),iAtm))
                 / (this->normBeckeW(Raw3Dg.gridPtCart(ipts))) ;
//    Build the Vxc for the ipts grid point (Vxc will be ready at the end of the two loop, to be finalized)
        this->buildVxc((Raw3Dg.gridPtCart(ipts)),weight);
        } //end loop over Raw grid points
    } // end loop natoms
//  Finishing the Vxc using the TF factor and the integration prefactor over a solid sphere
    (*this->vXA())    =  val * (*this->vXA());
    this->totalEx     =  val * CxEn * (this->totalEx);
    (*this->vCorA())  =  4.0 * math.pi * (*this->vCorA());
    if(!this->isClosedShell && this->Ref_ != TCS) (*this->vCorB())  =  4.0 * math.pi * (*this->vCorB());
    this->totalEcorr  =  4.0 * math.pi * (this->totalEcorr);
    if(!this->isClosedShell && this->Ref_ != TCS) (*this->vXA()) =  std::pow(2.0,(1.0/3.0)) * (*this->vXA()) ;  // For open shell averything has to be scaled by 2^(1/3)
    if(!this->isClosedShell && this->Ref_ != TCS) (*this->vXB()) =  std::pow(2.0,(1.0/3.0)) * val * (*this->vXB()) ;
//    if(!this->isClosedShell && this->Ref_ != TCS) this->totalEx = std::pow(2.0,(1.0/3.0)) * this->totalEx;
    prettyPrint(this->fileio_->out,(*this->vXA()),"LDA Vx alpha");
    prettyPrint(this->fileio_->out,(*this->vCorA()),"Vc alpha");
    if(!this->isClosedShell && this->Ref_ != TCS) prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vx beta");
    if(!this->isClosedShell && this->Ref_ != TCS) prettyPrint(this->fileio_->out,(*this->vCorB())," Vc beta");
    this->fileio_->out << "Total LDA Ex =" << this->totalEx << " Total VWN Corr= " << this->totalEcorr <<endl;
//    cout << "Total LDA Ex =" << this->totalEx << " Total VWN Corr= " << this->totalEcorr <<endl;
}; //End

template<>
void SingleSlater<double>::formCor(double rho, double spindensity){
// Parameter for the fit (according Eq 4.4 Vosko Can. J. Phys. 1980
   double A1  = 0.0621814; // In the text page 1207 (A^P)
   double A_p   = A1/2.0; // to Hartree
   double A_f   = A1/4.0; // to Hartree
   double A_a   = -(1.0/(6.0*math.pi*math.pi)) ;// to hartree already
   double b_p  = 0.0;
   double b_f  = 0.0;
   double b_a  = 0.0;
   double c_p  = 0.0;
   double c_f  = 0.0;
   double c_a  = 0.0;
   double x0_p = 0.0;
   double x0_f = 0.0;
   double x0_a = 0.0;
   double eps_p = 0.0;
   double eps_f = 0.0;
   double alpha = 0.0;
   double mu_p = 0.0;
   double mu_f = 0.0;
   double over3 = 1.0/3.0;
   double df2_deta2_0 = 4.0*over3*over3*(-1.0+std::pow((2.0),(1.0/3.0)));
   double beta = 0.0;
   double delta_eps_etha = 0.0;
   double delta_eps_1    = 0.0;
   double S1    = 0.0;
   double S2    = 0.0;
   double S3    = 0.0;
   double S4    = 0.0;
   double M1    = 0.0;
   double M2    = 0.0;
   double M3_A    = 0.0;
   double M3_B    = 0.0;
   double rs_db_drs    = 0.0;
   double rs_da_drs    = 0.0;
   double spindensity_4 = std::pow(spindensity,4.0);
   double spindensity_3 = std::pow(spindensity,3.0);
//   VWN5
   if (this->cor_type == 1){
     b_f  =  7.06042;  // Caption Table 5
     c_f  = 18.0578;   // Caption Table 5
     x0_f = -0.32500;  // Caption Table 5
     b_p  =  3.72744;  // Caption Table 5
     c_p  = 12.9352;   // Caption Table 5
     x0_p = -0.10498;   // Caption Table 5
     b_a  =  1.13107;   // intext page
     c_a  = 13.0045;    // intext page
     x0_a = -0.00475840; // intext page
   }else if(this->cor_type == 2){
//  VWN3
     b_p  =  13.0720;  // into text page 1207
     c_p  =  42.7198;  // into text page 1207
     x0_p =  -0.409286; // into text page 1207
   }
// Closed Shell
   if(this->isClosedShell && this->Ref_ != TCS) {
     this->eps_corr = 0.0;
     this->mu_corr  = 0.0;
     this->eps_corr =  EvepsVWN(0,A_p,b_p,c_p,x0_p,rho);
     this->mu_corr  = -over3*EvepsVWN(1,A_p,b_p,c_p,x0_p,rho);
     this->mu_corr += this->eps_corr ;
   }else{
//   Energy density Eq.2.4/2.2 of ref
/*
     this->eps_corr  = 0.0;
     this->mu_corr   = 0.0;
     this->mu_corr_B = 0.0;
     eps_p = EvepsVWN(0,A_p,b_p,c_p,x0_p,rho);
     eps_f = EvepsVWN(0,A_f,b_f,c_f,x0_f,rho);
     alpha = EvepsVWN(0,A_a,b_a,c_a,x0_a,rho);
     delta_eps_1 = eps_f - eps_p;
     beta  = df2_deta2_0 * delta_eps_1 / alpha;
     beta  += -1.0;
     delta_eps_etha = alpha;
     delta_eps_etha *= f_spindens(0,spindensity)/df2_deta2_0;
     delta_eps_etha *= (1.0 + beta*spindensity_4);
     this->eps_corr  = eps_p + delta_eps_etha ;
//   build the potential
     mu_p = EvepsVWN(0,A_p,b_p,c_p,x0_p,rho) - over3*EvepsVWN(1,A_p,b_p,c_p,x0_p,rho);
     rs_da_drs = EvepsVWN(1,A_a,b_a,c_a,x0_a,rho);
     rs_db_drs = EvepsVWN(1,A_f,b_f,c_f,x0_f,rho) - EvepsVWN(1,A_p,b_p,c_p,x0_p,rho);
     rs_db_drs *= df2_deta2_0/alpha;
     rs_db_drs += -df2_deta2_0*delta_eps_1*rs_da_drs/(alpha*alpha);// complete 
     S1 = -over3*rs_da_drs*(1.0 + beta*spindensity_4);
     S2 = -alpha*over3*rs_db_drs*spindensity_4;
     S3 = 4.0*beta*f_spindens(0,spindensity)*spindensity_3;
     S4 = (1.0 + beta*spindensity_4)*df_spindens(spindensity);
     M1 = f_spindens(0,spindensity)/df2_deta2_0;
     M2 = alpha/df2_deta2_0;
     M3_A = 1.0 - spindensity; 
     M3_B = 1.0 + spindensity;
    
     this->mu_corr  = mu_p + delta_eps_etha; 
     this->mu_corr += M1*(S1+S2);
     this->mu_corr_B = this->mu_corr;
     this->mu_corr   += M2*(S3+S4)*M3_A; 
     this->mu_corr_B += -M2*(S3+S4)*M3_B; 
*/
     }  //Open Shell
};

} // Namespace ChronusQ
