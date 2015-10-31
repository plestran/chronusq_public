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
    int nAtom    = this->molecule_->nAtoms(); // Number of Atoms
    int nRad     = 100;       // Number of Radial grid points for each center
    int nAng     = 302;       // Number of Angular grid points for each center 
                              //  (only certain values are allowed - see grid.h)

    int npts     = nRad*nAng; // Total Number of grid point for each center

    double weight= 0.0;                            

    //TF LDA Prefactor (for Vx)
    double CxVx = -(std::pow((3.0/math.pi),(1.0/3.0)));    
    double CxEn =  (3.0/4.0);      //TF LDA Prefactor to finish the X-Energy
    double val = 4.0*math.pi*CxVx;
    this->totalEx = 0.0;                                  // Total Exchange Energy
    this->totalEcorr = 0.0;                               // Total Correlation Energy
/*  
 *  Generate grids 
 *
 *    Raw grid, it has to be centered and integrated over each center and 
 *    centered over each atom
 */
    GaussChebyshev1stGridInf Rad(nRad,0.0,1.0);   // Radial Grid
    LebedevGrid GridLeb(nAng);                    // Angular Grid
    GridLeb.genGrid();                            // Generate Angular Grid
    this->vXA()->setZero();   // Set to zero every occurence of the SCF
    this->vCorA()->setZero(); // Set to zero every occurence of the SCF
    if(!this->isClosedShell && this->Ref_ != TCS) this->vXB()->setZero();
    if(!this->isClosedShell && this->Ref_ != TCS) this->vCorB()->setZero();
    // Loop over atomic centers
    for(int iAtm = 0; iAtm < nAtom; iAtm++){

      // The Radial grid is generated and scaled for each atom
      Rad.genGrid();
      // Scale the grid according the Atomic Bragg-Slater Radius 
      Rad.scalePts((elements[this->molecule_->index(iAtm)].sradius)) ;  
      // Final Raw (not centered) 3D grid (Radial times Angular grid)
      TwoDGrid Raw3Dg(npts,&Rad,&GridLeb);                              

      //Center the Grid at iAtom
      Raw3Dg.centerGrid(
        (*this->molecule_->cart())(0,iAtm),
        (*this->molecule_->cart())(1,iAtm),
        (*this->molecule_->cart())(2,iAtm)
      );

      // Loop over grid points
      for(int ipts = 0; ipts < npts; ipts++){

        // Evaluate each Becke fuzzy call weight, normalize it and muliply by 
        //   the Raw grid weight at that point
        weight = Raw3Dg.getweightsGrid(ipts)  
                 * (this->formBeckeW((Raw3Dg.gridPtCart(ipts)),iAtm))
                 / (this->normBeckeW(Raw3Dg.gridPtCart(ipts))) ;
        // Build the Vxc for the ipts grid point 
        //  ** Vxc will be ready at the end of the two loop, to be finalized ** 
        this->buildVxc((Raw3Dg.gridPtCart(ipts)),weight);

      } // loop ipts
    } // loop natoms

    //  Finishing the Vxc using the TF factor and the integration 
    //    prefactor over a solid sphere
    (*this->vXA())    =  val * (*this->vXA());
    this->totalEx     =  val * CxEn * (this->totalEx);
    (*this->vCorA())  =  4.0 * math.pi * (*this->vCorA());
    if(!this->isClosedShell && this->Ref_ != TCS) (*this->vCorB())  =  4.0 * math.pi * (*this->vCorB());
    this->totalEcorr  =  4.0 * math.pi * (this->totalEcorr);
    // For open shell averything has to be scaled by 2^(1/3)
    if(!this->isClosedShell && this->Ref_ != TCS){
      (*this->vXA()) *= std::pow(2.0,(1.0/3.0));  
      (*this->vXB()) *= std::pow(2.0,(1.0/3.0)) * val;
    }

    if(this->printLevel_ >= 3) {
      prettyPrint(this->fileio_->out,(*this->vXA()),"LDA Vx alpha");
      prettyPrint(this->fileio_->out,(*this->vCorA()),"Vc Vc alpha");
      if(!this->isClosedShell && this->Ref_ != TCS) 
        prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vx beta");

      this->fileio_->out << "Total LDA Ex ="    << this->totalEx 
                         << " Total VWN Corr= " << this->totalEcorr << endl;
    }
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
   double over3 = 1.0/3.0;
   double delta_eps_1    = 0.0;
   double S1    = 0.0;
   double S2    = 0.0;
   double M3_A    = 0.0;
   double M3_B    = 0.0;
   double rs      = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/3.0));

   double alpha = 0.0;
   double mu_p = 0.0;
   double mu_f = 0.0;
   double beta = 0.0;
   double S3    = 0.0;
   double S4    = 0.0;
   double M1    = 0.0;
   double M2    = 0.0;
   double db_dr = 0.0; 
   double delta_eps_etha = 0.0;
   double spindensity_4 = std::pow(spindensity,4.0);
   double spindensity_3 = std::pow(spindensity,3.0);


//   VWN5
   if (this->CorrKernel_ == VWN5){
     b_f  =  7.06042;  // Caption Table 5
     c_f  = 18.0578;   // Caption Table 5
     x0_f = -0.32500;  // Caption Table 5
     b_p  =  3.72744;  // Caption Table 5
     c_p  = 12.9352;   // Caption Table 5
     x0_p = -0.10498;   // Caption Table 5
     b_a  =  1.13107;   // intext page 1209
     c_a  = 13.0045;    // intext page 1209
     x0_a = -0.00475840; // intext page 1209
   }else if(this->CorrKernel_ == VWN3){
//  VWN3
     b_p  =  13.0720;   // into text page 1207
     c_p  =  42.7198;   // into text page 1207
     x0_p =  -0.409286; // into text page 1207
     b_f  =  20.1231;   // into text pagr 1207
     c_f  =  101.578;   // into text pagr 1207
     x0_f = -0.743294;  // into text pagr 1207
     b_a  =  1.13107;   // intext page 1209
     c_a  = 13.0045;    // intext page 1209
     x0_a = -0.00475840; // intext page 1209
   }
/*  // Debug
    cout << "**********" <<endl;
    double rho1;
    rho1 = 0.238732414637843;   //rs=1
    cout << "EpsP " <<   EvepsVWN(0,A_p,b_p,c_p,x0_p,rho1) << endl; 
    cout << "EpsF " <<   EvepsVWN(0,A_f,b_f,c_f,x0_f,rho1) << endl; 
    cout << "dEpsP " <<  EvepsVWN(2,A_p,b_p,c_p,x0_p,rho1) << endl; 
    cout << "dEpsF " <<  EvepsVWN(2,A_f,b_f,c_f,x0_f,rho1) << endl; 
    cout << "**********" <<endl;
*/
// Closed Shell
   if(this->isClosedShell && this->Ref_ != TCS) {
     this->eps_corr = 0.0;
     this->mu_corr  = 0.0;
     this->eps_corr =  EvepsVWN(0,A_p,b_p,c_p,x0_p,rho);
//     this->mu_corr  = -over3*EvepsVWN(1,A_p,b_p,c_p,x0_p,rho);
     this->mu_corr  = -over3*rs*EvepsVWN(2,A_p,b_p,c_p,x0_p,rho);
     this->mu_corr += this->eps_corr ;
   }else{
//   Open Shell Case
     if(this->CorrKernel_ == VWN3){
//     Used Linear Interpolation between parg and ferr 
//     Eq 2.4 and its analytic derivative for VWN3
       this->eps_corr  = 0.0;
       this->mu_corr   = 0.0;
       this->mu_corr_B = 0.0;
       eps_p = EvepsVWN(0,A_p,b_p,c_p,x0_p,rho);
       eps_f = EvepsVWN(0,A_f,b_f,c_f,x0_f,rho);
       delta_eps_1 = eps_f - eps_p;
       this->eps_corr  = eps_p + delta_eps_1*f_spindens(0,spindensity);
       S1 =  -rs*over3*EvepsVWN(2,A_p,b_p,c_p,x0_p,rho);
       S2 =  -rs*over3*f_spindens(0,spindensity)*(EvepsVWN(2,A_f,b_f,c_f,x0_f,rho) - EvepsVWN(2,A_p,b_p,c_p,x0_p,rho));
//     S1 =  -over3*EvepsVWN(1,A_p,b_p,c_p,x0_p,rho);
//     S2 =  -over3*f_spindens(0,spindensity)*(EvepsVWN(1,A_f,b_f,c_f,x0_f,rho) - EvepsVWN(1,A_p,b_p,c_p,x0_p,rho));
       M3_A =   1.0 - spindensity; 
       M3_B = -(1.0 + spindensity);
       this->mu_corr   = S1 + S2 + this->eps_corr;
       this->mu_corr_B = S1 + S2 + this->eps_corr;     
       this->mu_corr   +=  delta_eps_1*M3_A*df_spindens(spindensity);
       this->mu_corr_B +=  delta_eps_1*M3_B*df_spindens(spindensity);
     }else if(this->CorrKernel_ == VWN5){
//     Used improved Interpolation between parg and ferr 
//     Eq 3.2  and 3.3 and its analytic derivative for VWN5

     alpha = EvepsVWN(0,A_a,b_a,c_a,x0_a,rho);
     eps_p = EvepsVWN(0,A_p,b_p,c_p,x0_p,rho);
     eps_f = EvepsVWN(0,A_f,b_f,c_f,x0_f,rho);
     delta_eps_1 = eps_f - eps_p;
     beta  = this->df2_spindens(0.0) * delta_eps_1 / alpha;
     beta  += -1.0;
     delta_eps_etha = alpha;
     delta_eps_etha *= (f_spindens(0,spindensity)/this->df2_spindens(0.0));
     delta_eps_etha *= (1.0 + beta*spindensity_4);
     this->eps_corr  = eps_p + delta_eps_etha ;
//   build the potential

//   dbeta/dr
     db_dr = -delta_eps_1 * EvepsVWN(2,A_a,b_a,c_a,x0_a,rho);
     db_dr += (EvepsVWN(2,A_f,b_f,c_f,x0_f,rho) - EvepsVWN(2,A_p,b_p,c_p,x0_p,rho)) * alpha;
     db_dr *= this->df2_spindens(0.0);
     db_dr /= alpha*alpha;
//   S1 = da/dr * (f(zeta))*(1+zeta^4*beta)/ df2_spindens(0.0)
     S1 = this->f_spindens(0,spindensity);
     S1 *= (1.0 + beta*spindensity_4);
     S1 *= EvepsVWN(2,A_a,b_a,c_a,x0_a,rho);
     S1 /= this->df2_spindens(0.0);
//   S2 = d eps_p/ dr
     S2  = EvepsVWN(2,A_p,b_p,c_p,x0_p,rho);
//   S3 = df(zeta)/dr * alpha* (1+zeta^4*beta)/df2_spindens(0.0)
     S3  = alpha;
     S3 *= (1.0 + beta*spindensity_4);
     S3 *= this->df_spindens(spindensity); 
     S3 /= this->df2_spindens(0.0);
//   M1 alpha * f(zeta)/this->df2_spindens(0.0)
     M1  = alpha;
     M1 *= this->f_spindens(0,spindensity); 
     M1 /= this->df2_spindens(0.0);
//   M2  zeta^4 * dbeta/dr
     M2  = spindensity_4 * db_dr;
//   dzeta/drho_x 
     M3_A =   1.0 - spindensity; 
     M3_B = -(1.0 + spindensity);
     M3_A *= spindensity_3 * beta * 4.0;
     M3_B *= spindensity_3 * beta * 4.0;
     M3_A += M2;
     M3_B += M2;
     M3_A *= M1;
     M3_B *= M1;
     M3_A +=  S3*(1.0 - spindensity);   
     M3_B += -S3*(1.0 + spindensity);   
     this->mu_corr   = -rs*over3*(S1 + S2);
     this->mu_corr_B = -rs*over3*(S1 + S2);
     
     this->mu_corr     += M3_A;
     this->mu_corr_B   += M3_B;


     this->mu_corr     += this->eps_corr;
     this->mu_corr_B   += this->eps_corr;


     }
  }  //Open Shell
};

} // Namespace ChronusQ
