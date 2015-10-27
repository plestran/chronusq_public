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
    this->totalEx = 0.0;           // Total Exchange Energy
    this->totalEcorr = 0.0;        // Total Correlation Energy

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
    this->totalEcorr  =  4.0 * math.pi * (this->totalEcorr);

    // For open shell averything has to be scaled by 2^(1/3)
    if(!this->isClosedShell && this->Ref_ != TCS){
      (*this->vXA()) *= std::pow(2.0,(1.0/3.0));  
      (*this->vXB()) *= std::pow(2.0,(1.0/3.0)) * val;
      this->totalEx  *= std::pow(2.0,(1.0/3.0));
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
   double A  = A1/2.0; // to Hartree
   double b;
   double c;
   double x0;
//   VWN5
   if (this->CorrKernel_ == VWN5){
     b  = 3.72744;  // Caption Table 5
     c  = 12.9352;  // Caption Table 5
     x0 = -0.10498; // Caption Table 5
   }else if(this->CorrKernel_ == VWN3){
//  VWN3
     b  = 13.0720;  // into text page 1207
     c  = 42.7198;  // into text page 1207
     x0 = -0.409286; // into text page 1207
   }
//  Derivatives therms
   double b1 = (b*x0 - c)/(c*x0); 
   double b2 = (x0 - b)/(c*x0); 
   double b3 = (-1.0)/(c*x0); 
   double Q =std::pow((4.0*c - b*b),(1.0/2.0)); 
   double r_s = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/3.0));
   double r_s_sqrt = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/6.0));
   double r_s_32 = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/2.0));
   double X = r_s + b*(r_s_sqrt) + c; 
   double X_x0 = x0*x0 + b*x0 + c; 
   this->eps_corr = 0.0;
   this->mu_corr  = 0.0;
   this->eps_corr = A *
          ( std::log(r_s/X) + 
            2.0*b*(std::atan(Q/(2.0*r_s_sqrt + b)))/Q  -
            b*x0*(1.0/X_x0) * ( 
                 (std::log( (std::pow((r_s_sqrt-x0),2.0))/X )) +
                 (2.0*(b + 2.0*x0)*(1.0/Q)*(std::atan(Q/(2.0*r_s_sqrt + b))) ) 
                 ) );
   this->mu_corr = -A * ( (1.0 + b1*r_s_sqrt)/(1.0 + b1*r_s_sqrt + b2*r_s + b3*r_s_32)) / 3.0 ;
   this->mu_corr += this->eps_corr ;
//   cout << "r_s " << r_s << " r_s_sqrt " << r_s_sqrt << " r_s_32 " << r_s_32 <<endl;
//   cout << " eps "<< this->eps_corr*1000.0 << " mu " << this->mu_corr <<endl;
};


template<>
double SingleSlater<double>::spindens(double rho_A, double rho_B){
      double spindens;
      double spinrho;
      spindens = (rho_A + rho_B)/ (rho_A + rho_B);
      spinrho = -2.0;
      spinrho += std::pow((1.0+spindens),(4.0/3.0)); 
      spinrho += std::pow((1.0-spindens),(4.0/3.0)); 
      spinrho = spinrho/(2.0) ;
      spinrho = spinrho/(-1.0+std::pow((2.0),(1.0/3.0))); 
      return spinrho;
};  //end

} // Namespace ChronusQ
