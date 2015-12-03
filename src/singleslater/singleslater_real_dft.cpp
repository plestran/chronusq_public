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

template<>
double SingleSlater<double>::formBeckeW(cartGP gridPt, int iAtm){
//     Generate Frisch (not-normalized yet) Weights (if frischW) according to the partition schems in
//     (Chem. Phys. Let., 257, 213-223 (1996)) using Eq. 11 and 14
//     Note these Weights have to be normailzed 
//     Generate Becke not-normalized Weights (if becke) according to the partition schems in
//     (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii Fuzzi Cells
//     Note these Weights have to be normailzed (see normBeckeW) 
       int nAtom = this->molecule_->nAtoms();
       double WW = 1.0;
       double tmp;
       double muij;   ///< elliptical coordinate (ri -rj / Rij)
       cartGP rj;     ///< Cartisian position of Atom j
       cartGP ri;     ///< Cartisian position of Atom i
       ri.set<0>((*this->molecule_->cart())(0,iAtm) );
       ri.set<1>((*this->molecule_->cart())(1,iAtm) );
       ri.set<2>((*this->molecule_->cart())(2,iAtm) );
       for(auto jAtm = 0; jAtm < nAtom; jAtm++){
         if(jAtm != iAtm){
           muij = 0.0;
//       Vector rj (Atoms (j.ne.i) position)
           rj.set<0>((*this->molecule_->cart())(0,jAtm));
           rj.set<1>((*this->molecule_->cart())(1,jAtm));
           rj.set<2>((*this->molecule_->cart())(2,jAtm));
//       Coordinate of the Grid point in elleptical (Eq. 11) 
           muij = (boost::geometry::distance(gridPt,ri) - 
                   boost::geometry::distance(gridPt,rj))/
                   (*this->molecule_->rIJ())(iAtm,jAtm) ;
//       Do the product over all atoms i .ne. j (Eq. 13 using definition Eq. 21 with k=3)
           if (this->weightScheme_ == FRISCH) 
             WW *= 0.5*(1.0-this->twodgrid_->frischpol(muij,0.64));
           else if (this->weightScheme_ == BECKE)  
             WW *= 0.5 * 
                   (1.0-this->twodgrid_->voronoii(
                          this->twodgrid_->voronoii(
                            this->twodgrid_->voronoii(muij))));
           }
         }
       return WW;
}; //End formBeckeW


template<>
double SingleSlater<double>::normBeckeW(cartGP gridPt){
//     Normalization of Becke/Frisch Weights
//     (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii Fuzzi Cells
//     Eq. 22
       int   nAtom = this->molecule_->nAtoms();
       double norm = 0.0;
       for(auto iAtm = 0; iAtm < nAtom; iAtm++){
         norm += this->formBeckeW(gridPt,iAtm);
         }
       return norm ;
}; //End normBeckeW

template<>
double SingleSlater<double>::f_spindens(int iop, double spindensity){
      double f_spindensity;
      double thrs = 1.11e-16;
      double fact = (-2.0+std::pow((2.0),(4.0/3.0)));
      if (iop == 0) {
      f_spindensity = 0.0;
      if ((1.0+spindensity) >= thrs)   f_spindensity += std::pow((1.0+spindensity),(4.0/3.0)); 
      if ((1.0-spindensity) >= thrs)   f_spindensity += std::pow((1.0-spindensity),(4.0/3.0)); 
      f_spindensity += -2.0;
      f_spindensity /= fact; 
      }else if(iop ==1){
      f_spindensity  = std::pow((1.0+spindensity),(4.0/3.0)); 
      f_spindensity += std::pow((1.0-spindensity),(4.0/3.0)); 
      f_spindensity /= (2.0) ;
      }
      return f_spindensity;
};  //end


template<>
double SingleSlater<double>::df_spindens(double spindensity){
      double df_spindensity;
      double thrs = 1.11e-16;
      double fact = (-2.0+std::pow((2.0),(4.0/3.0)));
      df_spindensity = 0.0;
      if ((1.0+spindensity) >= thrs) df_spindensity +=  std::pow((1.0+spindensity),(1.0/3.0)); 
      if ((1.0-spindensity) >= thrs) df_spindensity -= std::pow((1.0-spindensity),(1.0/3.0)); 
      df_spindensity *= (4.0/3.0);
      df_spindensity /= fact; 
      return df_spindensity;
};  //end


template<>
double SingleSlater<double>::df2_spindens(double spindensity){
      double df2_spindensity;
      double thrs = 1.11e-16;
      double fact = (-1.0+std::pow((2.0),(1.0/3.0)));
      df2_spindensity = 0.0;
      if ((1.0+spindensity) >= thrs) df2_spindensity +=  std::pow((1.0+spindensity),(-2.0/3.0)); 
      if ((1.0-spindensity) >= thrs) df2_spindensity +=  std::pow((1.0-spindensity),(-2.0/3.0)); 
      df2_spindensity *= (2.0/9.0);
      df2_spindensity /= fact; 
      return df2_spindensity;
};  //end


template<> 
double SingleSlater<double>::spindens(double rho_A, double rho_B) {
return (rho_A - rho_B)/ (rho_A + rho_B);
};  // 



template<>
void SingleSlater<double>::genSparseBasisMap(){
  this->basisset_->radcut(this->epsScreen, this->maxiter, this->epsConv);
/*
  RealMatrix overlapR(this->nBasis_,this->nBasis_);        ///< Overlap at grid point
  overlapR.setZero();
*/
/*
  this->nRadDFTGridPts_ = 50;
  this->nAngDFTGridPts_ = 434;
  this->screenVxc = false;
*/
//  this->epsScreen = 1.0e-16;
  this->ngpts = this->nRadDFTGridPts_*this->nAngDFTGridPts_;  // Number of grid point for each center
  OneDGrid * Rad ;                              // Pointer for Radial Grid
  LebedevGrid GridLeb(this->nAngDFTGridPts_);   // Angular Grid
  int nDer = 0 ; //Order of differenzation
  if (this->dftGrid_ == GAUSSCHEB)  
    Rad = new GaussChebyshev1stGridInf(this->nRadDFTGridPts_,0.0,1.0);   
  else if (this->dftGrid_ == EULERMACL) 
    Rad = new  EulerMaclaurinGrid(this->nRadDFTGridPts_,0.0,1.0);   

//Generare Angular Grid
  GridLeb.genGrid();                            
  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    this->sparseMap_.push_back(RealSparseMatrix(this->nTCS_*this->nBasis_,this->ngpts));
    this->sparseWeights_.push_back(RealSparseMatrix(this->ngpts,1));
    this->sparseDoRho_.push_back(RealSparseMatrix(this->ngpts,1));
//  Derivative
    this->sparsedmudX_.push_back(RealSparseMatrix(this->nTCS_*this->nBasis_,this->ngpts));
    this->sparsedmudY_.push_back(RealSparseMatrix(this->nTCS_*this->nBasis_,this->ngpts));
    this->sparsedmudZ_.push_back(RealSparseMatrix(this->nTCS_*this->nBasis_,this->ngpts));
    RealSparseMatrix *Map        = &this->sparseMap_[iAtm];
    RealSparseMatrix *WeightsMap = &this->sparseWeights_[iAtm];
    RealSparseMatrix *DoRhoMap   = &this->sparseDoRho_[iAtm];
    RealSparseMatrix *MapdX      = &this->sparsedmudX_[iAtm];
    RealSparseMatrix *MapdY      = &this->sparsedmudY_[iAtm];
    RealSparseMatrix *MapdZ      = &this->sparsedmudZ_[iAtm];
    double val;
    // Generate grids
    Rad->genGrid(); 
    Rad->atomGrid((elements[this->molecule_->index(iAtm)].sradius)) ;  
    TwoDGrid Raw3Dg(this->ngpts,Rad,&GridLeb);             
    //Center the Grid at iAtom
    Raw3Dg.centerGrid(
      (*this->molecule_->cart())(0,iAtm),
      (*this->molecule_->cart())(1,iAtm),
      (*this->molecule_->cart())(2,iAtm)
    );
    for (auto ipts =0; ipts < this->ngpts; ipts++){
      
      cartGP pt = Raw3Dg.gridPtCart(ipts); 
      auto mapRad_ = this->basisset_->MapGridBasis(pt);
      // Evaluate each Becke fuzzy call weight, normalize it and muliply by 
      //   the Raw grid weight at that point
      auto bweight = (this->formBeckeW(pt,iAtm)) 
                     / (this->normBeckeW(pt)) ;
      if (this->screenVxc) {
        if (mapRad_[0] || (bweight < this->epsScreen)) continue;
        }
         WeightsMap->insert(ipts,0) = Raw3Dg.getweightsGrid(ipts) * bweight;
         //skyp all point
         for (int s1=0; s1 < this->basisset_->nShell(); s1++){
           if (this->screenVxc) {
             if (!mapRad_[s1+1]) continue;
             }
           int bf1_s = this->basisset_->mapSh2Bf(s1);
           auto shSize = this->basisset_->shells(s1).size(); 
           libint2::Shell shTmp = this->basisset_->shells(s1);
           double * s1Eval = this->basisset_->basisDEval(nDer,shTmp,&pt);
           double * ds1EvalX = s1Eval + shSize;
           double * ds1EvalY = ds1EvalX + shSize;
           double * ds1EvalZ = ds1EvalY + shSize;
           for (auto mu =0; mu < shSize; mu++){ 
             val = s1Eval[mu];
             if (std::abs(val) > this->epsScreen){
               Map->insert(bf1_s+mu,ipts) = val;
               }
             if (nDer ==1) {
               val = ds1EvalX[mu];
               if (std::abs(val) > this->epsScreen){
                 MapdX->insert(bf1_s+mu,ipts) = val;
                 }
               val = ds1EvalY[mu];
               if (std::abs(val) > this->epsScreen){
                 MapdY->insert(bf1_s+mu,ipts) = val;
                 }
               val = ds1EvalZ[mu];
               if (std::abs(val) > this->epsScreen){
                 MapdZ->insert(bf1_s+mu,ipts) = val;
                 }
               }
             } //loop over basis (within a given ishell)
         } //loop over shells
//AP
//            overlapR += (*WeightsMap).coeff(ipts,0)*(*Map).col(ipts)*(*MapdX).col(ipts).transpose();
//AP
         if (this->screenVxc && Map->col(ipts).norm() > this->epsScreen){
         DoRhoMap->insert(ipts,0) = 2;}
    } //loop over pts
//     cout << "non Zero " << this->sparseMap_[iAtm].nonZeros() << " " << this->ngpts <<endl; 
  } // loop over atoms
/*
         overlapR *= 4.0*math.pi;
         prettyPrint(cout,overlapR,"dipolex ");
*/
};// End genSparseBasisMap

template<>
double SingleSlater<double>::EvepsVWN(int iop, double A_x, double b_x, double c_x, double x0_x, double rho){
//    From Reference Vosko en Al., Can. J. Phys., 58, 1200 (1980). VWN3 and VWN5 interpolation formula   
//    IOP 0 -> Eq 4.4 
//    IOP 1 Eq. 4.3 (finishing the derivate of eps , rs factor already included)
//    IOP 2 Analitic Derv of Eq 4.4 (note this one has to be moltiplied outside by rs to get the final needed term)
    double val      = 0.0;
    double b1       = (b_x*x0_x - c_x)/(c_x*x0_x); 
    double b2       = (x0_x - b_x)/(c_x*x0_x); 
    double b3       = (-1.0)/(c_x*x0_x); 
    double Q        = std::pow((4.0*c_x - b_x*b_x),(1.0/2.0)); 
    double r_s      = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/3.0));
    double r_s_sqrt = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/6.0));
    double r_s_32   = std::pow(((3.0)/(4.0*math.pi*rho)),(1.0/2.0));
    double X        = r_s + b_x*(r_s_sqrt) + c_x; 
    double X_x0     = x0_x*x0_x + b_x*x0_x + c_x; 
    double x        = r_s_sqrt;
    double S1 = 0.0;
    double S2 = 0.0;
    double S3 = 0.0;
    if (iop == 0){
     val = A_x *
         ( std::log(r_s/X) + 
           2.0*b_x*(std::atan(Q/(2.0*r_s_sqrt + b_x)))/Q  -
           b_x*x0_x*(1.0/X_x0) * ( 
                                 (std::log( (std::pow((r_s_sqrt-x0_x),2.0))/X )) +
                                 (2.0*(b_x + 2.0*x0_x)*(1.0/Q)*(std::atan(Q/(2.0*r_s_sqrt + b_x))) ) 
                                 ) 
         );
     }else if (iop == 1){
//           A_x * ( (1.0 + b1*r_s_sqrt)/(1.0 + b1*r_s_sqrt + b2*r_s + b3*r_s_32)) / 3.0 ;
            val = A_x* ( (1.0 + b1*r_s_sqrt)/(1.0 + b1*r_s_sqrt + b2*r_s + b3*r_s_32));

     }else if (iop == 2){
            S1 = r_s_sqrt - x0_x;  //dxx0
            S2 = b_x * r_s_sqrt * x0_x; // bxx0
            S3 = X ; //c +bx+r
            val  = A_x*(c_x*S1 - S2);
            val /= (r_s*S3*S1);
     }
     return val;
};  //end

/*
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
   double S5    = 0.0;
   double M1    = 0.0;
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
//   S4  zeta^4 * dbeta/dr
     S4  = spindensity_4 * db_dr;
//   S5  4 * beta * zeta^3
     S5  = 4.0 * beta * spindensity_3;
//   dzeta/drho_x 
     M3_A =   1.0 - spindensity; 
     M3_B = -(1.0 + spindensity);

     this->mu_corr   = -rs*over3*(S1 + S2 + M1*S4);
     this->mu_corr_B = this->mu_corr;
     
     this->mu_corr     += M3_A*(M1*S5 + S3);
     this->mu_corr_B   += M3_B*(M1*S5 + S3);


     this->mu_corr     += this->eps_corr;
     this->mu_corr_B   += this->eps_corr;


     }
  }  //Open Shell
}; //End formCor
*/

template<>
std::array<double,3> SingleSlater<double>::formVExSlater (double rho, double spindensity){
//  epsmu ={energydens_c, potential_c_alpha, potential_c_beta}
  std::array<double,3> epsmu = {0.0,0.0,0.0};
  double d1over3    = 1.0/3.0;
  double d4over3  = 4.0/3.0;
  double eps0     = std::pow(rho,d1over3);      
  if(!this->isClosedShell && this->Ref_ != TCS) {
    epsmu[0]  =  eps0 * this->f_spindens(1,spindensity); 
    epsmu[1]  =  d4over3*eps0*std::pow((1.0+spindensity),d1over3);
    epsmu[2]  =  d4over3*eps0*std::pow((1.0-spindensity),d1over3);
/*
    double twoat1o3 = std::pow(2.0,(1.0/3.0)); 
    epsmu[0]  = eps0 + (d1over3*eps0 - eps0)*this->f_spindens(0,spindensity);
    epsmu[1]  = d1over3*epsmu[0];
    epsmu[1] += (twoat1o3*eps0 - eps0)*df_spindens(spindensity); 
    epsmu[2]  = epsmu[1] -( 1.0 + spindensity);
    epsmu[1] += 1.0 - spindensity;    
*/
  }else {
    epsmu[0]  = eps0 ;
    epsmu[1]  = d4over3*eps0;
  }
  return epsmu;
}; //End formVexSlater

template<>
double SingleSlater<double>::gB88 (int nDer, double x){
  double beta = 0.0042;
  double Cx   = 0.930525736349100;  // (3/2)*((3/(4*pi))^(1/3))
  double gx;
//  Eq A4 Pople, J. Chem. Phys. 5612, (1992) nder =0
  if (nDer == 0 ){
    gx  = -beta*x*x;
    gx /= (1.0 + 6.0 * x * boost::math::asinh(x)) ;
    gx += -Cx;
  }else if (nDer == 1){
//  Eq A8 Pople, J. Chem. Phys. 5612, (1992) nder =1
    gx  = x / (std::sqrt(x*x+1.0));
    gx += -boost::math::asinh(x) ;
    gx *= 6.0 * beta * beta * x * x;
    gx += -2.0*beta*x;
    gx /= (1.0 + 6.0 * boost::math::asinh(x))*(1.0 + 6.0 * x * boost::math::asinh(x));
  }
 return gx; 
};  //End Form g function for B88 Exchange

template<>
std::array<double,3> SingleSlater<double>::formVExB88 (double rhoA, double rhoB, 
double gammaAA, double gammaBB){
//  Becke Exchange : Eq 8 From Becke, Phys. Rev A., 3098 (1988)
//  Becke Exchange : Implemented as Eq A3 From Pople, J. Chem. Phys. 5612, (1992)
  std::array<double,3> epsmu = {0.0,0.0,0.0};
  double rhoA4ov3 = std::pow(rhoA,(4.0/3.0));
  double rhoB4ov3 = std::pow(rhoB,(4.0/3.0));
  double xA   = gammaAA / rhoA4ov3; 
  double xB   = gammaBB / rhoA4ov3; 
  epsmu[0]    = rhoA4ov3*this->gB88(0,xA);
  epsmu[0]   += rhoB4ov3*this->gB88(0,xB);

};  //End Form B88 Exchange

template<>
std::array<double,3> SingleSlater<double>::formVCLYP (double rhoA, double rhoB, 
double drhoA, double drhoB){
};  //End Form LYP Correlation

template<>
std::array<double,3> SingleSlater<double>::formVCVWN (double rho, double spindensity){
//  epsmu ={energydens_c, potential_c_alpha, potential_c_beta}
  std::array<double,3> epsmu = {0.0,0.0,0.0};
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
   double S5    = 0.0;
   double M1    = 0.0;
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
     epsmu[0] =  EvepsVWN(0,A_p,b_p,c_p,x0_p,rho);
//     epsmu[1]  = -over3*EvepsVWN(1,A_p,b_p,c_p,x0_p,rho);
     epsmu[1]  = -over3*rs*EvepsVWN(2,A_p,b_p,c_p,x0_p,rho);
     epsmu[1] += epsmu[0] ;
   }else{
//   Open Shell Case
     if(this->CorrKernel_ == VWN3){
//     Used Linear Interpolation between parg and ferr 
//     Eq 2.4 and its analytic derivative for VWN3
       eps_p = EvepsVWN(0,A_p,b_p,c_p,x0_p,rho);
       eps_f = EvepsVWN(0,A_f,b_f,c_f,x0_f,rho);
       delta_eps_1 = eps_f - eps_p;
       epsmu[0]  = eps_p + delta_eps_1*f_spindens(0,spindensity);
       S1 =  -rs*over3*EvepsVWN(2,A_p,b_p,c_p,x0_p,rho);
       S2 =  -rs*over3*f_spindens(0,spindensity)*(EvepsVWN(2,A_f,b_f,c_f,x0_f,rho) - EvepsVWN(2,A_p,b_p,c_p,x0_p,rho));
//     S1 =  -over3*EvepsVWN(1,A_p,b_p,c_p,x0_p,rho);
//     S2 =  -over3*f_spindens(0,spindensity)*(EvepsVWN(1,A_f,b_f,c_f,x0_f,rho) - EvepsVWN(1,A_p,b_p,c_p,x0_p,rho));
       M3_A =   1.0 - spindensity; 
       M3_B = -(1.0 + spindensity);
       epsmu[1]   = S1 + S2 + epsmu[0];
       epsmu[2] = S1 + S2 + epsmu[0];     
       epsmu[1]   +=  delta_eps_1*M3_A*df_spindens(spindensity);
       epsmu[2] +=  delta_eps_1*M3_B*df_spindens(spindensity);
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
     epsmu[0]  = eps_p + delta_eps_etha ;
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
//   S4  zeta^4 * dbeta/dr
     S4  = spindensity_4 * db_dr;
//   S5  4 * beta * zeta^3
     S5  = 4.0 * beta * spindensity_3;
//   dzeta/drho_x 
     M3_A =   1.0 - spindensity; 
     M3_B = -(1.0 + spindensity);

     epsmu[1] = -rs*over3*(S1 + S2 + M1*S4);
     epsmu[2] = epsmu[1];
     
     epsmu[1] += M3_A*(M1*S5 + S3);
     epsmu[2] += M3_B*(M1*S5 + S3);


     epsmu[1]   += epsmu[0];
     epsmu[2]   += epsmu[0];


     }
  }  //Open Shell
  return epsmu;
}; //END VWN

template<>
std::array<double,3> SingleSlater<double>::formVCGGA (double rhoA, double rhoB, 
     double drhoA, double drhoB){

    std::array<double,3> corrEpsMu;
//    if (this->CorrKernel_ == VWN3 || this->CorrKernel_ == VWN5) {
       corrEpsMu = this->formVCLYP(rhoA, rhoB, drhoA, drhoB);
    return corrEpsMu;
}; //END GENERIC FORMCOR

template<>
std::array<double,3> SingleSlater<double>::formVC (double rho, double spindensity){

    std::array<double,3> corrEpsMu;
    if (this->CorrKernel_ == VWN3 || this->CorrKernel_ == VWN5) {
       corrEpsMu = this->formVCVWN(rho, spindensity);}
    return corrEpsMu;
}; //END GENERIC FORMCOR

template<>
std::array<double,3> SingleSlater<double>::formVEx (double rho, double spindensity){

    std::array<double,3> exEpsMu;
//  Put sime Slater definition
       exEpsMu = this->formVExSlater(rho, spindensity);
    return exEpsMu;
}; //END GENERIC FORMVex

template<>
std::array<double,3> SingleSlater<double>::formVExGGA (double rhoA, double rhoB, 
     double drhoA, double drhoB){

    std::array<double,3> corrEpsMu;
//    if (this->CorrKernel_ == VWN3 || this->CorrKernel_ == VWN5) {
       corrEpsMu = this->formVExB88(rhoA, rhoB, drhoA, drhoB);
    return corrEpsMu;
}; //END GENERIC FORMCOR

/*
template<>
void SingleSlater<double>::evalVXC(cartGP gridPt, double weight, std::vector<bool> mapRad_,
       double & energyX, double & energyC, RealMatrix * VXA, RealMatrix * VXB, RealMatrix * VCA, 
       RealMatrix * VCB){
//  Build the Vxc therm at each Grid Points  

////T
//   auto startVxc = std::chrono::high_resolution_clock::now();  
////T

   double rhor = 0.0;
   double rhor_B = 0.0;
   bool   RHF  = this->Ref_ == RHF;
   bool   doTCS  = this->Ref_ == TCS;
//   double shMax;

////T   
//   auto start_2 = std::chrono::high_resolution_clock::now();  // Timing to allocate and set to zero S(r)
////T

   std::unique_ptr<RealMatrix>  PA_;        ///< Overlap at grid point
   std::unique_ptr<RealMatrix>  overlapR_;        ///< Overlap at grid point
   overlapR_ = std::unique_ptr<RealMatrix>(
     new RealMatrix(this->nBasis_,this->nBasis_));
   overlapR_->setZero();
////T
//   auto finish_2 = std::chrono::high_resolution_clock::now();  
//   this->duration_2 += finish_2 - start_2;
////T
// Loops over shells
////T
////   auto start_4 = std::chrono::high_resolution_clock::now();  // Timing to allocate and set to zero S(r)
   for(auto s1=0l, s12=0l; s1 < this->basisset_->nShell(); s1++){
      if (mapRad_[s1+1]){
        int bf1_s = this->basisset_->mapSh2Bf(s1);
        int n1    = this->basisset_->shells(s1).size();
        for(int s2=0; s2 <= s1; s2++, s12++){
          if (mapRad_[s2+1]){
            int bf2_s   = this->basisset_->mapSh2Bf(s2);
            int n2      = this->basisset_->shells(s2).size();
////T          
//            auto start_7 = std::chrono::high_resolution_clock::now();
////T

//            if (this->screenVxc){
///              if(this->isClosedShell || this->Ref_ !=TCS){
//                shMax = (*this->basisset_->shBlkNormAlpha)(s1,s2);
//                } else {
//                shMax = std::max((*this->basisset_->shBlkNormAlpha)(s1,s2),
//                         (*this->basisset_->shBlkNormBeta)(s1,s2));
//              }
//             if(shMax < (this->controls_->thresholdSchawrtz/ngpts) ) continue;
//             if(shMax < (this->epsScreen/ngpts) ) continue;
//            }
//             double * bfun = new double [n2];
//            libint2::Shell s2sh = this->basisset_->shells(s2); 
//            bfun = 
//              this->basisset_->basisEval(s2sh,&gridPt);
//            RealMap fBuff(bfun,n2,1);
//            fBuff += this->densityA()->block(bf1_s,bf2_s,n1,n2) * fBuff;
            auto pointProd = 
              this->basisset_->basisProdEval(
                this->basisset_->shells(s1),
                this->basisset_->shells(s2),
                &gridPt
              );
////T
//            auto finish_7 = std::chrono::high_resolution_clock::now();  
//            this->duration_7 += finish_7 - start_7;
////            auto start_8 = std::chrono::high_resolution_clock::now();
////T
            RealMap fBuff(pointProd,n1,n2);
            overlapR_->block(bf1_s,bf2_s,n1,n2) = fBuff; 
            delete[] pointProd;
////T
////            auto finish_8 = std::chrono::high_resolution_clock::now();  
////            this->duration_8 += finish_8 - start_8;
////T
          }
        }
      }
    }  
////T
////   auto finish_4 = std::chrono::high_resolution_clock::now();  
////   this->duration_4 += finish_4 - start_4;
////T
////T
//   auto start_3 = std::chrono::high_resolution_clock::now();  // Timing S contraction
////T
     (*overlapR_) = overlapR_->selfadjointView<Lower>();;
   if(this->isClosedShell && this->Ref_ != TCS) {
    rhor = overlapR_->frobInner(this->densityA()->conjugate());
////T
//   auto finish_3 = std::chrono::high_resolution_clock::now();  
//   this->duration_3 += finish_3 - start_3;
////T
    if (this->screenVxc ) {
      if(rhor    <= 0.0 ) {
        if((std::abs(rhor)) <= 1.0e10) rhor = 0.0;
        }
      }
//  LDA Slater Exchange
    if (this->ExchKernel_ != NOEXCH) {
      (*VXA)  += weight*(*overlapR_)*(std::pow(rhor,(1.0/3.0)));
      energyX += weight*(std::pow(rhor,(4.0/3.0)));
    }
//  VWN Correlation
    if (this->CorrKernel_ != NOCORR) {
    if (rhor > 1.0e-20) {                       //this if statement prevent numerical instability with zero guesses
//      this->formVWNPara(2.38732414637843e-04);
//      this->formVWNPara(2.98415518297304e-05);
//      this->formVWNPara(2.3873e-07);
      this->formCor(rhor,0.0);
      (*VCA)  += weight*(*overlapR_)*this->mu_corr;
      energyC += weight*rhor*this->eps_corr;
     }
    }
   }
   if(!this->isClosedShell && this->Ref_ != TCS) {
     rhor   = overlapR_->frobInner(this->densityA()->conjugate());
     rhor_B = overlapR_->frobInner(this->densityB()->conjugate());
//   Avoid numerical noise 
     if (this->screenVxc ) {
       if(rhor    <= 0.0 ) {
          if((std::abs(rhor)) <= 1.0e10) rhor = 0.0;
          }
       if(rhor_B    <= 0.0 ) {
          if((std::abs(rhor_B)) <= 1.0e10) rhor_B = 0.0;
          }
      }
    if (this->ExchKernel_ != NOEXCH){
     (*VXA) += weight*(*overlapR_)*(std::pow(rhor,(1.0/3.0)));
     (*VXB) += weight*(*overlapR_)*(std::pow(rhor_B,(1.0/3.0)));
     if((rhor+rhor_B) > 1.0e-20) this->totalEx   += 
                                   weight*(std::pow((rhor+rhor_B),(4.0/3.0)))*
                                     (this->f_spindens(1,this->spindens(rhor,rhor_B)));
      }
    if (this->CorrKernel_ != NOCORR) {
    if (rhor+rhor_B > 1.0e-20) {                       //this if statement prevent numerical instability with zero guesses
     this->formCor((rhor+rhor_B),(this->spindens(rhor,rhor_B)));
      (*VCA)    += weight*(*overlapR_)*this->mu_corr;
      (*VCB)    += weight*(*overlapR_)*this->mu_corr_B;
      this->totalEcorr    += weight*(rhor+rhor_B)*this->eps_corr;
      }
     }
    }

////T
//    auto finish_Vxc = std::chrono::high_resolution_clock::now();  
//    this->duration_1 += finish_Vxc - startVxc;
////T
//  }

}; //END
*/

template<>
// Cleaned version to handle parallelism (no global variable)
void SingleSlater<double>::evalVXC(cartGP gridPt, double weight, std::vector<bool> mapRad_,
       double & energyX, double & energyC, RealMatrix * VXA, RealMatrix * VXB, RealMatrix * VCA, 
       RealMatrix * VCB) {

   double rhor  = 0.0;  // Total density at point
   double rhorA = 0.0;  // alpha density at point
   double rhorB = 0.0;  // beta  density at point
   bool   RHF  = this->Ref_ == RHF;
   bool   doTCS  = this->Ref_ == TCS;
   double shMax;
   RealMatrix overlapR_(this->nBasis_,this->nBasis_);        ///< Overlap at grid point
   overlapR_.setZero();
   std::array<double,3>  epsMuCor = {0.0,0.0,0.0}; ///< {energydens_corr, potential_corr_alpha, potential_corr_B}
   std::array<double,3>  epsMuExc = {0.0,0.0,0.0}; ///< {energydend_exchange, potential_exchenge_alpha, potential_exchange_beta}

// Loops over shells
   for(auto s1=0l, s12=0l; s1 < this->basisset_->nShell(); s1++){
      if (mapRad_[s1+1]){
        int bf1_s = this->basisset_->mapSh2Bf(s1);
        int n1    = this->basisset_->shells(s1).size();
        for(int s2=0; s2 <= s1; s2++, s12++){
          if (mapRad_[s2+1]){
            int bf2_s   = this->basisset_->mapSh2Bf(s2);
            int n2      = this->basisset_->shells(s2).size();
/*
            if (this->screenVxc){
              if(this->isClosedShell || this->Ref_ !=TCS){
                shMax = (*this->basisset_->shBlkNormAlpha)(s1,s2);
                } else {
                shMax = std::max((*this->basisset_->shBlkNormAlpha)(s1,s2),
                         (*this->basisset_->shBlkNormBeta)(s1,s2));
              }
//             if(shMax < (this->controls_->thresholdSchawrtz/ngpts) ) continue;
             if(shMax < (this->epsScreen/ngpts) ) continue;
            }
*/
            auto pointProd = 
              this->basisset_->basisProdEval(
                this->basisset_->shells(s1),
                this->basisset_->shells(s2),
                &gridPt
              );
            RealMap fBuff(pointProd,n1,n2);
            overlapR_.block(bf1_s,bf2_s,n1,n2) = fBuff; 
            delete[] pointProd;
          }
        }
      }
    }  
//   Overlap at r is ready
    overlapR_ = overlapR_.selfadjointView<Lower>();
//    if (ipts > 1500 && ipts < 1520){ 
//    cout << "Old PTS = " << ipts <<endl;
//    prettyPrint(cout,(overlapR_)," S(ri) Old");
//    }
//  Handle the total density at r for RKS or UKS
    if(!this->isClosedShell && this->Ref_ != TCS) {
      rhorA = overlapR_.frobInner(this->densityA()->conjugate());
      rhorB = overlapR_.frobInner(this->densityB()->conjugate());
      rhor = rhorA + rhorB;
    } else {
      rhor = overlapR_.frobInner(this->densityA()->conjugate()) ;
    }
//  Handle numerical instability if screening on
    if (this->screenVxc ) {
      if(rhor    <= 0.0 ) {
        if((std::abs(rhor)) <= 1.0e10) rhor = 0.0;
      }
    }
//this if statement prevent numerical instability with zero guesses
    if (rhor > 1.0e-20) {     
//  Exchange
      if (this->ExchKernel_ != NOEXCH) {
        if(!this->isClosedShell && this->Ref_ != TCS){
            epsMuExc = this->formVEx(rhor,this->spindens(rhorA,rhorB));
            (*VXB)  += weight*overlapR_*epsMuExc[2];
        } else {
            epsMuExc = this->formVEx(rhor,0.0);
        }
        (*VXA)  += weight*overlapR_*epsMuExc[1];
        energyX += weight*rhor*epsMuExc[0];
      }
//  Correlation
      if (this->CorrKernel_ != NOCORR) {
        if(!this->isClosedShell && this->Ref_ != TCS){
          epsMuCor = this->formVC(rhor,(this->spindens(rhorA,rhorB)));
          (*VCB)  += weight*overlapR_*epsMuCor[2];
        } else{
          epsMuCor = this->formVC(rhor,0.0);
        }
        (*VCA)  += weight*overlapR_*epsMuCor[1];
        energyC += weight*rhor*epsMuCor[0];
      }

    }
}; //END

template<>
// Cleaned version to handle parallelism (no global variable)
void SingleSlater<double>::evalVXC_store(int iAtm, int ipts, double & energyX, 
       double & energyC, RealMatrix * VXA, RealMatrix * VXB, RealMatrix * VCA, 
       RealMatrix * VCB, RealMatrix *STmp) {

   RealSparseMatrix *Map        = &this->sparseMap_[iAtm];
   RealSparseMatrix *WeightsMap = &this->sparseWeights_[iAtm];
   RealSparseMatrix *MapdX      = &this->sparsedmudX_[iAtm];
   RealSparseMatrix *MapdY      = &this->sparsedmudY_[iAtm];
   RealSparseMatrix *MapdZ      = &this->sparsedmudZ_[iAtm];

//   if (this->screenVxc && Map->col(ipts).norm() < this->epsScreen)
//     return;

// Eventually Decleare outside as STmp;
   RealMatrix dSTmp(this->nBasis_,this->nBasis_);        ///< d(Overlap) at grid point
   dSTmp.setZero();
   RealMatrix dSTmpR(this->nBasis_,this->nBasis_);        ///< d(Overlap) at grid point
   dSTmpR.setZero();

   double rhor  = 0.0;  // Total density at point
   double rhorA = 0.0;  // alpha density at point
   double rhorB = 0.0;  // beta  density at point
   double gammaAA = 0.0;  // alpha  density at point der
   double gammaBB = 0.0;  // beta  density at point  der
   double tmpval = 0.0;  // beta  density at point  der
   int    nDer = 0   ;    // Order of Der
   bool   RHF  = this->Ref_ == RHF;
   bool   doTCS  = this->Ref_ == TCS;
// RealMatrix overlapR_(this->nBasis_,this->nBasis_);        ///< Overlap at grid point
// overlapR_.setZero();
// STmp->setZero();
   std::array<double,3>  epsMuCor = {0.0,0.0,0.0}; ///< {energydens_corr, potential_corr_alpha, potential_corr_B}
   std::array<double,3>  epsMuExc = {0.0,0.0,0.0}; ///< {energydend_exchange, potential_exchenge_alpha, potential_exchange_beta}

//   Build Overlap
//  overlapR_ = Map->col(ipts)*Map->col(ipts).transpose();
    (*STmp) = Map->col(ipts)*Map->col(ipts).transpose();
//  Handle the total density at r for RKS or UKS
    if(!this->isClosedShell && this->Ref_ != TCS) {
//    rhorA = overlapR_.frobInner(this->densityA()->conjugate());
//    rhorB = overlapR_.frobInner(this->densityB()->conjugate());
      rhorA = STmp->frobInner(this->densityA()->conjugate());
      rhorB = STmp->frobInner(this->densityB()->conjugate());
      rhor = rhorA + rhorB;
      if (nDer == 1 ){
      dSTmp   = 2.0*Map->col(ipts)*MapdX->col(ipts).transpose();
      tmpval  = dSTmp.frobInner(this->densityA()->conjugate());
      gammaAA = tmpval*tmpval;
      tmpval  = dSTmp.frobInner(this->densityB()->conjugate());
      gammaBB = tmpval*tmpval;
      dSTmp   = 2.0*Map->col(ipts)*MapdY->col(ipts).transpose();
      tmpval  = dSTmp.frobInner(this->densityA()->conjugate());
      gammaAA += tmpval*tmpval;
      tmpval  = dSTmp.frobInner(this->densityB()->conjugate());
      gammaBB += tmpval*tmpval;
      dSTmp   = 2.0*Map->col(ipts)*MapdZ->col(ipts).transpose();
      tmpval  = dSTmp.frobInner(this->densityA()->conjugate());
      gammaAA += tmpval*tmpval;
      tmpval  = dSTmp.frobInner(this->densityB()->conjugate());
      gammaBB += tmpval*tmpval;
      }
    } else {
      rhor = STmp->frobInner(this->densityA()->conjugate()) ;
//    rhor = overlapR_.frobInner(this->densityA()->conjugate()) ;
      if (nDer == 1 ){
      dSTmp    = 2.0*(Map->col(ipts)*MapdX->col(ipts).transpose());
      tmpval   = dSTmp.frobInner(this->densityA()->conjugate());
      gammaAA += tmpval*tmpval;
      dSTmp    = 2.0*Map->col(ipts)*MapdY->col(ipts).transpose();
      tmpval   = dSTmp.frobInner(this->densityA()->conjugate());
      gammaAA += tmpval*tmpval;
      dSTmp    = 2.0*Map->col(ipts)*MapdZ->col(ipts).transpose();
      tmpval   = dSTmp.frobInner(this->densityA()->conjugate());
      gammaAA += tmpval*tmpval;
      dSTmpR   = Map->col(ipts)*MapdY->col(ipts).transpose();
      dSTmpR  += MapdY->col(ipts)*Map->col(ipts).transpose();
//      dSTmp   = 2.0*(Map->col(ipts)*MapdX->col(ipts).transpose());
//      dSTmp  += 2.0*(Map->col(ipts)*MapdY->col(ipts).transpose());
//      dSTmp  += 2.0*(Map->col(ipts)*MapdZ->col(ipts).transpose());
//      tmpval  = dSTmp.frobInner(this->densityA()->conjugate());
//      gammaBB = tmpval*tmpval;
//      if ( std::abs(gammaAA-gammaBB) > 1.0e-8 ) cout << abs(gammaAA-gammaBB) <<endl;
//      if ( std::abs(gammaAA) > 1.0e-8 ) cout << gammaAA <<endl;
      }
    }
//  Handle numerical instability if screening on
    if (this->screenVxc ) {
//    check if are noise
      if(rhor    <= 0.0 ) {
        if((std::abs(rhor)) <= 1.0e-10) {
          return;
        }else{ 
          CErr("Numerical noise in the density");
        }
//    skyp points based on small density
      }else if(rhor < this->epsScreen){
        return;
      }
    }
//this if statement prevent numerical instability with zero guesses
    if (rhor > 1.0e-20) {     
//  Exchange
      if (this->ExchKernel_ != NOEXCH) {
        if(!this->isClosedShell && this->Ref_ != TCS){
            epsMuExc = this->formVEx(rhor,this->spindens(rhorA,rhorB));
//          (*VXB)  += ((*WeightsMap).coeff(ipts,0))*overlapR_*epsMuExc[2];
            (*VXB)  += ((*WeightsMap).coeff(ipts,0))*(*STmp)*epsMuExc[2];
        } else {
            epsMuExc = this->formVEx(rhor,0.0);
        }
//      (*VXA)  += ((*WeightsMap).coeff(ipts,0))*overlapR_*epsMuExc[1];
        (*VXA)  += ((*WeightsMap).coeff(ipts,0))*(*STmp)*epsMuExc[1];
        energyX += ((*WeightsMap).coeff(ipts,0))*rhor*epsMuExc[0];
      }
//  Correlation
      if (this->CorrKernel_ != NOCORR) {
        if(!this->isClosedShell && this->Ref_ != TCS){
          epsMuCor = this->formVC(rhor,(this->spindens(rhorA,rhorB)));
//        (*VCB)  += ((*WeightsMap).coeff(ipts,0))*overlapR_*epsMuCor[2];
          (*VCB)  += ((*WeightsMap).coeff(ipts,0))*(*STmp)*epsMuCor[2];
        } else{
          epsMuCor = this->formVC(rhor,0.0);
        }
//      (*VCA)  += ((*WeightsMap).coeff(ipts,0))*overlapR_*epsMuCor[1];
        (*VCA)  += ((*WeightsMap).coeff(ipts,0))*(*STmp)*epsMuCor[1];
//Debug       (*VCA)  += ((*WeightsMap).coeff(ipts,0))*(dSTmpR);
//Debug      this->nElectrons_ += ((*WeightsMap).coeff(ipts,0))*STmp->frobInner(this->densityA()->conjugate()) ;
        energyC += ((*WeightsMap).coeff(ipts,0))*rhor*epsMuCor[0];
      }

    }
}; //END

//----------------------------//
// form the Vxc matrix        //
//----------------------------//

/*
template<>
void SingleSlater<double>::formVXC(){
////Timing
    this->basisset_->duration_1 = std::chrono::seconds(0) ;
    this->basisset_->duration_2 = std::chrono::seconds(0) ;
    this->basisset_->duration_3 = std::chrono::seconds(0) ;
    this->basisset_->duration_4 = std::chrono::seconds(0) ;
    this->basisset_->duration_5 = std::chrono::seconds(0) ;
//    this->duration_1 = std::chrono::seconds(0) ;
//    this->duration_2 = std::chrono::seconds(0) ;
//    this->duration_3 = std::chrono::seconds(0) ;
//    this->duration_4 = std::chrono::seconds(0) ;
//    this->duration_5 = std::chrono::seconds(0) ;
//    this->duration_6 = std::chrono::seconds(0) ;
//    this->duration_7 = std::chrono::seconds(0) ;
//    this->duration_8 = std::chrono::seconds(0) ;
////T


    int nAtom   = this->molecule_->nAtoms(); // Number of Atoms
    // Total Number of grid point for each center
    this->ngpts = this->nRadDFTGridPts_*this->nAngDFTGridPts_; 

    int nPtsPerThread = this->ngpts / omp_get_max_threads();

    //double weight  = 0.0;                            
    //double bweight = 0.0;
    //TF LDA Prefactor (for Vx)
    double CxVx = -(std::pow((3.0/math.pi),(1.0/3.0)));    
    double CxEn =  (3.0/4.0);      //TF LDA Prefactor to finish the X-Energy
    double val = 4.0*math.pi*CxVx;
    this->totalEx = 0.0;    // Total Exchange Energy
    this->totalEcorr = 0.0; // Total Correlation Energy
    //bool nodens;
    std::vector<bool> tmpnull(this->basisset_->nShell()+1);
    OneDGrid * Rad ;
//   
// *  Generate grids 
// *
// *    Raw grid, it has to be centered and integrated over each center and 
// *    centered over each atom
//

//  Evaluate average cutoff radia for shells given epsScreen - if screenVxc ON
    if (this->screenVxc ) {
      this->basisset_->radcut(this->epsScreen, this->maxiter, this->epsConv);
    } else {
      std::fill(tmpnull.begin(),tmpnull.end(),true);
    }
//  Select Radial Grid
    if (this->dftGrid_ == GAUSSCHEB)  
      Rad = new GaussChebyshev1stGridInf(this->nRadDFTGridPts_,0.0,1.0);   
    else if (this->dftGrid_ == EULERMACL) 
      Rad = new  EulerMaclaurinGrid(this->nRadDFTGridPts_,0.0,1.0);   

    LebedevGrid GridLeb(this->nAngDFTGridPts_);   // Angular Grid
    GridLeb.genGrid();                            // Generate Angular Grid
    this->vXA()->setZero();   // Set to zero every occurence of the SCF
    this->vCorA()->setZero(); // Set to zero every occurence of the SCF
    if(!this->isClosedShell && this->Ref_ != TCS) this->vXB()->setZero();
    if(!this->isClosedShell && this->Ref_ != TCS) this->vCorB()->setZero();
    // Loop over atomic centers
    for(int iAtm = 0; iAtm < nAtom; iAtm++){
      Rad->genGrid(); 
      // The Radial grid is generated and scaled for each atom
      Rad->atomGrid((elements[this->molecule_->index(iAtm)].sradius)) ;  
      TwoDGrid Raw3Dg(this->ngpts,Rad,&GridLeb);             
      //Center the Grid at iAtom
      Raw3Dg.centerGrid(
        (*this->molecule_->cart())(0,iAtm),
        (*this->molecule_->cart())(1,iAtm),
        (*this->molecule_->cart())(2,iAtm)
      );
      // Loop over grid points
      for(int ipts = 0; ipts < this->ngpts; ipts++){
        //cout << ipts << " " << ipts/nPtsPerThread << endl;
////T   
//   auto start_5 = std::chrono::high_resolution_clock::now();  // Timing weights
////T
        bool nodens = false;
        // Evaluate each Becke fuzzy call weight, normalize it and muliply by 
        //   the Raw grid weight at that point
        auto bweight = (this->formBeckeW((Raw3Dg.gridPtCart(ipts)),iAtm)) 
                     / (this->normBeckeW(Raw3Dg.gridPtCart(ipts))) ;
        auto weight = Raw3Dg.getweightsGrid(ipts) * bweight;
        
////T
//   auto finish_5 = std::chrono::high_resolution_clock::now();  
//   this->duration_5 += finish_5 - start_5;
////T
        // Build the Vxc for the ipts grid point 
        //  ** Vxc will be ready at the end of the two loop, to be finalized ** 
        if (this->screenVxc ) {
          auto mapRad_ = this->basisset_->MapGridBasis(Raw3Dg.gridPtCart(ipts));
          if (mapRad_[0] || (bweight < this->epsScreen)) 
            nodens = true;
          if(!nodens) 
            this->evalVXC((Raw3Dg.gridPtCart(ipts)),weight,mapRad_,this->totalEx, this->totalEcorr, 
                           this->vXA_.get(),this->vXB_.get(),this->vCorA_.get(),this->vCorB_.get() );
        } else {
          this->evalVXC((Raw3Dg.gridPtCart(ipts)),weight,tmpnull,this->totalEx, this->totalEcorr, 
                         this->vXA_.get(),this->vXB_.get(),this->vCorA_.get(),this->vCorB_.get() );
        }

      } // loop ipts
    } // loop natoms
////T   
//   auto start_6 = std::chrono::high_resolution_clock::now();  // Timing Digestion VXC
////T

    //  Finishing the Vxc using the TF factor and the integration 
    //    prefactor over a solid sphere
    (*this->vXA())    =  val * (*this->vXA());
    this->totalEx     =  val * CxEn * (this->totalEx);
    (*this->vCorA())  =  4.0 * math.pi * (*this->vCorA());
    if(!this->isClosedShell && this->Ref_ != TCS) 
      (*this->vCorB())  =  4.0 * math.pi * (*this->vCorB());
    this->totalEcorr  =  4.0 * math.pi * (this->totalEcorr);
    // For open shell averything has to be scaled by 2^(1/3)
    if(!this->isClosedShell && this->Ref_ != TCS){
      (*this->vXA()) *= std::pow(2.0,(1.0/3.0));  
      (*this->vXB()) *= std::pow(2.0,(1.0/3.0)) * val;
    }
////T
//   auto finish_6 = std::chrono::high_resolution_clock::now();  
//   this->duration_6 += finish_6 - start_6;
////T

    if(this->printLevel_ >= 3) {
      prettyPrint(this->fileio_->out,(*this->vXA()),"LDA Vx alpha");
      prettyPrint(this->fileio_->out,(*this->vCorA()),"Vc Vc alpha");
      if(!this->isClosedShell && this->Ref_ != TCS) 
        prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vx beta");

      this->fileio_->out << "Total LDA Ex ="    << this->totalEx 
                         << " Total VWN Corr= " << this->totalEcorr << endl;
//    this->fileio_->out << "Weights Evaluation       Total Time " << this->duration_5.count() <<endl;
//    this->fileio_->out << "Overlap Alloc + set Zero Total Time " << this->duration_2.count() <<endl;
//    this->fileio_->out << "Overlap ProdEval         Total Time " << this->duration_7.count() <<endl;
//    this->fileio_->out << "Overlap BuildDend        Total Time " << this->duration_8.count() <<endl;
    this->fileio_->out << "Overlap Creation Part1(a)   Total Time " 
                       << this->basisset_->duration_1.count() << endl;
    this->fileio_->out << "Overlap Creation Part1(b)   Total Time " 
                       << this->basisset_->duration_4.count() << endl;
    this->fileio_->out << "Overlap Creation Part1(c)   Total Time " 
                       << this->basisset_->duration_5.count() << endl;
    this->fileio_->out << "Overlap Creation Part2   Total Time " 
                       << this->basisset_->duration_2.count() << endl;
    this->fileio_->out << "Overlap Creation Part3   Total Time " 
                       << this->basisset_->duration_3.count() << endl;
//    this->fileio_->out << "Overlap Creation         Total Time " << this->duration_4.count() <<endl;
//    this->fileio_->out << "Overlap Contraction      Total Time " << this->duration_3.count() <<endl;
//    this->fileio_->out << "Form (Vx + Vc)           Total Time " << this->duration_1.count() <<endl;
//    this->fileio_->out << "Vxc Digestion            Total Time " << this->duration_6.count() <<endl;
//  CErr("DIE DIE DIE");
    }

//  Cleaning
    delete Rad;
}; //End
*/

//----------------------------//
// form the Vxc matrix        //
//----------------------------//

template<>
void SingleSlater<double>::formVXC(){
////Timing
//    this->basisset_->duration_1 = std::chrono::seconds(0) ;
//    this->basisset_->duration_2 = std::chrono::seconds(0) ;
//    this->basisset_->duration_3 = std::chrono::seconds(0) ;
//    this->basisset_->duration_4 = std::chrono::seconds(0) ;
//    this->basisset_->duration_5 = std::chrono::seconds(0) ;
//    this->duration_1 = std::chrono::seconds(0) ;
//    this->duration_2 = std::chrono::seconds(0) ;
//    this->duration_3 = std::chrono::seconds(0) ;
//    this->duration_4 = std::chrono::seconds(0) ;
//    this->duration_5 = std::chrono::seconds(0) ;
//    this->duration_6 = std::chrono::seconds(0) ;
//    this->duration_7 = std::chrono::seconds(0) ;
//    this->duration_8 = std::chrono::seconds(0) ;
////T
//    this->screenVxc  = false;

    int nAtom   = this->molecule_->nAtoms();                    // Number of Atoms
    this->ngpts = this->nRadDFTGridPts_*this->nAngDFTGridPts_;  // Number of grid point for each center
    int nPtsPerThread = this->ngpts / omp_get_max_threads();    //  Number of Threads
//  Get the number of tread to dimension this last two
//  std::array<double,1> tmpEnergyEx  = {0.0} ;
//  std::array<double,1> tmpEnergyCor = {0.0} ;
    std::vector<double> tmpEnergyEx(omp_get_max_threads())  ;
    std::vector<double> tmpEnergyCor(omp_get_max_threads()) ;
    std::vector<int> tmpnpts(omp_get_max_threads()) ;

    int nRHF;
    if(this->isClosedShell || this->Ref_ == TCS) nRHF = 1;
    else    nRHF = 2;
    std::vector<std::vector<RealMatrix>> 
      tmpVX(nRHF,std::vector<RealMatrix>(omp_get_max_threads(),
              RealMatrix::Zero(this->nBasis_,this->nBasis_)
      )
    );
    std::vector<std::vector<RealMatrix>> 
      tmpVC(nRHF,std::vector<RealMatrix>(omp_get_max_threads(),
              RealMatrix::Zero(this->nBasis_,this->nBasis_)
      )
    );

    double CxVx  = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));  //TF LDA Prefactor (for Vx)  
    double val = 4.0*math.pi*CxVx;                                  // to take into account Ang Int
    this->totalEx = 0.0;    // Zero out Total Exchange Energy
    this->totalEcorr = 0.0; // Zero out Total Correlation Energy
    this->vXA()->setZero();   // Set to zero every occurence of the SCF
    this->vCorA()->setZero(); // Set to zero every occurence of the SCF
    if(!this->isClosedShell && this->Ref_ != TCS) {
      this->vXB()->setZero();
      this->vCorB()->setZero();
    }
    //bool nodens;
    std::vector<bool> tmpnull(this->basisset_->nShell()+1);
    OneDGrid * Rad ;                              // Pointer for Radial Grid
    LebedevGrid GridLeb(this->nAngDFTGridPts_);   // Angular Grid
/*  
 *  Generate grids 
 *
 *    Raw grid, it has to be centered and integrated over each center and 
 *    centered over each atom
 */

//  Evaluate average cutoff radia for shells given epsScreen - if screenVxc ON
    if (this->screenVxc ) {
      this->basisset_->radcut(this->epsScreen, this->maxiter, this->epsConv);
    } else {
      std::fill(tmpnull.begin(),tmpnull.end(),true);
    }
//  Select Radial Grid
    if (this->dftGrid_ == GAUSSCHEB)  
      Rad = new GaussChebyshev1stGridInf(this->nRadDFTGridPts_,0.0,1.0);   
    else if (this->dftGrid_ == EULERMACL) 
      Rad = new  EulerMaclaurinGrid(this->nRadDFTGridPts_,0.0,1.0);   
//  Generare Angular Grid
    GridLeb.genGrid();                            

//  Generate Sparse Matrix
    auto batch_dft = [&] (int thread_id,int iAtm, TwoDGrid &Raw3Dg) {
      auto loopSt = nPtsPerThread * thread_id;
      auto loopEn = nPtsPerThread * (thread_id + 1);
      if (thread_id == (omp_get_max_threads() - 1))
        loopEn = this->ngpts;
      for(int ipts = loopSt; ipts < loopEn; ipts++){
//      printf("%d_%d_%d_%d\n", thread_id, ipts/nPtsPerThread,  ipts, iAtm);
//      if(ipts/nPtsPerThread != thread_id) continue;
        tmpnpts[thread_id]++;
        bool nodens = false;
        // Evaluate each Becke fuzzy call weight, normalize it and muliply by 
        //   the Raw grid weight at that point
        auto bweight = (this->formBeckeW((Raw3Dg.gridPtCart(ipts)),iAtm)) 
                     / (this->normBeckeW(Raw3Dg.gridPtCart(ipts))) ;
        auto weight = Raw3Dg.getweightsGrid(ipts) * bweight;
        
        // Build the Vxc for the ipts grid point 
        //  ** Vxc will be ready at the end of the two loop, to be finalized ** 
        if (this->screenVxc ) {
          auto mapRad_ = this->basisset_->MapGridBasis(Raw3Dg.gridPtCart(ipts));
          if (mapRad_[0] || (bweight < this->epsScreen)) 
            nodens = true;
          if(!nodens) 
            this->evalVXC((Raw3Dg.gridPtCart(ipts)),weight,mapRad_,
              tmpEnergyEx[thread_id],tmpEnergyCor[thread_id],&tmpVX[0][thread_id],
              &tmpVX[1][thread_id],&tmpVC[0][thread_id],&tmpVC[1][thread_id]);
        } else {
            this->evalVXC((Raw3Dg.gridPtCart(ipts)),weight,tmpnull,
              tmpEnergyEx[thread_id],tmpEnergyCor[thread_id],&tmpVX[0][thread_id],
              &tmpVX[1][thread_id],&tmpVC[0][thread_id],&tmpVC[1][thread_id]);
        }

      } // loop ipts
    }; // end batch_dft

   for(int iAtm = 0; iAtm < nAtom; iAtm++){
      std::fill(tmpEnergyEx.begin(),tmpEnergyEx.end(),0.0);
      std::fill(tmpEnergyCor.begin(),tmpEnergyCor.end(),0.0);
      std::fill(tmpnpts.begin(),tmpnpts.end(),0);
      Rad->genGrid(); 
      // The Radial grid is generated and scaled for each atom
      Rad->atomGrid((elements[this->molecule_->index(iAtm)].sradius)) ;  
      TwoDGrid Raw3Dg(this->ngpts,Rad,&GridLeb);             
      //Center the Grid at iAtom
      Raw3Dg.centerGrid(
        (*this->molecule_->cart())(0,iAtm),
        (*this->molecule_->cart())(1,iAtm),
        (*this->molecule_->cart())(2,iAtm)
      );
//    if (this->screenVxc ) this->genSparseBasisMap(Raw3Dg,iAtm); 
//    else this->genSparseBasisMap(Raw3Dg,iAtm);
//     CErr("Final");
    
   #ifdef _OPENMP
     #pragma omp parallel
     {
       int thread_id = omp_get_thread_num();
       tmpVX[0][thread_id].setZero();  
       tmpVC[0][thread_id].setZero();  
       if(!this->isClosedShell && this->Ref_ != TCS) {
         tmpVX[1][thread_id].setZero();  
         tmpVC[1][thread_id].setZero();  
       }
       batch_dft(thread_id,iAtm,Raw3Dg);
     }
   #else
     tmpVX[0][0].setZero();  
     tmpVC[0][0].setZero();  
     if(!this->isClosedShell && this->Ref_ != TCS) {
       tmpVX[1][0].setZero();  
       tmpVC[1][0].setZero();  
     }
     batch_dft(0,iAtm,Raw3Dg);
   #endif
      for(auto iThread = 0; iThread < omp_get_max_threads(); iThread++) {
        (*this->vXA())   += tmpVX[0][iThread];
        (*this->vCorA()) += tmpVC[0][iThread];
        this->totalEx += tmpEnergyEx[iThread];
        this->totalEcorr += tmpEnergyCor[iThread];
        if(!this->isClosedShell && this->Ref_ != TCS) {
          (*this->vXB())   += tmpVX[1][iThread];
          (*this->vCorB()) += tmpVC[1][iThread];
        }
      }
    }; //loop atoms

    //  Finishing the Vxc using the TF factor and the integration 
    //    prefactor over a solid sphere
    (*this->vXA())      =  val * (*this->vXA());
    (*this->vCorA())    =  4.0 * math.pi * (*this->vCorA());
    this->totalEx       =  val * this->totalEx;
    this->totalEcorr    =  4.0 * math.pi * (this->totalEcorr);
    if(!this->isClosedShell && this->Ref_ != TCS) {
        (*this->vCorB())  =  4.0 * math.pi * (*this->vCorB());
        (*this->vXB())    =  val * (*this->vXB());
      }

    if(this->printLevel_ >= 3) {
      prettyPrint(this->fileio_->out,(*this->vXA()),"LDA Vx alpha");
      prettyPrint(this->fileio_->out,(*this->vCorA()),"Vc Vc alpha");
      if(!this->isClosedShell && this->Ref_ != TCS) 
        prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vx beta");

      this->fileio_->out << "Total LDA Ex ="    << this->totalEx 
                         << " Total VWN Corr= " << this->totalEcorr << endl;
/*
    this->fileio_->out << "Overlap Creation Part1(a)   Total Time " 
                       << this->basisset_->duration_1.count() << endl;
    this->fileio_->out << "Overlap Creation Part1(b)   Total Time " 
                       << this->basisset_->duration_4.count() << endl;
    this->fileio_->out << "Overlap Creation Part1(c)   Total Time " 
                       << this->basisset_->duration_5.count() << endl;
    this->fileio_->out << "Overlap Creation Part2   Total Time " 
                       << this->basisset_->duration_2.count() << endl;
    this->fileio_->out << "Overlap Creation Part3   Total Time " 
                       << this->basisset_->duration_3.count() << endl;
*/
//    this->fileio_->out << "Overlap Creation         Total Time " << this->duration_4.count() <<endl;
//    this->fileio_->out << "Overlap Contraction      Total Time " << this->duration_3.count() <<endl;
//    this->fileio_->out << "Form (Vx + Vc)           Total Time " << this->duration_1.count() <<endl;
//    this->fileio_->out << "Vxc Digestion            Total Time " << this->duration_6.count() <<endl;
    }

//  Cleaning
    delete Rad;
}; //End


template<>
void SingleSlater<double>::formVXC_store(){
//    this->nElectrons_= 0.0;
    int nAtom   = this->molecule_->nAtoms();                    // Number of Atoms
    this->ngpts = this->nRadDFTGridPts_*this->nAngDFTGridPts_;  // Number of grid point for each center
    int nPtsPerThread = this->ngpts / omp_get_max_threads();    //  Number of Threads
    std::vector<double> tmpEnergyEx(omp_get_max_threads())  ;
    std::vector<double> tmpEnergyCor(omp_get_max_threads()) ;
    std::vector<int> tmpnpts(omp_get_max_threads()) ;
    int nRHF;
    if(this->isClosedShell || this->Ref_ == TCS) nRHF = 1;
    else    nRHF = 2;
    std::vector<std::vector<RealMatrix>> 
      tmpVX(nRHF,std::vector<RealMatrix>(omp_get_max_threads(),
              RealMatrix::Zero(this->nBasis_,this->nBasis_)
      )
    );
    std::vector<std::vector<RealMatrix>> 
      tmpVC(nRHF,std::vector<RealMatrix>(omp_get_max_threads(),
              RealMatrix::Zero(this->nBasis_,this->nBasis_)
      )
    );
    std::vector<RealMatrix> overlapR_(omp_get_max_threads(),RealMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));

    double CxVx  = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));  //TF LDA Prefactor (for Vx)  
    double val = 4.0*math.pi*CxVx;                                  // to take into account Ang Int
    this->totalEx    = 0.0;   // Zero out Total Exchange Energy
    this->totalEcorr = 0.0;   // Zero out Total Correlation Energy
    this->vXA()->setZero();   // Set to zero every occurence of the SCF
    this->vCorA()->setZero(); // Set to zero every occurence of the SCF
    if(!this->isClosedShell && this->Ref_ != TCS) {
      this->vXB()->setZero();
      this->vCorB()->setZero();
    }

/*Parellel
//  Generate Sparse Matrix
    auto batch_dft = [&] (int thread_id,int iAtm, TwoDGrid &Raw3Dg) {
      auto loopSt = nPtsPerThread * thread_id;
      auto loopEn = nPtsPerThread * (thread_id + 1);
      if (thread_id == (omp_get_max_threads() - 1))
        loopEn = this->ngpts;
      for(int ipts = loopSt; ipts < loopEn; ipts++){
//      printf("%d_%d_%d_%d\n", thread_id, ipts/nPtsPerThread,  ipts, iAtm);
//      if(ipts/nPtsPerThread != thread_id) continue;
        tmpnpts[thread_id]++;
        bool nodens = false;
        // Evaluate each Becke fuzzy call weight, normalize it and muliply by 
        //   the Raw grid weight at that point
        auto bweight = (this->formBeckeW((Raw3Dg.gridPtCart(ipts)),iAtm)) 
                     / (this->normBeckeW(Raw3Dg.gridPtCart(ipts))) ;
        auto weight = Raw3Dg.getweightsGrid(ipts) * bweight;
        
        // Build the Vxc for the ipts grid point 
        //  ** Vxc will be ready at the end of the two loop, to be finalized ** 
        if (this->screenVxc ) {
          auto mapRad_ = this->basisset_->MapGridBasis(Raw3Dg.gridPtCart(ipts));
          if (mapRad_[0] || (bweight < this->epsScreen)) 
            nodens = true;
          if(!nodens) 
            this->evalVXC((Raw3Dg.gridPtCart(ipts)),weight,mapRad_,
              tmpEnergyEx[thread_id],tmpEnergyCor[thread_id],&tmpVX[0][thread_id],
              &tmpVX[1][thread_id],&tmpVC[0][thread_id],&tmpVC[1][thread_id]);
        } else {
            this->evalVXC((Raw3Dg.gridPtCart(ipts)),weight,tmpnull,
              tmpEnergyEx[thread_id],tmpEnergyCor[thread_id],&tmpVX[0][thread_id],
              &tmpVX[1][thread_id],&tmpVC[0][thread_id],&tmpVC[1][thread_id]);
        }

      } // loop ipts
    }; // end batch_dft
*/
///for(int iAtm = 0; iAtm < nAtom; iAtm++){
///   for(int ipts = 0; ipts < this->ngpts; ipts++){

/*    Parallel
      std::fill(tmpEnergyEx.begin(),tmpEnergyEx.end(),0.0);
      std::fill(tmpEnergyCor.begin(),tmpEnergyCor.end(),0.0);
      std::fill(tmpnpts.begin(),tmpnpts.end(),0);
   #ifdef _OPENMP
     #pragma omp parallel
     {
       int thread_id = omp_get_thread_num();
       tmpVX[0][thread_id].setZero();  
       tmpVC[0][thread_id].setZero();  
       if(!this->isClosedShell && this->Ref_ != TCS) {
         tmpVX[1][thread_id].setZero();  
         tmpVC[1][thread_id].setZero();  
       }
       batch_dft(thread_id,iAtm,Raw3Dg);
     }
   #else
     tmpVX[0][0].setZero();  
     tmpVC[0][0].setZero();  
     if(!this->isClosedShell && this->Ref_ != TCS) {
       tmpVX[1][0].setZero();  
       tmpVC[1][0].setZero();  
     }
     batch_dft(0,iAtm,Raw3Dg);
   #endif
      for(auto iThread = 0; iThread < omp_get_max_threads(); iThread++) {
        (*this->vXA())   += tmpVX[0][iThread];
        (*this->vCorA()) += tmpVC[0][iThread];
        this->totalEx += tmpEnergyEx[iThread];
        this->totalEcorr += tmpEnergyCor[iThread];
        if(!this->isClosedShell && this->Ref_ != TCS) {
          (*this->vXB())   += tmpVX[1][iThread];
          (*this->vCorB()) += tmpVC[1][iThread];
        }
        this->evalVXC_store(iAtm,ipts,this->totalEx,tmpEnergyCor[thread_id],
              &tmpVX[0][thread_id],&tmpVX[1][thread_id],&tmpVC[0][thread_id],
              &tmpVC[1][thread_id]);
      }
*/
  // Build the Vxc for the ipts grid point 
  //  ** Vxc will be ready at the end of the two loop, to be finalized ** 

///     this->evalVXC_store(iAtm,ipts,this->totalEx,this->totalEcorr,
///           (this->vXA()),(this->vXB()),(this->vCorA()),(this->vCorB()),
///           &overlapR_);
///   }; //loop over gridpts
/// }; //loop atoms

    
    std::vector<std::chrono::duration<double>> thread_timers(omp_get_max_threads());
    auto batch_dft = [&] (int thread_id,int iAtm) {
      auto loopSt = nPtsPerThread * thread_id;
      auto loopEn = nPtsPerThread * (thread_id + 1);
      RealSparseMatrix *DoRhoMap   = &this->sparseDoRho_[iAtm];
//      auto start = std::chrono::high_resolution_clock::now();
//    for(int ipts = loopSt; ipts < loopEn; ipts++){
      for(auto ipts = 0; ipts < this->ngpts; ipts++) {
        if(ipts % omp_get_max_threads() != thread_id) continue;
        if( this->screenVxc && ((*DoRhoMap).coeff(ipts,0) < 1) ) continue;
        this->evalVXC_store(iAtm,ipts,tmpEnergyEx[thread_id],tmpEnergyCor[thread_id],
              &tmpVX[0][thread_id],&tmpVX[1][thread_id],&tmpVC[0][thread_id],
              &tmpVC[1][thread_id],&overlapR_[thread_id]);
      } // loop ipts
//      auto finish = std::chrono::high_resolution_clock::now();
//      thread_timers[thread_id] = finish - start;
    }; // batch_dft

    for(int iAtm = 0; iAtm < nAtom; iAtm++){
      std::fill(tmpEnergyEx.begin(),tmpEnergyEx.end(),0.0);
      std::fill(tmpEnergyCor.begin(),tmpEnergyCor.end(),0.0);
      std::fill(tmpnpts.begin(),tmpnpts.end(),0);
    #ifdef _OPENMP
      #pragma omp parallel
      {
        int thread_id = omp_get_thread_num();
        tmpVX[0][thread_id].setZero();  
        tmpVC[0][thread_id].setZero();  
        if(!this->isClosedShell && this->Ref_ != TCS) {
          tmpVX[1][thread_id].setZero();  
          tmpVC[1][thread_id].setZero();  
        }
        batch_dft(thread_id,iAtm);
      }
    #else
      tmpVX[0][0].setZero();  
      tmpVC[0][0].setZero();  
      if(!this->isClosedShell && this->Ref_ != TCS) {
        tmpVX[1][0].setZero();  
        tmpVC[1][0].setZero();  
      }
      batch_dft(0,iAtm);
    #endif
      for(auto iThread = 0; iThread < omp_get_max_threads(); iThread++) {
        (*this->vXA())   += tmpVX[0][iThread];
        (*this->vCorA()) += tmpVC[0][iThread];
        this->totalEx += tmpEnergyEx[iThread];
        this->totalEcorr += tmpEnergyCor[iThread];
        if(!this->isClosedShell && this->Ref_ != TCS) {
          (*this->vXB())   += tmpVX[1][iThread];
          (*this->vCorB()) += tmpVC[1][iThread];
        }
      }
    }; // loop over atoms

/* DBWY Thread Timings
    cout << "Thread Timings" << endl;
    for(auto i = thread_timers.begin(); i != thread_timers.end(); i++)
      cout << "   " << i->count() << endl;
*/

    //  Finishing the Vxc using the TF factor and the integration 
    //    prefactor over a solid sphere
    (*this->vXA())      =  val * (*this->vXA());
    (*this->vCorA())    =  4.0 * math.pi * (*this->vCorA());
    this->totalEx       =  val * this->totalEx;
    this->totalEcorr    =  4.0 * math.pi * (this->totalEcorr);
    if(!this->isClosedShell && this->Ref_ != TCS) {
        (*this->vCorB())  =  4.0 * math.pi * (*this->vCorB());
        (*this->vXB())    =  val * (*this->vXB());
      }

    if(this->printLevel_ >= 3) {
      prettyPrint(this->fileio_->out,(*this->vXA()),"LDA Vx alpha");
      prettyPrint(this->fileio_->out,(*this->vCorA()),"Vc Vc alpha");
      if(!this->isClosedShell && this->Ref_ != TCS) 
        prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vx beta");

      this->fileio_->out << "Total LDA Ex ="    << this->totalEx 
                         << " Total VWN Corr= " << this->totalEcorr << endl;
    }
/*
      cout << "nElectrons " << 10.0-4.0*math.pi*this->nElectrons_ << endl;
      prettyPrint(cout,(*this->vCorA()),"Vc Vc alpha");
      cout << "Max " << (*this->vCorA()).lpNorm<Infinity>() << endl;
    CErr();
*/
}; //End

} // Namespace ChronusQ
