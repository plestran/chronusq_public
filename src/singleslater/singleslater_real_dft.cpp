/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
#include<dft.h>
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
//Populate cutoff radii vector over shell
  this->basisset_->radcut(this->epsScreen, this->maxiter, this->epsConv);
  this->ngpts = this->nRadDFTGridPts_*this->nAngDFTGridPts_;  // Number of grid point for each center
  OneDGrid * Rad ;                              // Pointer for Radial Grid
  LebedevGrid GridLeb(this->nAngDFTGridPts_);   // Angular Grid
  int nDer = 0 ; //Order of differenzation
  if (this->isGGA) nDer = 1;
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
void SingleSlater<double>::genSparseRcrosP(){
//Populate cutoff radii vector over shell
  this->screenVxc = false;
  this->basisset_->radcut(this->epsScreen, this->maxiter, this->epsConv);
  this->ngpts = this->nRadDFTGridPts_*this->nAngDFTGridPts_;  // Number of grid point for each center
  OneDGrid * Rad ;                              // Pointer for Radial Grid
  LebedevGrid GridLeb(this->nAngDFTGridPts_);   // Angular Grid
  int nDer = 1 ; //Order of differenzation
//  if (this->isGGA) nDer = 1;
  if (this->dftGrid_ == GAUSSCHEB)  
    Rad = new GaussChebyshev1stGridInf(this->nRadDFTGridPts_,0.0,1.0);   
  else if (this->dftGrid_ == EULERMACL) 
    Rad = new  EulerMaclaurinGrid(this->nRadDFTGridPts_,0.0,1.0);   
//Generare Angular Grid
  GridLeb.genGrid();                            
  std::unique_ptr<RealMatrix>  rdotpX;        ///< r cross p - X at grid point
  std::unique_ptr<RealMatrix>  rdotpY;        ///< r cross p - Y at grid point
  std::unique_ptr<RealMatrix>  rdotpZ;        ///< r cross p - Z at grid point
  rdotpX = std::unique_ptr<RealMatrix>(
    new RealMatrix(this->nBasis_,this->nBasis_));
  rdotpY = std::unique_ptr<RealMatrix>(
    new RealMatrix(this->nBasis_,this->nBasis_));
  rdotpZ = std::unique_ptr<RealMatrix>(
    new RealMatrix(this->nBasis_,this->nBasis_));
  rdotpX->setZero();
  rdotpY->setZero();
  rdotpZ->setZero();
  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    this->sparseMap_.push_back(RealSparseMatrix(this->nTCS_*this->nBasis_,this->ngpts));
    this->sparseWeights_.push_back(RealSparseMatrix(this->ngpts,1));
    this->sparseDoRho_.push_back(RealSparseMatrix(this->ngpts,1));
//  Derivative
    this->sparsedmudX_.push_back(RealSparseMatrix(this->nTCS_*this->nBasis_,this->ngpts));
    this->sparsedmudY_.push_back(RealSparseMatrix(this->nTCS_*this->nBasis_,this->ngpts));
    this->sparsedmudZ_.push_back(RealSparseMatrix(this->nTCS_*this->nBasis_,this->ngpts));
//  dipole like
    RealSparseMatrix *Map        = &this->sparseMap_[iAtm];
    RealSparseMatrix *WeightsMap = &this->sparseWeights_[iAtm];
    RealSparseMatrix *DoRhoMap   = &this->sparseDoRho_[iAtm];
    RealSparseMatrix *MapdX      = &this->sparsedmudX_[iAtm];
    RealSparseMatrix *MapdY      = &this->sparsedmudY_[iAtm];
    RealSparseMatrix *MapdZ      = &this->sparsedmudZ_[iAtm];
    double val;
    double x;
    double y;
    double z;
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
            x = bg::get<0>(pt) ;
            y = bg::get<1>(pt) ;
            z = bg::get<2>(pt) ;
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
// dipole like
               }
             } //loop over basis (within a given ishell)
         } //loop over shells
//AP
//            (*rdotpX) += (*WeightsMap).coeff(ipts,0)*(*MapmudotX).col(ipts)*(*Map).col(ipts).transpose();
            (*rdotpX) +=  (*WeightsMap).coeff(ipts,0)*(
                          ( (*Map).col(ipts)*(*MapdZ).col(ipts).transpose()*y) -
                          ( (*Map).col(ipts)*(*MapdY).col(ipts).transpose()*z) );

            (*rdotpY) +=  (*WeightsMap).coeff(ipts,0)*(
                          ( (*Map).col(ipts)*(*MapdX).col(ipts).transpose()*z) -
                          ( (*Map).col(ipts)*(*MapdZ).col(ipts).transpose()*x) );

            (*rdotpZ) +=  (*WeightsMap).coeff(ipts,0)*(
                          ( (*Map).col(ipts)*(*MapdY).col(ipts).transpose()*x) -
                          ( (*Map).col(ipts)*(*MapdX).col(ipts).transpose()*y) );
//AP
         if (this->screenVxc && Map->col(ipts).norm() > this->epsScreen){
         DoRhoMap->insert(ipts,0) = 2;}
    } //loop over pts
//     cout << "non Zero " << this->sparseMap_[iAtm].nonZeros() << " " << this->ngpts <<endl; 
  } // loop over atoms

         (*rdotpX) *= 4.0*math.pi;
         (*rdotpY) *= 4.0*math.pi;
         (*rdotpZ) *= 4.0*math.pi;
         prettyPrint(cout,(*rdotpX),"Numeric <dipole vel> - x comp");
         prettyPrint(cout,(*rdotpY),"Numeric <dipole vel> - y comp");
         prettyPrint(cout,(*rdotpZ),"Numeric <dipole vel> - z comp");

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
std::array<double,6> SingleSlater<double>::formVExSlater (double rho, double spindensity){
//  epsmu ={energydens_c, potential_c_alpha, potential_c_beta}
  std::array<double,6> epsmu = {0.0,0.0,0.0,0.0,0.0,0.0};
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
//  if(!this->isClosedShell && this->Ref_ != TCS) Cx *= std::pow(2.0,(1.0/3.0));
  
//  double Cx  = (3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));  //TF LDA Prefactor (for Vx)  
  double gx;
//  x = 0.0;
//  Eq A4 Pople, J. Chem. Phys. 5612, (1992) nder =0
  if (nDer == 0 ){
    gx  = -beta*x*x;
    gx /= (1.0 + 6.0 * beta * x * boost::math::asinh(x)) ;
    gx += -Cx;
  }else if (nDer == 1){
//  Eq A8 Pople, J. Chem. Phys. 5612, (1992) nder =1
    gx  = x / (std::sqrt(x*x+1.0));
    gx += -boost::math::asinh(x) ;
    gx *= 6.0 * beta * beta * x * x;
    gx += -2.0*beta*x;
    gx /= (1.0 + 6.0 * beta * x * boost::math::asinh(x))*(1.0 + 6.0 * x * beta* boost::math::asinh(x));
  }
 return gx; 
};  //End Form g function for B88 Exchange

template<>
std::array<double,6> SingleSlater<double>::formVExB88 (double rhoA, double rhoB, 
double gammaAA, double gammaBB){
//  Becke Exchange : Eq 8 From Becke, Phys. Rev A., 3098 (1988)
//  Becke Exchange : Implemented as Eq A3 From Pople, J. Chem. Phys. 5612, (1992)
  std::array<double,6> epsmu = {0.0,0.0,0.0,0.0,0.0,0.0};
  double rhoA1ov3 = std::pow(rhoA,(1.0/3.0));
  double rhoA4ov3 = std::pow(rhoA,(4.0/3.0));
  double rhoB1ov3 = std::pow(rhoB,(1.0/3.0));
  double rhoB4ov3 = std::pow(rhoB,(4.0/3.0));
// Note that in Eq A5 xA   = gammaAA / rhoA4ov3; 
// but actually they meants xA = sqrt(gammaAA) /rhoA4ov3
// and also eq A6 rho is rho^(1/3) instead of rho^(4/3)

  double xA   = std::sqrt(gammaAA) / rhoA4ov3; 
  double xB   ;
  epsmu[0]    = rhoA4ov3*this->gB88(0,xA);
  epsmu[1]    =this->gB88(0,xA) - xA*this->gB88(1,xA);
  epsmu[1]   *= 4.0*rhoA1ov3/3.0;
  epsmu[3]    = 0.5*(this->gB88(1,xA))/(std::sqrt(gammaAA));
//Work  epsmu[3]    = 0.5*(this->gB88(1,xA))/(gammaAA);
//  cout       << "-----------------------------" <<endl;
//  cout       << rhoA   << " " <<gammaAA << " "  << (gammaAA / rhoA4ov3) <<endl;
//  cout       << rhoA4ov3 <<" " << xA << " "<<epsmu[3] <<endl;
//  cout       << "-----------------------------" <<endl;
  if(!this->isClosedShell && this->Ref_ != TCS) {
//Paper  xB   = gammaBB / rhoA4ov3; 
  xB   = std::sqrt(gammaBB) / rhoB4ov3; 
  epsmu[0]   += rhoB4ov3*this->gB88(0,xB);
  epsmu[2]    =this->gB88(0,xB) - xB*this->gB88(1,xB);
  epsmu[2]   *= 4.0*rhoB1ov3/3.0;
  epsmu[4]    = 0.5*(this->gB88(1,xB))/(std::sqrt(gammaBB));
//Work  epsmu[4]    = 0.5*(this->gB88(1,xB))/(gammaBB);
  }
/*
  if (xA < 1.0e-10){
  cout << "---" <<endl;
  cout << xA                    <<  " " << this->gB88(0,xA) << " "<<  rhoA <<endl;
  cout << xA*this->gB88(1,xA)   <<  " " << 4.0*rhoA4ov3/3.0 << " "<<  epsmu[1] <<endl;
  cout << "---" <<endl;
  }
*/
  if (this->isClosedShell && this->Ref_ != TCS) epsmu[0]  *= 2.0; 
  return epsmu;
};  //End Form B88 Exchange

//LYP Correlation from Preuss et Al Chem. Phys. Lett. 157, 200 (1989) LYP1
//                     Implemented as in Pople et Al, J. Chem. Phys.m 98, 5612 (1993) LYP2
template<>
double SingleSlater<double>::omegaLYP(int iop, double c, double d, double rho){ //function used in LYP Correlation (Eq. A26, A32 Ref.LYP2)
  double omega;
  double rho1over3 = std::pow(rho,(1.0/3.0));
    omega  =  std::exp(-c / rho1over3);
    omega /=  1.0 + (d / rho1over3);
    omega *= std::pow(rho,(-11.0/3.0));
    if (iop == 1) {
    omega *= (11.0*rho1over3 -c - ( d / (1 + d/rho1over3) )  );
    omega /= -3.0*rho*rho1over3 ;
    }
    return omega;
};//End omegaLYP

template<>
double SingleSlater<double>::deltaLYP(int iop, double c, double d, double rho){ //function used in LYP Correlation (Eq. A27, A33 Ref.LYP2)
  double delta;
  double rho1over3 = std::pow(rho,(1.0/3.0));
  delta  = d / rho1over3;
  delta /= (1.0 + d / rho1over3);
  delta += c / rho1over3;
  if (iop == 1){
//EqA33
  delta   /= -3.0*rho;
  delta   += (d *d / std::pow(rho,(5.0/3.0))) /( 3.0 * (1.0 + d / rho1over3) * (1.0 + d / rho1over3) );
  }
  return delta; 
};//End deltaLYP

template<>
double SingleSlater<double>::derLYP(int iop, double a, double b, double c, double d, double rhoA, double rhoB){ //function used in LYP Correlation (A23, A24,A25 ref.LYP2)
  double del;
  double delta = this->deltaLYP(0, c, d, (rhoA+rhoB));
  if (iop == 0 ){
  // Eq A23(delLYP/delgammaAA) and eq A25(delLYP/delgammaBB) (remind for A25 call the function inverting rhoA with rhoB)
     del  = 11.0 - delta;
     del *= rhoA/(rhoA+rhoB);
     del += 1.0 - 3.0*delta;
     del *= rhoA*rhoB/9.0;
     del -= rhoB*rhoB;
     del *= -a * b * this->omegaLYP(0, c, d, (rhoA+rhoB));
  } else if (iop ==1) {
//  Eq A24 (delLYP/delgammaAB)
     del  = 47.0 - 7.0 * delta;
     del *= rhoA * rhoB / 9.0;
     del -= 4.0 * (rhoA+rhoB) * (rhoA+rhoB) / 3.0;
     del *= -a * b * this->omegaLYP(0, c, d, (rhoA+rhoB)); 
  }
   return del;
};//End derLYP

template<>
double SingleSlater<double>::der2LYP(int iop, double a, double b, double c, double d, double rhoA, double rhoB, double der1){ //function second der used in LYP Correlation (A29, A30,A31 ref.LYP2)
  double del;
  double rhoT = rhoA+rhoB;
  double delta = this->deltaLYP(0, c, d, rhoT);
  if (iop == 0 ) {
// Eq A29 (del^2 LYP / (del rho_X del gammaXX))
// note pass der1 as dLYP/dgammaAA - debugged 
  del =   - (rhoA*rhoB/9.0) 
             *(  ( 3.0 + rhoA/rhoT ) * deltaLYP(1, c, d, rhoT) + rhoB*(delta - 11.0)/(rhoT*rhoT)  ) 
           + (rhoB/9.0)
             *( 1.0 - 3.0*delta -rhoA*(delta - 11)/rhoT );
  } else if (iop == 1) {
// Eq A30 (del^2 LYP / (del rho_X del gammaXY))
// note pass der1 as dLYP/dgammaAB 
//   debugged
   del  =  - 8.0*rhoT/3.0;
   del +=  - deltaLYP(1, c, d, rhoT)*(7.0*rhoA*rhoB/9.0);
   del +=  rhoB*(47.0-7.0*delta)/9.0;
  } else if (iop == 2 ){
// Eq A31 (del^2 LYP / (del rho_X del gammaYY)) (debugged alread)
// note pass der1 as dLYP/dgammaBB - debugged 
  del =   - (rhoA*rhoB/9.0) 
             *(  ( 3.0 + rhoB/rhoT ) * deltaLYP(1, c, d, rhoT) - rhoB*(delta - 11.0)/(rhoT*rhoT)  ) 
           + (rhoB/9.0)
             *( 1.0 - 3.0*delta -rhoB*(delta - 11)/rhoT )
           - 2.0 * rhoA;
   }
// finishing (note just pass different der1 to get the right term)
   del *= -a * b * this->omegaLYP(0, c, d, rhoT);
   del += der1*this->omegaLYP(1, c, d, rhoT)/this->omegaLYP(0, c, d, rhoT);
   return del;
};//End der2LYP

template<>
std::array<double,6> SingleSlater<double>::formVCLYP (double rhoA, double rhoB, 
     double gammaAA, double gammaBB, double gammaAB){
    double Cfact = 36.4623989787648;  //-(2.0^(11.0/3.0))*(3.0/10.0)*((3.0*pi*pi)^(2.0/3.0))
    double a = 0.04918;
    double b = 0.132;
    double c = 0.2533;
    double d = 0.349;
    double rhoA8over3 = std::pow(rhoA,(8.0/3.0));
    double rhoB8over3 = std::pow(rhoB,(8.0/3.0));
    std::array<double,6> epsmu = {0.0,0.0,0.0,0.0,0.0,0.0};

//  Eq. A23  dLYP/dgammaAA (debugged)
    epsmu[3]  = this->derLYP(0, a, b, c, d, rhoA, rhoB);

//  Eq. A23* dLYP/dgammaBB (debugged)
    epsmu[4]  = this->derLYP(0, a, b, c, d, rhoB, rhoA);

//  Eq. A24  dLYP/dgammaAB (debugged)
    epsmu[5]  = this->derLYP(1, a, b, c, d, rhoA, rhoB);

//  Eq. A28  dLYP/dRhoA (debugged)
    epsmu[1]   = - 4.0 * a * rhoA * rhoB / ( (rhoA+rhoB)*(1.0 + d / std::pow((rhoA+rhoB),(1.0/3.0)) ) );
    epsmu[1]  *= ( (1.0/rhoA)
                  -(1.0/(rhoA+rhoB)) 
                  +((d/3.0) *(std::pow((rhoA+rhoB),(-4.0/3.0)))/ (1.0 + d* std::pow((rhoA+rhoB),(-1.0/3.0))))
                 );
    epsmu[1]  += - Cfact * a * b 
             *( 
              (this->omegaLYP(1, c, d, (rhoA+rhoB)) * rhoA * rhoB * (rhoA8over3 + rhoB8over3))
              +(this->omegaLYP(0, c, d, (rhoA+rhoB)) * rhoB * ( (11.0*rhoA8over3/3.0) + rhoB8over3))
             );
    epsmu[1]  += gammaAA*der2LYP(0, a, b, c, d, rhoA, rhoB, epsmu[3] );
    epsmu[1]  += gammaAB*der2LYP(1, a, b, c, d, rhoA, rhoB, epsmu[5] );
    epsmu[1]  += gammaBB*der2LYP(2, a, b, c, d, rhoA, rhoB, epsmu[4] );
//  Eq. A28* dLYP/dRhoB (debugged)
    epsmu[2]   = - 4.0 * a * rhoA * rhoB / ( (rhoA+rhoB)*(1.0 + d / std::pow((rhoA+rhoB),(1.0/3.0)) ) );
    epsmu[2]  *= ( (1.0/rhoB)
                  -(1.0/(rhoA+rhoB)) 
                  +((d/3.0) *(std::pow((rhoA+rhoB),(-4.0/3.0)))/ (1.0 + d* std::pow((rhoA+rhoB),(-1.0/3.0))))
                 );
    epsmu[2]  += - Cfact * a * b  
                 *( 
                   (this->omegaLYP(1, c, d, (rhoA+rhoB)) * rhoA * rhoB * (rhoA8over3 + rhoB8over3))
                  +(this->omegaLYP(0, c, d, (rhoA+rhoB)) * rhoA * ( (11.0*rhoB8over3/3.0) + rhoA8over3))
                 );
    epsmu[2]  += gammaBB*der2LYP(0, a, b, c, d, rhoB, rhoA, epsmu[4] );
    epsmu[2]  += gammaAB*der2LYP(1, a, b, c, d, rhoB, rhoA, epsmu[5] );
    epsmu[2]  += gammaAA*der2LYP(2, a, b, c, d, rhoB, rhoA, epsmu[3] );
//  Eq. A22  LYP  (debugged)
    epsmu[0]  = - 4.0 * a * rhoA * rhoB / ( (rhoA+rhoB)*(1.0 + d / std::pow((rhoA+rhoB),(1.0/3.0)) ) );
    epsmu[0] += - Cfact * a * b * this->omegaLYP(0, c, d, (rhoA+rhoB)) 
               * rhoA * rhoB * (rhoA8over3 + rhoB8over3); 
    epsmu[0] += epsmu[3] * gammaAA;
    epsmu[0] += epsmu[4] * gammaBB;
    epsmu[0] += epsmu[5] * gammaAB;


//AP Debug 
/*
    rhoA = 0.0;
    rhoB = 0.1;
    double der1 = this->derLYP(1, a, b, c, d, rhoA, rhoB);
    cout << rhoA << " " << this->der2LYP(1,0.04918,0.132,0.2533, 0.349, rhoA, rhoB,der1) << " " << this->der2LYP(1,0.04918,0.132,0.2533, 0.349, rhoB, rhoA,der1) <<  endl;
    rhoA = 0.1;
    der1 = this->derLYP(1, a, b, c, d, rhoA, rhoB);
    cout << rhoA << " " << this->der2LYP(1,0.04918,0.132,0.2533, 0.349, rhoA, rhoB,der1) << " " << this->der2LYP(1,0.04918,0.132,0.2533, 0.349, rhoB, rhoA,der1) <<  endl;
    rhoA = 1.0;
    der1 = this->derLYP(1, a, b, c, d, rhoA, rhoB);
    cout << rhoA << " " << this->der2LYP(1,0.04918,0.132,0.2533, 0.349, rhoA, rhoB,der1) << " " << this->der2LYP(1,0.04918,0.132,0.2533, 0.349, rhoB, rhoA,der1)  << endl;
    rhoA = 2.333;
    der1 = this->derLYP(1, a, b, c, d, rhoA, rhoB);
    cout << rhoA << " " << this->der2LYP(1,0.04918,0.132,0.2533, 0.349, rhoA, rhoB, der1) << " " << this->der2LYP(1,0.04918,0.132,0.2533, 0.349, rhoB, rhoA,der1) << endl;
    CErr();
*/
    return epsmu;
};  //End Form LYP Correlation

template<>
std::array<double,6> SingleSlater<double>::formVCVWN (double rho, double spindensity){
//  epsmu ={energydens_c, potential_c_alpha, potential_c_beta}
  std::array<double,6> epsmu = {0.0,0.0,0.0,0.0,0.0,0.0};
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
   if(this->CorrKernel_ == VWN5) {
     b_f  =  7.06042;  // Caption Table 5
     c_f  = 18.0578;   // Caption Table 5
     x0_f = -0.32500;  // Caption Table 5
     b_p  =  3.72744;  // Caption Table 5
     c_p  = 12.9352;   // Caption Table 5
     x0_p = -0.10498;   // Caption Table 5
     b_a  =  1.13107;   // intext page 1209
     c_a  = 13.0045;    // intext page 1209
     x0_a = -0.00475840; // intext page 1209
   } else if(this->CorrKernel_ == VWN3) {
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
   } else {
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
std::array<double,6> SingleSlater<double>::formVCGGA (double rhoA, double rhoB, 
     double gammaAA, double gammaBB, double gammaAB){

    std::array<double,6> corrEpsMu;
    if (this->CorrKernel_ == VWN3 || this->CorrKernel_ == VWN5) {
       corrEpsMu = this->formVCVWN((rhoA+rhoB),(this->spindens(rhoA,rhoB)));
    } else if (this->CorrKernel_ == LYP) {
       corrEpsMu = this->formVCLYP(rhoA, rhoB, gammaAA, gammaBB, gammaAB);
/*AP
       rhoB = 0.11;
       gammaAA = 0.15;
       gammaBB = 0.25;
       gammaAB = 1.0;
       rhoA    = 1.0;
       corrEpsMu = this->formVCLYP(rhoA, rhoB, gammaAA, gammaBB, gammaAB);
       cout << "rhoA " << rhoA << " rhoB "<< rhoB <<endl;
       cout << "gammaAA " << gammaAA << " gammaBB "<< gammaBB << " gammaAB "<< gammaAB <<endl; 
       cout << "tot " <<corrEpsMu[2] <<endl;
       gammaAB = 0.10;
       rhoA    = 0.10;
       corrEpsMu = this->formVCLYP(rhoA, rhoB, gammaAA, gammaBB, gammaAB);
       cout << "rhoA " << rhoA << " rhoB "<< rhoB <<endl;
       cout << "gammaAA " << gammaAA << " gammaBB "<< gammaBB << " gammaAB "<< gammaAB <<endl; 
       cout << "tot " <<corrEpsMu[2] <<endl;
       gammaAB = 2.33;
       rhoA    = 2.33;
       corrEpsMu = this->formVCLYP(rhoA, rhoB, gammaAA, gammaBB, gammaAB);
       cout << "rhoA " << rhoA << " rhoB "<< rhoB <<endl;
       cout << "gammaAA " << gammaAA << " gammaBB "<< gammaBB << " gammaAB "<< gammaAB <<endl; 
       cout << "tot " <<corrEpsMu[2] <<endl;
       CErr();
*/
    }
    return corrEpsMu;
}; //END GENERIC FORMCOR

template<>
std::array<double,6> SingleSlater<double>::formVC (double rho, double spindensity){

    std::array<double,6> corrEpsMu;
    if (this->CorrKernel_ == VWN3 || this->CorrKernel_ == VWN5) {
       corrEpsMu = this->formVCVWN(rho, spindensity);}
    return corrEpsMu;
}; //END GENERIC FORMCOR

template<>
std::array<double,6> SingleSlater<double>::formVEx (double rho, double spindensity){

    std::array<double,6> exEpsMu;
//  Put sime Slater definition
       exEpsMu = this->formVExSlater(rho, spindensity);
    return exEpsMu;
}; //END GENERIC FORMVex

template<>
std::array<double,6> SingleSlater<double>::formVExGGA (double rhoA, double rhoB, 
    double drhoA, double drhoB){
   
    std::array<double,6> exEpsMu;
    if (this->isGGA) {
    if (this->ExchKernel_ == B88) {
      exEpsMu = this->formVExB88(rhoA, rhoB, drhoA, drhoB);
      }
//       cout << exEpsMu[3] << endl;
    }else {
      exEpsMu = this->formVExSlater((rhoA+rhoB),(this->spindens(rhoA,rhoB)));
    }
    return exEpsMu;
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
//            fBuff += this->onePDMA()->block(bf1_s,bf2_s,n1,n2) * fBuff;
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
    rhor = overlapR_->frobInner(this->onePDMA()->conjugate());
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
     rhor   = overlapR_->frobInner(this->onePDMA()->conjugate());
     rhor_B = overlapR_->frobInner(this->onePDMB()->conjugate());
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
   std::array<double,6>  epsMuCor = {0.0,0.0,0.0,0.0,0.0,0.0}; ///< {energydens_corr, potential_corr_alpha, potential_corr_B}
   std::array<double,6>  epsMuExc = {0.0,0.0,0.0,0.0,0.0,0.0}; ///< {energydend_exchange, potential_exchenge_alpha, potential_exchange_beta}

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
      rhorA = overlapR_.frobInner(this->onePDMA()->conjugate());
      rhorB = overlapR_.frobInner(this->onePDMB()->conjugate());
      rhor = rhorA + rhorB;
    } else {
      rhor = overlapR_.frobInner(this->onePDMA()->conjugate()) ;
    }
//  Handle numerical instability if screening on
    if (this->screenVxc ) {
      if(rhor    <= 0.0 ) {
        if((std::abs(rhor)) <= 1.0e-10) rhor = 0.0;
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
       RealMatrix * VCB, RealMatrix *STmp, RealMatrix *dSTmpX, RealMatrix *dSTmpY, 
       RealMatrix *dSTmpZ) {

   RealSparseMatrix *Map        = &this->sparseMap_[iAtm];
   RealSparseMatrix *WeightsMap = &this->sparseWeights_[iAtm];
   RealSparseMatrix *MapdX      = &this->sparsedmudX_[iAtm];
   RealSparseMatrix *MapdY      = &this->sparsedmudY_[iAtm];
   RealSparseMatrix *MapdZ      = &this->sparsedmudZ_[iAtm];

//   if (this->screenVxc && Map->col(ipts).norm() < this->epsScreen)
//     return;

// Eventually Decleare outside as STmp;
/*
   RealMatrix dSTmpX(this->nBasis_,this->nBasis_);        ///< d(Overlap) at grid point
   dSTmpX.setZero();
   RealMatrix dSTmpY(this->nBasis_,this->nBasis_);        ///< d(Overlap) at grid point
   dSTmpY.setZero();
   RealMatrix dSTmpZ(this->nBasis_,this->nBasis_);        ///< d(Overlap) at grid point
   dSTmpZ.setZero();
*/
   
   double rhor  = 0.0;  // Total density at point
   double rhorA = 0.0;  // alpha density at point
   double rhorB = 0.0;  // beta  density at point
   double gammaAA = 0.0;  // Gradient Invariant: Del alpha dot Del alpha density at point (GGA)
   double gammaBB = 0.0;  // Gradient Invariant: Del beta  dot Del beta  density at point (GGA)
   double gammaAB = 0.0;  // Grdient  Invariant: Del alpha dot Del beta  density at point (GGA)
   std::array<double,3>  drhoA = {0.0,0.0,0.0}; ///< array pf density gradient components
   std::array<double,3>  drhoB = {0.0,0.0,0.0}; ///< array pf density gradient components
   int    nDer   = 0;    // Order of Der
  if (this->isGGA) nDer = 1;
// cout << "nDer" << nDer << endl;
   bool   RHF  = this->Ref_ == RHF;
   bool   doTCS  = this->Ref_ == TCS;
// RealMatrix overlapR_(this->nBasis_,this->nBasis_);        ///< Overlap at grid point
// overlapR_.setZero();
// STmp->setZero();
   std::array<double,6>  epsMuCor = {0.0,0.0,0.0,0.0,0.0,0.0}; ///< {energydens_Cor, potential_Cor_A, potential_Cor_B, potential_Cor_gammaAA_GGA, potential_Cor_gammaBB_GGA, potential_Cor_gammaAB_GGA}
   std::array<double,6>  epsMuExc = {0.0,0.0,0.0,0.0,0.0,0.0}; ///< {energydend_Exc, potential_Exc_A, potential_Exc_B, potential_Exc_gammaAA_GGA,potential_Exc_gammaBB_GGA,potential_Exc_gammaAB_GGA}

//   Build Overlap
    (*STmp) = Map->col(ipts)*Map->col(ipts).transpose();
//  Build Gradient components GGA
    if (nDer == 1){
      (*dSTmpX)   = Map->col(ipts)*MapdX->col(ipts).transpose();
      (*dSTmpY)   = Map->col(ipts)*MapdY->col(ipts).transpose();
      (*dSTmpZ)   = Map->col(ipts)*MapdZ->col(ipts).transpose();
    }
//  Handle the total density at r for RKS or UKS
    if(!this->isClosedShell && this->Ref_ != TCS) {
//    rhorA = overlapR_.frobInner(this->onePDMA()->conjugate());
//    rhorB = overlapR_.frobInner(this->onePDMB()->conjugate());
      rhorA = STmp->frobInner(this->onePDMA()->conjugate());
      rhorB = STmp->frobInner(this->onePDMB()->conjugate());
      rhor = rhorA + rhorB;
      if (nDer == 1 ){
        drhoA[0]   = 2.0*dSTmpX->frobInner(this->onePDMA()->conjugate());
        drhoA[1]   = 2.0*dSTmpY->frobInner(this->onePDMA()->conjugate());
        drhoA[2]   = 2.0*dSTmpZ->frobInner(this->onePDMA()->conjugate());
        gammaAA    = (drhoA[0]*drhoA[0] + drhoA[1]*drhoA[1] + drhoA[2]*drhoA[2]);
        drhoB[0]   = 2.0*dSTmpX->frobInner(this->onePDMB()->conjugate());
        drhoB[1]   = 2.0*dSTmpY->frobInner(this->onePDMB()->conjugate());
        drhoB[2]   = 2.0*dSTmpZ->frobInner(this->onePDMB()->conjugate());
        gammaBB    = (drhoB[0]*drhoB[0] + drhoB[1]*drhoB[1] + drhoB[2]*drhoB[2]);
        gammaAB    = (drhoA[0]*drhoB[0] + drhoA[1]*drhoB[1] + drhoA[2]*drhoB[2]);
      }
    } else {
//      rhor    = STmp->frobInner(this->onePDMA()->conjugate()) ;
      if (this->isGGA) {
        rhorA    = STmp->frobInner(this->onePDMA()->conjugate()/2.0) ;
        rhor = 2.0*rhorA;
        rhorB = rhorA;
      }else{ 
        rhorA    = STmp->frobInner(this->onePDMA()->conjugate()) ;
        rhor = rhorA;
      }
      if (nDer == 1 ){
//        drhoA[0]   = 2.0*dSTmpX.frobInner(this->onePDMA()->conjugate());
        drhoA[0]   = 2.0*(dSTmpX->frobInner(this->onePDMA()->conjugate()/2.0));
        drhoA[1]   = 2.0*(dSTmpY->frobInner(this->onePDMA()->conjugate()/2.0));
        drhoA[2]   = 2.0*(dSTmpZ->frobInner(this->onePDMA()->conjugate()/2.0));
        gammaAA    = (drhoA[0]*drhoA[0] + drhoA[1]*drhoA[1] + drhoA[2]*drhoA[2]);
        gammaBB    = gammaAA;
        gammaAB    = gammaAA;
      }
    }
//    rhor = overlapR_.frobInner(this->onePDMA()->conjugate()) ;
    if (nDer ==  1 ){ 
      (*dSTmpX)   += MapdX->col(ipts)*Map->col(ipts).transpose();
      (*dSTmpY)   += MapdY->col(ipts)*Map->col(ipts).transpose();
      (*dSTmpZ)   += MapdZ->col(ipts)*Map->col(ipts).transpose();
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
        epsMuExc = this->formVExGGA(rhorA,rhorB,gammaAA,gammaBB);
        (*VXA)  +=   ((*WeightsMap).coeff(ipts,0))*(*STmp)*epsMuExc[1];
        if (this->isGGA && this->ExchKernel_ == B88) {
          (*VXA)  += 2.0*drhoA[0]*((*WeightsMap).coeff(ipts,0))*(*dSTmpX)*epsMuExc[3];
          (*VXA)  += 2.0*drhoA[1]*((*WeightsMap).coeff(ipts,0))*(*dSTmpY)*epsMuExc[3];
          (*VXA)  += 2.0*drhoA[2]*((*WeightsMap).coeff(ipts,0))*(*dSTmpZ)*epsMuExc[3];
          energyX += ((*WeightsMap).coeff(ipts,0))*epsMuExc[0];
        } else { 
          energyX += ((*WeightsMap).coeff(ipts,0))*rhor*epsMuExc[0];
        }
        if(!this->isClosedShell && this->Ref_ != TCS){
          (*VXB)  += ((*WeightsMap).coeff(ipts,0))*(*STmp)*epsMuExc[2];
          if (this->isGGA && this->ExchKernel_ == B88) {
            (*VXB)  += 2.0*drhoB[0]*((*WeightsMap).coeff(ipts,0))*(*dSTmpX)*epsMuExc[4];
            (*VXB)  += 2.0*drhoB[1]*((*WeightsMap).coeff(ipts,0))*(*dSTmpY)*epsMuExc[4];
            (*VXB)  += 2.0*drhoB[2]*((*WeightsMap).coeff(ipts,0))*(*dSTmpZ)*epsMuExc[4];
          }
        }
      } //End Exch 

//  Correlation
      if (this->CorrKernel_ != NOCORR) {
         epsMuCor = this->formVCGGA(rhorA,rhorB,gammaAA,gammaBB,gammaAB);
        (*VCA)  += ((*WeightsMap).coeff(ipts,0))*(*STmp)*epsMuCor[1];
        if (this->isGGA && this->CorrKernel_ == LYP) {
//        Open Shell alpha LYP 
          if(!this->isClosedShell && this->Ref_ != TCS) { 
            (*VCA)  += 2.0*drhoA[0]*((*WeightsMap).coeff(ipts,0))*(*dSTmpX)*epsMuCor[3];
            (*VCA)  +=     drhoB[0]*((*WeightsMap).coeff(ipts,0))*(*dSTmpX)*epsMuCor[5];
            (*VCA)  += 2.0*drhoA[1]*((*WeightsMap).coeff(ipts,0))*(*dSTmpY)*epsMuCor[3];
            (*VCA)  +=     drhoB[1]*((*WeightsMap).coeff(ipts,0))*(*dSTmpY)*epsMuCor[5];
            (*VCA)  += 2.0*drhoA[2]*((*WeightsMap).coeff(ipts,0))*(*dSTmpZ)*epsMuCor[3];
            (*VCA)  +=     drhoB[2]*((*WeightsMap).coeff(ipts,0))*(*dSTmpZ)*epsMuCor[5];
          }else{
//        Closed Shell
            (*VCA)  += drhoA[0]*((*WeightsMap).coeff(ipts,0))*(*dSTmpX)*(2.0*epsMuCor[3]+epsMuCor[5]);
            (*VCA)  += drhoA[1]*((*WeightsMap).coeff(ipts,0))*(*dSTmpY)*(2.0*epsMuCor[3]+epsMuCor[5]);
            (*VCA)  += drhoA[2]*((*WeightsMap).coeff(ipts,0))*(*dSTmpZ)*(2.0*epsMuCor[3]+epsMuCor[5]);
          }
          energyC += ((*WeightsMap).coeff(ipts,0))*epsMuCor[0];
        } else { 
          energyC += ((*WeightsMap).coeff(ipts,0))*rhor*epsMuCor[0];
        }
        if(!this->isClosedShell && this->Ref_ != TCS) { 
          (*VCB)  += ((*WeightsMap).coeff(ipts,0))*(*STmp)*epsMuCor[2];
          if (this->isGGA && this->CorrKernel_ == LYP) {
            (*VCB)  += 2.0*drhoB[0]*((*WeightsMap).coeff(ipts,0))*(*dSTmpX)*epsMuCor[4];
            (*VCB)  +=     drhoA[0]*((*WeightsMap).coeff(ipts,0))*(*dSTmpX)*epsMuCor[5];
            (*VCB)  += 2.0*drhoB[1]*((*WeightsMap).coeff(ipts,0))*(*dSTmpY)*epsMuCor[4];
            (*VCB)  +=     drhoA[1]*((*WeightsMap).coeff(ipts,0))*(*dSTmpY)*epsMuCor[5];
            (*VCB)  += 2.0*drhoB[2]*((*WeightsMap).coeff(ipts,0))*(*dSTmpZ)*epsMuCor[4];
            (*VCB)  +=     drhoA[2]*((*WeightsMap).coeff(ipts,0))*(*dSTmpZ)*epsMuCor[5];
          }
        }
//        energyC += ((*WeightsMap).coeff(ipts,0))*rhor*epsMuCor[0];
      }  //End Corr
    } //End VXC
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

    int nPtsPerThread = this->ngpts / nthreads;

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

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    int nAtom   = this->molecule_->nAtoms();                    // Number of Atoms
    this->ngpts = this->nRadDFTGridPts_*this->nAngDFTGridPts_;  // Number of grid point for each center
    int nPtsPerThread = this->ngpts / nthreads;    //  Number of Threads
//  Get the number of tread to dimension this last two
//  std::array<double,1> tmpEnergyEx  = {0.0} ;
//  std::array<double,1> tmpEnergyCor = {0.0} ;
    std::vector<double> tmpEnergyEx(nthreads)  ;
    std::vector<double> tmpEnergyCor(nthreads) ;
    std::vector<int> tmpnpts(nthreads) ;

    int nRHF;
    if(this->isClosedShell || this->Ref_ == TCS) nRHF = 1;
    else    nRHF = 2;
    std::vector<std::vector<RealMatrix>> 
      tmpVX(nRHF,std::vector<RealMatrix>(nthreads,
              RealMatrix::Zero(this->nBasis_,this->nBasis_)
      )
    );
    std::vector<std::vector<RealMatrix>> 
      tmpVC(nRHF,std::vector<RealMatrix>(nthreads,
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
      if (thread_id == (nthreads - 1))
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
          auto GP = Raw3Dg.gridPtCart(ipts);
          auto mapRad_ = this->basisset_->MapGridBasis(GP);
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
      for(auto iThread = 0; iThread < nthreads; iThread++) {
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
        prettyPrint(this->fileio_->out,(*this->vCorB()),"Vc Vc beta");

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
#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    int nAtom   = this->molecule_->nAtoms();                    // Number of Atoms
    this->ngpts = this->nRadDFTGridPts_*this->nAngDFTGridPts_;  // Number of grid point for each center
    int nPtsPerThread = this->ngpts / nthreads;    //  Number of Threads
    std::vector<double> tmpEnergyEx(nthreads)  ;
    std::vector<double> tmpEnergyCor(nthreads) ;
    std::vector<int> tmpnpts(nthreads) ;
    int nRHF;
    int nDerMatrix = 0;                                       // Number of Matrix to be allocated for GGA
    if(this->isClosedShell || this->Ref_ == TCS) nRHF = 1;
    else    nRHF = 2;
    std::vector<std::vector<RealMatrix>> 
      tmpVX(nRHF,std::vector<RealMatrix>(nthreads,
              RealMatrix::Zero(this->nBasis_,this->nBasis_)
      )
    );
    std::vector<std::vector<RealMatrix>> 
      tmpVC(nRHF,std::vector<RealMatrix>(nthreads,
              RealMatrix::Zero(this->nBasis_,this->nBasis_)
      )
    );
    std::vector<RealMatrix> overlapR_(nthreads,RealMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));

    double CxVx  = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));  //TF LDA Prefactor (for Vx)  
    double val ;                                  // to take into account Ang Int
    if (this->isGGA) {
      val = 4.0*math.pi;
      nDerMatrix = nthreads; 
    } else {
      val = 4.0*math.pi*CxVx;                                  // to take into account Ang Int
    }
//GGA
    std::vector<RealMatrix> dXOverlapR_(nDerMatrix,
         RealMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
    std::vector<RealMatrix> dYOverlapR_(nDerMatrix,
         RealMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));
    std::vector<RealMatrix> dZOverlapR_(nDerMatrix,
         RealMatrix(this->nTCS_*this->nBasis_,this->nTCS_*this->nBasis_));

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
      if (thread_id == (nthreads - 1))
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
      for(auto iThread = 0; iThread < nthreads; iThread++) {
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

    
    std::vector<std::chrono::duration<double>> thread_timers(nthreads);
    auto batch_dft = [&] (int thread_id,int iAtm) {
      auto loopSt = nPtsPerThread * thread_id;
      auto loopEn = nPtsPerThread * (thread_id + 1);
      RealSparseMatrix *DoRhoMap   = &this->sparseDoRho_[iAtm];
//      auto start = std::chrono::high_resolution_clock::now();
//    for(int ipts = loopSt; ipts < loopEn; ipts++){
      for(auto ipts = 0; ipts < this->ngpts; ipts++) {
        if(ipts % nthreads != thread_id) continue;
        if( this->screenVxc && ((*DoRhoMap).coeff(ipts,0) < 1) ) continue;
        this->evalVXC_store(iAtm,ipts,tmpEnergyEx[thread_id],tmpEnergyCor[thread_id],
              &tmpVX[0][thread_id],&tmpVX[1][thread_id],&tmpVC[0][thread_id],
              &tmpVC[1][thread_id],&overlapR_[thread_id],&dXOverlapR_[thread_id],
              &dYOverlapR_[thread_id],&dZOverlapR_[thread_id]);
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
      for(auto iThread = 0; iThread < nthreads; iThread++) {
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
      if(!this->isClosedShell && this->Ref_ != TCS) { 
        prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vx beta");
        prettyPrint(this->fileio_->out,(*this->vCorB()),"Vc Vc beta");
        }
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





/// START NEW CODE
template<>
void SingleSlater<dcomplex>::formVXC_new(){;};

template<>
void SingleSlater<double>::formVXC_new(){
//Timing

  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point finish;
  std::chrono::duration<double> duration_formVxc;
  std::chrono::duration<double> T1(0.0);
  std::chrono::duration<double> T2(0.0);
  std::chrono::duration<double> T3(0.0);
  std::chrono::duration<double> T4(0.0);
  std::chrono::duration<double> T5(0.0);
  std::chrono::duration<double> T6(0.0);
  
  std::vector<std::chrono::duration<double>> 
    TF(this->dftFunctionals_.size(),std::chrono::duration<double>(0.0));
/*
  if(this->printLevel_ >= 3) {
    start = std::chrono::high_resolution_clock::now();
  }
*/
  bool isGGA   = true;
  RealMatrix SCRATCH2(this->nBasis_,this->nBasis_);
  VectorXd   SCRATCH1(this->nBasis_);
  RealMatrix SCRATCH2X(this->nBasis_,this->nBasis_);
  RealMatrix SCRATCH2Y(this->nBasis_,this->nBasis_);
  RealMatrix SCRATCH2Z(this->nBasis_,this->nBasis_);
  VectorXd   SCRATCH1X(this->nBasis_);
  VectorXd   SCRATCH1Y(this->nBasis_);
  VectorXd   SCRATCH1Z(this->nBasis_);
  int NDer = 0;
  if(isGGA) NDer = 1; 

/*
  std::array<double,3>  drhoT = {0.0,0.0,0.0}; ///< array TOTAL density gradient components
  std::array<double,3>  drhoS = {0.0,0.0,0.0}; ///< array SPIN  density gradient components
  std::array<double,3>  drhoA = {0.0,0.0,0.0}; ///< array ALPHA  density gradient components
  std::array<double,3>  drhoB = {0.0,0.0,0.0}; ///< array BETA  density gradient components
  RealVecMap GradRhoT(&drhoT[0],3);
  RealVecMap GradRhoS(&drhoS[0],3);
  RealVecMap GradRhoA(&drhoA[0],3);
  RealVecMap GradRhoB(&drhoB[0],3);
*/

  auto valVxc = [&](ChronusQ::IntegrationPoint pt, 
  KernelIntegrand<double> &result) -> void {

    auto Newstart = std::chrono::high_resolution_clock::now();
    SCRATCH1.setZero();
    SCRATCH1X.setZero();
    SCRATCH1Y.setZero();
    SCRATCH1Z.setZero();

    std::array<double,3>  drhoT = {0.0,0.0,0.0}; ///< array TOTAL density gradient components
    std::array<double,3>  drhoS = {0.0,0.0,0.0}; ///< array SPIN  density gradient components
    std::array<double,3>  drhoA = {0.0,0.0,0.0}; ///< array ALPHA  density gradient components
    std::array<double,3>  drhoB = {0.0,0.0,0.0}; ///< array BETA  density gradient components
    RealVecMap GradRhoT(&drhoT[0],3);
    RealVecMap GradRhoS(&drhoS[0],3);
    RealVecMap GradRhoA(&drhoA[0],3);
    RealVecMap GradRhoB(&drhoB[0],3);
    cartGP GP = pt.pt;
    double rhoA;
    double rhoB;
    double gammaAA;
    double gammaBB;
    double gammaAB;
    auto shMap = this->basisset_->MapGridBasis(GP); 

    auto Newend = std::chrono::high_resolution_clock::now();
    T1 += Newend - Newstart;
    if(shMap[0]) {return;}
    Newstart = std::chrono::high_resolution_clock::now();
    for(auto iShell = 0; iShell < this->basisset_->nShell(); iShell++){
      if(!shMap[iShell+1]) continue;


      int shSize= this->basisset_->shells(iShell).size();
//    T6 += Newend - Newstart;
      double * buff = this->basisset_->basisDEval(NDer,this->basisset_->shells(iShell),&pt.pt);

      RealMap bMap(buff,shSize,1);
//    Newstart = std::chrono::high_resolution_clock::now();
      SCRATCH1.block(this->basisset_->mapSh2Bf(iShell),0,shSize,1) = bMap;
      if(NDer>0){
        double * ds1EvalX = buff + shSize;
        double * ds1EvalY = ds1EvalX + shSize;
        double * ds1EvalZ = ds1EvalY + shSize;
        RealMap bMapX(ds1EvalX,shSize,1);
        SCRATCH1X.block(this->basisset_->mapSh2Bf(iShell),0,shSize,1) = bMapX;
        RealMap bMapY(ds1EvalY,shSize,1);
        SCRATCH1Y.block(this->basisset_->mapSh2Bf(iShell),0,shSize,1) = bMapY;
        RealMap bMapZ(ds1EvalZ,shSize,1);
        SCRATCH1Z.block(this->basisset_->mapSh2Bf(iShell),0,shSize,1) = bMapZ;
      }

      //delete [] buff;
    };
    Newend = std::chrono::high_resolution_clock::now();
    T2 += Newend - Newstart;

    if(SCRATCH1.norm() < 1e-8) return;
    Newstart = std::chrono::high_resolution_clock::now();
    SCRATCH2.noalias() = SCRATCH1 * SCRATCH1.transpose();
    double rhoT = this->template computeProperty<double,TOTAL>(SCRATCH2);
    double rhoS = this->template computeProperty<double,MZ>(SCRATCH2);
    rhoA = 0.5 * (rhoT + rhoS);
    rhoB = 0.5 * (rhoT - rhoS);

    if(NDer>0){
      //Closed Shell GGA
      SCRATCH2X.noalias() = SCRATCH1 * SCRATCH1X.transpose();
      SCRATCH2Y.noalias() = SCRATCH1 * SCRATCH1Y.transpose();
      SCRATCH2Z.noalias() = SCRATCH1 * SCRATCH1Z.transpose();
      drhoT[0] = 2.0*this->template computeProperty<double,TOTAL>(SCRATCH2X); 
      drhoT[1] = 2.0*this->template computeProperty<double,TOTAL>(SCRATCH2Y); 
      drhoT[2] = 2.0*this->template computeProperty<double,TOTAL>(SCRATCH2Z); 
      drhoS[0] = 2.0*this->template computeProperty<double,MZ>(SCRATCH2X); 
      drhoS[1] = 2.0*this->template computeProperty<double,MZ>(SCRATCH2Y); 
      drhoS[2] = 2.0*this->template computeProperty<double,MZ>(SCRATCH2Z); 
      GradRhoA.noalias() = 0.5 * (GradRhoT + GradRhoS);
      GradRhoB.noalias() = 0.5 * (GradRhoT - GradRhoS);
      gammaAA = GradRhoA.dot(GradRhoA);           
      gammaBB = GradRhoB.dot(GradRhoB);           
      gammaAB = GradRhoA.dot(GradRhoB);           
      SCRATCH2X.noalias() += SCRATCH1X * SCRATCH1.transpose();
      SCRATCH2Y.noalias() += SCRATCH1Y * SCRATCH1.transpose();
      SCRATCH2Z.noalias() += SCRATCH1Z * SCRATCH1.transpose();
    }
    Newend = std::chrono::high_resolution_clock::now();
    T3 += Newend - Newstart;
    
    if  (rhoT < 1.0e-08) return;
    DFTFunctional::DFTInfo kernelXC;
    for(auto i = 0; i < this->dftFunctionals_.size(); i++){
      Newstart = std::chrono::high_resolution_clock::now();
      if (NDer > 0) {
        kernelXC += this->dftFunctionals_[i]->eval(
            rhoA,rhoB,gammaAA,gammaAB,gammaBB);
      } else {
        kernelXC += this->dftFunctionals_[i]->eval(rhoA, rhoB);
      }
      Newend = std::chrono::high_resolution_clock::now();
      TF[i] += Newend - Newstart;
    } // loop over kernels
    Newstart = std::chrono::high_resolution_clock::now();
    if(NDer > 0) {
      result.VXCA.real() += pt.weight * SCRATCH2X  
        * (  2.0 * GradRhoA[0]*kernelXC.ddgammaAA 
           + GradRhoB[0]* kernelXC.ddgammaAB);  
      result.VXCA.real() += pt.weight * SCRATCH2Y  
        * (  2.0 * GradRhoA[1]*kernelXC.ddgammaAA 
           + GradRhoB[1]* kernelXC.ddgammaAB);  
      result.VXCA.real() += pt.weight * SCRATCH2Z  
        * (  2.0 * GradRhoA[2]*kernelXC.ddgammaAA 
           + GradRhoB[2]* kernelXC.ddgammaAB);  
      result.Energy += pt.weight * kernelXC.eps;
      if(!this->isClosedShell && this->nTCS_ != 2) {
        result.VXCB.real() += pt.weight * SCRATCH2X  
          * (  2.0 * GradRhoB[0]*kernelXC.ddgammaBB 
             + GradRhoA[0]* kernelXC.ddgammaAB);  
        result.VXCB.real() += pt.weight * SCRATCH2Y  
          * (  2.0 * GradRhoB[1]*kernelXC.ddgammaBB 
             + GradRhoA[1]* kernelXC.ddgammaAB);  
        result.VXCB.real() += pt.weight * SCRATCH2Z  
          * (  2.0 * GradRhoB[2]*kernelXC.ddgammaBB 
             + GradRhoA[2]* kernelXC.ddgammaAB);  
      }
    } else {
      result.Energy += pt.weight * (rhoA+rhoB) * kernelXC.eps;
    }
    result.VXCA.real()   += pt.weight * SCRATCH2 * kernelXC.ddrhoA; 
    if(!this->isClosedShell && this->nTCS_ != 2) 
      result.VXCB.real()   += pt.weight * SCRATCH2 * kernelXC.ddrhoB; 
//  if(rhoT < 1e-5) cout << "****** " << rhoT << endl;
    Newend = std::chrono::high_resolution_clock::now();
    T5 += Newend - Newstart;
  };

  ChronusQ::AtomicGrid AGrid(100,302,ChronusQ::GRID_TYPE::EULERMAC,
      ChronusQ::GRID_TYPE::LEBEDEV,ChronusQ::ATOMIC_PARTITION::BECKE,
      this->molecule_->cartArray(),0,1.0,false);
   
  KernelIntegrand<double> res(this->vXA_->cols());
  this->basisset_->radcut(1.0e-08, 50, 1.0e-7);
  this->totalEx    = 0.0;
  this->vXA()->setZero();   // Set to zero every occurence of the SCF
  for(auto iAtm = 0; iAtm < this->molecule_->nAtoms(); iAtm++){
    AGrid.center() = iAtm;
    AGrid.scalingFactor()=0.5 *
      elements[this->molecule_->index(iAtm)].sradius/phys.bohr;
    AGrid.integrate<KernelIntegrand<double>>(valVxc,res);
  };
  (*this->vXA_) = 4*math.pi*res.VXCA;
  this->totalEx = 4*math.pi*res.Energy;
  if(!this->isClosedShell && this->nTCS_ != 2){
    (*this->vXB_) = 4*math.pi*res.VXCB;
  }
  cout << "T1 = " << T1.count() << endl;
  cout << "T2 = " << T2.count() << endl;
  cout << "T3 = " << T3.count() << endl;
  cout << "T4 = " << T4.count() << endl;
  cout << "T5 = " << T5.count() << endl;
  cout << "T6 = " << T6.count() << endl;
  for(auto i : TF) cout << "TF " << i.count() << endl;
  if(this->printLevel_ >= 3) {
    finish = std::chrono::high_resolution_clock::now();
    duration_formVxc = finish - start;
    prettyPrint(this->fileio_->out,(*this->vXA()),"LDA Vxc alpha");
    if(!this->isClosedShell && this->nTCS_ != 2){
      prettyPrint(this->fileio_->out,(*this->vXB()),"LDA Vxc beta");
    }
    this->fileio_->out << "VXC Energy= " <<  this->totalEx << endl, 
    this->fileio_->out << endl << "CPU time for VXC integral:  "
                       << duration_formVxc.count() << " seconds." 
                       << endl;
  }
}

} // Namespace ChronusQ
