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
/***********************
 * Form Becke Weights *
 ***********************/
template<typename T>
double SingleSlater<T>::formBeckeW(cartGP gridPt, int iAtm){
//     Generate Becke Weights according to the partition schems in
//     (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii Fuzzi Cells
//     Note these Weights have to be normailzed (see normBeckeW) 
       int nAtom = this->molecule_->nAtoms();
       double WW = 1.0;
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
           muij = (boost::geometry::distance(gridPt,ri) - boost::geometry::distance(gridPt,rj))/(*this->molecule_->rIJ())(iAtm,jAtm) ;
//       Do the product over all atoms i .ne. j (Eq. 13 using definition Eq. 21 with k=3)
           WW *= 0.5*(1.0-this->twodgrid_->voronoii(this->twodgrid_->voronoii(this->twodgrid_->voronoii(muij))));
           }
         }
       return WW;
}


template<typename T>
double SingleSlater<T>::normBeckeW(cartGP gridPt){
//     Normalization of Becke Weights
//     (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii Fuzzi Cells
//     Eq. 22
       int   nAtom = this->molecule_->nAtoms();
       double norm = 0.0;
       for(auto iAtm = 0; iAtm < nAtom; iAtm++){
         norm += this->formBeckeW(gridPt,iAtm);
         }
       return norm ;
};

template<typename T>
double SingleSlater<T>::radsphe(int iAtm, double thr){
       double rad = this->basisset_->radcut(iAtm,1.0e-10);
       return rad ;
};

template<typename T>
void SingleSlater<T>::buildVxc(cartGP gridPt, double weight){
//  Build the Vxc therm at each Grid Points  
   double *pointProd; 
   double rhor = 0.0;
   double rhor_B = 0.0;
// bool   do_corr = false;                   //Activate Correlation
// bool   do_ex = true;                      //Activate Exchange
   std::unique_ptr<RealMatrix>  overlapR_;        ///< Overlap at grid point
   overlapR_ = std::unique_ptr<RealMatrix>(
     new RealMatrix(this->nBasis_,this->nBasis_));
   overlapR_->setZero();
// Loops over shells

///Timing
   
   auto start_dens = std::chrono::high_resolution_clock::now();  
   for(auto s1=0l, s12=0l; s1 < this->basisset_->nShell(); s1++){
      int bf1_s = this->basisset_->mapSh2Bf(s1);
      int n1    = this->basisset_->shells(s1).size();
      for(int s2=0; s2 <= s1; s2++, s12++){
        int bf2_s   = this->basisset_->mapSh2Bf(s2);
        int n2      = this->basisset_->shells(s2).size();
        auto center = this->basisset_->shells(s1).O;
        double *Buff = new double [n1*n2];
   auto start_5 = std::chrono::high_resolution_clock::now();  
        RealMap fBuff(Buff,n1,n2);
        fBuff.setZero();
   auto finish_5 = std::chrono::high_resolution_clock::now();  
   this->duration_5 += finish_5 - start_5;
   auto start_1 = std::chrono::high_resolution_clock::now();  
        pointProd = 
          this->basisset_->basisProdEval(
            this->basisset_->shells(s1),
            this->basisset_->shells(s2),
            &gridPt
          );
   auto finish_1 = std::chrono::high_resolution_clock::now();  
   this->duration_1 += finish_1 - start_1;
   auto start_2 = std::chrono::high_resolution_clock::now();  
        Buff = this->twodgrid_->BuildDensity(Buff,pointProd,n1,n2);
        overlapR_->block(bf1_s,bf2_s,n1,n2) = fBuff; 
   auto finish_2 = std::chrono::high_resolution_clock::now();  
   this->duration_2 += finish_2 - start_2;
        }
     }
    auto finish_dens = std::chrono::high_resolution_clock::now();  
    this->duration_dens += finish_dens - start_dens;
     (*overlapR_) = overlapR_->selfadjointView<Lower>();;
   if(this->isClosedShell && this->Ref_ != TCS) {
    rhor = overlapR_->frobInner(this->densityA()->conjugate());
//
//  LDA Slater Exchange
    if (this->ExchKernel_ != NOEXCH) {
    (*this->vXA())  += weight*(*overlapR_)*(std::pow(rhor,(1.0/3.0)));
    this->totalEx   += weight*(std::pow(rhor,(4.0/3.0)));
    }
//  VWN Correlation
    if (this->CorrKernel_ != NOCORR) {
    if (rhor > 1.0e-20) {                       //this if statement prevent numerical instability with zero guesses
//      this->formVWNPara(2.38732414637843e-04);
//      this->formVWNPara(2.98415518297304e-05);
//      this->formVWNPara(2.3873e-07);
      this->formCor(rhor,0.0);
      (*this->vCorA())    += weight*(*overlapR_)*this->mu_corr;
      this->totalEcorr    += weight*rhor*this->eps_corr;
     }
    }
   }
   if(!this->isClosedShell && this->Ref_ != TCS) {
     rhor   = overlapR_->frobInner(this->densityA()->conjugate());
     rhor_B = overlapR_->frobInner(this->densityB()->conjugate());
    if (this->ExchKernel_ != NOEXCH){
     (*this->vXA()) += weight*(*overlapR_)*(std::pow(rhor,(1.0/3.0)));
     (*this->vXB()) += weight*(*overlapR_)*(std::pow(rhor_B,(1.0/3.0)));
     if((rhor+rhor_B) > 1.0e-20) this->totalEx   += 
                                   weight*(std::pow((rhor+rhor_B),(4.0/3.0)))*
                                     (this->f_spindens(1,this->spindens(rhor,rhor_B)));
      }
    if (this->CorrKernel_ != NOCORR) {
    if (rhor+rhor_B > 1.0e-20) {                       //this if statement prevent numerical instability with zero guesses
     this->formCor((rhor+rhor_B),(this->spindens(rhor,rhor_B)));
      (*this->vCorA())    += weight*(*overlapR_)*this->mu_corr;
      (*this->vCorB())    += weight*(*overlapR_)*this->mu_corr_B;
      this->totalEcorr    += weight*(rhor+rhor_B)*this->eps_corr;
      }
     }
    }
//  }
};  //End


template<typename T>
double SingleSlater<T>::EvepsVWN(int iop, double A_x, double b_x, double c_x, double x0_x, double rho){
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

template<typename T>
double SingleSlater<T>::f_spindens(int iop, double spindensity){
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


template<typename T>
double SingleSlater<T>::df_spindens(double spindensity){
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


template<typename T>
double SingleSlater<T>::df2_spindens(double spindensity){
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


  template<typename T> 
  double SingleSlater<T>::spindens(double rho_A, double rho_B) {
  return (rho_A - rho_B)/ (rho_A + rho_B);
  };  // 
