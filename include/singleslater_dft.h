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
//       Coordinate of the Grid point in elleptical 
           muij = (boost::geometry::distance(gridPt,ri) - boost::geometry::distance(gridPt,rj))/(*this->molecule_->rIJ())(iAtm,jAtm) ;
//       Do the product over all atoms i .ne. j
           WW *= 0.5*(1.0-this->twodgrid_->voronoii(this->twodgrid_->voronoii(this->twodgrid_->voronoii(muij))));
           }
         }
       return WW;
}


template<typename T>
double SingleSlater<T>::normBeckeW(cartGP gridPt){
//     Normalization of Becke Weights
       int   nAtom = this->molecule_->nAtoms();
       double norm = 0.0;
       for(auto iAtm = 0; iAtm < nAtom; iAtm++){
         norm += this->formBeckeW(gridPt,iAtm);
         }
       return norm ;
};

template<typename T>
void SingleSlater<T>::buildVxc(cartGP gridPt, double weight){
//  Build the Vxc therm at each Grid Points  
   double *pointProd; 
   double rhor;
   std::unique_ptr<RealMatrix>  overlapR_;        ///< Overlap at grid point
   overlapR_ = std::unique_ptr<RealMatrix>(
     new RealMatrix(this->nBasis_,this->nBasis_));
   rhor = 0.0;
   overlapR_->setZero();
// Loops over shells
   for(auto s1=0l, s12=0l; s1 < this->basisset_->nShell(); s1++){
      int bf1_s = this->basisset_->mapSh2Bf(s1);
      int n1    = this->basisset_->shells(s1).size();
      for(int s2=0; s2 <= s1; s2++, s12++){
        int bf2_s   = this->basisset_->mapSh2Bf(s2);
        int n2      = this->basisset_->shells(s2).size();
        auto center = this->basisset_->shells(s1).O;
        double *Buff = new double [n1*n2];
        RealMap fBuff(Buff,n1,n2);
        fBuff.setZero();
        pointProd = this->basisset_->basisProdEval(this->basisset_->shells(s1),this->basisset_->shells(s2),&gridPt);
        Buff = this->twodgrid_->BuildDensity(Buff,pointProd,n1,n2);
        overlapR_->block(bf1_s,bf2_s,n1,n2) = fBuff; 
        }
     }
     (*overlapR_) = overlapR_->selfadjointView<Lower>();;
     rhor = overlapR_->frobInner(this->densityA()->conjugate());
//   Slater LDA (not yet ready. completated outside)       
     (*this->vXCA()) += weight*(*overlapR_)*(std::pow(rhor,(1.0/3.0)));
     if(!this->isClosedShell && this->Ref_ != TCS) (*this->vXCB()) += weight*(*overlapR_)*(std::pow(rhor,(1.0/3.0)));
};
