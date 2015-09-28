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
#include <grid.h>
#include <cerr.h>

namespace ChronusQ{

  
//// TEST FUNCTIONS  ///

  double OneDGrid::f_val(double rad){
  // Test Function to be integrated by One-dimensional grid
  // INT[0,1] r^2 * exp(-r^2);
/*
//         return (std::pow(rad,2.0))*(std::exp(-((std::pow(rad,2.0)))));
//         return (std::pow(rad,2.0))*(std::exp(-((std::pow(rad,2.0)))));
//         return std::exp(-(std::pow(rad,2.0)));
*/
   double a0=0.9996651;
   double val;
// 1s H
   val = rad*(std::exp(-rad/a0))/ ((std::pow(a0,1.5))*(std::sqrt(math.pi)));
   return val*val;
         } 

  double OneDGrid::f2_val(double elevation,double azimut){
  // Test Function to be integrated by One-dimensional grid over a solid angle
         return (15.0)*(std::pow(sin(elevation),4.0))/(32.0*math.pi);
         } 

/*/
  double * TwoDGrid::ftestVal(cartGP *pt){
  // Test Function to be integrated by Two-dimensional grid
     sph3GP  ptSPH;
     double val;
     double * gEval;   
     bg::transform(*pt,ptSPH);
//     val = this->ftest(bg::get<2>(ptSPH),bg::get<1>(ptSPH),bg::get<0>(ptSPH)); 
     cout << "Test1 " <<endl;
     val = this->ftest(bg::get<2>(ptSPH),bg::get<1>(ptSPH),bg::get<0>(ptSPH)); 
//     cout << "Test2 " <<endl;
//     *gEval = val;
//     cout << "Test3 " <<endl;
//     cout << " value " << *gEval << endl;
     return &val;
   } 

  double TwoDGrid::ftest(double rad, double elevation, double azimut) {
//  return std::sin(elevation)*(std::pow(rad,2.0))*(std::exp( -std::pow(rad,2.0)-std::pow((elevation-azimut),2.0)) );
 // return (std::pow(rad,2.0))*(std::exp( -std::pow(rad,2.0)-std::pow((elevation-azimut),2.0)) );
   double a0=0.9996651;
   double val;
//  H 1s
//   val = rad*(std::exp(-rad/a0))/ ((std::pow(a0,1.5))*(std::sqrt(math.pi)));
//  H 2p_0
//     cout << "2p_0 " << endl;
     val = rad*rad*std::sin(elevation)*(std::exp(-rad/(2.0*a0))) / (8.0*(std::pow(a0,2.5))*(std::sqrt(math.pi)));
   return val*val;
//         return std::exp(-(std::pow(rad,2.0)));
//    return (std::pow(rad,2.0))*(std::exp(-((std::pow(rad,2.0)))));
}

  double TwoDGrid::foxy(cartGP pt, cartGP O,double a1, double a2, double a3, double d1, double d2, double d3, double lx, double ly, double lz) {
     double x = bg::get<0>(pt)-bg::get<0>(O);
     double y = bg::get<1>(pt)-bg::get<1>(O);
     double z = bg::get<2>(pt)-bg::get<2>(O);
//     cout << "x "<<x <<" y "<< y << " z " << z <<endl;
     double fun = 0.0;
     double rSq;
     rSq = (x*x + y*y + z*z);
      fun  += d1*std::exp(-a1*rSq);
      fun  += d2*std::exp(-a2*rSq);
      fun  += d3*std::exp(-a3*rSq);
      fun *= std::pow(x,lx);
      fun *= std::pow(y,ly);
      fun *= std::pow(z,lz);
     return fun*fun*rSq;
  }

// END TEST FUNCTIONS
*/


///  ONE GRID GENERAL ///

double OneDGrid::integrate(){
  // Integrate a test function for a one dimensional grid radial part
  // intas2GPt_ is a logical to integrad a 2D angular gris as OneD Grid
   double sum = 0.0;
   std::cout << "Number of One-grid points= "<< this->nPts_  <<std::endl;
   for(int i = 0; i < this->nPts_; i++){
     if(this->intas2GPt_){
       sum += (this->f2_val(bg::get<0>(this->grid2GPts_[i]),bg::get<1>(this->grid2GPts_[i])))*(this->weights_[i]);
       }else{
       sum += (this->f_val(this->gridPts_[i]))*(this->weights_[i]);
       }
     }
   if(this->intas2GPt_){
     return 4.0*math.pi*sum;
     }else{
     cout << "before norm 1d "  << sum <<endl;
     return sum*this->norm_;
//       return sum;
     }
 }

void OneDGrid::printGrid(){
  for(int i = 0; i < this->nPts_; i++){
//  Printing for mathematica
    if(this->intas2GPt_){
      cout << "{ 1.0, "<<bg::get<1>(this->grid2GPts_[i])<<", " <<bg::get<0>(this->grid2GPts_[i]) <<"}, "<< endl;
      }else{
      cout << this->gridPts_[i] << endl;
      }    
    }
 }

///  TWO GRID GENERAL ///

double * TwoDGrid::Buffintegrate(double * Sum,double * Buff,int n1, int n2, double fact){
  //  Integration over batches : Overlap at each point (numerical)
  ConstRealMap fBuff(Buff,n1,n2); 
  RealMap Sout(Sum,n1,n2); 
  Sout += fBuff*fact;  
  return Sum;
} //End

double * TwoDGrid::BuildDensity(double * Sum,double * Buff,int n1, int n2){
  //  Integration over batches : Density at each point (numerical)
  ConstRealMap fBuff(Buff,n1,n2); 
  RealMap Sout(Sum,n1,n2); 
  Sout = fBuff;  
  return Sum;
} //End

/*
  RealMatrix * TwoDGrid::integrateO(){

  // Integrated Over a TWODGrid (Radial x Angular) the basisProdEval (shells(s1),shells(s2))
  // function, returning the pointer of the whole (Nbasis,Nbasis) Matrix
    int    nBase = basisSet_->nBasis();
    int    Ngridr = Gr_->npts();
    int    NLeb   = Gs_->npts();
    double fact;
    sph3GP ptSPH;
    cartGP ptCar;
    RealMatrix *Integral = new RealMatrix(nBase,nBase);  ///< (NBase,Nbase) Integral ove Grid Point
    std::cout <<" --- Numerical Quadrature ---- " <<std::endl;
    std::cout << "Number of Radial-grid points      = "<< Ngridr  <<std::endl;
    std::cout << "Number of Solid Angle-grid points = "<< NLeb  <<std::endl;
    for(auto s1=0l, s12=0l; s1 < basisSet_->nShell(); s1++){
      int bf1_s = basisSet_->mapSh2Bf(s1);
      int n1    = basisSet_->shells(s1).size();
      for(int s2=0; s2 <= s1; s2++, s12++){
        int bf2_s   = basisSet_->mapSh2Bf(s2);
        int n2      = basisSet_->shells(s2).size();
        auto center = basisSet_->shells(s1).O;
        double *pointProd; 
        double *SumInt = new double [n1*n2];
        double val;
        RealMap BlockInt(SumInt,n1,n2);
        BlockInt.setZero();
        for(int i = 0; i < Ngridr; i++)
        for(int j = 0; j < NLeb; j++){
        ptSPH = this->gridPt(i,j);
        bg::transform(ptSPH,ptCar);
        ptCar.set<0>(bg::get<0>(ptCar) + center[0]);
        ptCar.set<1>(bg::get<1>(ptCar) + center[1]);
        ptCar.set<2>(bg::get<2>(ptCar) + center[2]);
        fact = (Gs_->weights()[j])*(Gr_->weights()[i])*(std::pow(Gr_->gridPts()[i],2.0));
        pointProd = basisSet_->basisProdEval(basisSet_->shells(s1),basisSet_->shells(s2),&ptCar);
        SumInt=this->Buffintegrate(SumInt,pointProd,n1,n2,fact);
        }
      Integral->block(bf1_s,bf2_s,n1,n2) = 4*math.pi*BlockInt;
      delete [] SumInt;
      }
    }
    (*Integral) = Integral->selfadjointView<Lower>(); 
//    cout << (*Integral)  << endl;
    return Integral;
 
}  //End
*/

   RealMatrix * TwoDGrid::integrateAtoms(){
// Integrated Over a TWODGrid (Radial x Angular) the basisProdEval (shells(s1),shells(s2))
// function, returning the pointer of the whole (Nbasis,Nbasis) Matrix
   int    nBase = basisSet_->nBasis();
   int    Ngridpts = (Gr_->npts()*Gs_->npts()*this->molecule_->nAtoms());
   cartGP ptCar;
   RealMatrix *Integral = new RealMatrix(nBase,nBase);  ///< (NBase,Nbase) Integral ove Grid Point
   std::cout <<" --- Numerical Quadrature to build overlap ---- " <<std::endl;
   std::cout << "Total Number of grid points = "<< Ngridpts  <<std::endl;
// Loop Over Shells To build the overlap at each grid point
   for(auto s1=0l, s12=0l; s1 < basisSet_->nShell(); s1++){
    int bf1_s = basisSet_->mapSh2Bf(s1);
    int n1    = basisSet_->shells(s1).size();
    for(int s2=0; s2 <= s1; s2++, s12++){
      int bf2_s   = basisSet_->mapSh2Bf(s2);
      int n2      = basisSet_->shells(s2).size();
      auto center = basisSet_->shells(s1).O;
      double *pointProd; 
      double *SumInt = new double [n1*n2];
      double val;
      RealMap BlockInt(SumInt,n1,n2);
      BlockInt.setZero();
//    Loop over grid points
      for(int ipts = 0; ipts < Ngridpts; ipts++){
        ptCar = this->gridPtCart(ipts);
        pointProd = basisSet_->basisProdEval(basisSet_->shells(s1),basisSet_->shells(s2),&ptCar);
        SumInt=this->Buffintegrate(SumInt,pointProd,n1,n2,getweightsGrid(ipts));
        }
      Integral->block(bf1_s,bf2_s,n1,n2) = (4.0*math.pi*BlockInt);
      delete [] SumInt;
      }
    }
    (*Integral) = Integral->selfadjointView<Lower>(); 
    return Integral;
}  //End

double TwoDGrid::integrateDensity(){
//  Build the density at each Grid Points end 
//  return the LDA XC 
   int    nBase = basisSet_->nBasis();
   int    Ngridpts = (Gr_->npts()*Gs_->npts()*this->molecule_->nAtoms());
   double sum = 0.0;
   double Cx = -(3.0/4.0)*(std::pow((3.0/math.pi),(1.0/3.0)));   //TF LDA Prefactor
   double val;
   double *pointProd; 
   double rhor;
   cartGP ptCar;
   RealMatrix *OveratR = new RealMatrix(nBase,nBase);  ///< (NBase,Nbase) Integral ove Grid Point
   std::cout <<" --- Numerical Quadrature for LDA ---- " <<std::endl;
   std::cout << "Number Radial "<< Gr_->npts() << " Number of Angular " << Gs_->npts() <<std::endl;
   std::cout << "Total Number of grid points = "<< Ngridpts  <<std::endl;
// Loop Over Grid Points
   for(int ipts = 0; ipts < Ngridpts; ipts++){
     ptCar = this->gridPtCart(ipts);
     rhor = 0.0;
//   Evaluate the density at each grid points (rhor)
//   Loops over shells
     for(auto s1=0l, s12=0l; s1 < basisSet_->nShell(); s1++){
        int bf1_s = basisSet_->mapSh2Bf(s1);
        int n1    = basisSet_->shells(s1).size();
        for(int s2=0; s2 <= s1; s2++, s12++){
          int bf2_s   = basisSet_->mapSh2Bf(s2);
          int n2      = basisSet_->shells(s2).size();
          auto center = basisSet_->shells(s1).O;
          double *Buff = new double [n1*n2];
          RealMap fBuff(Buff,n1,n2);
          fBuff.setZero();
          pointProd = basisSet_->basisProdEval(basisSet_->shells(s1),basisSet_->shells(s2),&ptCar);
          Buff = this->BuildDensity(Buff,pointProd,n1,n2);
          OveratR->block(bf1_s,bf2_s,n1,n2) = fBuff;
          }
       }
       (*OveratR) = OveratR->selfadjointView<Lower>(); 
//     Ask David what is better
//     rhor = ((*OveratR)*(this->singleSlater_->densityA()->conjugate())).trace();
       rhor = ((*OveratR).frobInner(this->singleSlater_->densityA()->conjugate()));
//     Grid points weights
       val = 4.0*math.pi*getweightsGrid(ipts);
//     Slater LDA        
       sum  +=  val*(std::pow(rhor,(4.0/3.0)));
//     Uncomment to get the Number of Electron
//     sum  +=  val*rhor;
    }
    return Cx*sum;
  
}  //End


  double TwoDGrid::integrate(){
//  OLD
   double sum = 0.0;
   for(int i = 0; i < Gr_->npts(); i++){
    for(int j = 0; j < Gs_->npts(); j++){
//           sum += this->fEVal[i*Gs_->npts() + j]*(Gs_->weights()[j])*(Gr_->weights()[i])*(std::pow(Gr_->gridPts()[i],2.0));
      }
    }
    return 4.0*math.pi*sum;
  }

/* OLD
  double TwoDGrid::integrate(){
//  auto basisset     	= std::unique_ptr<BasisSet>(new BasisSet());
  // Integrate a test function for a one dimensional grid radial part
   double sum;
   cartGP pt(0.01,0.02,0.03);
   sph3GP ptSph;
   bg::transform(pt,ptSph);
   cout << this->basisSet_->shells(2) << endl;
 //  double *f = basisset->basisEval(2,basisset->shells(2).O,&ptSph);
     std::cout << "Number of Radial-grid points= "<< Gr_->npts()  <<std::endl;
     std::cout << "Number of Solid Angle-grid points= "<< Gs_->npts()  <<std::endl;
     for(int i = 0; i < Gr_->npts(); i++){
      for(int j = 0; j < Gs_->npts(); j++){
            
//          ptSph.set<0>(bg::get<0>(Gs_->grid2GPts()[j])); 
//          ptSph.set<1>(bg::get<1>(Gs_->grid2GPts()[j])); 
//          ptSph.set<2>(Gr_->gridPts()[i]); 
//          cout << bg::get<0>(ptSph)-bg::get<0>(Gs_->grid2GPts()[j]);
//          cout << bg::get<1>(ptSph)-bg::get<1>(Gs_->grid2GPts()[j]);
//          cout << bg::get<2>(ptSph)-(Gr_->gridPts()[i]);
//         double  *val = basisset->basisEval(2,basisset->shells(2).O,&ptSph);
       for(auto k = 0; k < 3; k++) { 
//          sum += *(val+k)*(Gs_->weights()[j])*(Gr_->weights()[i]);
//            cout << *(val+k) <<endl;
         }
        }
      }
        return 4.0*(math.pi)*sum*(Gr_->norm());
  }
*/

void TwoDGrid::transformPts(){
};

double TwoDGrid::voronoii( double mu){
       double p;
       p = (1.5)*(mu) - 0.5*(std::pow(mu,3.0));
       return p;
};

double TwoDGrid::step_fun( double mu){
       double p = 0.0;
       if (mu < 0 ){p = 1.0;}
       return p;
};

double TwoDGrid::BeckeW(cartGP GridPt, int iAtm){
//     Generate Becke Weights according to the partition schems in
//     (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii Fuzzi Cells
//     Note these Weights have to be normailzed (see NormBeckeW) 
       int nAtom = this->molecule_->nAtoms();
       double WW = 1.0;
       double muij;   /// elliptical coordinate (ri -rj / Rij)
       cartGP rj;  ///< Cartisian position of Atom j
       cartGP ri;  ///< Cartisian position of Atom i
       ri.set<0>((*this->molecule_->cart())(0,iAtm) );
       ri.set<1>((*this->molecule_->cart())(1,iAtm) );
       ri.set<2>((*this->molecule_->cart())(2,iAtm) );
       for(auto jAtm = 0; jAtm < nAtom; jAtm++){
         if (jAtm != iAtm){
           muij = 0.0;
//       Vector rj (Atoms (j.ne.i) position)
           rj.set<0>((*this->molecule_->cart())(0,jAtm));
           rj.set<1>((*this->molecule_->cart())(1,jAtm));
           rj.set<2>((*this->molecule_->cart())(2,jAtm));
//       Coordinate of the Grid point in elleptical 
           muij = (boost::geometry::distance(GridPt,ri) - boost::geometry::distance(GridPt,rj))/(*this->molecule_->rIJ())(iAtm,jAtm) ;
//       Do the product over all atoms i .ne. j
           WW *= 0.5*(1.0-voronoii(voronoii(voronoii(muij))));
           }
         }
       return WW;
};

double TwoDGrid::NormBeckeW(cartGP GridPt){
//     Normalization of Becke Weights
       int nAtom = this->molecule_->nAtoms();
       double norm = 0.0;
       for(auto iAtm = 0; iAtm < nAtom; iAtm++){
         norm += BeckeW(GridPt,iAtm);
         }
       return norm ;
};

void TwoDGrid::genGrid(){
//   Given a general origin center 3D grid (r_p,Omega_p), it will
//   generate nAtoms 3D grid center on the each atom position
//   by transforming the grid in spherical into cartesian first and
//   adding sequentially each atom cartian coordinates
//   in the end we will have NAtoms time 3D grids (NRad times N Ang)
//   all the grid points cartesian components will be collected in the
//   GridCar_(NtotGrid,3)
//
//   This routine also will build the final weight for each grip points
//   by multiplying the actual gridweight (due to the single center integration scheme)
//   by the Becke (J. Chem. Phys., 88 (4),2457 (1988)) partition scheme for the Atomic
//   weights (based on Voronoii Cells). The final weights will be storered into 
//   this->weightsGrid_

     int nAtom = this->molecule_->nAtoms();
     sph3GP ptSPH; ///< Temp spherical point to store the 3D Grid point (not yet translated over atoms centers)
     cartGP ptCarGrid; /// Several Temp Cartesian Points to perform the translation, cell wieghts funtion
     int ipts  = 0;
     int    Ngridr = Gr_->npts();
     int    NLeb   = Gs_->npts();
     double Cartx = 0.0;
     double Carty = 0.0;
     double Cartz = 0.0;
//   Loop over 3D grid points
     for(int i = 0; i < Ngridr; i++)
     for(int j = 0; j < NLeb; j++){
       double norm = 0;    ///< Voronoi weights normalization factor 
       ptSPH = this->gridPt(i,j);
       bg::transform(ptSPH,ptCarGrid);
      ///Loop over NAtoms
      for(auto iAtm = 0; iAtm < nAtom; iAtm++){
//     Center each 3D over each Atom centers
       Cartx = (bg::get<0>(ptCarGrid) + (*this->molecule_->cart())(0,iAtm) );
       Carty = (bg::get<1>(ptCarGrid) + (*this->molecule_->cart())(1,iAtm) );
       Cartz = (bg::get<2>(ptCarGrid) + (*this->molecule_->cart())(2,iAtm) );
       this->SetgridPtCart(ipts,Cartx, Carty, Cartz);
//     store all in the GridCar_(NtotGrid,3) thanks to this function
//     Start to Evaluate WA (over each center/atom)       
//     Final Weight W_Ang * W_Rad * W_Atom * _ r^2 :
//     W_Atom = BeckeW/NormBeckeW for each Atom i given a grid poin (ipts) 
       this->weightsGrid_[ipts] = (Gs_->weights()[j])
                                * (Gr_->weights()[i])
                                * (std::pow(Gr_->gridPts()[i],2.0))
                                * ((this->BeckeW((this->gridPtCart(ipts)),iAtm))/(this->NormBeckeW(gridPtCart(ipts))) );
       ipts ++;
      }
    }
}; //End

void TwoDGrid::printGrid(){
//  Call to print Grid point to be poletted (Mathematica Format)
    int    Ngridpts = (Gr_->npts()*Gs_->npts()*this->molecule_->nAtoms());
    cartGP ptCar;
    for(int ipts = 0; ipts < Ngridpts; ipts++){
       ptCar = this->gridPtCart(ipts);
       cout << "{" <<bg::get<0>(ptCar) << ", "<<bg::get<1>(ptCar)<<", " <<bg::get<2>(ptCar) <<"}, "<< endl;
      }
}; // End


// Specific Grid Functions Declaration //

// Function Gauss-Chebyshev 1st kind 
void GaussChebyshev1stGrid::genGrid(){
     // Gauss-Chebyshev 1st kind grid
     // Int {a,b} f(x) = Int {-1,1} g(x')/sqrt(1-x'^2) ~ Sum [1, NGrid] weights[i] * g(zeta[i])
     // Where g(zeta[i]) = [(b-a)/2] * sqrt(1-zeta[i]^2) f( [(b-a)/2]*zeta[i] + [(b+a)/2] )
     // and   zeta[i] = gridPts_[i] = cos ( [(2*i-1)*math.pi/2*NGrid])
     // weights[i] = math.pi/NGrid
     // Note I have included the "sqrt(1-zeta[i]^2)" of the transformation in the actual weights[i] .
     // Note in c++ i starts from 0, so has been shifted i = i+1 
  for(int i = 0; i < this->nPts_; i++) {
    this->gridPts_[i] = cos(( (2.0*(i+1)-1.0)/(2.0*this->nPts_))*math.pi);
    this->weights_[i] = (math.pi/this->nPts_)*(sqrt(1-(this->gridPts_[i]*this->gridPts_[i]))) ;
//    the weights are only (math.pi/this->nPts_), the second term is including the integrand transformation factor in Eq.25.4.38 Abramowitz Handbook     
//    std::cout << i <<" " <<this->gridPts_[i] << " "<< this->weights_[i]  << std::endl;
    }
  this->transformPts(); 
}

void GaussChebyshev1stGrid::transformPts(){
  for(int i = 0; i < this->nPts_; i++)
    this->gridPts_[i] = 
     (this->range_[1] - this->range_[0]) / 2.0 * this->gridPts_[i] +
     (this->range_[1] + this->range_[0]) / 2.0;
    } 
  

// Function Gauss-Chebyshev 1st kind (From 0 to Inf)
void GaussChebyshev1stGridInf::genGrid(){
     // Gauss-Chebyshev 1st kind grid
     // Int {a,b} f(x) = Int {-1,1} g(x')/sqrt(1-x'^2) ~ Sum [1, NGrid] weights[i] * g(zeta[i])
     // weights[i] = math.pi/NGrid
     // Note in c++ i starts from 0, so has been shifted i = i+1 
  for(int i = 0; i < this->nPts_; i++) {
    this->gridPts_[i] = cos(( (2.0*(i+1)-1.0)/(2.0*this->nPts_))*math.pi);
    this->weights_[i] = (math.pi/this->nPts_);
    }
  this->transformPts(); 
 }

  
void GaussChebyshev1stGridInf::transformPts(){
//   Hydrogen
//   double ralpha= 0.529;
//   Nitrogen
//     double ralpha= 0.65/2.0;
//   Oxygen
//     double ralpha= 0.60/2.0;
//   Lithium
   double ralpha= 1.45/2.0;
     double toau = (1.0)/phys.bohr;
     double val;
     double dmu;
     double den;
     for(int i = 0; i < this->nPts_; i++){
       dmu = (std::pow((1-this->gridPts_[i]),2.0)) / (2.0*toau*ralpha);
       val = toau * ralpha *  (1+this->gridPts_[i]) / (1-this->gridPts_[i]);
       den = std::sqrt(1.0-(std::pow(this->gridPts_[i],2.0)));
       this->weights_[i] = this->weights_[i]*den/dmu;
       this->gridPts_[i] = val;
       }
      cout << " Transformed Becke " << endl;
 }

 
// Function definition for Lebedev
   void LebedevGrid::genGrid(){
     double one = 1.0;
     double overradtwo = std::sqrt(0.5);
     double overradthree = std::sqrt(1.0/3.0);
     double A1;
     double A2;
     double A3;
     double C1;
     double B1;
     double q1;
     double l1;
// We are using the values in Lebedev75 Zh. vychisl Mat mat Fiz 15, 1, 48-54, 1975
     if(this->nPts_ == 6){
       A1 = 1.0/6.0;
       gen6_A1(0,one,A1);
      }else if(this->nPts_== 14){
       A1= 0.6666666666666667e-1;
       A3= 0.7500000000000000e-1;
       gen6_A1(0,one,A1);
       gen8_A3(6,overradthree,A3);
    }else if(this->nPts_ == 26){
       A1 = 0.4761904761904762e-1;
       A2 = 0.3809523809523810e-1;
       A3 = 0.3214285714285714e-1; 
       gen6_A1(0,one,A1);
       gen12_A2(6,overradtwo,A2);
       gen8_A3(18,overradthree,A3); 
    }else if(this->nPts_ == 38){

// Values of the nodes and weights of ninth to seventeenth order gauss-markov quadrature formulae invariant under the octahedron group with inversion
// V.I. Lebedev Vol 15, pg.44-51, 1975
// USSR Computational Mathematics and Mathematical Physics 
// http://dx.doi.org//10.1016/0041-5553(75)90133-0

// Lebedev N=38; n=9, eta 0.877 Lebedev 1976 ZVMMF_15_48 table 9.1;
      A1 = double(1)/double(105);
      A3 = double(9)/double(280);
      q1 = 0.4597008433809831;
      C1 = double(1)/double(35);
      gen6_A1(0,one,A1);
      gen8_A3(6,overradthree,A3);
      gen24_Cn(14,q1,C1);
    }else if(this->nPts_ == 50){
// Lebedev N=50; n= 11, eta 0.960 Lebedev 1976 ZVMMF_15_48 table 11.1;
      A1 = double(4)/double(315);
      A2 = double(64)/double(2835);
      A3 = double(27)/double(1280);
      B1 = double(std::pow(11.0,4.0))/double(725760);
      l1 = 0.301511344578;
      gen6_A1(0,one,A1);
      gen12_A2(6,overradtwo,A2);
      gen8_A3(18,overradthree,A3);
      gen24_Bn(26,l1,B1);
    }else if(this->nPts_ == 110){
// Lebedev N=110; n=17, eta 0.982 Lebedev 1976 ZVMMF_15_48 table 17.1;
// Note the commented ones are from the original paper ...
      A1 = 0.00382827049494;
//      A3 = 0.00988550016044;
      A3 = 0.009793737512487512;
//      B1 = 0.00844068048232;
      B1 = 0.008211737283191111;
      l1 = 0.185115635345;
      gen6_A1(0,one,A1); 
      gen8_A3(6,overradthree,A3); 
      gen24_Bn(14,l1,B1);
      double B2 = 0.00959547133607;
//      double l2 = 0.383386152638; 
      double l2 = 0.3956894730559419; 
      gen24_Bn(38,l2,B2);
      double B3 = 0.00994281489118;
      double l3 = 0.690421048382;
      gen24_Bn(62,l3,B3);
//      C1 = 4.0 * (std::pow(17.0,3.0)) / 2027025.0;
//      in the paper see few lines before the table (the same for table 17.1,17.2)
      C1 = 0.00969499636166;
      q1 = 0.478369028812;
      gen24_Cn(86,q1,C1);
      }else if(this->nPts_ == 194){
// Quadratures on a sphere 
// V.I. Lebedev Vol 16, pg.10-24, 1976
// USSR Computational Mathematics and Mathematical Physics 
// http://dx.doi.org/10.1016/0041-5553(76)90100-2

// Lebedev N=194; n=23, eta 0.990  
// page 18;

      A1 = 0.001782340447244611; 
      gen6_A1(0,one,A1); 
      A2 = 0.005716905949977102;
      gen12_A2(6,overradtwo,A2);
      A3 = 0.005573383178848738;
      gen8_A3(18,overradthree,A3); 
      l1 = 0.4446933178717437;
      B1 = 0.005518771467273614;
      gen24_Bn(26,l1,B1);
      double l2 = 0.2892465627575439;
      double B2 = 0.005158237711805383;
      gen24_Bn(50,l2,B2);
      double l3 = 0.6712973442695226;
      double B3 = 0.005608704082587997;
      gen24_Bn(74,l3,B3);
      double l4 = 0.1299335447650067;
      double B4 = 0.004106777028169394;
      gen24_Bn(98,l4,B4);
      q1 = 0.3457702197611283;
      C1 = 0.005051846064614808;
      gen24_Cn(122,q1,C1);
      double u1 = 0.1590417105383530;
      double r1 = 0.8360360154824589;
      double D1 = 0.005530248916233094;
      gen48_Dn(146, u1, r1, D1);
      }else{
      CErr("Number of points not available in Lebedev quadrature");
      }
  };


    void LebedevGrid::transformPts(){
};

void LebedevGrid::gen6_A1(int num, double a, double v){
//  v is A1 in Lebedev Tables. 
    cartGP tmpCart;
    tmpCart.set<0>(a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(0.0);
    this->weights_[num+0] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+0]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(0.0);
    this->weights_[num+1] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+1]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(a);
    tmpCart.set<2>(0.0);
    this->weights_[num+2] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+2]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(0.0);
    this->weights_[num+3] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+3]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(a);
    this->weights_[num+4] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+4]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(-a);
    this->weights_[num+5] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+5]));
}

void LebedevGrid::gen12_A2(int num, double a, double v){
//  v is A2 in Lebedev Tables. 
    
    cartGP tmpCart;

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(a);
    tmpCart.set<2>(a);
    this->weights_[num+0] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+0]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(a);
    this->weights_[num+1] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+1]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(a);
    tmpCart.set<2>(-a);
    this->weights_[num+2] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+2]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(-a);
    this->weights_[num+3] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+3]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(a);
    this->weights_[num+4] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+4]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(a);
    this->weights_[num+5] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+5]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(-a);
    this->weights_[num+6] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+6]));
 
    tmpCart.set<0>(-a);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>( -a);
    this->weights_[num+7] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+7]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(a);
    tmpCart.set<2>(0.0);
    this->weights_[num+8] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+8]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>( a);
    tmpCart.set<2>(0.0);
    this->weights_[num+9] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+9]));

    tmpCart.set<0>( a);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(0.0);
    this->weights_[num+10] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+10]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(0.0);
    this->weights_[num+11] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+11]));

}

void LebedevGrid::gen8_A3(int num, double a, double v){
//  v is A3 in Lebedev Tables. 
    
    cartGP tmpCart;
    tmpCart.set<0>(a);
    tmpCart.set<1>(a);
    tmpCart.set<2>(a);
    this->weights_[num+0] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+0]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(a);
    tmpCart.set<2>(a);
    this->weights_[num+1] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+1]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(a);
    this->weights_[num+2] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+2]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(a);
    this->weights_[num+3] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+3]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(a);
    tmpCart.set<2>(-a);
    this->weights_[num+4] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+4]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(a);
    tmpCart.set<2>(-a);
    this->weights_[num+5] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+5]));

    tmpCart.set<0>(a);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(-a);
    this->weights_[num+6] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+6]));

    tmpCart.set<0>(-a);
    tmpCart.set<1>(-a);
    tmpCart.set<2>(-a);
    this->weights_[num+7] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+7]));
}

void LebedevGrid::gen24_Cn(int num, double q, double v){
//  v is Cn in Lebedev
    double p= std::sqrt(1.0 - q * q);
    cartGP tmpCart;

    tmpCart.set<0>(q);
    tmpCart.set<1>(p);
    tmpCart.set<2>(0.0);
    this->weights_[num+0] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+0]));

    tmpCart.set<0>(-q);
    tmpCart.set<1>(p);
    tmpCart.set<2>(0.0);
    this->weights_[num+1] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+1]));
    
    tmpCart.set<0>(q);
    tmpCart.set<1>(-p);
    tmpCart.set<2>(0.0);
    this->weights_[num+2] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+2]));

    tmpCart.set<0>(-q);
    tmpCart.set<1>(-p);
    tmpCart.set<2>(0.0);
    this->weights_[num+3] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+3]));
 
    tmpCart.set<0>(p);
    tmpCart.set<1>(q);
    tmpCart.set<2>(0.0);
    this->weights_[num+4] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+4]));

    tmpCart.set<0>(-p);
    tmpCart.set<1>(q);
    tmpCart.set<2>(0.0);
    this->weights_[num+5] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+5]));

    tmpCart.set<0>(p);
    tmpCart.set<1>(-q);
    tmpCart.set<2>(0.0);
    this->weights_[num+6] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+6]));

    tmpCart.set<0>(-p);
    tmpCart.set<1>(-q);
    tmpCart.set<2>(0.0);
    this->weights_[num+7] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+7]));

    tmpCart.set<0>(q);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(p);
    this->weights_[num+8] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+8]));

    tmpCart.set<0>(-q);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(p);
    this->weights_[num+9] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+9]));
//x[10]
    tmpCart.set<0>(q);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(-p);
    this->weights_[num+10] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+10]));

    tmpCart.set<0>(-q);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(-p);
    this->weights_[num+11] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+11]));
//check
    tmpCart.set<0>(p);
    tmpCart.set<1>(0);
    tmpCart.set<2>(q);
    this->weights_[num+12] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+12]));

    tmpCart.set<0>(-p);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(q);
    this->weights_[num+13] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+13]));

    tmpCart.set<0>(p);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(-q);
    this->weights_[num+14] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+14]));

    tmpCart.set<0>(-p);
    tmpCart.set<1>(0.0);
    tmpCart.set<2>(-q);
    this->weights_[num+15] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+15]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(q);
    tmpCart.set<2>(p);
    this->weights_[num+16] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+16]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-q);
    tmpCart.set<2>(p);
    this->weights_[num+17] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+17]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(q);
    tmpCart.set<2>(-p);
    this->weights_[num+18] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+18]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-q);
    tmpCart.set<2>(-p);
    this->weights_[num+19] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+19]));
//20
    tmpCart.set<0>(0.0);
    tmpCart.set<1>(p);
    tmpCart.set<2>(q);
    this->weights_[num+20] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+20]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-p);
    tmpCart.set<2>(q);
    this->weights_[num+21] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+21]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(p);
    tmpCart.set<2>(-q);
    this->weights_[num+22] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+22]));

    tmpCart.set<0>(0.0);
    tmpCart.set<1>(-p);
    tmpCart.set<2>(-q);
    this->weights_[num+23] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+23]));

}

void LebedevGrid::gen24_Bn(int num, double l, double v){
    double m;
    cartGP tmpCart;

    m = std::sqrt ( 1.0 - 2.0 * l * l ); 
    tmpCart.set<0>(l);
    tmpCart.set<1>(l);
    tmpCart.set<2>(m);
    this->weights_[num+0] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+0]));

    tmpCart.set<0>(-l);
    tmpCart.set<1>( l);
    tmpCart.set<2>( m);
    this->weights_[num+1] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+1]));
    
    tmpCart.set<0>(l);
    tmpCart.set<1>(-l);
    tmpCart.set<2>( m);
    this->weights_[num+2] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+2]));

    tmpCart.set<0>(-l);
    tmpCart.set<1>(-l);
    tmpCart.set<2>(m);
    this->weights_[num+3] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+3]));
 
    tmpCart.set<0>(l);
    tmpCart.set<1>(l);
    tmpCart.set<2>(-m);
    this->weights_[num+4] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+4]));

    tmpCart.set<0>(-l);
    tmpCart.set<1>( l);
    tmpCart.set<2>(-m);
    this->weights_[num+5] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+5]));

    tmpCart.set<0>(l);
    tmpCart.set<1>(-l);
    tmpCart.set<2>(-m);
    this->weights_[num+6] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+6]));

    tmpCart.set<0>(-l);
    tmpCart.set<1>(-l);
    tmpCart.set<2>(-m);
    this->weights_[num+7] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+7]));

    tmpCart.set<0>(l);
    tmpCart.set<1>(m);
    tmpCart.set<2>(l);
    this->weights_[num+8] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+8]));

    tmpCart.set<0>(-l);
    tmpCart.set<1>( m);
    tmpCart.set<2>( l);
    this->weights_[num+9] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+9]));
//x[10]
    tmpCart.set<0>( l);
    tmpCart.set<1>(-m);
    tmpCart.set<2>( l);
    this->weights_[num+10] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+10]));

    tmpCart.set<0>(-l);
    tmpCart.set<1>(-m);
    tmpCart.set<2>( l);
    this->weights_[num+11] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+11]));
//check
    tmpCart.set<0>(l);
    tmpCart.set<1>(m);
    tmpCart.set<2>(-l);
    this->weights_[num+12] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+12]));

    tmpCart.set<0>(-l);
    tmpCart.set<1>(m);
    tmpCart.set<2>(-l);
    this->weights_[num+13] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+13]));

    tmpCart.set<0>(l);
    tmpCart.set<1>(-m);
    tmpCart.set<2>(-l);
    this->weights_[num+14] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+14]));

    tmpCart.set<0>(-l);
    tmpCart.set<1>(-m);
    tmpCart.set<2>(-l);
    this->weights_[num+15] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+15]));

    tmpCart.set<0>(m);
    tmpCart.set<1>(l);
    tmpCart.set<2>(l);
    this->weights_[num+16] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+16]));

    tmpCart.set<0>(-m);
    tmpCart.set<1>( l);
    tmpCart.set<2>( l);
    this->weights_[num+17] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+17]));

    tmpCart.set<0>( m);
    tmpCart.set<1>(-l);
    tmpCart.set<2>( l);
    this->weights_[num+18] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+18]));

    tmpCart.set<0>(-m);
    tmpCart.set<1>(-l);
    tmpCart.set<2>(l);
    this->weights_[num+19] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+19]));
//20
    tmpCart.set<0>(m);
    tmpCart.set<1>(l);
    tmpCart.set<2>(-l);
    this->weights_[num+20] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+20]));

    tmpCart.set<0>(-m);
    tmpCart.set<1>( l);
    tmpCart.set<2>(-l);
    this->weights_[num+21] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+21]));

    tmpCart.set<0>(m);
    tmpCart.set<1>(-l);
    tmpCart.set<2>(-l);
    this->weights_[num+22] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+22]));

    tmpCart.set<0>(-m);
    tmpCart.set<1>(-l);
    tmpCart.set<2>(-l);
    this->weights_[num+23] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+23]));

}

void LebedevGrid::gen48_Dn(int num, double u, double r, double v){
//  v is Dn in Lebedev
//  Because We are on unit sphere the 3rd coefficient w is:
    double w = std::sqrt(1.0 - (u * u) - (r * r));
    cartGP tmpCart;

    tmpCart.set<0>(r);
    tmpCart.set<1>(u);
    tmpCart.set<2>(w);
    this->weights_[num+0] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+0]));

    tmpCart.set<0>(-r);
    tmpCart.set<1>(u);
    tmpCart.set<2>(w);
    this->weights_[num+1] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+1]));
    
    tmpCart.set<0>(r);
    tmpCart.set<1>(-u);
    tmpCart.set<2>(w);
    this->weights_[num+2] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+2]));

    tmpCart.set<0>(-r);
    tmpCart.set<1>(-u);
    tmpCart.set<2>(w);
    this->weights_[num+3] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+3]));
 
    tmpCart.set<0>(u);
    tmpCart.set<1>(r);
    tmpCart.set<2>(w);
    this->weights_[num+4] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+4]));

    tmpCart.set<0>(-u);
    tmpCart.set<1>(r);
    tmpCart.set<2>(w);
    this->weights_[num+5] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+5]));

    tmpCart.set<0>(u);
    tmpCart.set<1>(-r);
    tmpCart.set<2>(w);
    this->weights_[num+6] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+6]));

    tmpCart.set<0>(-u);
    tmpCart.set<1>(-r);
    tmpCart.set<2>(w);
    this->weights_[num+7] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+7]));

    tmpCart.set<0>(r);
    tmpCart.set<1>(w);
    tmpCart.set<2>(u);
    this->weights_[num+8] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+8]));

    tmpCart.set<0>(-r);
    tmpCart.set<1>(w);
    tmpCart.set<2>(u);
    this->weights_[num+9] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+9]));
//x[10]
    tmpCart.set<0>(r);
    tmpCart.set<1>(w);
    tmpCart.set<2>(-u);
    this->weights_[num+10] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+10]));

    tmpCart.set<0>(-r);
    tmpCart.set<1>(w);
    tmpCart.set<2>(-u);
    this->weights_[num+11] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+11]));
//check
    tmpCart.set<0>(u);
    tmpCart.set<1>(w);
    tmpCart.set<2>(r);
    this->weights_[num+12] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+12]));

    tmpCart.set<0>(-u);
    tmpCart.set<1>(w);
    tmpCart.set<2>(r);
    this->weights_[num+13] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+13]));

    tmpCart.set<0>(u);
    tmpCart.set<1>(w);
    tmpCart.set<2>(-r);
    this->weights_[num+14] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+14]));

    tmpCart.set<0>(-u);
    tmpCart.set<1>(w);
    tmpCart.set<2>(-r);
    this->weights_[num+15] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+15]));

    tmpCart.set<0>(w);
    tmpCart.set<1>(r);
    tmpCart.set<2>(u);
    this->weights_[num+16] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+16]));

    tmpCart.set<0>(w);
    tmpCart.set<1>(-r);
    tmpCart.set<2>(u);
    this->weights_[num+17] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+17]));

    tmpCart.set<0>(w);
    tmpCart.set<1>(r);
    tmpCart.set<2>(-u);
    this->weights_[num+18] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+18]));

    tmpCart.set<0>(w);
    tmpCart.set<1>(-r);
    tmpCart.set<2>(-u);
    this->weights_[num+19] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+19]));
//20
    tmpCart.set<0>(w);
    tmpCart.set<1>(u);
    tmpCart.set<2>(r);
    this->weights_[num+20] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+20]));

    tmpCart.set<0>(w);
    tmpCart.set<1>(-u);
    tmpCart.set<2>(r);
    this->weights_[num+21] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+21]));

    tmpCart.set<0>(w);
    tmpCart.set<1>(u);
    tmpCart.set<2>(-r);
    this->weights_[num+22] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+22]));

    tmpCart.set<0>(w);
    tmpCart.set<1>(-u);
    tmpCart.set<2>(-r);
    this->weights_[num+23] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+23]));


    tmpCart.set<0>(r);
    tmpCart.set<1>(u);
    tmpCart.set<2>(-w);
    this->weights_[num+24] = v;
    bg::transform(tmpCart,(this->grid2GPts_[num+24]));

    tmpCart.set<0>(-r);
    tmpCart.set<1>(u);
    tmpCart.set<2>(-w);
    this->weights_[num+25] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+25]));
    
    tmpCart.set<0>(r);
    tmpCart.set<1>(-u);
    tmpCart.set<2>(-w);
    this->weights_[num+26] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+26]));

    tmpCart.set<0>(-r);
    tmpCart.set<1>(-u);
    tmpCart.set<2>(-w);
    this->weights_[num+27] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+27]));
 
    tmpCart.set<0>(u);
    tmpCart.set<1>(r);
    tmpCart.set<2>(-w);
    this->weights_[num+28] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+28]));

    tmpCart.set<0>(-u);
    tmpCart.set<1>(r);
    tmpCart.set<2>(-w);
    this->weights_[num+29] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+29]));

    tmpCart.set<0>(u);
    tmpCart.set<1>(-r);
    tmpCart.set<2>(-w);
    this->weights_[num+30] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+30]));

    tmpCart.set<0>(-u);
    tmpCart.set<1>(-r);
    tmpCart.set<2>(-w);
    this->weights_[num+31] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+31]));

    tmpCart.set<0>(r);
    tmpCart.set<1>(-w);
    tmpCart.set<2>(u);
    this->weights_[num+32] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+32]));

    tmpCart.set<0>(-r);
    tmpCart.set<1>(-w);
    tmpCart.set<2>(u);
    this->weights_[num+33] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+33]));
//x[10]
    tmpCart.set<0>(r);
    tmpCart.set<1>(-w);
    tmpCart.set<2>(-u);
    this->weights_[num+34] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+34]));

    tmpCart.set<0>(-r);
    tmpCart.set<1>(-w);
    tmpCart.set<2>(-u);
    this->weights_[num+35] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+35]));
//check
    tmpCart.set<0>(u);
    tmpCart.set<1>(-w);
    tmpCart.set<2>(r);
    this->weights_[num+36] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+36]));

    tmpCart.set<0>(-u);
    tmpCart.set<1>(-w);
    tmpCart.set<2>(r);
    this->weights_[num+37] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+37]));

    tmpCart.set<0>(u);
    tmpCart.set<1>(-w);
    tmpCart.set<2>(-r);
    this->weights_[num+38] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+38]));

    tmpCart.set<0>(-u);
    tmpCart.set<1>(-w);
    tmpCart.set<2>(-r);
    this->weights_[num+39] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+39]));

    tmpCart.set<0>(-w);
    tmpCart.set<1>(r);
    tmpCart.set<2>(u);
    this->weights_[num+40] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+40]));

    tmpCart.set<0>(-w);
    tmpCart.set<1>(-r);
    tmpCart.set<2>(u);
    this->weights_[num+41] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+41]));

    tmpCart.set<0>(-w);
    tmpCart.set<1>(r);
    tmpCart.set<2>(-u);
    this->weights_[num+42] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+42]));

    tmpCart.set<0>(-w);
    tmpCart.set<1>(-r);
    tmpCart.set<2>(-u);
    this->weights_[num+43] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+43]));
//20
    tmpCart.set<0>(-w);
    tmpCart.set<1>(u);
    tmpCart.set<2>(r);
    this->weights_[num+44] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+44]));

    tmpCart.set<0>(-w);
    tmpCart.set<1>(-u);
    tmpCart.set<2>(r);
    this->weights_[num+45] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+45]));

    tmpCart.set<0>(-w);
    tmpCart.set<1>(u);
    tmpCart.set<2>(-r);
    this->weights_[num+46] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+46]));

    tmpCart.set<0>(-w);
    tmpCart.set<1>(-u);
    tmpCart.set<2>(-r);
    this->weights_[num+47] =  v;
    bg::transform(tmpCart,(this->grid2GPts_[num+47]));


}

}; // namespace ChronusQ
