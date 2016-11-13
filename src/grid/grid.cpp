/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
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
////   double a0=0.9996651;
////   double val;
// 1s H
////   val = rad*(std::exp(-rad/a0))/ ((std::pow(a0,1.5))*(std::sqrt(math.pi)));
////   return val*val;
         return std::exp(-(std::pow(rad,2.0)));
         } 

  double OneDGrid::f2_val(double elevation,double azimut){
  // Test Function to be integrated by One-dimensional grid over a solid angle
         return (15.0)*(std::pow(sin(elevation),4.0))/(32.0*math.pi);
         } 

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

double * TwoDGrid::BuildDensity(double * Sum,double * Buff,int n1, int n2){
  ConstRealMap fBuff(Buff,n1,n2); 
  RealMap Sout(Sum,n1,n2); 
  Sout = fBuff;  
  return Sum;
} //End

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

void TwoDGrid::transformPts(){
};

double TwoDGrid::voronoii( double mu){
//     Generate Becke Weights according to the partition schems in
//     (J. Chem. Phys., 88 (4),2457 (1988)) using Voronoii Fuzzi Cells
//     Eq. 19
       double p;
       p = (1.5)*(mu) - 0.5*(std::pow(mu,3.0));
       return p;
};

double TwoDGrid::frischpol( double mu, double alpha){
//     Eq. 14 of Chem. Phys. Let. 257 , 213-223 (1996).
//     used in FrischWeight
       double p;
       double mua = mu /alpha;
       if (mu <= -alpha) { 
          p = -1.0;
       }else if (mu >= alpha){
          p = 1.0;
       }else {
          p  = 35.0*(mua) ;
          p += -35.0*(std::pow(mua,3.0));
          p +=  21.0*(std::pow(mua,5.0));
          p +=  -5.0*(std::pow(mua,7.0));
          p /=  16.0;
       }
       return p;
};

double TwoDGrid::step_fun( double mu){
       double p = 0.0;
       if (mu < 0 ){p = 1.0;}
       return p;
};

void TwoDGrid::genGrid(){
};

void TwoDGrid::centerGrid(double cx, double cy, double cz){
//   Given a general 3D grid (r_p,Omega_p), it will
//   generate a centerd 3D grid center on iAtm position
//   by transforming the grid in spherical into cartesian first and
//   adding atom cartian coordinates (cx, cy, cz)
//   we will have a  3D grid (NRad times N Ang)
//   all the grid points cartesian components will be collected in the
//   GridCar_(NtotGrid,3) and the weight will be still raw.
//

     sph3GP ptSPH; ///< Temp spherical point to store the 3D Grid point (not yet translated over atoms centers)
     cartGP ptCarGrid; /// Several Temp Cartesian Points to perform the translation, cell wieghts funtion
     int ipts  = 0;
     int    nRad   = Gr_->npts();
     int    nAng   = Gs_->npts();
     double Cartx = 0.0;
     double Carty = 0.0;
     double Cartz = 0.0;
//   Loop over 3D grid points
     for(int i = 0; i < nRad; i++)
     for(int j = 0; j < nAng; j++){
       ptSPH = this->gridPt(i,j);
       bg::transform(ptSPH,ptCarGrid);
//     Center each 3D over each Atom centers
       Cartx = (bg::get<0>(ptCarGrid) + cx );
       Carty = (bg::get<1>(ptCarGrid) + cy );
       Cartz = (bg::get<2>(ptCarGrid) + cz );
       this->SetgridPtCart(ipts,Cartx, Carty, Cartz);
//     store all in the GridCar_(NtotGrid,3) thanks to this function
       this->weightsGrid_[ipts] = (Gs_->weights()[j])
                                * (Gr_->weights()[i])
                                * (std::pow(Gr_->gridPts()[i],2.0));
       ipts ++;
      }
}; //End

void TwoDGrid::printGrid(){
//  Call to print Grid point to be poletted (Mathematica Format)
    int    npts  = Gr_->npts()*Gs_->npts();
    cout << std::fixed;
    for(int i = 0; i < npts; i++){
    cout << "{ "<<this->GridCarX_[i] << ", " <<this->GridCarY_[i] << ", "<<this->GridCarZ_[i] <<"}, "<< endl;
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
  

void GaussChebyshev1stGrid::atomGrid(double srad){
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
//  this->transformPts(); 
 }

void GaussChebyshev1stGridInf::atomGrid(double sradius) {
//   Given the Slater Radius of the given atom, all the radial point and weight are scaled
//   according to Becke Fuzzi Cell Prescription 
     sradius = 0.5*sradius;
     double toau = (1.0)/phys.bohr;
     double val;
     double dmu;
     double den;
     for(int i = 0; i < this->nPts_; i++){
       dmu = (std::pow((1-this->gridPts_[i]),2.0)) / (2.0*toau*sradius);
       val = toau * sradius *  (1+this->gridPts_[i]) / (1-this->gridPts_[i]);
       den = std::sqrt(1.0-(std::pow(this->gridPts_[i],2.0)));
       this->weights_[i] = this->weights_[i]*den/dmu;
       this->gridPts_[i] = val;
       }
}  //End
  
void GaussChebyshev1stGridInf::transformPts(){
//     if (INuc == 7){
//   Hydrogen
//   double ralpha= 0.529;
//   Nitrogen
//     double ralpha= 0.65/2.0;
//   Oxygen
     double ralpha= 0.60/2.0;
//   Lithium
//   double ralpha= 1.45/2.0;
//     }
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
//      cout << " Transformed Becke " << endl;
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
      }else if(this->nPts_ == 302){
// Spherical Quadrature Formula Exact to Orders 25-29 
// V.I. Lebedev Vol 18, pg.99-107, 1977
// Siberian Mathematical Journal 
// http://dx.doi.org/10.1007/BF00966954
// Lebedev N=302; n=29, eta 0.993  
// page 7;
      A1 = 0.0008545911725128148;
      gen6_A1(0,one,A1);
      A3 = 0.003599119285025571;
      gen8_A3(6,overradthree,A3);
      B1 = 0.003650045807677255;
      l1 = 0.7011766416089545;
      gen24_Bn(14,l1,B1);
      double B2 = 0.003604822601419882;
      double l2 = 0.6566329410219612;
      gen24_Bn(38,l2,B2);
      double B3 = 0.003576729661743367;
      double l3 = 0.4729054132581005;
      gen24_Bn(62,l3,B3);
      double B4 = 0.003449788424305883;
      double l4 = 0.3515640345570105;
      gen24_Bn(86,l4,B4);
      double B5 = 0.003108953122413675;
      double l5 = 0.2219645236294178;
      gen24_Bn(110,l5,B5);
      double B6 = 0.002352101413689164;
      double l6 = 0.09618308522614784;
      gen24_Bn(134,l6,B6);
      C1 = 0.003600820932216460;
      q1 = 0.5718955891878961;
      gen24_Cn(158,q1,C1);
      double C2 = 0.002982344963171804;
      double q2 = 0.2644152887060663;
      gen24_Cn(182,q2,C2);
      double D1 = 0.003571540554273387;
//    Note u1 correspond to r1 in the paper
//    Note r1 correspond to s1 in the paper
      double u1 = 0.2510034751770465;
      double r1 = 0.8000727494073952;
      gen48_Dn(206,u1,r1,D1);
//    Note u2 corresponds to w2 in the paper
//    Note r2 corresponds to s2 in the paper
      double u2 = 0.1233548532583327;
      double r2 = 0.4127724083168531;
      double D2 = 0.003392312205006170;
      gen48_Dn(254,u2,r2,D2);
      }else if(this->nPts_ == 434){
      A1 = 0.0005265897968224436;
      gen6_A1(0,one,A1);
      A2 = 0.002548219972002607;
      gen12_A2(6,overradtwo,A2);
      A3 = 0.002512317418927307;
      gen8_A3(18,overradthree,A3); 
      gen24_Bn(26,0.6909346307509111,0.002530403801186355);
      gen24_Bn(50,0.1774836054609158,0.002014279020918528);
      gen24_Bn(74,0.4914342637784746,0.002501725168402936);
      gen24_Bn(98,0.6456664707424256,0.002513267174597564);
      gen24_Bn(122,0.2861289010307638,0.002302694782227416);
      gen24_Bn(146,0.07568084367178018,0.001462495621594614);
      gen24_Bn(170,0.3927259763368002,0.002445373437312980);
      gen24_Cn(194,0.8818132877794288,0.002417442375638981);
      gen24_Cn(218,0.9776428111182649,0.001910951282179532);
      gen48_Dn(242,0.2054823696403044,0.8689460322872412,0.002416930044324775);
      gen48_Dn(290,0.5905157048925271,0.7999278543857286,0.002512236854563495);
      gen48_Dn(338,0.5550152361076807,0.7717462626915901,0.002496644054553086);
      gen48_Dn(386,0.9371809858553722,0.3344363145343455,0.002236607760437849);
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

void LebedevGrid::atomGrid(double sradius){
}

// Euler Maclaurin one dimensional grid functions
void EulerMaclaurinGrid::genGrid(){
 }                                      

void EulerMaclaurinGrid::transformPts(){

}                                      

// Euler Maclaurin one dimensional grid functions
void EulerMaclaurinGrid::atomGrid(double sradius){
//   Eq. 24 and 25 from "Modern Density Functional Theroy: A tool
//   for chemistry. Theor. Comp. Chem., Vol 2, 169-219 (1995)
//   Note that i = i+1 (do 0 initialization of i)
     double rd3 = sradius*sradius*sradius;
     for(int i = 0; i < this->nPts_; i++) {
     this->gridPts_[i]  = sradius*(i+1.0)*(i+1.0);
     this->gridPts_[i] /= std::pow((this->nPts_+1.0-i-1.0),2.0);
//     this->weights_[i]  = 2.0*(this->nPts_+1.0)*rd3*std::pow((i+1),5.0);
//     this->weights_[i] /= std::pow((this->nPts_+1.0-i-1.0),7.0);
     this->weights_[i]  = 2.0*(this->nPts_+1.0)*sradius*(i+1);
     this->weights_[i] /= std::pow((this->nPts_+1.0-i-1.0),3.0);
      } 
 }                                      


}; // namespace ChronusQ
