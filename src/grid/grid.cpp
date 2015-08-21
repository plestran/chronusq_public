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


//// TEST FUNCTIONS
  double OneDGrid::f_val(double rad){
  // Test Function to be integrated by One-dimensional grid
  // INT[0,1] r^2 * exp(-r^2);
//         return (std::pow(rad,2.0))*(std::exp(-((std::pow(rad,2.0)))));
//         return (std::pow(rad,2.0))*(std::exp(-((std::pow(rad,2.0)))));
//         return std::exp(-(std::pow(rad,2.0)));

   double a0=0.9996651;
   double val;
// 1s H
   val = rad*(std::exp(-rad/a0))/ ((std::pow(a0,1.5))*(std::sqrt(math.pi)));


   return val*val;
         //return std::pow(math.pi,-1.5)*(std::pow(rad,2.0))*(std::exp(-((std::pow(rad,2.0)))));
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

double OneDGrid::integrate(){
  // Integrate a test function for a one dimensional grid radial part
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

  double TwoDGrid::integrate(){
  // Integrate a test function for a one dimensional grid radial part
   double sum = 0.0;
//     std::cout << "Number of Radial-grid points= "<< Gr_->npts()  <<std::endl;
//     std::cout << "Number of Solid Angle-grid points= "<< Gs_->npts()  <<std::endl;
     for(int i = 0; i < Gr_->npts(); i++){
      for(int j = 0; j < Gs_->npts(); j++){
           sum += this->fEVal[i*Gs_->npts() + j]*(Gs_->weights()[j])*(Gr_->weights()[i])*(std::pow(Gr_->gridPts()[i],2.0));
        }
      }
//           cout << "before norm "<< sum <<endl;
//        return 4.0*(math.pi)*sum*(Gr_->norm());
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

void TwoDGrid::printGrid(){
  for(int i = 0; i < Gr_->npts(); i++){
    for(int j = 0; j < Gs_->npts(); j++){
      cout << "{" << Gr_->gridPts()[i] << ", "<<bg::get<1>(Gs_->grid2GPts()[j])<<", " <<bg::get<0>(Gs_->grid2GPts()[j]) <<"}, "<< endl;
       }
    }
  };

void TwoDGrid::genGrid(){
  };

// Specific Grid Functions Declaration

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
   double ralpha= 0.529;
//   Oxygen
//     double ralpha= 0.60/2.0;
//   Lithium
//   double ralpha= 1.45/2.0;
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
// Lebedev N=38; eta 0.877 Lebedev 1976 ZVMMF_15_48 table 9.1;
      A1 = double(1)/double(105);
      A3 = double(9)/double(280);
      q1 = 0.4597008433809831;
      C1 = double(1)/double(35);
      gen6_A1(0,one,A1);
      gen8_A3(6,overradthree,A3);
      gen24_Cn(14,q1,C1);
    }else if(this->nPts_ == 50){
// Lebedev N=50; eta 0.96 Lebedev 1976 ZVMMF_15_48 table 11.1;
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
// Lebedev N=110; eta 0.982 Lebedev 1976 ZVMMF_15_48 table 11.1;
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
      C1 = 0.00969499636166;
      q1 = 0.478369028812;
      gen24_Cn(86,q1,C1);
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

}; // namespace ChronusQ
