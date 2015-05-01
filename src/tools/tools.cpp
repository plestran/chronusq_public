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
#include <tools.h>
#include <atoms.h>
//--------------------------------------//
// factorial function:  t! = 1*2*3...*t //
//--------------------------------------//
namespace ChronusQ { 
double factorial(int t){
  int i;
  double tmp = 1.0;
  if (t<0 ) CErr("Factorial (t!) only defined on domain t in [0,inf)");
  if (t==0) return 1.0;
  else {
    for(i=1;i<=t;i++) tmp *= i;
    return tmp;
  };
};
//-----------------------------------------------//
// double factorial:  (2t-1)!! = 1*3*5*..*(2t-1) //
//-----------------------------------------------//
double doubleFact(int t){
  int i;
  double tmp = 1.0;
  if (t<0 ) CErr("Double factorial (t!!) only defined on domain t in [0,inf)");
  if (t==0) return 1.0;
  else  {
    for(i=1;i<=t;i++) tmp *= (2*i-1);
    return tmp;
  };
};
//--------------------------//
// (base)^order:  0**0 == 1 //
//--------------------------//
double powerInt(double base, int order){
  int i;
  double tmp=1.0;
  if(abs(base)<math.small&&order!=0 )  return 0.0;
  if(order==0)  return 1.0;
  else{
    for(i=0;i<order;i++) tmp *= base;
    return tmp;
  };
};
//---------------------------------------------------------//
// polynomial coeff:   (x+a)^l= sum(coeff * x^i * a^(l-i)) //
//---------------------------------------------------------//
double polyCoeff(int l, int i){
  if (l>= i) return factorial(l)/( factorial(i)*factorial(l-i) );
  else throw 10002; // FIXME need a CErr for this, don't understand the error
};
//-----------------------------------------------//
// (x+pA)^a*(x+pB)^b = sum[coeff * x^k]; k<=a+b  //
//-----------------------------------------------//
double kCoeff(int k, int a, int b, double pA, double pB){
  if(a<0||b<0||k>a+b) throw 10003; // FIXME need a CErr for this, don't understand the error
  int from, to;
  int it,i,j;
  double tmp=0;
  if( k<=(2*a-k) ) to  = k; else to  =2*a-k;
  if( k<=(2*b-k) ) from=-k; else from=k-2*b;
  for(it=from;it<=to;it=it+2){
    i=(k+it)/2;
    j=(k-it)/2;
    tmp += ( powerInt(pA,a-i)*polyCoeff(a,i)*powerInt(pB,b-j)*polyCoeff(b,j) );
  };
  return tmp;
};
//------------------------------//
// make an uppercase copy of s  //
//------------------------------//
void strupr(char *a) {
  while (*a != '\0') {
    if (islower (*a)) *a = toupper (*a);
    ++a;
  }
};
//-----------------------------//
// make a lowercase copy of s  //
//-----------------------------//
void strlwr(char *a) {
  while (*a != '\0') {
    if (isupper (*a)) *a = tolower (*a);
    ++a;
  }
};
//-----------------------------//
//  make give string to lower  //
//  ---------------------------//
std::string stringlower(std::string str){
  for (int i =0;i<str.size();i++){
    str[i]=tolower(str[i]);
  }; 
  return str;
};

//-----------------------------//
// return uppercase of the str //
// ----------------------------//
std::string stringupper(std::string str){
  std::string ucStr(str);
  for (auto i=0;i<str.size();i++){
    ucStr[i]=toupper(str[i]);
  }
  return ucStr;
}
//-----------------------------------------------------------------------------//
// convert atomic symbol and mass to an index in the table of atom in atoms.h  //
//-----------------------------------------------------------------------------//
int HashAtom(std::string element, int massNumber) { 
  element = stringlower(element);
  std::string currentAtom;
  
  for (auto i=0;i<atom.size();i++){
    currentAtom=atom[i].symbol;
    currentAtom=stringlower(currentAtom);
    if (massNumber>0){
      if (!element.compare(currentAtom)&&massNumber==atom[i].massNumber){
        return i;
        break;
      };
    };
    
    if (!element.compare(currentAtom)&&!atom[i].stable.compare("Y")){
      return i;
      break;
    };
  };
  return -1; 
}; 
//-------------------------------------------------------------------------//
// convert shell symbol or angular momentum to number of AOs in the shell  //
//-------------------------------------------------------------------------//
int HashNAOs(int L) { 
  return (L+1)*(L+2)/2;
};

int HashNAOs(std::string symbol) { 
  symbol = stringupper(symbol);
  if(!symbol.compare("S")) return 1; 
  if(!symbol.compare("P")) return 3; 
  if(!symbol.compare("D")) return 6; 
  if(!symbol.compare("F")) return 10; 
  if(!symbol.compare("G")) return 15;
  if(!symbol.compare("H")) return 21;
  return -1; 
}; 
//-------------------------------------------//
// convert shell symbol to angular momentum  //
//-------------------------------------------//
int HashL(std::string symbol) {
  if(!symbol.compare("S")) return 0; 
  if(!symbol.compare("P")) return 1; 
  if(!symbol.compare("D")) return 2; 
  if(!symbol.compare("F")) return 3; 
  if(!symbol.compare("G")) return 4; 
  if(!symbol.compare("H")) return 5; 
  if(!symbol.compare("I")) return 6; 
  if(!symbol.compare("J")) return 7; 
  if(!symbol.compare("K")) return 8; 
  if(!symbol.compare("L")) return 9; 
  if(!symbol.compare("M")) return 10; 
  if(!symbol.compare("N")) return 11; 
  if(!symbol.compare("O")) return 12; 
  if(!symbol.compare("P")) return 13; 
  if(!symbol.compare("Q")) return 14; 
  if(!symbol.compare("R")) return 15; 
  if(!symbol.compare("T")) return 16; 
  if(!symbol.compare("U")) return 17; 
  if(!symbol.compare("V")) return 18; 
  if(!symbol.compare("W")) return 19; 
  if(!symbol.compare("X")) return 20; 
  if(!symbol.compare("Y")) return 21; 
  if(!symbol.compare("Z")) return 22; 
  return -1; 
};
//----------------------------------------------------//
// convert angular momentum to AO index in a shell    //
// for L = 1: (0,1,2) -> (x,y,z)                      //
// for L = 2: (0,1,2,3,4,5) -> (x^2,y^2,z^2,xy,xz,yz) //
//----------------------------------------------------//
int HashIAO(int L,int *l) {
  if(L==0) return 0;
  else if(L==1) {
    if(l[0]==1) return 0;
    else if(l[1]==1) return 1;
    else if(l[2]==1) return 2;
  } else if(L==2) {
    if(l[0]==2) return 0;
    else if(l[1]==2) return 1;
    else if(l[2]==2) return 2;
    else if(l[0]==1) {
      if(l[1]==1) return 3;
      else return 4;
    } else return 5;
  };
};

} // namespace ChronusQ
