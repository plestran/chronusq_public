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
#include "tools.h"
//--------------------------------------//
// factorial function:  t! = 1*2*3...*t //
//--------------------------------------//
namespace ChronusQ { 
double factorial(int t){
  int i;
  double tmp = 1.0;
  if (t<0 ) throw 10000;
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
  if (t<0) throw 10001;
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
  else throw 10002;
};
//-----------------------------------------------//
// (x+pA)^a*(x+pB)^b = sum[coeff * x^k]; k<=a+b  //
//-----------------------------------------------//
double kCoeff(int k, int a, int b, double pA, double pB){
  if(a<0||b<0||k>a+b) throw 10003;
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
//-----------------------------------------------------------------------------//
// convert atomic symbol and mass to an index in the table of atom in atoms.h  //
//-----------------------------------------------------------------------------//
int HashAtom(char *symbol, int massNumber) { 
  strlwr(symbol);
  if(!strcmp(symbol,"h")&&(massNumber==0||massNumber==1)) return 0; 
  if(!strcmp(symbol,"d")||(!strcmp(symbol,"h")&&massNumber==2)) return 1; 
  if(!strcmp(symbol,"t")||(!strcmp(symbol,"h")&&massNumber==3)) return 2; 
  if(!strcmp(symbol,"he")&&(massNumber==0||massNumber==3)) return 3; 
  if(!strcmp(symbol,"he")&&massNumber==4) return 4; 
  if(!strcmp(symbol,"li")&&(massNumber==0||massNumber==6)) return 5; 
  if(!strcmp(symbol,"li")&&massNumber==7) return 6; 
  if(!strcmp(symbol,"c")) return 11;
  return -1; 
}; 
//-------------------------------------------------------------------------//
// convert shell symbol or angular momentum to number of AOs in the shell  //
//-------------------------------------------------------------------------//
int HashNAOs(int L) { 
  return (L+1)*(L+2)/2;
};

int HashNAOs(char *symbol) { 
  strupr(symbol);
  if(!strcmp(symbol,"S")) return 1; 
  if(!strcmp(symbol,"P")) return 3; 
  if(!strcmp(symbol,"D")) return 6; 
  if(!strcmp(symbol,"F")) return 10; 
  if(!strcmp(symbol,"G")) return 15;
  if(!strcmp(symbol,"H")) return 21;
  return -1; 
}; 
//-------------------------------------------//
// convert shell symbol to angular momentum  //
//-------------------------------------------//
int HashL(char *symbol) {
  strupr(symbol);
  if(!strcmp(symbol,"S")) return 0; 
  if(!strcmp(symbol,"P")) return 1; 
  if(!strcmp(symbol,"D")) return 2; 
  if(!strcmp(symbol,"F")) return 3; 
  if(!strcmp(symbol,"G")) return 4;
  if(!strcmp(symbol,"H")) return 5;
  if(!strcmp(symbol,"I")) return 6;
  if(!strcmp(symbol,"J")) return 7;
  if(!strcmp(symbol,"K")) return 8;
  if(!strcmp(symbol,"L")) return 9;
  if(!strcmp(symbol,"M")) return 10;
  if(!strcmp(symbol,"N")) return 11;
  if(!strcmp(symbol,"O")) return 12;
  if(!strcmp(symbol,"P")) return 13;
  if(!strcmp(symbol,"Q")) return 14;
  if(!strcmp(symbol,"R")) return 15;
  if(!strcmp(symbol,"T")) return 16;
  if(!strcmp(symbol,"U")) return 17;
  if(!strcmp(symbol,"V")) return 18;
  if(!strcmp(symbol,"W")) return 19;
  if(!strcmp(symbol,"X")) return 20;
  if(!strcmp(symbol,"Y")) return 21;
  if(!strcmp(symbol,"Z")) return 22;
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
