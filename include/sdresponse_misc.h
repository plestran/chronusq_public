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


  void formAOTDen(const TVecMap &, TMatrix &, TMatrix &);
  void formMOTDen(TVecMap &, const TMatrix &, const TMatrix &);
  void placeVOOV(const TVecMap &,TMatrix &, TMatrix&);
  void placeVVOO(const TVecMap &,TMatrix &);
  void retrvVOOV(TVecMap &,const TMatrix &,const TMatrix&);
  void retrvVVOO(TVecMap &,const TMatrix &);
  void initRMu();
  void scaleDagPPRPA(bool,TVecMap &,TVecMap &,TVecMap *AX=NULL); 
  void initMeth();
  void checkValid();

  /*************************/
  /* MPI Related Routines  */
  /*************************/
  void mpiSend(int,int tag=tagSDResponse);
  void mpiRecv(int,int tag=tagSDResponse);

  // Proof of concept routines
  void incorePPRPA();
  void incorePPRPAnew();
  void incoreCIS();
  void incoreRPA();
  void formRM();
  TMatrix formRM2(TMatrix &XMO);
