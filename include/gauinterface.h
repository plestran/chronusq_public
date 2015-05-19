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
#ifndef INCLUDED_GAUINTERFACE
#define INCLUDED_GAUINTERFACE
const int LEN_GAU_STR = 64*sizeof(char);
class GauMatEl{
  bool         doPrint_;

  ifstream     infile_;

  std::string  fname_;
  std::string  labFil_;
  std::string  gVers_;
  std::string  jobTitle_;
  
  int          *iAn_;
  int          *iAtTyp_;
  int          *iBfAtm_;
  int          *iBfTyp_;
  int          iVers_;
  int          nLab_;
  int          nAtoms_;
  int          nBasis_; 
  int          nBsUse_;
  int          iCharge_;
  int          multip_;
  int          nE_;
  int          len12L_;
  int          len4L_;
  int          iOpCl_;
  int          iCGU_;  
  int          nFC_;
  int          nFV_;
  int          iTran_;
  int          iDum_;
  int          nInitRem_;
  int          nShellAO_;
  int          nPrimAO_;
  int          nShellDB_; 
  int          nPrimDB_;
  int          nBTot_;

  double       *atmChg_;
  double       *cart_;
  double       *atmWgt_;

  void readGauRec1_();
  void readGauRec2_();
  void readGauRec3_();
  void readGauRec4_();
  void readGauRec5_();
  void readGauRec6_();
  void readGauRec7_();
  void readGauRec8_();
  void readGauRec9_();
  void readGauRec10_();
  void readGauRec11_();

  void init_();

public:
  GauMatEl(){this->init_();};
  ~GauMatEl(){
    delete[] iAn_;
    delete[] iAtTyp_;
    delete[] iBfAtm_;
    delete[] iBfTyp_;
    delete[] atmChg_;
    delete[] cart_;
    delete[] atmWgt_;
  };
  
  void readInitRecs(); 
}

#endif
