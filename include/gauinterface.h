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
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

const int LEN_GAU_STR = 64*sizeof(char);
class GauMatEl{
  bool          doPrint_;

  std::ifstream infile_;

  std::string   fname_;
  std::string   labFil_;   // Matrix element file type
  std::string   gVers_;    // Gaussian version which generated the file
  std::string   jobTitle_; // Title card of Gaussian job
  
  /*
   * ICGU = KLM
   *
   * K = 1/2 for spin aligned vs GHF
   * L = 1/2 for real vs complex
   * M = 1/2 for RHF/GHF vs UHF (1 vs 2 spin blocks)
   *
   */
  int           *iAn_;     // Atomic numbers
  int           *iAtTyp_;
  int           *iBfAtm_;
  int           *iBfTyp_;
  int           iVers_;    // Version number of file format
  int           nLab_;     // Number on general data records (??)
  int           nAtoms_;   // Number of atoms
  int           nBasis_;   // Number of basis functions
  int           nBsUse_;   // Number of linerly independant basis functions
  int           iCharge_;  // Molecular charge
  int           multip_;   // Spin multiplicity
  int           nE_;       // Number of electrons
  int           len12L_;   // Number of bytes for integer labels for 1d/2d
  int           len4L_;    // Number of bytes for integer labels for 4d
  int           iOpCl_;    // Open/Closed shell flag
  int           iCGU_;     // Complex and / or GHF flag
  int           nFC_;
  int           nFV_;
  int           iTran_;
  int           iDum_;
  int           nInitRem_;
  int           nShellAO_;
  int           nPrimAO_;
  int           nShellDB_; 
  int           nPrimDB_;
  int           nBTot_;

  double        *atmChg_;
  double        *cart_;
  double        *atmWgt_;

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

  void init_(std::string &);

public:
  GauMatEl(std::string name=""){
    this->init_(name);
  }
  ~GauMatEl(){
    delete[] iAn_;
    delete[] iAtTyp_;
    delete[] iBfAtm_;
    delete[] iBfTyp_;
    delete[] atmChg_;
    delete[] cart_;
    delete[] atmWgt_;
    infile_.close();
  };
  
  void readInitRecs(); 
};

#endif
