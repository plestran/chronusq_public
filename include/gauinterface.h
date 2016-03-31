/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
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
#ifndef INCLUDED_GAUINTERFACE
#define INCLUDED_GAUINTERFACE
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>

const int LEN_GAU_STR = 64*sizeof(char);
class GauMatEl{
  bool          doPrint_;
  bool          haveHeader_;
  bool          readInit_;

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

  int lenRec1_;
  int lenRec2_;
  int lenRec3_;
  int lenRec4_;
  int lenRec5_;
  int lenRec6_;
  int lenRec7_;
  int lenRec8_;
  int lenRec9_;
  int lenRec10_;
  int lenRec11_;
  int lenInit_;

  void init_(std::string &);

  std::vector<std::string> GauHeader;
  void initHeader_();

public:
  enum {
    dipole,
    quadrupole,
    octupole,
    overlap,
    corehama,
    corehamb,
    kinetic,
    orthbasis,
    epsa,
    epsb,
    moa,
    mob,
    dena,
    denb,
    scfdena,
    scfdenb,
    focka,
    fockb
  };
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
  inline std::string   fname(){   return this->fname_;};
  inline std::string   labFil(){  return this->labFil_;};   
  inline std::string   gVers(){   return this->gVers_;};    
  inline std::string   jobTitle(){return this->jobTitle_;}; 
  inline int *         iAn(){     return this->iAn_;};     
  inline int *         iAtTyp(){  return this->iAtTyp_;};
  inline int *         iBfAtm(){  return this->iBfAtm_;};
  inline int *         iBfTyp(){  return this->iBfTyp_;};
  inline int           iVers(){   return this->iVers_;};    
  inline int           nLab(){    return this->nLab_;};     
  inline int           nAtoms(){  return this->nAtoms_;};   
  inline int           nBasis(){  return this->nBasis_;};   
  inline int           nBsUse(){  return this->nBsUse_;};   
  inline int           iCharge(){ return this->iCharge_;};  
  inline int           multip(){  return this->multip_;};   
  inline int           nE(){return this->nE_;};       
  inline int           len12L(){return this->len12L_;};   
  inline int           len4L(){return this->len4L_;};    
  inline int           iOpCl(){return this->iOpCl_;};    
  inline int           iCGU(){return this->iCGU_;};     
  inline int           nFC(){return this->nFC_;};
  inline int           nFV(){return this->nFV_;};
  inline int           iTran(){return this->iTran_;};
  inline int           iDum(){return this->iDum_;};
  inline int           nInitRem(){return this->nInitRem_;};
  inline int           nShellAO(){return this->nShellAO_;};
  inline int           nPrimAO(){return this->nPrimAO_;};
  inline int           nShellDB(){return this->nShellDB_;}; 
  inline int           nPrimDB(){return this->nPrimDB_;};
  inline int           nBTot(){return this->nBTot_;};

  inline double *      atmChg(){return this->atmChg_;};
  inline double *      cart(){return this->cart_;};
  inline double *      atmWgt(){return this->atmWgt_;};

  double * readRec(int);
};

class GauJob{
  std::ofstream gauInput_;
  std::string   gauFName_;
  std::string   basisName_;
  bool          doOpt_;
  double*       cart_;
  int           charge_;
  int           multip_;
  std::vector<int> atoms_;
  std::unique_ptr<GauMatEl>     matEl_;

  void genInput();
  

public:
  GauJob(bool doOpt, std::string basis, double* cart, std::vector<int> &atoms, int charge, int multip, std::string name="TEMP"){

    this->doOpt_      = doOpt;
    this->gauFName_   = name;
    this->basisName_  = basis;
    this->cart_ = cart;
    this->charge_ = charge;
    this->multip_ = multip;
    this->atoms_ = atoms;

    this->genInput();
  };

  ~GauJob(){this->gauInput_.close();};

   void run();
   
};

#endif
