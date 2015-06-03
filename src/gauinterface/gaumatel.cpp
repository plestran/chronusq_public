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
#include <gauinterface.h>

void GauMatEl::init_(std::string& name){
  this->fname_ = name;
//if(name.compare("")) return;

  infile_.open(name, std::ios::binary | std::ios::in);

/*  if(!infile_){
    std::cout << "Fatal Error: Could not open MatEl file" << std::endl;
    exit(1);
  }
*/

  labFil_   = '\0';
  gVers_    = '\0';
  jobTitle_ = '\0';

  labFil_.resize(  LEN_GAU_STR,' ');
  gVers_.resize(   LEN_GAU_STR,' ');
  jobTitle_.resize(LEN_GAU_STR,' ');

  this->iAn_    = NULL;
  this->iAtTyp_ = NULL;
  this->iBfAtm_ = NULL;
  this->iBfTyp_ = NULL;
  this->atmChg_ = NULL;
  this->cart_   = NULL;
  this->atmWgt_ = NULL;

  this->doPrint_ = false;
  this->readInitRecs();
  this->initHeader_();
}

void GauMatEl::readGauRec1_(){
  infile_.read((char*)&this->labFil_[0], LEN_GAU_STR ); 
  infile_.read((char*)&this->iVers_,     sizeof(int) ); 
  infile_.read((char*)&this->nLab_,      sizeof(int) ); 
  infile_.read((char*)&this->gVers_[0],  LEN_GAU_STR ); 

  if(doPrint_) {
    std::cout << std::setw(15) << std::left << "Record 1:"  << std::endl;
    std::cout << std::setw(15) << std::left << "  LabFil: " << 
                 this->labFil_ << std::endl;
    std::cout << std::setw(15) << std::left << "  IVers:  " << 
                 this->iVers_  << std::endl;
    std::cout << std::setw(15) << std::left << "  NLab:   " << 
                 this->nLab_   << std::endl;
    std::cout << std::setw(15) << std::left << "  GVers:  " << 
                 this->gVers_  << std::endl;
  } 
//else std::cout << "Record 1 Sucessfully Read" << std::endl;
  this->lenRec1_ = 2*LEN_GAU_STR + 2*sizeof(int);

}
void GauMatEl::readGauRec2_(){
  infile_.read((char*)&this->jobTitle_[0], LEN_GAU_STR ); 
  infile_.read((char*)&this->nAtoms_,      sizeof(int) ); 
  infile_.read((char*)&this->nBasis_,      sizeof(int) ); 
  infile_.read((char*)&this->nBsUse_,      sizeof(int) ); 
  infile_.read((char*)&this->iCharge_,     sizeof(int) ); 
  infile_.read((char*)&this->multip_,      sizeof(int) ); 
  infile_.read((char*)&this->nE_,          sizeof(int) ); 
  infile_.read((char*)&this->len12L_,      sizeof(int) ); 
  infile_.read((char*)&this->len4L_,       sizeof(int) ); 
  infile_.read((char*)&this->iOpCl_,       sizeof(int) ); 
  infile_.read((char*)&this->iCGU_,        sizeof(int) ); 
  /*
   * ICGU = KLM
   *
   * K = 1/2 for spin aligned vs GHF
   * L = 1/2 for real vs complex
   * M = 1/2 for RHF/GHF vs UHF (1 vs 2 spin blocks)
   *
   */
  
  if(doPrint_) {
    std::cout << std::setw(15) << std::left << "Record 2:"   << std::endl;
    std::cout << std::setw(15) << std::left << "  Title:   " << 
                 this->jobTitle_            << std::endl;
    std::cout << std::setw(15) << std::left << "  this->nAtoms_:  " << 
                 this->nAtoms_ << std::endl; 
    std::cout << std::setw(15) << std::left << "  this->nBasis_:  " << 
                 this->nBasis_ << std::endl; 
    std::cout << std::setw(15) << std::left << "  NBsUse:  " << 
                 this->nBsUse_ << std::endl; 
    std::cout << std::setw(15) << std::left << "  ICharge: " << 
                 this->iCharge_<< std::endl; 
    std::cout << std::setw(15) << std::left << "  Multip:  " << 
                 this->multip_ << std::endl; 
    std::cout << std::setw(15) << std::left << "  NE:      " <<
                 this->nE_     << std::endl; 
    std::cout << std::setw(15) << std::left << "  Len12L:  " << 
                 this->len12L_ << std::endl; 
    std::cout << std::setw(15) << std::left << "  Len4L:   " << 
                 this->len4L_  << std::endl; 
    std::cout << std::setw(15) << std::left << "  IOpCl:   " << 
                 this->iOpCl_  << std::endl; 
    std::cout << std::setw(15) << std::left << "  ICGU:    " << 
                 this->iCGU_   << std::endl; 
  }
//else std::cout << "Record 2 Sucessfully Read" << std::endl;
  this->lenRec2_ = LEN_GAU_STR + 10*sizeof(int);

}
void GauMatEl::readGauRec3_(){
  this->iAn_ = new int[this->nAtoms_+this->nAtoms_%2];
  infile_.read((char*)this->iAn_,(this->nAtoms_+this->nAtoms_%2)*
                sizeof(int));
  
  if(doPrint_) {
    std::cout << std::setw(15) << std::left << "Record 3:" << std::endl;
    std::cout << std::setw(15) << std::left << "  IAn: "   << std::endl;
    for(auto i = 0; i < this->nAtoms_+this->nAtoms_%2; i++) 
      std::cout << std::setw(15) << std::left << " " << 
                   *(this->iAn_+i) << std::endl;
  }
//else std::cout << "Record 3 Sucessfully Read" << std::endl;
  this->lenRec3_ = (this->nAtoms_+this->nAtoms_%2)*sizeof(int);
}
void GauMatEl::readGauRec4_(){
  this->iAtTyp_ = new int[this->nAtoms_+this->nAtoms_%2];
  infile_.read((char*)this->iAtTyp_,(this->nAtoms_+this->nAtoms_%2)*
              sizeof(int));
  
  if(doPrint_) {
    std::cout << std::setw(15) << std::left << "Record 4:" << std::endl;
    std::cout << std::setw(15) << std::left << "  IAtTyp: " << std::endl;
    for(auto i = 0; i < this->nAtoms_+this->nAtoms_%2; i++) 
      std::cout << std::setw(15) << std::left << " " << 
                   *(this->iAtTyp_+i) << std::endl;
  } 
//else std::cout << "Record 4 Sucessfully Read" << std::endl;
  this->lenRec4_ = (this->nAtoms_+this->nAtoms_%2)*sizeof(int);
}
void GauMatEl::readGauRec5_(){
  this->atmChg_ = new double[this->nAtoms_];
  infile_.read((char*)this->atmChg_,this->nAtoms_*sizeof(double));

  if(doPrint_){
    std::cout << std::setw(15) << std::left << "Record 5:" << std::endl;
    std::cout << std::setw(15) << std::left << "  AtmChg_: " << std::endl;
    for(auto i = 0; i < this->nAtoms_; i++) 
      std::cout << std::setw(15) << std::fixed << std::left << " " << *(this->atmChg_+i) << std::endl;
  } 
//else std::cout << "Record 5 Sucessfully Read" << std::endl;
  this->lenRec5_ = (this->nAtoms_)*sizeof(double);
}
void GauMatEl::readGauRec6_(){
  this->cart_ = new double[3*this->nAtoms_];
  infile_.read((char*)cart_,3*this->nAtoms_*sizeof(double));

  if(doPrint_){
    std::cout << std::setw(15) << std::left << "Record 6:" << std::endl;
    std::cout << std::setw(15) << std::left << "  C: " << std::endl;
    for(auto i = 0;    i < this->nAtoms_; i++){ 
      for(auto ixyz = 0; ixyz < 3;   ixyz++){
        std::cout << std::setw(15) << std::fixed << std::right << *(cart_+ixyz+i*3);
      }
      std::cout << std::endl;
    }
  } 
//else std::cout << "Record 6 Sucessfully Read" << std::endl;
  this->lenRec6_ = (3*this->nAtoms_)*sizeof(double);
}
void GauMatEl::readGauRec7_(){
  this->iBfAtm_ = new int[this->nBasis_];
  this->iBfTyp_ = new int[this->nBasis_];
  infile_.read((char*)this->iBfAtm_,(this->nBasis_)*sizeof(int));
  infile_.read((char*)this->iBfTyp_,(this->nBasis_)*sizeof(int));
  
  if(doPrint_){
    std::cout << std::setw(15) << std::left << "Record 7:" << std::endl;
    std::cout << std::setw(15) << std::left << "  IBfAtm_: " << std::endl;
    for(auto i = 0; i < this->nBasis_; i++) 
      std::cout << std::setw(15) << std::left << " " << *(this->iBfAtm_+i) << std::endl;
    std::cout << std::setw(15) << std::left << "  IBfTyp_: " << std::endl;
    for(auto i = 0; i < this->nBasis_; i++) 
      std::cout << std::setw(15) << std::left << " " << *(this->iBfTyp_+i) << std::endl;
  }
//else std::cout << "Record 7 Sucessfully Read" << std::endl;
  this->lenRec7_ = (2*this->nBasis_)*sizeof(int);
}
void GauMatEl::readGauRec8_(){
  this->atmWgt_ = new double[this->nAtoms_];
  infile_.read((char*)this->atmWgt_,this->nAtoms_*sizeof(double));

  if(doPrint_) {
    std::cout << std::setw(15) << std::left << "Record 8:" << std::endl;
    std::cout << std::setw(15) << std::left << "  AtmWgt_: " << std::endl;
    for(auto i = 0; i < this->nAtoms_; i++) 
      std::cout << std::setw(15) << std::fixed << std::left << " " << *(this->atmWgt_+i) <<std::endl;
  }
//else std::cout << "Record 8 Sucessfully Read" << std::endl;
  this->lenRec8_ = (this->nAtoms_)*sizeof(double);
}
void GauMatEl::readGauRec9_(){
  infile_.read((char*)&this->nFC_,  sizeof(int) );
  infile_.read((char*)&this->nFV_,  sizeof(int) );
  infile_.read((char*)&this->iTran_,sizeof(int) );
  infile_.read((char*)&this->iDum_, sizeof(int) );
  
  if(doPrint_) {
    std::cout << std::setw(15) << std::left << "Record 9:" << std::endl;
    std::cout << std::setw(15) << std::left << "  NFC: "   << this->nFC_   << std::endl;
    std::cout << std::setw(15) << std::left << "  NFV: "   << this->nFV_   << std::endl;
    std::cout << std::setw(15) << std::left << "  ITran: "   << this->iTran_   << std::endl;
    std::cout << std::setw(15) << std::left << "  IDum: "   << this->iDum_   << std::endl;
  } 
//else std::cout << "Record 9 Sucessfully Read" << std::endl;
  this->lenRec9_ = 4*sizeof(int);
}
void GauMatEl::readGauRec10_(){
  int IX;
  infile_.read((char*)&this->nInitRem_,sizeof(int));
  infile_.read((char*)&IX,sizeof(int)); // Weird offset!

  if(doPrint_) {
    std::cout << std::setw(15) << std::left << std::endl;
    std::cout << std::setw(15) << std::left << "Record 10:" << std::endl;
    std::cout << std::setw(15) << std::left << "  Number of Remaining Recs: "   << this->nInitRem_   << std::endl;
  } 
//else std::cout << "Record 10 Sucessfully Read" << std::endl;
  this->lenRec10_ = 2*sizeof(int);
}
void GauMatEl::readGauRec11_(){
  int NRecs = 5;
  int IX;


  infile_.read((char*)&this->nShellAO_,sizeof(int) );
  infile_.read((char*)&this->nPrimAO_, sizeof(int) );
  infile_.read((char*)&this->nShellDB_,sizeof(int) );
  infile_.read((char*)&this->nPrimDB_, sizeof(int) );
  infile_.read((char*)&this->nBTot_,   sizeof(int) );
  for(auto i = 0; i < this->nInitRem_-NRecs; i++)
    infile_.read((char*)&IX,sizeof(int));

  
  if(doPrint_) {
    std::cout << std::setw(15) << std::left << std::endl;
    std::cout << std::setw(15) << std::left << "Record 11:" << std::endl;
    std::cout << std::setw(15) << std::left << "  NShellAO: "   << this->nShellAO_   << std::endl;
    std::cout << std::setw(15) << std::left << "  NPrimAO: "   << this->nPrimAO_   << std::endl;
    std::cout << std::setw(15) << std::left << "  NShellDB: "   << this->nShellDB_   << std::endl;
    std::cout << std::setw(15) << std::left << "  NPrimDB: "   << this->nPrimDB_   << std::endl;
    std::cout << std::setw(15) << std::left << "  NBTot: "   << this->nBTot_   << std::endl;
  }
//else std::cout << "Record 11 Sucessfully Read" << std::endl;
  this->lenRec11_ = this->nInitRem_*sizeof(int);
}
void GauMatEl::readInitRecs(){
  readGauRec1_();
  readGauRec2_();
  readGauRec3_();
  readGauRec4_();
  readGauRec5_();
  readGauRec6_();
  readGauRec7_();
  readGauRec8_();
  readGauRec9_();
  readGauRec10_();
  readGauRec11_();
  this->readInit_=true;
  this->lenInit_ = 
        this->lenRec1_ +  this->lenRec2_ + this->lenRec3_ +
        this->lenRec4_ +  this->lenRec5_ + this->lenRec6_ +
        this->lenRec7_ +  this->lenRec8_ + this->lenRec9_ +
        this->lenRec10_ + this->lenRec11_;
}
void GauMatEl::initHeader_(){
  GauHeader = {
    "DIPOLE INTEGRALS",
    "QUADRUPOLE INTEGRALS",
    "OCTUPOLE INTEGRALS",
    "OVERLAP",
    "CORE HAMILTONIAN ALPHA",
    "CORE HAMILTONIAN BETA",
    "KINETIC",
    "ORTHOGONAL BASIS",
    "ALPHA ORBITAL ENERGIES",
    "BETA ORBITAL ENERGIES",
    "ALPHA MO COEFFICIENTS",
    "BETA MO COEFFICIENTS",
    "ALPHA DENSITY MATRIX",
    "BETA DENSITY MATRIX",
    "ALPHA SCF DENSITY MATRIX",
    "BETA SCF DENSITY MATRIX",
    "ALPHA FOCK MATRIX",
    "BETA FOCK MATRIX" 
  };

  for(auto i = 0; i < GauHeader.size(); i++)
    GauHeader[i].resize(LEN_GAU_STR,' ');
}
double * GauMatEl::readRec(int rec){
  bool found = false;

  std::string FileTerm = "END";
  FileTerm.resize(LEN_GAU_STR,' ');
  infile_.seekg(this->lenInit_,infile_.beg);

  while(!infile_.eof()){
    int      *ID = NULL;
    double   *DX = NULL;
    std::string Label(LEN_GAU_STR,'\0');
    int NI, NR, NTot, NPerRec, N1, N2, N3, N4, N5, ISym;

    infile_.read(&Label[0],      LEN_GAU_STR );
    if(!Label.compare(FileTerm)) break;

    infile_.read((char*)&NI,     sizeof(int) );
    infile_.read((char*)&NR,     sizeof(int) );
    infile_.read((char*)&NTot,   sizeof(int) );
    infile_.read((char*)&NPerRec,sizeof(int) );
    infile_.read((char*)&N1,     sizeof(int) );
    infile_.read((char*)&N2,     sizeof(int) );
    infile_.read((char*)&N3,     sizeof(int) );
    infile_.read((char*)&N4,     sizeof(int) );
    infile_.read((char*)&N5,     sizeof(int) );
    infile_.read((char*)&ISym,   sizeof(int) );
    int NRec = (NTot+NPerRec-1)/NPerRec;
    ID = new int[NRec*NI*NPerRec];
    DX = new double[NRec*NR*NPerRec];
    for(auto iRec = 0; iRec < NRec; iRec++){
      if(NI != 0) {
        infile_.read((char*)(ID+iRec*NI*NPerRec),NI*NPerRec*sizeof(int));
      }
      if(NR >  0) {
        infile_.read((char*)(DX+iRec*NR*NPerRec),NR*NPerRec*sizeof(double));
      }
    }
    if(!Label.compare(GauHeader[rec])){
      double *data = new double[NTot];
      for(auto i = 0; i < NTot; i++){
        data[i] = DX[i];
      }
     
      found = true;
      delete [] ID;
      delete [] DX;
      return data;
      break; 
    }
    
    delete [] ID;
    delete [] DX;
  }
//if(found) std::cout << "FOUND " << GauHeader[rec] << std::endl;

}
