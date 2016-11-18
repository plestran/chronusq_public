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
#include <basisset.h>
namespace ChronusQ{
/**
 *  Generate shell index -> starting basis function map
 */
void BasisSet::makeMapSh2Bf(){
  this->mapSh2Bf_.clear();
  auto n = 0;
  for(auto shell : this->shells_){
     this->mapSh2Bf_.push_back(n);
     n += shell.size();
  } // loop shell
  this->haveMapSh2Bf = true;
} // BasisSet::makeMapSh2Bf

/**
 *  Generate shell index -> atomic center map
 */
void BasisSet::makeMapSh2Cen(Molecule *mol){
  this->mapSh2Cen_.clear();
  for(auto shell : this->shells_){
    for(auto iAtom = 0; iAtom < mol->nAtoms(); iAtom++){
      std::array<double,3> center = {{ (*mol->cart())(0,iAtom),
                                       (*mol->cart())(1,iAtom),
                                       (*mol->cart())(2,iAtom) }};
      if(shell.O == center){
        this->mapSh2Cen_.push_back(iAtom+1);
        break;
      }
    } // loop iAtom 
  } // loop shell
  this->haveMapSh2Cen = true;
} // BasisSet::makeMapSh2Cen

/**
 *  Generate atomic center index -> starting basis function map
 */
void BasisSet::makeMapCen2Bf(Molecule *mol){
  if(!this->haveMapSh2Bf ) this->makeMapSh2Bf();
  if(!this->haveMapSh2Cen) this->makeMapSh2Cen(mol);
  this->mapCen2Bf_.clear();

  for(auto iAtm = 0; iAtm < mol->nAtoms(); iAtm++){
    auto nSize = 0;
    for(auto iShell = 0; iShell < this->nShell_; iShell++){
      if((iAtm+1) == this->mapSh2Cen_[iShell]) 
        nSize += this->shells_[iShell].size();
    } // loop iShell
    auto iSt = -1;
    for(auto iShell = 0; iShell < this->nShell_; iShell++){
      if((iAtm+1) == this->mapSh2Cen_[iShell]){
       iSt = this->mapSh2Bf_[iShell];
       break;
      }
    } // loop iShell
    if(iSt == -1) 
      CErr("Could not find Center in Basis definition",this->fileio_->out);
    this->mapCen2Bf_.push_back({{iSt,nSize}});
  } // loop iAtm
  this->haveMapCen2Bf = true;
} // BasisSet::makeMapCen2Bf

void BasisSet::makeBasisMap(){
  // Pople Sets
  this->basisMap[PopleSTO3G]    = std::string(BASIS_PATH) + "/sto3g.gbs";
  this->basisMap[PopleSTO6G]    = std::string(BASIS_PATH) + "/sto6g.gbs";
  this->basisMap[Pople321G]     = std::string(BASIS_PATH) + "/3-21g.gbs";
  this->basisMap[Pople431G]     = std::string(BASIS_PATH) + "/4-31g.gbs";
  this->basisMap[Pople631G]     = std::string(BASIS_PATH) + "/6-31g.gbs";
  this->basisMap[Pople631ppGs]  = std::string(BASIS_PATH) + "/6-31++g*.gbs";
  this->basisMap[Pople6311pGs]  = std::string(BASIS_PATH) + "/6-311+g*.gbs";
  this->basisMap[Pople6311pGss] = std::string(BASIS_PATH) + "/6-311+g**.gbs";
  this->basisMap[Pople6311pG2dp]= std::string(BASIS_PATH) + "/6-311+g_2d_p.gbs";

  // Dunning Sets
  this->basisMap[ccpVDZ]        = std::string(BASIS_PATH) + "/cc-pvdz.gbs";
  this->basisMap[ccpVTZ]        = std::string(BASIS_PATH) + "/cc-pvtz.gbs";

  // Turbomol sets
  this->basisMap[def2SVP]       = std::string(BASIS_PATH) + "/def2-svp.gbs";
  this->basisMap[def2SVPD]      = std::string(BASIS_PATH) + "/def2-svpd.gbs";
  this->basisMap[def2TZVP]      = std::string(BASIS_PATH) + "/def2-tzvp.gbs";

  // Sapporo Sets
    
  // Standard Sets
  this->basisMap[SapporoDZP]    = 
    std::string(BASIS_PATH) + "/sapporo-dzp_no.gbs";
  this->basisMap[SapporoTZP]    = 
    std::string(BASIS_PATH) + "/sapporo-tzp_no.gbs";
  this->basisMap[SapporoQZP]    = 
    std::string(BASIS_PATH) + "/sapporo-qzp_no.gbs";

  this->basisMap[SapporoDZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-dzp_sp.gbs";
  this->basisMap[SapporoTZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-tzp_sp.gbs";
  this->basisMap[SapporoQZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-qzp_sp.gbs";

  this->basisMap[SapporoDZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-dzp_all.gbs";
  this->basisMap[SapporoTZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-tzp_all.gbs";
  this->basisMap[SapporoQZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-qzp_all.gbs";

  // 2012 Sets
  this->basisMap[Sapporo2012DZP]    = 
    std::string(BASIS_PATH) + "/sapporo-dzp-2012_no.gbs";
  this->basisMap[Sapporo2012TZP]    = 
    std::string(BASIS_PATH) + "/sapporo-tzp-2012_no.gbs";
  this->basisMap[Sapporo2012QZP]    = 
    std::string(BASIS_PATH) + "/sapporo-qzp-2012_no.gbs";

  this->basisMap[Sapporo2012DZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-dzp-2012_sp.gbs";
  this->basisMap[Sapporo2012TZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-tzp-2012_sp.gbs";
  this->basisMap[Sapporo2012QZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-qzp-2012_sp.gbs";

  this->basisMap[Sapporo2012DZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-dzp-2012_all.gbs";
  this->basisMap[Sapporo2012TZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-tzp-2012_all.gbs";
  this->basisMap[Sapporo2012QZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-qzp-2012_all.gbs";

  // DKH3 Sets
  this->basisMap[SapporoDKH3DZP]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-dzp_no.gbs";
  this->basisMap[SapporoDKH3TZP]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-tzp_no.gbs";
  this->basisMap[SapporoDKH3QZP]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-qzp_no.gbs";

  this->basisMap[SapporoDKH3DZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-dzp_sp.gbs";
  this->basisMap[SapporoDKH3TZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-tzp_sp.gbs";
  this->basisMap[SapporoDKH3QZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-qzp_sp.gbs";

  this->basisMap[SapporoDKH3DZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-dzp_all.gbs";
  this->basisMap[SapporoDKH3TZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-tzp_all.gbs";
  this->basisMap[SapporoDKH3QZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-qzp_all.gbs";

  // 2012 DKH3 Sets
  this->basisMap[Sapporo2012DKH3DZP]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-dzp-2012_no.gbs";
  this->basisMap[Sapporo2012DKH3TZP]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-tzp-2012_no.gbs";
  this->basisMap[Sapporo2012DKH3QZP]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-qzp-2012_no.gbs";

  this->basisMap[Sapporo2012DKH3DZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-dzp-2012_sp.gbs";
  this->basisMap[Sapporo2012DKH3TZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-tzp-2012_sp.gbs";
  this->basisMap[Sapporo2012DKH3QZPsp]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-qzp-2012_sp.gbs";

  this->basisMap[Sapporo2012DKH3DZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-dzp-2012_all.gbs";
  this->basisMap[Sapporo2012DKH3TZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-tzp-2012_all.gbs";
  this->basisMap[Sapporo2012DKH3QZPall]    = 
    std::string(BASIS_PATH) + "/sapporo-dkh3-qzp-2012_all.gbs";

  // Pople Sets
  this->basisKey["STO-3G"]        = PopleSTO3G;
  this->basisKey["STO-6G"]        = PopleSTO6G;
  this->basisKey["3-21G"]         = Pople321G;
  this->basisKey["4-31G"]         = Pople431G;
  this->basisKey["6-31G"]         = Pople631G;

  this->basisKey["6-31++G*"]      = Pople631ppGs;
  this->basisKey["6-31++G(D)"]    = Pople631ppGs;

  this->basisKey["6-311+G*"]      = Pople6311pGs;
  this->basisKey["6-311+G(D)"]    = Pople6311pGs;

  this->basisKey["6-311+G**"]     = Pople6311pGss;
  this->basisKey["6-311+G(D,P)"]  = Pople6311pGss;

  this->basisKey["6-311+G(2D,P)"] = Pople6311pG2dp;

  // Dunning Sets
  this->basisKey["CC-PVDZ"]       = ccpVDZ;
  this->basisKey["CC-PVTZ"]       = ccpVTZ;

  // Turbomol Sets
  this->basisKey["DEF2-SVP"]      = def2SVP;
  this->basisKey["DEF2-SVPD"]     = def2SVPD;
  this->basisKey["DEF2-TZVP"]     = def2TZVP;

  // Sapporo Sets
    
  // Standard Sets
  this->basisKey["SAPPORO-DZP"]      = SapporoDZP;
  this->basisKey["SAPPORO-TZP"]      = SapporoTZP;
  this->basisKey["SAPPORO-QZP"]      = SapporoQZP;

  this->basisKey["SAPPORO-DZP-SP"]   = SapporoDZPsp;
  this->basisKey["SAPPORO-TZP-SP"]   = SapporoTZPsp;
  this->basisKey["SAPPORO-QZP-SP"]   = SapporoQZPsp;

  this->basisKey["SAPPORO-DZP-ALL"]  = SapporoDZPall;
  this->basisKey["SAPPORO-TZP-ALL"]  = SapporoTZPall;
  this->basisKey["SAPPORO-QZP-ALL"]  = SapporoQZPall;

  // 2012 Sets
  this->basisKey["SAPPORO-2012-DZP"]      = Sapporo2012DZP;
  this->basisKey["SAPPORO-2012-TZP"]      = Sapporo2012TZP;
  this->basisKey["SAPPORO-2012-QZP"]      = Sapporo2012QZP;

  this->basisKey["SAPPORO-2012-DZP-SP"]   = Sapporo2012DZPsp;
  this->basisKey["SAPPORO-2012-TZP-SP"]   = Sapporo2012TZPsp;
  this->basisKey["SAPPORO-2012-QZP-SP"]   = Sapporo2012QZPsp;

  this->basisKey["SAPPORO-2012-DZP-ALL"]  = Sapporo2012DZPall;
  this->basisKey["SAPPORO-2012-TZP-ALL"]  = Sapporo2012TZPall;
  this->basisKey["SAPPORO-2012-QZP-ALL"]  = Sapporo2012QZPall;

  // DKH3 Sets
  this->basisKey["SAPPORO-DKH3-DZP"]      = SapporoDKH3DZP;
  this->basisKey["SAPPORO-DKH3-TZP"]      = SapporoDKH3TZP;
  this->basisKey["SAPPORO-DKH3-QZP"]      = SapporoDKH3QZP;

  this->basisKey["SAPPORO-DKH3-DZP-SP"]   = SapporoDKH3DZPsp;
  this->basisKey["SAPPORO-DKH3-TZP-SP"]   = SapporoDKH3TZPsp;
  this->basisKey["SAPPORO-DKH3-QZP-SP"]   = SapporoDKH3QZPsp;

  this->basisKey["SAPPORO-DKH3-DZP-ALL"]  = SapporoDKH3DZPall;
  this->basisKey["SAPPORO-DKH3-TZP-ALL"]  = SapporoDKH3TZPall;
  this->basisKey["SAPPORO-DKH3-QZP-ALL"]  = SapporoDKH3QZPall;

  // 2012 DKH3 Sets
  this->basisKey["SAPPORO-2012-DKH3-DZP"]      = Sapporo2012DKH3DZP;
  this->basisKey["SAPPORO-2012-DKH3-TZP"]      = Sapporo2012DKH3TZP;
  this->basisKey["SAPPORO-2012-DKH3-QZP"]      = Sapporo2012DKH3QZP;

  this->basisKey["SAPPORO-2012-DKH3-DZP-SP"]   = Sapporo2012DKH3DZPsp;
  this->basisKey["SAPPORO-2012-DKH3-TZP-SP"]   = Sapporo2012DKH3TZPsp;
  this->basisKey["SAPPORO-2012-DKH3-QZP-SP"]   = Sapporo2012DKH3QZPsp;

  this->basisKey["SAPPORO-2012-DKH3-DZP-ALL"]  = Sapporo2012DKH3DZPall;
  this->basisKey["SAPPORO-2012-DKH3-TZP-ALL"]  = Sapporo2012DKH3TZPall;
  this->basisKey["SAPPORO-2012-DKH3-QZP-ALL"]  = Sapporo2012DKH3QZPall;

  // Sapporo Sets
}; // BasisSet::makeBasisMap

void BasisSet::makeMapPrim2Bf(){
  this->mapPrim2Bf_ = 
   std::unique_ptr<RealMatrix>(new RealMatrix(this->nBasis_,this->nPrimitive_));
  double * memAddress = this->mapPrim2Bf_->data();
  for (auto iSh = 0; iSh < this->nShell_; iSh++){
    int nPrim = this->shells_[iSh].contr[0].coeff.size();
    int nBf   = this->shells_[iSh].size();
    for (auto iP = 0; iP < nPrim; iP++)
    for (auto iBf = 0; iBf < nBf; iBf++){
      memAddress[(iP*nBf + iBf)*this->nBasis_ + iBf] = 
        this->unNormCons_[iSh][iP];
      }
    memAddress += this->nBasis_ * (nPrim * nBf) + nBf;
  }
  //  prettyPrint(cout,*this->mapPrim2Bf_,"MAP");

//
/*
for(auto iSh = 0, iBf = 0, iPrim = 0;
      iSh < this->nShell_; 
      iBf += this->shells_[iSh].size(), 
        iPrim += this->shells_[iSh].contr[0].coeff.size() * this->shells_[iSh].size(), 
        iSh++){
  
    int nPrim = this->shells_[iSh].contr[0].coeff.size();
    int nBf   = this->shells_[iSh].size();
    //RealMap PrimCoeff(&this->shells_[iSh].contr[0].coeff[0],1,nPrim);
    RealMap PrimCoeff(&this->unNormCons_[iSh][0],1,nPrim);
    for(auto jBf = iBf, jPrim = iPrim; jBf < iBf + nBf; jBf++, jPrim += nPrim){
      this->mapPrim2Bf_->block(jBf,jPrim,1,nPrim) = PrimCoeff;
    }
  }
*/
//  prettyPrint(cout,*this->mapPrim2Bf_,"MAP");
};

}; // namespace ChronusQ
