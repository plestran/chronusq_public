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
template <typename T>
void RealTime<T>::alloc() {
  ssPropagator_ = std::unique_ptr<SingleSlater<dcomplex>>(
    new SingleSlater<dcomplex>(groundState_));
  ssPropagator_->isNotPrimary();

//prettyPrintSmart(cout,*groundState_->onePDMScalar(),"GSPS");
//prettyPrintSmart(cout,*groundState_->fockScalar(),"GSFS");
//prettyPrintSmart(cout,*ssPropagator_->onePDMScalar(),"PS");
//prettyPrintSmart(cout,*ssPropagator_->fockScalar(),"FS");

//CErr();
  auto NB  = ssPropagator_->basisset()->nBasis();
  auto NBT = NB * ssPropagator_->nTCS();

  bool needsFOSav = iRstScheme_ == ExplicitMagnus2 or 
                    iScheme_    == ExpMagnus2;

  NBTSqScratch_  = memManager_->template malloc<dcomplex>(NBT*NBT);
  NBTSqScratch2_ = memManager_->template malloc<dcomplex>(NBT*NBT);
  NBTSqScratch3_ = memManager_->template malloc<dcomplex>(NBT*NBT);



  UTransScalar_ = memManager_->template malloc<dcomplex>(NB*NB);
  POScalarSav_  = std::unique_ptr<ComplexMap>(
    new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

  POSav_.emplace_back(POScalarSav_.get());

  if(needsFOSav) {
    FOScalarSav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    FOSav_.emplace_back(FOScalarSav_.get());
  }

  if(ssPropagator_->nTCS() == 2 or !ssPropagator_->isClosedShell){
    UTransMz_= memManager_->template malloc<dcomplex>(NB*NB);
    POMzSav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    POSav_.emplace_back(POMzSav_.get());

    if(needsFOSav) {
      FOMzSav_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav_.emplace_back(FOMzSav_.get());
    }
  }
  if(ssPropagator_->nTCS() == 2){ 
    UTransMx_= memManager_->template malloc<dcomplex>(NB*NB);
    UTransMy_= memManager_->template malloc<dcomplex>(NB*NB);
    POMySav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
    POMxSav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    POSav_.emplace_back(POMySav_.get());
    POSav_.emplace_back(POMxSav_.get());

    if(needsFOSav) {
      FOMySav_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
      FOMxSav_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav_.emplace_back(FOMySav_.get());
      FOSav_.emplace_back(FOMxSav_.get());
    }
  }

  groundState_->transformOrthoMO();
}

