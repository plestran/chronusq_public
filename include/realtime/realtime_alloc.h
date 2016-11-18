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
                    iRstScheme_ == ExplicitMagnus3 or 
                    iScheme_    == ExpMagnus2 or
                    iScheme_    == ExpMagnus3;

  bool needsFOSav2 = iRstScheme_ == ExplicitMagnus3 or 
                     iScheme_    == ExpMagnus3;
  bool needsFOSav3 = iRstScheme_ == ExplicitMagnus3 or 
                     iScheme_    == ExpMagnus3;
  bool needsFOSav4 = iRstScheme_ == ExplicitMagnus3 or 
                     iScheme_    == ExpMagnus3;
  bool needsFOSav5 = iRstScheme_ == ExplicitMagnus3 or 
                     iScheme_    == ExpMagnus3;


  // Allocate Scratch Space
  NBTSqScratch_  = memManager_->template malloc<dcomplex>(NBT*NBT);
  NBTSqScratch2_ = memManager_->template malloc<dcomplex>(NBT*NBT);
  NBTSqScratch3_ = memManager_->template malloc<dcomplex>(NBT*NBT);



  // U**H Scalar
  UTransScalar_ = memManager_->template malloc<dcomplex>(NB*NB);

  // Saved PO Scalar
  POScalarSav_  = std::unique_ptr<ComplexMap>(
    new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

  POSav_.emplace_back(POScalarSav_.get());

  // Saved FO Scalar (1)
  if(needsFOSav) {
    FOScalarSav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    FOSav_.emplace_back(FOScalarSav_.get());
  }

  // Saved FO Scalar (2)
  if(needsFOSav2) {
    FOScalarSav2_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    FOSav2_.emplace_back(FOScalarSav2_.get());
  }

  // Saved FO Scalar (3)
  if(needsFOSav3) {
    FOScalarSav3_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    FOSav3_.emplace_back(FOScalarSav3_.get());
  }

  // Saved FO Scalar (4)
  if(needsFOSav4) {
    FOScalarSav4_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    FOSav4_.emplace_back(FOScalarSav4_.get());
  }

  // Saved FO Scalar (5)
  if(needsFOSav5) {
    FOScalarSav5_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    FOSav5_.emplace_back(FOScalarSav5_.get());
  }

  if(ssPropagator_->nTCS() == 2 or !ssPropagator_->isClosedShell){
    // U**H Mz
    UTransMz_= memManager_->template malloc<dcomplex>(NB*NB);
  
    // Saved PO Mz
    POMzSav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    POSav_.emplace_back(POMzSav_.get());

    // Saved FO Mz (1)
    if(needsFOSav) {
      FOMzSav_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav_.emplace_back(FOMzSav_.get());
    }

    // Saved FO Mz (2)
    if(needsFOSav2) {
      FOMzSav2_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav2_.emplace_back(FOMzSav2_.get());
    }

    // Saved FO Mz (3)
    if(needsFOSav3) {
      FOMzSav3_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav3_.emplace_back(FOMzSav3_.get());
    }


    // Saved FO Mz (4)
    if(needsFOSav4) {
      FOMzSav4_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav4_.emplace_back(FOMzSav4_.get());
    }


    // Saved FO Mz (5)
    if(needsFOSav5) {
      FOMzSav5_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav5_.emplace_back(FOMzSav5_.get());
    }


  }

  if(ssPropagator_->nTCS() == 2){ 
    // U**H Mx/My
    UTransMx_= memManager_->template malloc<dcomplex>(NB*NB);
    UTransMy_= memManager_->template malloc<dcomplex>(NB*NB);

    // Saved PO Mx/My
    POMySav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
    POMxSav_  = std::unique_ptr<ComplexMap>(
      new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

    POSav_.emplace_back(POMySav_.get());
    POSav_.emplace_back(POMxSav_.get());

    // Saved FO Mx/My (1)
    if(needsFOSav) {
      FOMySav_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
      FOMxSav_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav_.emplace_back(FOMySav_.get());
      FOSav_.emplace_back(FOMxSav_.get());
    }

    // Saved FO Mx/My (2)
    if(needsFOSav2) {
      FOMySav2_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
      FOMxSav2_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav2_.emplace_back(FOMySav2_.get());
      FOSav2_.emplace_back(FOMxSav2_.get());
    }

    // Saved FO Mx/My (3)
    if(needsFOSav3) {
      FOMySav3_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
      FOMxSav3_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav3_.emplace_back(FOMySav3_.get());
      FOSav3_.emplace_back(FOMxSav3_.get());
    }


    // Saved FO Mx/My (4)
    if(needsFOSav4) {
      FOMySav4_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
      FOMxSav4_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav4_.emplace_back(FOMySav4_.get());
      FOSav4_.emplace_back(FOMxSav4_.get());
    }


    // Saved FO Mx/My (5)
    if(needsFOSav5) {
      FOMySav5_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));
      FOMxSav5_  = std::unique_ptr<ComplexMap>(
        new ComplexMap(memManager_->template malloc<dcomplex>(NB*NB),NB,NB));

      FOSav5_.emplace_back(FOMySav5_.get());
      FOSav5_.emplace_back(FOMxSav5_.get());
    }


  }

  groundState_->transformOrthoMO();
}

