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
#include <fileio.h>
using ChronusQ::FileIO;

namespace ChronusQ {
FileIO::FileIO(const std::string nm_input) {
  if(nm_input.empty()) 
    CErr("Fatal: Input File Required");

  this->doRestart = false;
  
  this->name_in = nm_input + ".inp";
  this->name_out = nm_input + ".out";
  this->name_scr = nm_input + ".scr";
  this->name_restart = nm_input + ".bin";

  this->operatorGroupPath = "/OPERATORS";
  this->SCFGroupPath      = "/SCF";

  this->overlapPath  = this->operatorGroupPath + "/OVERLAP";
  this->kineticPath  = this->operatorGroupPath + "/KINETIC";
  this->nucReplPath  = this->operatorGroupPath + "/NUCLEAR_REPULSION";
  this->coreHamPath  = this->operatorGroupPath + "/CORE_HAMILTONIAN";
  this->dipolePath   = this->operatorGroupPath + "/DIPOLE";
  this->quadpolePath = this->operatorGroupPath + "/QUADRUPOLE";
  this->octupolePath = this->operatorGroupPath + "/OCTUPOLE";

  this->alphaSCFDenPath = this->SCFGroupPath + "/ALPHA_DENSITY";
  this->betaSCFDenPath  = this->SCFGroupPath + "/BETA_DENSITY";
  this->alphaMOPath     = this->SCFGroupPath + "/ALPHA_MO";
  this->betaMOPath      = this->SCFGroupPath + "/BETA>MO";

  this->in.open(name_in,ios::in);
  this->out.open(name_out,ios::out);

  this->iniCompType();
};

void FileIO::iniH5Files(){
  if(doRestart) {
    this->scr = std::unique_ptr<H5::H5File>(
      new H5::H5File(this->name_scr,H5F_ACC_RDWR)
    );
    this->restart = std::unique_ptr<H5::H5File>(
      new H5::H5File(this->name_restart,H5F_ACC_RDWR)
    );
  } else {
    this->scr = std::unique_ptr<H5::H5File>(
      new H5::H5File(this->name_scr,H5F_ACC_TRUNC)
    );
    this->restart = std::unique_ptr<H5::H5File>(
      new H5::H5File(this->name_restart,H5F_ACC_TRUNC)
    );
  }
}

void FileIO::iniStdGroups(){
  if(doRestart){
    this->Operators = std::unique_ptr<H5::Group>(
      new H5::Group(this->restart->openGroup(operatorGroupPath))
    );
    this->SCF = std::unique_ptr<H5::Group>(
      new H5::Group(this->restart->openGroup(SCFGroupPath))
    );
  } else {
    this->Operators = std::unique_ptr<H5::Group>(
      new H5::Group(this->restart->createGroup(operatorGroupPath))
    );
    this->SCF = std::unique_ptr<H5::Group>(
      new H5::Group(this->restart->createGroup(SCFGroupPath))
    );
  }
}

void FileIO::iniStdOpFiles(int nBasis){
  hsize_t NBSq[] = {nBasis,nBasis};
  hsize_t dipoleDim[] = {3,nBasis,nBasis};
  hsize_t quadpoleDim[] = {6,nBasis,nBasis};
  hsize_t octpoleDim[] = {10,nBasis,nBasis};

  H5::DataSpace NBSqDataSpace(2,NBSq);
  H5::DataSpace DipoleDataSpace(3,dipoleDim);
  H5::DataSpace QuadrupoleDataSpace(3,quadpoleDim);
  H5::DataSpace OctupoleDataSpace(3,octpoleDim);

  if(this->doRestart) {
    this->overlap = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->overlapPath))
    );
 
    this->kinetic = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->kineticPath))
    );
 
    this->nucRepl = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->nucReplPath))
    );
 
    this->coreHam = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->coreHamPath))
    );
 
    this->dipole = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->dipolePath))
    );
 
    this->quadpole = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->quadpolePath))
    );
 
    this->octupole = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->octupolePath))
    );

  } else { 

    this->overlap = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->overlapPath,H5::PredType::NATIVE_DOUBLE,NBSqDataSpace
        )
      )
    );
 
    this->kinetic = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->kineticPath,H5::PredType::NATIVE_DOUBLE,NBSqDataSpace
        )
      )
    );
 
    this->nucRepl = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->nucReplPath,H5::PredType::NATIVE_DOUBLE,NBSqDataSpace
        )
      )
    );
 
    this->coreHam = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->coreHamPath,H5::PredType::NATIVE_DOUBLE,NBSqDataSpace
        )
      )
    );
 
    this->dipole = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->dipolePath,H5::PredType::NATIVE_DOUBLE,DipoleDataSpace
        )
      )
    );
 
    this->quadpole = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->quadpolePath,H5::PredType::NATIVE_DOUBLE,QuadrupoleDataSpace
        )
      )
    );
 
    this->octupole = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->octupolePath,H5::PredType::NATIVE_DOUBLE,OctupoleDataSpace
        )
      )
    );
  }
}

/*
template<>
void FileIO::iniStdSCFFiles<double>(bool allocBeta, int nBasis){
*/
//void FileIO::iniStdSCFFiles(bool allocBeta, int nBasis){
void FileIO::iniStdSCFFilesDouble(bool allocBeta, int nBasis){
  hsize_t NBSq[] = {nBasis,nBasis};
  H5::DataSpace NBSqDataSpace(2,NBSq);

  if(doRestart) {
    this->alphaSCFDen = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->alphaSCFDenPath))
    );
 
    this->alphaMO = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->alphaMOPath))
    );
 
    if(allocBeta) {
      this->betaSCFDen = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(this->restart->openDataSet(this->betaSCFDenPath))
      );
  
      this->betaMO = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(this->restart->openDataSet(this->betaMOPath))
      );
    }
 
  } else {

    this->alphaSCFDen = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->alphaSCFDenPath,H5::PredType::NATIVE_DOUBLE,NBSqDataSpace
        )
      )
    );
 
    this->alphaMO = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->alphaMOPath,H5::PredType::NATIVE_DOUBLE,NBSqDataSpace
        )
      )
    );
 
    if(allocBeta) {
      this->betaSCFDen = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(
          this->restart->createDataSet(
            this->betaSCFDenPath,H5::PredType::NATIVE_DOUBLE,NBSqDataSpace
          )
        )
      );
  
      this->betaMO = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(
          this->restart->createDataSet(
            this->betaMOPath,H5::PredType::NATIVE_DOUBLE,NBSqDataSpace
          )
        )
      );
    }

  }

}

void FileIO::iniStdSCFFilesComplex(bool allocBeta, int nBasis){
  hsize_t NBSq[] = {nBasis,nBasis};
  H5::DataSpace NBSqDataSpace(2,NBSq);

  if(this->doRestart) {
    this->alphaSCFDen = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->alphaSCFDenPath))
    );
 
    this->alphaMO = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(this->restart->openDataSet(this->alphaMOPath))
    );
 
    if(allocBeta) {
      this->betaSCFDen = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(this->restart->openDataSet(this->betaSCFDenPath))
      );
      this->betaMO = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(this->restart->openDataSet(this->betaMOPath))
      );
    }

  } else {

    this->alphaSCFDen = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->alphaSCFDenPath,*this->complexType,NBSqDataSpace
        )
      )
    );
 
    this->alphaMO = std::unique_ptr<H5::DataSet>(
      new H5::DataSet(
        this->restart->createDataSet(
          this->alphaMOPath,*this->complexType,NBSqDataSpace
        )
      )
    );
 
    if(allocBeta) {
 
      this->betaSCFDen = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(
          this->restart->createDataSet(
            this->betaSCFDenPath,*this->complexType,NBSqDataSpace
          )
        )
      );
  
      this->betaMO = std::unique_ptr<H5::DataSet>(
        new H5::DataSet(
          this->restart->createDataSet(
            this->betaMOPath,*this->complexType,NBSqDataSpace
          )
        )
      );
 
    }

  }
}

void FileIO::iniCompType(){
  typedef struct {
    double re;
    double im;
  } complex_t;

  this->complexType = std::unique_ptr<H5::CompType>(
    new H5::CompType(sizeof(complex_t))
  );
  this->complexType->insertMember(
    "RE",HOFFSET(complex_t,re),H5::PredType::NATIVE_DOUBLE
  );
  this->complexType->insertMember(
    "IM",HOFFSET(complex_t,im),H5::PredType::NATIVE_DOUBLE
  );


  this->metaDataTypeDouble = std::unique_ptr<H5::CompType>(
    new H5::CompType(sizeof(metaData<double>))
  );
  this->metaDataTypeDouble->insertMember(
    "VALUE",HOFFSET(metaData<double>,val),H5::PredType::NATIVE_DOUBLE
  );
  this->metaDataTypeDouble->insertMember(
    "DESC",HOFFSET(metaData<double>,desc),H5::StrType(H5::PredType::C_S1,25)
  );

  this->metaDataTypeInt = std::unique_ptr<H5::CompType>(
    new H5::CompType(sizeof(metaData<int>))
  );
  this->metaDataTypeInt->insertMember(
    "VALUE",HOFFSET(metaData<int>,val),H5::PredType::NATIVE_INT
  );
  this->metaDataTypeInt->insertMember(
    "DESC",HOFFSET(metaData<int>,desc),H5::StrType(H5::PredType::C_S1,25)
  );


}

}; // namespace ChronusQ
