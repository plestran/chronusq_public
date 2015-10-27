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
#ifndef INCLUDED_FILEIO
#define INCLUDED_FILEIO
#define MAXBLOCK 1000
#include <global.h>
#include <cerr.h>
#include <tools.h>

/****************************/
/* Error Messages 1000-1999 */
/****************************/
namespace ChronusQ {
class FileIO {

  std::string  name_in;                 // name of the input file
  std::string  name_out;                // name of the output file
  std::string  name_scr;                // name of the scratch file
  std::string  name_restart;            // name of the restart file

  std::string  metaDataGroupPath;
  std::string  operatorGroupPath;
  std::string  SCFGroupPath;

  std::string  referenceMetaPath;
  std::string  scfMetaPath;

  // Operator Paths
  std::string  overlapPath  ;
  std::string  kineticPath  ;
  std::string  nucReplPath  ;
  std::string  coreHamPath  ;
  std::string  dipolePath   ;
  std::string  quadpolePath ;
  std::string  octupolePath ;

  // SCF Paths
  std::string  alphaSCFDenPath ;
  std::string  betaSCFDenPath  ;
  std::string  alphaMOPath     ;
  std::string  betaMOPath      ;


public:

  fstream in;                    // file handler of the input file
  fstream out;                   // file handler of the output file

  std::unique_ptr<H5::H5File> scr;
  std::unique_ptr<H5::H5File> restart;
  
  std::unique_ptr<H5::Group>  Meta;
  std::unique_ptr<H5::Group>  Operators;
  std::unique_ptr<H5::Group>  SCF;

  std::unique_ptr<H5::DataSet> referenceMeta;
  std::unique_ptr<H5::DataSet> scfMeta;

  std::unique_ptr<H5::DataSet> overlap;
  std::unique_ptr<H5::DataSet> kinetic;
  std::unique_ptr<H5::DataSet> nucRepl;
  std::unique_ptr<H5::DataSet> coreHam;

  std::unique_ptr<H5::DataSet> dipole;
  std::unique_ptr<H5::DataSet> quadpole;
  std::unique_ptr<H5::DataSet> octupole;

  std::unique_ptr<H5::DataSet> alphaSCFDen;
  std::unique_ptr<H5::DataSet> betaSCFDen;
  std::unique_ptr<H5::DataSet> alphaMO;
  std::unique_ptr<H5::DataSet> betaMO;

  template<typename T> struct metaData {
    T val;
    char desc[25];
  };
  std::unique_ptr<H5::CompType> metaDataTypeDouble;
  std::unique_ptr<H5::CompType> metaDataTypeInt;
  std::unique_ptr<H5::CompType> complexType;

  bool doRestart;

  

  // constructor and destructor
  FileIO(const std::string);
  ~FileIO() {
    if(in.is_open())  in.close();
    if(out.is_open()) out.close();
  }

  enum STDFILES {
    Overlap,
    Kinetic,
    NuclearRepulsion,
    CoreHamiltonian,
    Dipole,
    Quadrupole,
    Octupole,
    AlphaSCFDensity,
    BetaSCFDensity,
    AlphaMO,
    BetaMO
  };

  void iniH5Files();
  void iniCompType();
  void iniStdGroups();
  void iniStdOpFiles(int);
//template<typename T> void iniStdSCFFiles(bool,int);
//void iniStdSCFFiles(bool,int);
  void iniStdSCFFilesDouble(bool,int);
  void iniStdSCFFilesComplex(bool,int);

  // Python API
  void write(std::string);

};
} // namespace ChronusQ
#endif
