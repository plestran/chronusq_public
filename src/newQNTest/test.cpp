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
#include <response.h>
#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;
int main(int argc, char **argv){
  Molecule molecule;
  BasisSet basis;
  Controls controls;
  AOIntegrals aoints;
  MOIntegrals<double> moints;
  SingleSlater<double> singleSlater;
  Response<double> resp;
  FileIO fileio("test.inp","test.out");

  initCQ(argc,argv);
  controls.iniControls();
  fileio.iniH5Files();
  fileio.iniStdGroups();
  CQSetNumThreads(1);


  // Molecule Specification for Water
  molecule.setNAtoms(3);
  molecule.setCharge(0);
  molecule.setMultip(1);
  molecule.alloc(fileio.out);
  molecule.setIndex(0,HashAtom("O",0));
  molecule.setIndex(1,HashAtom("H",0));
  molecule.setIndex(2,HashAtom("H",0));
  molecule.setCart(0,0.000000000 ,-0.07579184359, 0.0);
  molecule.setCart(1,0.866811829 ,0.6014357793  ,0.0);
  molecule.setCart(2,-0.866811829, 0.6014357793 ,0.0);
  molecule.setNTotalE(10);
  basis.findBasisFile("sto3g");

/*
  // Molecule Specification ofr BH
  molecule.setNAtoms(2);
  molecule.setCharge(0);
  molecule.setMultip(1);
  molecule.alloc(fileio.out);
  molecule.setIndex(0,HashAtom("B",0));
  molecule.setIndex(1,HashAtom("H",0));
  molecule.setCart(0,0.0,0.0,0.0);
  molecule.setCart(1,0.0,0.0,1.232);
  molecule.setNTotalE(6);
  basis.findBasisFile("3-21G");
*/

  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(SingleSlater<double>::RHF);
  singleSlater.isClosedShell = true;

  basis.communicate(fileio);
  basis.parseGlobal();
  basis.constructLocal(&molecule);
  basis.makeMaps(1,&molecule);
  basis.renormShells();

  aoints.communicate(molecule,basis,fileio,controls);
  singleSlater.communicate(molecule,basis,aoints,fileio,controls);
  moints.communicate(molecule,basis,fileio,controls,aoints,singleSlater);

  aoints.initMeta();
  aoints.integralAlgorithm = AOIntegrals::INCORE;
  singleSlater.initMeta();
  singleSlater.genMethString();

  aoints.alloc();
  singleSlater.alloc();

  singleSlater.formGuess();
  singleSlater.formFock();
  singleSlater.computeEnergy();
  singleSlater.SCF();
  singleSlater.computeProperties();
  singleSlater.printProperties();
  
  moints.communicate(molecule,basis,fileio,controls,aoints,singleSlater);
  moints.initMeta();
  resp.communicate(singleSlater,moints,fileio);  
  resp.setMeth(RESPONSE_TYPE::RPA);
  //resp.doSA();
  resp.setNSek(3);
  //resp.doFull();
  resp.doResponse();
  finalizeCQ(); 
  return 0;
};
