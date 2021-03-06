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

  molecule.setCharge(0);
  molecule.setMultip(1);
  molecule.setNAtoms(2);

  molecule.alloc(fileio.out);
  molecule.setIndex(0,HashAtom("H",0));
  molecule.setIndex(1,HashAtom("H",0));
  molecule.setCart(0,0.0,0.0,0.0);
  molecule.setCart(1,0.0,0.0,1.4);
  molecule.setNTotalE(2);

  //molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(SingleSlater<double>::RHF);
  singleSlater.isClosedShell = true;

  basis.communicate(fileio);
  basis.findBasisFile("sto3g");
  basis.parseGlobal();
  basis.constructLocal(&molecule);
  basis.makeMaps(&molecule);
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

  GaussChebyshev1stGridInf grid(70,0,1);
  LebedevGrid ang(302);
  grid.genGrid();
  ang.genGrid();
//grid.atomGrid(elements[molecule.index(0)].sradius);
  grid.atomGrid(0);

  TwoDGrid tgrid(302*70,&grid,&ang);
  tgrid.centerGrid(0,0,0);
  
  double *x  = new double[302*70];
  for(auto i = 0; i < 70*302; i++){
    x[i] = bg::get<0>(tgrid.gridPtCart(i));
    cout << bg::get<0>(tgrid.gridPtCart(i)) << "\t" << tgrid.getweightsGrid(i) << endl;
  }
  RealMap X(x,302,70);
  RealMap W(tgrid.weightsGrid(),302,70);
  prettyPrint(cout,X,"X");
  prettyPrint(cout,W,"W");
  delete [] x;
  
  finalizeCQ(); 
  return 0;
};

