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
#define EIGEN_DONT_PARALLELIZE
#include <response.h>
#include <workers.h>
#include <pythonapi.h>


using namespace ChronusQ;


Derivatives computeNAC2pt1D(Molecule&,Molecule&,Molecule&,FileIO&,FileIO&, 
  FileIO&,Controls&,std::string&,RESPONSE_TYPE,int,double);


void cartDiff(std::string &XYZFName, std::string &basisName, 
  RESPONSE_TYPE respType, int charge, int nFreq){

  initCQ();
  Controls controls;
  CQSetNumThreads(3);
  controls.iniControls();


  std::string PREFIX = XYZFName;

  Molecule geom_0;
  FileIO fileio_0( PREFIX+"_0.inp", PREFIX+".out_0", PREFIX+".rst_0" );

  fileio_0.iniH5Files();

  fileio_0.iniStdGroups();

  std::ifstream XYZ0File(XYZFName);
  if(!XYZ0File.good()) CErr("Cannot find Specified XYZ File");



  int nAtoms(0);
  std::string inLine;
  getline(XYZ0File,inLine);
  nAtoms = std::atoi(inLine.c_str());
  geom_0.setNAtoms(nAtoms);
  geom_0.setMultip(1);
  geom_0.alloc(fileio_0.out);
  // Skip line
  getline(XYZ0File,inLine);

  int nElec = -charge;
  for(auto iAtm = 0; iAtm < nAtoms; iAtm++){
    getline(XYZ0File,inLine);

    std::istringstream iss(inLine);

    std::vector<std::string> tokens{
      std::istream_iterator<std::string>{iss},
      std::istream_iterator<std::string>{ }
    };

    auto INDX = HashAtom(tokens[0],0);
    nElec += atom[INDX].atomicNumber;

    geom_0.setIndex(iAtm,INDX);


    geom_0.setCart(iAtm,
      std::atof(tokens[1].c_str()),
      std::atof(tokens[2].c_str()),
      std::atof(tokens[3].c_str())
    );
  }

  geom_0.setCharge(charge);
  geom_0.setNTotalE(nElec);
  geom_0.convBohr();
  geom_0.computeNucRep();
  geom_0.computeRij();
  geom_0.computeI();

  double step = 0.00001;
  std::vector<Derivatives> derv;
  for(auto iAtm = 0, IX = 0; iAtm < nAtoms; iAtm++)
  for(auto iXYZ = 0; iXYZ < 3     ; iXYZ++, IX++) {
    FileIO fileio_p(PREFIX+"_"+std::to_string(IX)+"_p.inp",
      PREFIX+"_"+std::to_string(IX)+".out_p",
      PREFIX+"_"+std::to_string(IX)+".rst_p");
    FileIO fileio_m(PREFIX+"_"+std::to_string(IX)+"_m.inp",
      PREFIX+"_"+std::to_string(IX)+".out_m",
      PREFIX+"_"+std::to_string(IX)+".rst_m");
 
    fileio_p.iniH5Files();
    fileio_m.iniH5Files();
 
    fileio_p.iniStdGroups();
    fileio_m.iniStdGroups();

    Molecule geom_p, geom_m;

    geom_p.setNAtoms(nAtoms);
    geom_p.setMultip(1);
    geom_p.alloc(fileio_p.out);
    geom_m.setNAtoms(nAtoms);
    geom_m.setMultip(1);
    geom_m.alloc(fileio_m.out);

    geom_p.setCharge(charge); geom_p.setNTotalE(nElec);
    geom_m.setCharge(charge); geom_m.setNTotalE(nElec);

    for(auto jAtm = 0; jAtm < nAtoms; jAtm++){
      geom_p.index(jAtm) = geom_0.index(jAtm);
      geom_m.index(jAtm) = geom_0.index(jAtm);
    }
    (*geom_p.cart()) = (*geom_0.cart());
    (*geom_m.cart()) = (*geom_0.cart());

    (*geom_p.cart())(iXYZ,iAtm) += step;
    (*geom_m.cart())(iXYZ,iAtm) -= step;

    geom_p.computeNucRep();
    geom_p.computeRij();
    geom_p.computeI();
    geom_m.computeNucRep();
    geom_m.computeRij();
    geom_m.computeI();

    derv.push_back(computeNAC2pt1D(geom_0,geom_p,geom_m,fileio_0,fileio_p,
      fileio_m,controls,basisName,respType,nFreq,step));
  }

  cout << "Summary" << endl;
  cout << endl; 

  cout << "GS Gradient" << endl;
  cout << std::setw(15) << "Atom Number";
  cout << std::setw(15) << "Atomic Type";
  cout << std::setw(15) << "X";
  cout << std::setw(15) << "Y";
  cout << std::setw(15) << "Z";
  cout << endl;
  for(auto iAtm = 0, IX = 0; iAtm < nAtoms; iAtm++){
    cout << std::setw(15) << iAtm;
    cout << std::setw(15) << atom[geom_0.index(iAtm)].atomicNumber;
    for(auto iXYZ = 0; iXYZ < 3     ; iXYZ++, IX++) {
      cout << std::setw(15) << derv[IX].GS_GRAD;
    }
    cout << endl;
  }


  cout << endl << endl;
  cout << "ES Gradients" << endl;
  for(auto iSt = 0; iSt < nFreq; iSt++){
    cout << "  Excited State " << iSt + 1<< endl;
    cout << std::setw(15) << "Atom Number";
    cout << std::setw(15) << "Atomic Type";
    cout << std::setw(15) << "X";
    cout << std::setw(15) << "Y";
    cout << std::setw(15) << "Z";
    cout << endl;
    for(auto iAtm = 0, IX = 0; iAtm < nAtoms; iAtm++){
      cout << std::setw(15) << iAtm;
      cout << std::setw(15) << atom[geom_0.index(iAtm)].atomicNumber;
      for(auto iXYZ = 0; iXYZ < 3     ; iXYZ++, IX++) {
        cout << std::setw(15) << derv[IX].ES_GRAD(iSt) + derv[IX].GS_GRAD;
      }
      cout << endl;
    }
  }



  cout << endl << endl;
  cout << "ES -> GS Non Adiabatic Couplings" << endl;
  for(auto iSt = 0; iSt < nFreq; iSt++){
    cout << "  " << iSt + 1 << " - > 0" << endl;
    cout << std::setw(15) << "Atom Number";
    cout << std::setw(15) << "Atomic Type";
    cout << std::setw(15) << "X";
    cout << std::setw(15) << "Y";
    cout << std::setw(15) << "Z";
    cout << endl;
    for(auto iAtm = 0, IX = 0; iAtm < nAtoms; iAtm++){
      cout << std::setw(15) << iAtm;
      cout << std::setw(15) << atom[geom_0.index(iAtm)].atomicNumber;
      for(auto iXYZ = 0; iXYZ < 3     ; iXYZ++, IX++) {
        cout << std::setw(15) << derv[IX].ES_GS_NACME(iSt);
      }
      cout << endl;
    }
  }

  cout << endl << endl;
  cout << "ES -> ES Non Adiabatic Couplings" << endl;
  for(auto iSt = 0; iSt < nFreq; iSt++)
  for(auto jSt = 0; jSt < nFreq; jSt++){
    cout << "  " << iSt + 1 << " - > " << jSt + 1 << endl;
    cout << std::setw(15) << "Atom Number";
    cout << std::setw(15) << "Atomic Type";
    cout << std::setw(15) << "X";
    cout << std::setw(15) << "Y";
    cout << std::setw(15) << "Z";
    cout << endl;
    for(auto iAtm = 0, IX = 0; iAtm < nAtoms; iAtm++){
      cout << std::setw(15) << iAtm;
      cout << std::setw(15) << atom[geom_0.index(iAtm)].atomicNumber;
      for(auto iXYZ = 0; iXYZ < 3     ; iXYZ++, IX++) {
        cout << std::setw(15) << derv[IX].ES_ES_NACME(iSt,jSt);
      }
      cout << endl;
    }
  }

};
