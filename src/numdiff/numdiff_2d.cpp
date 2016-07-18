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
#define EIGEN_DONT_PARALLELIZE
#include <response.h>
#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;
std::string genFName(double,double,int,int,std::string);
Eigen::VectorXd checkPhase(double*,double*,int,int);
RealMatrix genSpx(BasisSet&,BasisSet&);

std::vector<Eigen::VectorXd> ES_GS_NACME_CIS(int,bool,int,int,double,
  double,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&);
std::vector<Eigen::VectorXd> GS_ES_NACME_CIS(int,bool,int,int,double,
  double,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&);
std::vector<RealMatrix> ES_ES_NACME_CIS(int,bool,int,int,double,
  double,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&);
std::vector<RealMatrix> ES_ES_NACME_PPTDA(int,bool,int,int,double,
  double,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&);




void twoDScan(std::vector<double>& scanX, 
              std::vector<double>& scanY,
              std::string &basisName,
              RESPONSE_TYPE respType,
              int charge, int nFreq,
              int px, int py,
              std::string &PREFIX){

  std::string BINPREFIX = "BIN";
  std::string XYZPREFIX = "XYZ";

  initCQ();
  Controls controls;
  CQSetNumThreads(3);
  controls.iniControls();

  bool debug = true;

  std::string GS_ENERGY_FILE_NAME = PREFIX + "_gsEnergy.dat";
  std::string ES_ENERGY_FILE_NAME = PREFIX + "_esEnergy.dat";
  std::string GS_GRADIENT_DX_FILE_NAME = PREFIX + "_gsGrad_IX_1.dat";
  std::string GS_GRADIENT_DY_FILE_NAME = PREFIX + "_gsGrad_IX_2.dat";
  std::string ES_GRADIENT_DX_FILE_NAME = PREFIX + "_esGrad_IX_1.dat";
  std::string ES_GRADIENT_DY_FILE_NAME = PREFIX + "_esGrad_IX_2.dat";
  std::string ES_GS_NACME_DX_FILE_NAME = PREFIX + "_esgsNAC_IX_1.dat";
  std::string ES_GS_NACME_DY_FILE_NAME = PREFIX + "_esgsNAC_IX_2.dat";
  std::string ES_ES_NACME_DX_FILE_NAME = PREFIX + "_esesNAC_IX_1.dat";
  std::string ES_ES_NACME_DY_FILE_NAME = PREFIX + "_esesNAC_IX_2.dat";

  std::ofstream GS_ENERGY_FILE      (GS_ENERGY_FILE_NAME      ); 
  std::ofstream ES_ENERGY_FILE      (ES_ENERGY_FILE_NAME      );
  std::ofstream GS_GRADIENT_DX_FILE (GS_GRADIENT_DX_FILE_NAME );
  std::ofstream GS_GRADIENT_DY_FILE (GS_GRADIENT_DY_FILE_NAME );
  std::ofstream ES_GRADIENT_DX_FILE (ES_GRADIENT_DX_FILE_NAME );
  std::ofstream ES_GRADIENT_DY_FILE (ES_GRADIENT_DY_FILE_NAME );
  std::ofstream ES_GS_NACME_DX_FILE (ES_GS_NACME_DX_FILE_NAME );
  std::ofstream ES_GS_NACME_DY_FILE (ES_GS_NACME_DY_FILE_NAME );
  std::ofstream ES_ES_NACME_DX_FILE (ES_ES_NACME_DX_FILE_NAME );
  std::ofstream ES_ES_NACME_DY_FILE (ES_ES_NACME_DY_FILE_NAME );

  GS_ENERGY_FILE << std::setprecision(8) << std::fixed;
  ES_ENERGY_FILE << std::setprecision(8) << std::fixed;
  GS_GRADIENT_DX_FILE  << std::setprecision(8)    << std::fixed;
  GS_GRADIENT_DY_FILE  << std::setprecision(8)    << std::fixed;
  ES_GRADIENT_DX_FILE  << std::setprecision(8)    << std::fixed;
  ES_GRADIENT_DY_FILE  << std::setprecision(8)    << std::fixed;
  ES_GS_NACME_DX_FILE  << std::setprecision(8)    << std::fixed;
  ES_GS_NACME_DY_FILE  << std::setprecision(8)    << std::fixed;
  ES_ES_NACME_DX_FILE  << std::setprecision(8)    << std::fixed;
  ES_ES_NACME_DY_FILE  << std::setprecision(8)    << std::fixed;


  for(auto IX = 1; IX < (scanX.size()-1); IX++)
  for(auto JY = 1; JY < (scanY.size()-1); JY++){
    // Output header
    cout << bannerTop << endl;
    cout << "Starting Numerical Differentiation for:" << endl;
    cout << "  IX = " << IX << "  JY = " << JY << endl;
    FileIO fileio_00(PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+"_00.inp",PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".out_00",
      PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".rst_00");
    FileIO fileio_p0(PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+"_p0.inp",PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".out_p0",
      PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".rst_p0");
    FileIO fileio_m0(PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+"_m0.inp",PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".out_m0",
      PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".rst_m0");
    FileIO fileio_0p(PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+"_0p.inp",PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".out_0p",
      PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".rst_0p");
    FileIO fileio_0m(PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+"_0m.inp",PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".out_0m",
      PREFIX+"_"+std::to_string(IX)+"_"+std::to_string(JY)+".rst_0m");

    fileio_00.iniH5Files();
    fileio_p0.iniH5Files();
    fileio_m0.iniH5Files();
    fileio_0p.iniH5Files();
    fileio_0m.iniH5Files();

    fileio_00.iniStdGroups();
    fileio_p0.iniStdGroups();
    fileio_m0.iniStdGroups();
    fileio_0p.iniStdGroups();
    fileio_0m.iniStdGroups();

    // Get coordinate information
    double X = scanX[IX];
    double Y = scanY[JY];


    double XP = scanX[IX+1];
    double XM = scanX[IX-1];
    double YP = scanY[JY+1];
    double YM = scanY[JY-1];


    // Output coordinate information
    cout << endl;
    cout << "  X      = " << X << endl;
    cout << "  X + DX = " << XP << endl;
    cout << "  X - DX = " << XM << endl;
    cout << endl;
    cout << "  Y      = " << Y << endl;
    cout << "  Y + DY = " << YP << endl;
    cout << "  Y - DY = " << YM << endl;


    // Generate file names baed on convention
    std::string fBaseName00 = genFName(X ,Y ,px,py,PREFIX);
    std::string fBaseNamep0 = genFName(XP,Y ,px,py,PREFIX);
    std::string fBaseNamem0 = genFName(XM,Y ,px,py,PREFIX);
    std::string fBaseName0p = genFName(X ,YP,px,py,PREFIX);
    std::string fBaseName0m = genFName(X ,YM,px,py,PREFIX);

    std::string fNameGeom_00 = XYZPREFIX+"/"+fBaseName00 + ".xyz";
    std::string fNameGeom_p0 = XYZPREFIX+"/"+fBaseNamep0 + ".xyz";
    std::string fNameGeom_0p = XYZPREFIX+"/"+fBaseName0p + ".xyz";
    std::string fNameGeom_m0 = XYZPREFIX+"/"+fBaseNamem0 + ".xyz";
    std::string fNameGeom_0m = XYZPREFIX+"/"+fBaseName0m + ".xyz";

    // Output file names
    cout << endl;
    cout << "  Base File Name (X,Y):    " << fBaseName00 << endl;
    cout << "  Base File Name (X+DX,Y): " << fBaseNamep0 << endl;
    cout << "  Base File Name (X-DX,Y): " << fBaseNamem0 << endl;
    cout << "  Base File Name (X,Y+DY): " << fBaseName0p << endl;
    cout << "  Base File Name (X,Y-DY): " << fBaseName0m << endl;

    cout << endl;
    cout << "  XYZ File Name (X,Y):    " << fNameGeom_00 << endl;
    cout << "  XYZ File Name (X+DX,Y): " << fNameGeom_p0 << endl;
    cout << "  XYZ File Name (X-DX,Y): " << fNameGeom_m0 << endl;
    cout << "  XYZ File Name (X,Y+DY): " << fNameGeom_0p << endl;
    cout << "  XYZ File Name (X,Y-DY): " << fNameGeom_0m << endl;

    

    // Create CQ Objects
    Molecule geom_00;
    Molecule geom_p0;
    Molecule geom_m0;
    Molecule geom_0p;
    Molecule geom_0m;

    BasisSet basis00;
    BasisSet basisp0;
    BasisSet basism0;
    BasisSet basis0p;
    BasisSet basis0m;

    AOIntegrals aoints_00;
    AOIntegrals aoints_p0;
    AOIntegrals aoints_m0;
    AOIntegrals aoints_0p;
    AOIntegrals aoints_0m;

    SingleSlater<double> ss_00;
    SingleSlater<double> ss_p0;
    SingleSlater<double> ss_m0;
    SingleSlater<double> ss_0p;
    SingleSlater<double> ss_0m;

    MOIntegrals<double> moints_00;
    MOIntegrals<double> moints_p0;
    MOIntegrals<double> moints_m0;
    MOIntegrals<double> moints_0p;
    MOIntegrals<double> moints_0m;

    Response<double> resp_00;
    Response<double> resp_p0;
    Response<double> resp_m0;
    Response<double> resp_0p;
    Response<double> resp_0m;


    // Set up Geometries
    geom_00.setMultip(1);
    geom_p0.setMultip(1);
    geom_m0.setMultip(1);
    geom_0p.setMultip(1);
    geom_0m.setMultip(1);

    std::ifstream XYZ00File(fNameGeom_00);
    std::ifstream XYZp0File(fNameGeom_p0);
    std::ifstream XYZm0File(fNameGeom_m0);
    std::ifstream XYZ0pFile(fNameGeom_0p);
    std::ifstream XYZ0mFile(fNameGeom_0m);

    if(!XYZ00File.good()) CErr("Cannot find XYZFile00");
    if(!XYZp0File.good()) CErr("Cannot find XYZFilep0");
    if(!XYZm0File.good()) CErr("Cannot find XYZFilem0");
    if(!XYZ0pFile.good()) CErr("Cannot find XYZFile0p");
    if(!XYZ0mFile.good()) CErr("Cannot find XYZFile0m");

    // Get number of atoms
    int nAtoms(0);
    std::string inLine_00;
    std::string inLine_p0;
    std::string inLine_m0;
    std::string inLine_0p;
    std::string inLine_0m;

    getline(XYZ00File,inLine_00);
    getline(XYZp0File,inLine_p0);
    getline(XYZm0File,inLine_m0);
    getline(XYZ0pFile,inLine_0p);
    getline(XYZ0mFile,inLine_0m);

    nAtoms = std::atoi(inLine_00.c_str());
    geom_00.setNAtoms(nAtoms);
    geom_p0.setNAtoms(nAtoms);
    geom_m0.setNAtoms(nAtoms);
    geom_0p.setNAtoms(nAtoms);
    geom_0m.setNAtoms(nAtoms);

    // Allocate Molecule objects
    geom_00.alloc(fileio_00.out);
    geom_p0.alloc(fileio_p0.out);
    geom_m0.alloc(fileio_m0.out);
    geom_0p.alloc(fileio_0p.out);
    geom_0m.alloc(fileio_0m.out);

    // Skip line
    getline(XYZ00File,inLine_00);
    getline(XYZp0File,inLine_p0);
    getline(XYZm0File,inLine_m0);
    getline(XYZ0pFile,inLine_0p);
    getline(XYZ0mFile,inLine_0m);

    int nElec = -charge;
    for(auto iAtm = 0; iAtm < nAtoms; iAtm++){
      getline(XYZ00File,inLine_00);
      getline(XYZp0File,inLine_p0);
      getline(XYZm0File,inLine_m0);
      getline(XYZ0pFile,inLine_0p);
      getline(XYZ0mFile,inLine_0m);

      std::istringstream iss_00(inLine_00);
      std::istringstream iss_p0(inLine_p0);
      std::istringstream iss_m0(inLine_m0);
      std::istringstream iss_0p(inLine_0p);
      std::istringstream iss_0m(inLine_0m);

      std::vector<std::string> tokens_00{
        std::istream_iterator<std::string>{iss_00},
        std::istream_iterator<std::string>{ }
      };
      std::vector<std::string> tokens_p0{
        std::istream_iterator<std::string>{iss_p0},
        std::istream_iterator<std::string>{ }
      };
      std::vector<std::string> tokens_m0{
        std::istream_iterator<std::string>{iss_m0},
        std::istream_iterator<std::string>{ }
      };
      std::vector<std::string> tokens_0p{
        std::istream_iterator<std::string>{iss_0p},
        std::istream_iterator<std::string>{ }
      };
      std::vector<std::string> tokens_0m{
        std::istream_iterator<std::string>{iss_0m},
        std::istream_iterator<std::string>{ }
      };

      auto INDX = HashAtom(tokens_00[0],0);
      nElec += atom[INDX].atomicNumber;

      geom_00.setIndex(iAtm,INDX);
      geom_p0.setIndex(iAtm,INDX);
      geom_m0.setIndex(iAtm,INDX);
      geom_0p.setIndex(iAtm,INDX);
      geom_0m.setIndex(iAtm,INDX);


      geom_00.setCart(iAtm,
        std::atof(tokens_00[1].c_str()),
        std::atof(tokens_00[2].c_str()),
        std::atof(tokens_00[3].c_str())
      );
      geom_p0.setCart(iAtm,
        std::atof(tokens_p0[1].c_str()),
        std::atof(tokens_p0[2].c_str()),
        std::atof(tokens_p0[3].c_str())
      );
      geom_m0.setCart(iAtm,
        std::atof(tokens_m0[1].c_str()),
        std::atof(tokens_m0[2].c_str()),
        std::atof(tokens_m0[3].c_str())
      );
      geom_0p.setCart(iAtm,
        std::atof(tokens_0p[1].c_str()),
        std::atof(tokens_0p[2].c_str()),
        std::atof(tokens_0p[3].c_str())
      );
      geom_0m.setCart(iAtm,
        std::atof(tokens_0m[1].c_str()),
        std::atof(tokens_0m[2].c_str()),
        std::atof(tokens_0m[3].c_str())
      );
    }

    geom_00.setCharge(charge);
    geom_p0.setCharge(charge);
    geom_m0.setCharge(charge);
    geom_0p.setCharge(charge);
    geom_0m.setCharge(charge);

    geom_00.setNTotalE(nElec);
    geom_p0.setNTotalE(nElec);
    geom_m0.setNTotalE(nElec);
    geom_0p.setNTotalE(nElec);
    geom_0m.setNTotalE(nElec);

    geom_00.convBohr();
    geom_p0.convBohr();
    geom_m0.convBohr();
    geom_0p.convBohr();
    geom_0m.convBohr();

    geom_00.computeNucRep();
    geom_p0.computeNucRep();
    geom_m0.computeNucRep();
    geom_0p.computeNucRep();
    geom_0m.computeNucRep();

    geom_00.computeRij();
    geom_p0.computeRij();
    geom_m0.computeRij();
    geom_0p.computeRij();
    geom_0m.computeRij();

    geom_00.computeI();
    geom_p0.computeI();
    geom_m0.computeI();
    geom_0p.computeI();
    geom_0m.computeI();

    if(debug){
      cout << "GEOMETRY (X,Y):" << endl;
      geom_00.printInfo();
      cout << "GEOMETRY (X+DX,Y):" << endl;
      geom_p0.printInfo();
      cout << "GEOMETRY (X-DX,Y):" << endl;
      geom_m0.printInfo();
      cout << "GEOMETRY (X,Y+DY):" << endl;
      geom_0p.printInfo();
      cout << "GEOMETRY (X,Y-DY):" << endl;
      geom_0m.printInfo();
    }

    basis00.findBasisFile(basisName);
    basisp0.findBasisFile(basisName);
    basism0.findBasisFile(basisName);
    basis0p.findBasisFile(basisName);
    basis0m.findBasisFile(basisName);

    basis00.communicate(fileio_00);
    basisp0.communicate(fileio_00);
    basism0.communicate(fileio_00);
    basis0p.communicate(fileio_00);
    basis0m.communicate(fileio_00);

    basis00.parseGlobal();
    basisp0.parseGlobal();
    basism0.parseGlobal();
    basis0p.parseGlobal();
    basis0m.parseGlobal();

    basis00.constructLocal(&geom_00);
    basisp0.constructLocal(&geom_p0);
    basism0.constructLocal(&geom_m0);
    basis0p.constructLocal(&geom_0p);
    basis0m.constructLocal(&geom_0m);
   
    basis00.makeMaps(&geom_00);
    basisp0.makeMaps(&geom_p0);
    basism0.makeMaps(&geom_m0);
    basis0p.makeMaps(&geom_0p);
    basis0m.makeMaps(&geom_0m);


    RealMatrix S_00_p0 = genSpx(basis00,basisp0);
    RealMatrix S_00_0p = genSpx(basis00,basis0p);
    RealMatrix S_00_m0 = genSpx(basis00,basism0);
    RealMatrix S_00_0m = genSpx(basis00,basis0m);

    RealMatrix S_00_00,S_p0_p0,S_m0_m0,S_0p_0p,S_0m_0m;

    if(debug){
      S_00_00 = genSpx(basis00,basis00);
      S_p0_p0 = genSpx(basisp0,basisp0);
      S_m0_m0 = genSpx(basism0,basism0);
      S_0p_0p = genSpx(basis0p,basis0p);
      S_0m_0m = genSpx(basis0m,basis0m);
    }

    aoints_00.communicate(geom_00,basis00,fileio_00,controls);
    aoints_p0.communicate(geom_p0,basisp0,fileio_p0,controls);
    aoints_m0.communicate(geom_m0,basism0,fileio_m0,controls);
    aoints_0p.communicate(geom_0p,basis0p,fileio_0p,controls);
    aoints_0m.communicate(geom_0m,basis0m,fileio_0m,controls);

    ss_00.communicate(geom_00,basis00,aoints_00,fileio_00,controls);
    ss_p0.communicate(geom_p0,basisp0,aoints_p0,fileio_p0,controls);
    ss_m0.communicate(geom_m0,basism0,aoints_m0,fileio_m0,controls);
    ss_0p.communicate(geom_0p,basis0p,aoints_0p,fileio_0p,controls);
    ss_0m.communicate(geom_0m,basis0m,aoints_0m,fileio_0m,controls);

    moints_00.communicate(geom_00,basis00,fileio_00,controls,
      aoints_00, ss_00);
    moints_p0.communicate(geom_p0,basisp0,fileio_p0,controls,
      aoints_p0, ss_p0);
    moints_m0.communicate(geom_m0,basism0,fileio_m0,controls,
      aoints_m0, ss_m0);
    moints_0p.communicate(geom_0p,basis0p,fileio_0p,controls,
      aoints_0p, ss_0p);
    moints_0m.communicate(geom_0m,basis0m,fileio_0m,controls,
      aoints_0m, ss_0m);
    
    aoints_00.initMeta();
    aoints_p0.initMeta();
    aoints_m0.initMeta();
    aoints_0p.initMeta();
    aoints_0m.initMeta();

    aoints_00.integralAlgorithm = AOIntegrals::INCORE;
    aoints_p0.integralAlgorithm = AOIntegrals::INCORE;
    aoints_m0.integralAlgorithm = AOIntegrals::INCORE;
    aoints_0p.integralAlgorithm = AOIntegrals::INCORE;
    aoints_0m.integralAlgorithm = AOIntegrals::INCORE;

    ss_00.setRef(SingleSlater<double>::RHF);
    ss_p0.setRef(SingleSlater<double>::RHF);
    ss_m0.setRef(SingleSlater<double>::RHF);
    ss_0p.setRef(SingleSlater<double>::RHF);
    ss_0m.setRef(SingleSlater<double>::RHF);

    ss_00.isClosedShell = true;
    ss_p0.isClosedShell = true;
    ss_m0.isClosedShell = true;
    ss_0p.isClosedShell = true;
    ss_0m.isClosedShell = true;

    ss_00.initMeta();
    ss_p0.initMeta();
    ss_m0.initMeta();
    ss_0p.initMeta();
    ss_0m.initMeta();

    ss_00.genMethString();
    ss_p0.genMethString();
    ss_m0.genMethString();
    ss_0p.genMethString();
    ss_0m.genMethString();



    aoints_00.alloc();
    aoints_p0.alloc();
    aoints_m0.alloc();
    aoints_0p.alloc();
    aoints_0m.alloc();

    ss_00.alloc();
    ss_p0.alloc();
    ss_m0.alloc();
    ss_0p.alloc();
    ss_0m.alloc();

    // SCF (X,Y)
    cout << "  Performing SCF (X,Y)" << endl;
    ss_00.formGuess();
    ss_00.formFock();
    ss_00.computeEnergy();
    ss_00.SCF();
    ss_00.computeProperties();
    ss_00.printProperties();

    // SCF (X+DX,Y)
    cout << "  Performing SCF (X+DX,Y)" << endl;
    ss_p0.formGuess();
    ss_p0.formFock();
    ss_p0.computeEnergy();
    ss_p0.SCF();
    ss_p0.computeProperties();
    ss_p0.printProperties();

    // SCF (X-DX,Y)
    cout << "  Performing SCF (X-DX,Y)" << endl;
    ss_m0.formGuess();
    ss_m0.formFock();
    ss_m0.computeEnergy();
    ss_m0.SCF();
    ss_m0.computeProperties();
    ss_m0.printProperties();

    // SCF (X,Y+DY)
    cout << "  Performing SCF (X,Y+DY)" << endl;
    ss_0p.formGuess();
    ss_0p.formFock();
    ss_0p.computeEnergy();
    ss_0p.SCF();
    ss_0p.computeProperties();
    ss_0p.printProperties();

    // SCF (X,Y-DY)
    cout << "  Performing SCF (X,Y-DY)" << endl;
    ss_0m.formGuess();
    ss_0m.formFock();
    ss_0m.computeEnergy();
    ss_0m.SCF();
    ss_0m.computeProperties();
    ss_0m.printProperties();


    if(debug) {
      cout << endl;
      cout << "  Checking | C - C' | Before Phase Check:" << endl;
      
      cout << "  | C(X,Y) - C(X+DX,Y) | = " 
           << diffNorm((*ss_00.moA()),(*ss_p0.moA())) 
           << endl;  
      cout << "  | C(X,Y) - C(X-DX,Y) | = " 
           << diffNorm((*ss_00.moA()),(*ss_m0.moA())) 
           << endl;  
      cout << "  | C(X,Y) - C(X,Y+DY) | = " 
           << diffNorm((*ss_00.moA()),(*ss_0p.moA())) 
           << endl;  
      cout << "  | C(X,Y) - C(X,Y-DY) | = " 
           << diffNorm((*ss_00.moA()),(*ss_0m.moA())) 
           << endl;  
    }



    cout << endl << "  Performing Phase Check on MOs" << endl;

    RealMatrix SMO_00_00,SMO_p0_p0,SMO_m0_m0,SMO_0p_0p,SMO_0m_0m;

    if(debug){
    // < p (X,Y) | S[(X,Y),(X,Y)] | q(X,Y) >
     SMO_00_00 = ss_00.moA()->transpose() * S_00_00 * (*ss_00.moA());
    // < p (X+DX,Y) | S[(X+DX,Y),(X+DX,Y)] | q(X+DX,Y)(* >
     SMO_p0_p0 = ss_p0.moA()->transpose() * S_p0_p0 * (*ss_p0.moA());
    // < p (X-DX,Y) | S[(X-DX,Y),(X-DX,Y)] | q(X-DX,Y)(* >
     SMO_m0_m0 = ss_m0.moA()->transpose() * S_m0_m0 * (*ss_m0.moA());
    // < p (X,Y+DY) | S[(X,Y+DY),(X,Y+DY)] | q(X,Y+DY)(* >
     SMO_0p_0p = ss_0p.moA()->transpose() * S_0p_0p * (*ss_0p.moA());
    // < p (X,Y-DY) | S[(X,Y-DY),(X,Y-DY)] | q(X,Y-DY)(* >
     SMO_0m_0m = ss_0m.moA()->transpose() * S_0m_0m * (*ss_0m.moA());
    }

    // < p (X,Y) | S[(X,Y),(X+DX,Y)] | q(X+DX,Y) >
    RealMatrix SMO_00_p0 = 
      ss_00.moA()->transpose() * S_00_p0 * (*ss_p0.moA());
    // < p (X,Y) | S[(X,Y),(X-DX,Y)] | q(X-DX,Y) >
    RealMatrix SMO_00_m0 = 
      ss_00.moA()->transpose() * S_00_m0 * (*ss_m0.moA());

    // < p (X,Y+DY) | S[(X,Y),(X,Y+DY)] | q(X,Y+DY) >
    RealMatrix SMO_00_0p = 
      ss_00.moA()->transpose() * S_00_0p * (*ss_0p.moA());
    // < p (X,Y-DY) | S[(X,Y),(X,Y-DY)] | q(X,Y-DY) >
    RealMatrix SMO_00_0m = 
      ss_00.moA()->transpose() * S_00_0m * (*ss_0m.moA());




/*
    Eigen::VectorXd Trans_mo_0p;
    Eigen::VectorXd Trans_mo_0m;
    Eigen::VectorXd Trans_mo_p0;
    Eigen::VectorXd Trans_mo_m0;
    Trans_mo_p0 = checkPhase(ss_00.moA()->data(),
      ss_p0.moA()->data(),basis00.nBasis(),basis00.nBasis());
    Trans_mo_m0 = checkPhase(ss_00.moA()->data(),
      ss_m0.moA()->data(),basis00.nBasis(),basis00.nBasis());
    Trans_mo_0p = checkPhase(ss_00.moA()->data(),
      ss_0p.moA()->data(),basis00.nBasis(),basis00.nBasis());
    Trans_mo_0m = checkPhase(ss_00.moA()->data(),
      ss_0m.moA()->data(),basis00.nBasis(),basis00.nBasis());
*/

    RealMatrix O_00_p0(SMO_00_p0);
    RealMatrix O_00_m0(SMO_00_m0);
    RealMatrix O_00_0p(SMO_00_0p);
    RealMatrix O_00_0m(SMO_00_0m);



    cout << "  (X,Y) -> (X+DX,Y) MO Mapping" << endl;
    for(auto iMO = 0; iMO < SMO_00_p0.cols(); iMO++){
      RealMatrix::Index maxRow;
      double maxVal = O_00_p0.cwiseAbs().col(iMO).maxCoeff(&maxRow);
      std::stringstream mapStr, ovlpStr; 
      mapStr  << "    " << iMO << " -> " << maxRow;
      ovlpStr << "< " << iMO << "(X,Y) | " << maxRow << "(X+DX,Y) >"; 
              

      cout << std::left << std::setw(10) << mapStr.str();
      cout << "\t";
      cout << std::left << std::setw(25) << ovlpStr.str() << " = " << maxVal;
      cout << endl;
    }

    cout << endl;
    cout << "  (X,Y) -> (X-DX,Y) MO Mapping" << endl;
    for(auto iMO = 0; iMO < SMO_00_m0.cols(); iMO++){
      RealMatrix::Index maxRow;
      double maxVal = O_00_m0.cwiseAbs().col(iMO).maxCoeff(&maxRow);
      std::stringstream mapStr, ovlpStr; 
      mapStr  << "    " << iMO << " -> " << maxRow;
      ovlpStr << "< " << iMO << "(X,Y) | " << maxRow << "(X-DX,Y) >"; 
              

      cout << std::left << std::setw(10) << mapStr.str();
      cout << "\t";
      cout << std::left << std::setw(25) << ovlpStr.str() << " = " << maxVal;
      cout << endl;
    }

    cout << endl;
    cout << "  (X,Y) -> (X,Y+DY) MO Mapping" << endl;
    for(auto iMO = 0; iMO < SMO_00_0p.cols(); iMO++){
      RealMatrix::Index maxRow;
      double maxVal = O_00_0p.cwiseAbs().col(iMO).maxCoeff(&maxRow);
      std::stringstream mapStr, ovlpStr; 
      mapStr  << "    " << iMO << " -> " << maxRow;
      ovlpStr << "< " << iMO << "(X,Y) | " << maxRow << "(X,Y+DY) >"; 
              

      cout << std::left << std::setw(10) << mapStr.str();
      cout << "\t";
      cout << std::left << std::setw(25) << ovlpStr.str() << " = " << maxVal;
      cout << endl;
    }

    cout << endl;
    cout << "  (X,Y) -> (X,Y-DY) MO Mapping" << endl;
    for(auto iMO = 0; iMO < SMO_00_0m.cols(); iMO++){
      RealMatrix::Index maxRow;
      double maxVal = O_00_0m.cwiseAbs().col(iMO).maxCoeff(&maxRow);
      std::stringstream mapStr, ovlpStr; 
      mapStr  << "    " << iMO << " -> " << maxRow;
      ovlpStr << "< " << iMO << "(X,Y) | " << maxRow << "(X,Y-DY) >"; 
              

      cout << std::left << std::setw(10) << mapStr.str();
      cout << "\t";
      cout << std::left << std::setw(25) << ovlpStr.str() << " = " << maxVal;
      cout << endl;
    }

    cout << endl;


    for(auto mu = 0; mu < SMO_00_p0.rows(); mu++)
    for(auto nu = 0; nu < SMO_00_p0.rows(); nu++){
      if(std::abs(O_00_p0(mu,nu)) < 1e-2) O_00_p0(mu,nu) = 0.0;
      else if(O_00_p0(mu,nu) > 0.0)       O_00_p0(mu,nu) = 1.0;
      else                                  O_00_p0(mu,nu) = -1.0;
      if(std::abs(O_00_m0(mu,nu)) < 1e-2) O_00_m0(mu,nu) = 0.0;
      else if(O_00_m0(mu,nu) > 0.0)       O_00_m0(mu,nu) = 1.0;
      else                                  O_00_m0(mu,nu) = -1.0;
      if(std::abs(O_00_0p(mu,nu)) < 1e-2) O_00_0p(mu,nu) = 0.0;
      else if(O_00_0p(mu,nu) > 0.0)       O_00_0p(mu,nu) = 1.0;
      else                                  O_00_0p(mu,nu) = -1.0;
      if(std::abs(O_00_0m(mu,nu)) < 1e-2) O_00_0m(mu,nu) = 0.0;
      else if(O_00_0m(mu,nu) > 0.0)       O_00_0m(mu,nu) = 1.0;
      else                                  O_00_0m(mu,nu) = -1.0;
    }

    RealMatrix TMP;
    TMP = (*ss_p0.moA()) * O_00_p0; 
    (*ss_p0.moA()) = TMP;
    TMP = (*ss_m0.moA()) * O_00_m0; 
    (*ss_m0.moA()) = TMP;
    TMP = (*ss_0p.moA()) * O_00_0p; 
    (*ss_0p.moA()) = TMP;
    TMP = (*ss_0m.moA()) * O_00_0m; 
    (*ss_0m.moA()) = TMP;
  

    if(debug) {
      cout << endl;
      cout << "  Checking | C - C' | After Phase Check:" << endl;
      
      cout << "  | C(X,Y) - C(X+DX,Y) | = " 
           << diffNorm((*ss_00.moA()),(*ss_p0.moA())) 
           << endl;  
      cout << "  | C(X,Y) - C(X-DX,Y) | = " 
           << diffNorm((*ss_00.moA()),(*ss_m0.moA())) 
           << endl;  
      cout << "  | C(X,Y) - C(X,Y+DY) | = " 
           << diffNorm((*ss_00.moA()),(*ss_0p.moA())) 
           << endl;  
      cout << "  | C(X,Y) - C(X,Y-DY) | = " 
           << diffNorm((*ss_00.moA()),(*ss_0m.moA())) 
           << endl;  
    }

    if(debug){
      // < p (X,Y) | S[(X,Y),(X+DX,Y)] | q(X+DX,Y) >
      SMO_00_p0 = 
        ss_00.moA()->transpose() * S_00_p0 * (*ss_p0.moA());
      // < p (X,Y) | S[(X,Y),(X-DX,Y)] | q(X-DX,Y) >
      SMO_00_m0 = 
        ss_00.moA()->transpose() * S_00_m0 * (*ss_m0.moA());
     
      // < p (X,Y+DY) | S[(X,Y),(X,Y+DY)] | q(X,Y+DY) >
      SMO_00_0p = 
        ss_00.moA()->transpose() * S_00_0p * (*ss_0p.moA());
      // < p (X,Y-DY) | S[(X,Y),(X,Y-DY)] | q(X,Y-DY) >
      SMO_00_0m = 
        ss_00.moA()->transpose() * S_00_0m * (*ss_0m.moA());

      // Quantify deviation from C**H * S * C = I
      cout << endl;
      cout << "  Checking | C(X)**H * S(X,X') * C(X') - I |:" 
           << endl;
     
      cout << "  < (X,Y) | (X,Y) >       = " 
           << diffNormI(SMO_00_00) 
           << endl;
      cout << "  < (X+DX,Y) | (X+DX,Y) > = " 
           << diffNormI(SMO_p0_p0) 
           << endl;
      cout << "  < (X-DX,Y) | (X-DX,Y) > = " 
           << diffNormI(SMO_m0_m0) 
           << endl;
      cout << "  < (X,Y+DY) | (X,Y+DY) > = " 
           << diffNormI(SMO_0p_0p) 
           << endl;
      cout << "  < (X,Y-DY) | (X,Y-DY) > = " 
           << diffNormI(SMO_0m_0m) 
           << endl;
     
     
     
      cout << "  < (X,Y) | (X+DX,Y) >    = " 
           << diffNormI(SMO_00_p0) 
           << endl;
      cout << "  < (X,Y) | (X-DX,Y) >    = " 
           << diffNormI(SMO_00_m0) 
           << endl;
      cout << "  < (X,Y) | (X,Y+DY) >    = " 
           << diffNormI(SMO_00_0p) 
           << endl;
      cout << "  < (X,Y) | (X,Y-DY) >    = " 
           << diffNormI(SMO_00_0m) 
           << endl;
    }

    moints_00.initMeta();
    moints_p0.initMeta();
    moints_m0.initMeta();
    moints_0p.initMeta();
    moints_0m.initMeta();

    resp_00.communicate(ss_00,moints_00,fileio_00);
    resp_p0.communicate(ss_p0,moints_p0,fileio_p0);
    resp_m0.communicate(ss_m0,moints_m0,fileio_m0);
    resp_0p.communicate(ss_0p,moints_0p,fileio_0p);
    resp_0m.communicate(ss_0m,moints_0m,fileio_0m);

    resp_00.setMeth(respType);
    resp_p0.setMeth(respType);
    resp_m0.setMeth(respType);
    resp_0p.setMeth(respType);
    resp_0m.setMeth(respType);

    resp_00.doSA();
    resp_p0.doSA();
    resp_m0.doSA();
    resp_0p.doSA();
    resp_0m.doSA();

    resp_00.setNSek(nFreq);
    resp_p0.setNSek(nFreq);
    resp_m0.setNSek(nFreq);
    resp_0p.setNSek(nFreq);
    resp_0m.setNSek(nFreq);

    resp_00.doFull();
    resp_p0.doFull();
    resp_m0.doFull();
    resp_0p.doFull();
    resp_0m.doFull();

    cout << endl;
    cout << "  Performing Response Calculations (";
    if(respType == RESPONSE_TYPE::CIS)       cout << "CIS" ;
    else if(respType == RESPONSE_TYPE::PPTDA) cout << "PPTDA" ;
    cout << ")" << endl;

    cout << "  Performing Response (X,Y)" << endl;
    resp_00.doResponse();
    cout << "  Performing Response (X+DX,Y)" << endl;
    resp_p0.doResponse();
    cout << "  Performing Response (X-DX,Y)" << endl;
    resp_m0.doResponse();
    cout << "  Performing Response (X,Y+DY)" << endl;
    resp_0p.doResponse();
    cout << "  Performing Response (X,Y-DY)" << endl;
    resp_0m.doResponse();





    // Start Differentiation things
    //


    
    double scf_00 = ss_00.totalEnergy;
    double scf_p0 = ss_p0.totalEnergy;
    double scf_m0 = ss_m0.totalEnergy;
    double scf_0p = ss_0p.totalEnergy;
    double scf_0m = ss_0m.totalEnergy;

    GS_ENERGY_FILE << X << "," << Y << "," << scf_00 << endl;

    if(debug) {
      cout << endl;
      cout << "  SCF ENERGIES:" << endl;
      cout << "  SCF Energy (X,Y):    " << std::setprecision(8) 
           << scf_00 << " Eh" << endl;
      cout << "  SCF Energy (X+DX,Y): " << std::setprecision(8) 
           << scf_p0 << " Eh" << endl;
      cout << "  SCF Energy (X-DX,Y): " << std::setprecision(8) 
           << scf_m0 << " Eh" << endl;
      cout << "  SCF Energy (X,Y+DY): " << std::setprecision(8) 
           << scf_0p << " Eh" << endl;
      cout << "  SCF Energy (X,Y-DY): " << std::setprecision(8) 
           << scf_0m << " Eh" << endl;
   }

    Eigen::VectorXd freq_00,freq_p0,freq_m0,freq_0p,freq_0m;

    if(this->respType_ == RESPONSE_TYPE::CIS){
      freq_00 = resp_00.template frequencies<SINGLETS>().head(nFreq);
      freq_p0 = resp_p0.template frequencies<SINGLETS>().head(nFreq);
      freq_m0 = resp_m0.template frequencies<SINGLETS>().head(nFreq);
      freq_0p = resp_0p.template frequencies<SINGLETS>().head(nFreq);
      freq_0m = resp_0m.template frequencies<SINGLETS>().head(nFreq);
    } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
      freq_00 = resp_00.template frequencies<A_PPTDA_SINGLETS>().head(nFreq);
      freq_p0 = resp_p0.template frequencies<A_PPTDA_SINGLETS>().head(nFreq);
      freq_m0 = resp_m0.template frequencies<A_PPTDA_SINGLETS>().head(nFreq);
      freq_0p = resp_0p.template frequencies<A_PPTDA_SINGLETS>().head(nFreq);
      freq_0m = resp_0m.template frequencies<A_PPTDA_SINGLETS>().head(nFreq);
    }

    ES_ENERGY_FILE << X << "," << Y << ",";
    for(auto iSt = 0; iSt < nFreq; iSt++){
      ES_ENERGY_FILE << freq_00(iSt);
      if(iSt != (nFreq-1)) ES_ENERGY_FILE << ",";
    }
    ES_ENERGY_FILE << endl;

    if(debug) {
      cout << endl << "  Excitation template frequencies:" << endl;
     
      for(auto iSt = 0; iSt < nFreq; iSt++){
        cout << "  W(" << iSt << ") (X,Y):    " 
             << std::setprecision(8) << freq_00[iSt] << " Eh" 
             << endl;
        cout << "  W(" << iSt << ") (X+DX,Y): " 
             << std::setprecision(8) << freq_p0[iSt] << " Eh" 
             << endl;
        cout << "  W(" << iSt << ") (X-DX,Y): " 
             << std::setprecision(8) << freq_m0[iSt] << " Eh" 
             << endl;
        cout << "  W(" << iSt << ") (X,Y+DY): " 
             << std::setprecision(8) << freq_0p[iSt] << " Eh" 
             << endl;
        cout << "  W(" << iSt << ") (X,Y-DY): " 
             << std::setprecision(8) << freq_0m[iSt] << " Eh" 
             << endl;


        cout << endl;
      }
    }





    RealMatrix T_00, T_p0, T_m0, T_0p, T_0m; 

    if(this->respType_ == RESPONSE_TYPE::CIS){
      T_00 = resp_00.template transDen<SINGLETS>().block(0,0,
        resp_00.template nMatDim<SINGLETS>(),nFreq);
      T_p0 = resp_p0.template transDen<SINGLETS>().block(0,0,
        resp_p0.template nMatDim<SINGLETS>(),nFreq);
      T_m0 = resp_m0.template transDen<SINGLETS>().block(0,0,
        resp_m0.template nMatDim<SINGLETS>(),nFreq);
      T_0p = resp_0p.template transDen<SINGLETS>().block(0,0,
        resp_0p.template nMatDim<SINGLETS>(),nFreq);
      T_0m = resp_0m.template transDen<SINGLETS>().block(0,0,
        resp_0m.template nMatDim<SINGLETS>(),nFreq);
    } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
      T_00 = resp_00.template transDen<A_PPTDA_SINGLETS>().block(0,0,
        resp_00.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
      T_p0 = resp_p0.template transDen<A_PPTDA_SINGLETS>().block(0,0,
        resp_p0.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
      T_m0 = resp_m0.template transDen<A_PPTDA_SINGLETS>().block(0,0,
        resp_m0.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
      T_0p = resp_0p.template transDen<A_PPTDA_SINGLETS>().block(0,0,
        resp_0p.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
      T_0m = resp_0m.template transDen<A_PPTDA_SINGLETS>().block(0,0,
        resp_0m.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
    }



    if(debug){
      cout << endl;
      cout << "  Checking | T - T' | Before Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " << diffNorm(T_00,T_p0) 
           << endl;  
      cout << "  | T(X,Y) - T(X-DX,Y) | = " << diffNorm(T_00,T_m0) 
           << endl;  
      cout << "  | T(X,Y) - T(X,Y+DY) | = " << diffNorm(T_00,T_0p) 
           << endl;  
      cout << "  | T(X,Y) - T(X,Y-DY) | = " << diffNorm(T_00,T_0m) 
           << endl;  
    }


    Eigen::VectorXd Trans_t_0p;
    Eigen::VectorXd Trans_t_0m;
    Eigen::VectorXd Trans_t_p0;
    Eigen::VectorXd Trans_t_m0;
    if(this->respType_ == RESPONSE_TYPE::CIS){
      Trans_t_p0 = checkPhase(T_00.data(),T_p0.data(),
        resp_p0.template nMatDim<SINGLETS>(),nFreq);
      Trans_t_m0 = checkPhase(T_00.data(),T_m0.data(),
        resp_m0.template nMatDim<SINGLETS>(),nFreq);
      Trans_t_0p = checkPhase(T_00.data(),T_0p.data(),
        resp_0p.template nMatDim<SINGLETS>(),nFreq);
      Trans_t_0m = checkPhase(T_00.data(),T_0m.data(),
        resp_0m.template nMatDim<SINGLETS>(),nFreq);
    } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
      Trans_t_p0 = checkPhase(T_00.data(),T_p0.data(),
        resp_p0.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
      Trans_t_m0 = checkPhase(T_00.data(),T_m0.data(),
        resp_m0.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
      Trans_t_0p = checkPhase(T_00.data(),T_0p.data(),
        resp_0p.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
      Trans_t_0m = checkPhase(T_00.data(),T_0m.data(),
        resp_0m.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
    }


    if(debug){
      cout << endl;
      cout << "  Checking | T - T' | After Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " << diffNorm(T_00,T_p0) 
           << endl;  
      cout << "  | T(X,Y) - T(X-DX,Y) | = " << diffNorm(T_00,T_m0) 
           << endl;  
      cout << "  | T(X,Y) - T(X,Y+DY) | = " << diffNorm(T_00,T_0p) 
           << endl;  
      cout << "  | T(X,Y) - T(X,Y-DY) | = " << diffNorm(T_00,T_0m) 
           << endl;  
    }

    if(debug) {
      cout << endl << " Checking < T | T > = 1 " << endl;
      for(auto iSt = 0; iSt < nFreq; iSt++){
        cout << "   State " << iSt << ":" << endl;
        cout << "    < T(X,Y) | T(X,Y) >       = " 
             << selfInner(T_00.col(iSt)) 
             << endl;
        cout << "    < T(X+DX,Y) | T(X+DX,Y) > = " 
             << selfInner(T_p0.col(iSt)) 
             << endl;
        cout << "    < T(X-DX,Y) | T(X-DX,Y) > = " 
             << selfInner(T_m0.col(iSt)) 
             << endl;
        cout << "    < T(X,Y+DY) | T(X,Y+DY) > = " 
             << selfInner(T_0p.col(iSt)) 
             << endl;
        cout << "    < T(X,Y-DY) | T(X,Y-DY) > = " 
             << selfInner(T_0m.col(iSt)) 
             << endl;
      }
    }








    // Calculate Derivatives
      
    // SCF energy gradient
    double gsdx = (scf_p0 - scf_m0) / (2*(XP-X));
    double gsdy = (scf_0p - scf_0m) / (2*(YP-Y));
    double gsnormd = std::sqrt(gsdx*gsdx + gsdy*gsdy);

    GS_GRADIENT_DX_FILE << X << "," << Y << "," << gsdx << endl;
    GS_GRADIENT_DY_FILE << X << "," << Y << "," << gsdy << endl;

    if(debug)
      cout << endl << "  GS Gradient = (" << gsdx << "," << gsdy <<
          ")" << endl;

    // Excitation frequency gradient
    Eigen::VectorXd freqDX = (freq_p0 - freq_m0)/(2*(XP-X));
    Eigen::VectorXd freqDY = (freq_0p - freq_0m)/(2*(YP-Y));

    Eigen::VectorXd freqNorm = freqDX.cwiseProduct(freqDX);
    freqNorm += freqDY.cwiseProduct(freqDY);

    for(auto iFreq = 0; iFreq < nFreq; iFreq++)
      freqNorm(iFreq) = std::sqrt(freqNorm(iFreq));


    ES_GRADIENT_DX_FILE << X << "," << Y << ",";
    ES_GRADIENT_DY_FILE << X << "," << Y << ",";
    for(auto iSt = 0; iSt < nFreq; iSt++){
      ES_GRADIENT_DX_FILE << freqDX(iSt);
      ES_GRADIENT_DY_FILE << freqDY(iSt);
      if(iSt != (nFreq-1)) {
        ES_GRADIENT_DX_FILE << ",";
        ES_GRADIENT_DY_FILE << ",";
      }
    }
    ES_GRADIENT_DX_FILE << endl;
    ES_GRADIENT_DY_FILE << endl;

    if(debug){
      cout << endl << "  ES Gradients:" << endl;
      for(auto iSt = 0; iSt < nFreq; iSt++){
        cout << "   W(" << iSt << ")' = (" << freqDX(iSt) 
             << "," << freqDY(iSt) << ")" << endl;
      }
    }





    // Assemble ground to excited state couplings
    //
    
    int NOCC = ss_00.nOccA();
    int NVIR = ss_00.nVirA();

    if(debug){
     // Compute Overlaps of wavefunctions at different geometries
     double OvLp_00_00 = 
       SMO_00_00.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
     double OvLp_p0_p0 = 
       SMO_p0_p0.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
     double OvLp_m0_m0 = 
       SMO_m0_m0.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
     double OvLp_0p_0p = 
       SMO_0p_0p.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
     double OvLp_0m_0m = 
       SMO_0m_0m.block(0,0,NOCC,NOCC).eigenvalues().prod().real();

     double OvLp_00_0p = 
       SMO_00_0p.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
     double OvLp_00_0m = 
       SMO_00_0m.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
     double OvLp_00_p0 = 
       SMO_00_p0.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
     double OvLp_00_m0 = 
       SMO_00_m0.block(0,0,NOCC,NOCC).eigenvalues().prod().real();


     cout << endl;
     cout << "  Checking Wavefunction Overlaps < X | X' >:" << endl;
     cout << "  < (X,Y) | (X,Y) >          = " << OvLp_00_00 << endl;
     cout << "  < (X+DX,Y) | (X+DX,Y) >    = " << OvLp_p0_p0 << endl;
     cout << "  < (X-DX,Y) | (X-DX,Y) >    = " << OvLp_m0_m0 << endl;
     cout << "  < (X,Y+DY) | (X,Y+DY) >    = " << OvLp_0p_0p << endl;
     cout << "  < (X,Y-DY) | (X,Y-DY) >    = " << OvLp_0m_0m << endl;

     cout << "  < (X,Y) | (X+DX,Y) >       = " << OvLp_00_p0 << endl;
     cout << "  < (X,Y) | (X-DX,Y) >       = " << OvLp_00_m0 << endl;
     cout << "  < (X,Y) | (X,Y+DY) >       = " << OvLp_00_0p << endl;
     cout << "  < (X,Y) | (X,Y-DY) >       = " << OvLp_00_0m << endl;
   }


   std::vector<Eigen::VectorXd> NAC_ES_GS;
   std::vector<RealMatrix>      NAC_ES_ES;

   if(respType == RESPONSE_TYPE::CIS){
     NAC_ES_GS = ES_GS_NACME_CIS(nFreq,true,NOCC,NVIR,XP-X,YP-Y,
       T_00,(*ss_00.moA()),(*ss_p0.moA()),(*ss_m0.moA()),
       (*ss_0p.moA()),(*ss_0m.moA()),S_00_p0,S_00_m0,S_00_0p,
       S_00_0m);
     NAC_ES_ES = ES_ES_NACME_CIS(nFreq,true,NOCC,NVIR,XP-X,YP-Y,
       T_00,T_p0,T_m0,T_0p,T_0m,(*ss_00.moA()),(*ss_p0.moA()),
       (*ss_m0.moA()),(*ss_0p.moA()),(*ss_0m.moA()),S_00_p0,S_00_m0,
       S_00_0p,S_00_0m);

     ES_GS_NACME_DX_FILE << X << "," << Y << ",";
     ES_GS_NACME_DY_FILE << X << "," << Y << ",";
     for(auto iSt = 0; iSt < nFreq; iSt++){
       ES_GS_NACME_DX_FILE << NAC_ES_GS[0](iSt);
       ES_GS_NACME_DY_FILE << NAC_ES_GS[1](iSt);
       if(iSt != (nFreq-1)) {
         ES_GS_NACME_DX_FILE << ",";
         ES_GS_NACME_DY_FILE << ",";
       }
     }
     ES_GS_NACME_DX_FILE << endl;
     ES_GS_NACME_DY_FILE << endl;


   } else if(respType == RESPONSE_TYPE::PPTDA) {
     NAC_ES_ES = ES_ES_NACME_PPTDA(nFreq,true,NOCC,NVIR,XP-X,YP-Y,
       T_00,T_p0,T_m0,T_0p,T_0m,(*ss_00.moA()),(*ss_p0.moA()),
       (*ss_m0.moA()),(*ss_0p.moA()),(*ss_0m.moA()),S_00_p0,S_00_m0,
       S_00_0p,S_00_0m);
   }


   ES_ES_NACME_DX_FILE << X << "," << Y << ",";
   ES_ES_NACME_DY_FILE << X << "," << Y << ",";
   for(auto iSt = 0; iSt < nFreq; iSt++)
   for(auto jSt = 0; jSt < nFreq; jSt++){
     ES_ES_NACME_DX_FILE << NAC_ES_ES[0](iSt,jSt);
     ES_ES_NACME_DY_FILE << NAC_ES_ES[1](iSt,jSt);
     if(iSt != (nFreq-1) || jSt != (nFreq-1)) {
       ES_ES_NACME_DX_FILE << ",";
       ES_ES_NACME_DY_FILE << ",";
     }
   }
   ES_ES_NACME_DX_FILE << endl;
   ES_ES_NACME_DY_FILE << endl;

   if(debug){
     std::vector<Eigen::VectorXd> NAC_GS_ES;

     if(respType == RESPONSE_TYPE::CIS) {
       NAC_GS_ES = GS_ES_NACME_CIS(nFreq,false,NOCC,NVIR,XP-X,YP-Y,
         T_00,T_p0,T_m0,T_0p,T_0m,(*ss_00.moA()),(*ss_p0.moA()),
         (*ss_m0.moA()),(*ss_0p.moA()),(*ss_0m.moA()),S_00_p0,
         S_00_m0,S_00_0p,S_00_0m);
      
       prettyPrint(cout,NAC_ES_GS[0]*0.529177,"ES->GS NACME DX");
       prettyPrint(cout,NAC_ES_GS[1]*0.529177,"ES->GS NACME DY");
       prettyPrint(cout,NAC_GS_ES[0]*0.529177,"GS->ES NACME DX");
       prettyPrint(cout,NAC_GS_ES[1]*0.529177,"GS->ES NACME DY");
       prettyPrint(cout,NAC_ES_GS[0] + NAC_GS_ES[0],"DIFF DX");
       prettyPrint(cout,NAC_ES_GS[1] + NAC_GS_ES[1],"DIFF DY");


     }

     prettyPrint(cout,NAC_ES_ES[0]*0.529177,"ES->ES DX");
     prettyPrint(cout,NAC_ES_ES[1]*0.529177,"ES->ES DY");

   }

  }
  finalizeCQ();
}



// Naming convention PREFIX_xptx_ypty.xyz
std::string genFName(double x, double y, int px, int py,
  std::string PREFIX) {
    std::stringstream ss1, ss2;

    ss1 << std::fixed << std::setprecision(px) << x;
    ss2 << std::fixed << std::setprecision(py) << y;

    std::string str1 = ss1.str();
    std::string str2 = ss2.str();

    auto ptPos = str1.find(".");
    str1.replace(str1.begin()+ptPos,str1.begin()+ptPos+1,"pt");

    ptPos = str2.find(".");
    str2.replace(str2.begin()+ptPos,str2.begin()+ptPos+1,"pt");

    std::string fName = PREFIX + "_" + str1 + "_" + str2;    

    return fName;
};

Eigen::VectorXd checkPhase(double *A, double *B,int n, int m){
  RealMap AMap(A,n,m);
  RealMap BMap(B,n,m);

  Eigen::VectorXd Trans(m);
  Trans.setZero();

  for(auto j = 0; j < m; j++){
    int sgn1(1), sgn2(1);
    for(auto i = 0; i < n; i++){
      if(std::abs(AMap(i,j)) > 1e-8) {
        if(std::abs(BMap(i,j)) < 1e-8) continue;
        else { // check B
          sgn1 = AMap(i,j) / std::abs(AMap(i,j));
          sgn2 = BMap(i,j) / std::abs(BMap(i,j));
          break;
        } // get signs
      } // check A
    } // loop i
    if(sgn1 != sgn2) {
      BMap.col(j) *= -1.0;
      Trans(j) = -1.0;
    } else {
      Trans(j) = 1.0;
    }
  } // loop j

  return Trans;
};

RealMatrix genSpx(BasisSet &obs1, BasisSet &obs2){
  if(obs1.nShell() != obs2.nShell()) {
    cout << "genSpx cannot take two basis sets of different size" 
         << endl;
    std::exit(1);
  }
  RealMatrix S(obs1.nBasis(),obs2.nBasis());

  libint2::OneBodyEngine engine(libint2::OneBodyEngine::overlap,
    obs1.maxPrim(),obs1.maxL(),0);

  S.setZero();
  for(auto iSh = 0; iSh < obs1.nShell(); iSh++)
  for(auto jSh = 0; jSh < obs2.nShell(); jSh++){
    auto s1 = obs1.mapSh2Bf(iSh);
    auto s2 = obs1.mapSh2Bf(jSh);
    auto n1 = obs1.shells(iSh).size();
    auto n2 = obs2.shells(jSh).size();
    const auto *buff = engine.compute(
      obs1.shells(iSh),obs2.shells(jSh)
    );

    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,
                                   Eigen::Dynamic,
                                   Eigen::RowMajor>
      > sMap(buff,n1,n2);
    S.block(s1,s2,n1,n2) = sMap;
  }

/*
  // Renormalize
  Eigen::VectorXd norms1(obs1.nbf());
  for(auto iSh = 0; iSh < obs1.size(); iSh++){
    auto s1 = shell2bf[iSh];
    auto n1 = obs1[iSh].size();
    const auto *buff = engine.compute(obs1[iSh],obs1[iSh]);
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,
                                   Eigen::Dynamic,
                                   Eigen::RowMajor>
      > sMap(buff,n1,n1);

    for(auto iBf = s1, iloc = 0ul; iBf < s1 + n1; iBf++, iloc++)
      norms1(iBf) = std::sqrt(sMap(iloc,iloc));
  }

  Eigen::VectorXd norms2(obs2.nbf());
  for(auto iSh = 0; iSh < obs2.size(); iSh++){
    auto s2 = shell2bf[iSh];
    auto n2 = obs2[iSh].size();
    const auto *buff = engine.compute(obs2[iSh],obs2[iSh]);
    Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,
                                   Eigen::Dynamic,
                                   Eigen::RowMajor>
      > sMap(buff,n2,n2);

    for(auto iBf = s2, iloc = 0ul; iBf < s2 + n2; iBf++, iloc++)
      norms2(iBf) = std::sqrt(sMap(iloc,iloc));
  }
  
  for(auto iBf = 0; iBf < obs1.nbf(); iBf++)
  for(auto jBf = 0; jBf < obs2.nbf(); jBf++) {
    S(iBf,jBf) = S(iBf,jBf) / (norms1(iBf)*norms2(jBf));
  }
*/

  return S;
};
std::vector<Eigen::VectorXd> ES_GS_NACME_CIS(int nFreq,bool renorm,
  int nocc, int nvir, double step1, double step2, RealMatrix &T_00, 
  RealMatrix &MO_00, RealMatrix &MO_p0, RealMatrix &MO_m0, 
  RealMatrix &MO_0p, RealMatrix &MO_0m, RealMatrix &S_00_p0, 
  RealMatrix &S_00_m0, RealMatrix &S_00_0p, RealMatrix &S_00_0m){

  std::vector<Eigen::VectorXd> NACME(2,Eigen::VectorXd(nFreq));

  for(auto i = 0; i < NACME.size(); i++) NACME[i].setZero();

  if(renorm) T_00 *= std::sqrt(0.5);
//T_00.transposeInPlace();
//prettyPrint(cout,T_00,"T");

  int nThreads = omp_get_max_threads();

  std::vector<RealMatrix> SWAPPED_00(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_0p(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_0m(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_p0(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_m0(nThreads,RealMatrix(MO_00));

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < nFreq; iSt++){
    cout << "Working on iSt = " << iSt << " ES->GS" << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto i = 0, ia = 0; i < nocc; i++)
      for(auto a = 0; a < nvir; a++, ia++){
        if(ia % nThreads != thread_id) continue;
        if(std::abs(T_00(ia,iSt)) < 1e-8) continue;

        SWAPPED_00[thread_id].setZero();
        Prod_00_0p[thread_id].setZero();
        Prod_00_0m[thread_id].setZero();
        Prod_00_p0[thread_id].setZero();
        Prod_00_m0[thread_id].setZero();
        
        SWAPPED_00[thread_id] = MO_00;
        SWAPPED_00[thread_id].col(i).swap(SWAPPED_00[thread_id].col(nocc+a)); 
        
        Prod_00_0p[thread_id] = SWAPPED_00[thread_id].transpose() * S_00_0p * MO_0p;
        Prod_00_0m[thread_id] = SWAPPED_00[thread_id].transpose() * S_00_0m * MO_0m;
        Prod_00_p0[thread_id] = SWAPPED_00[thread_id].transpose() * S_00_p0 * MO_p0;
        Prod_00_m0[thread_id] = SWAPPED_00[thread_id].transpose() * S_00_m0 * MO_m0;
     
        double OvLp_IASwap_00_0p =
          Prod_00_0p[thread_id].block(0,0,nocc,nocc).determinant();
        double OvLp_IASwap_00_0m =
          Prod_00_0m[thread_id].block(0,0,nocc,nocc).determinant();
        double OvLp_IASwap_00_p0 =
          Prod_00_p0[thread_id].block(0,0,nocc,nocc).determinant();
        double OvLp_IASwap_00_m0 =
          Prod_00_m0[thread_id].block(0,0,nocc,nocc).determinant();
     
        #pragma omp critical
        {
          NACME[0](iSt) += T_00(ia,iSt) *
            (OvLp_IASwap_00_p0 - OvLp_IASwap_00_m0) /
            (2*step1);
          NACME[1](iSt) += T_00(ia,iSt) *
            (OvLp_IASwap_00_0p - OvLp_IASwap_00_0m) /
            (2*step2);
        }
      } // loop ove ia
    }
  } // loop over states

  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "ES -> GS NACME took " << elapsed.count() << " s" << endl;
  NACME[0] = 2*NACME[0];
  NACME[1] = 2*NACME[1];
  return NACME;
};

std::vector<Eigen::VectorXd> GS_ES_NACME_CIS(int nFreq,bool renorm,
  int nocc, int nvir, double step1, double step2, RealMatrix &T_00, 
  RealMatrix &T_p0, RealMatrix &T_m0, RealMatrix &T_0p, 
  RealMatrix &T_0m, RealMatrix &MO_00, RealMatrix &MO_p0, 
  RealMatrix &MO_m0, RealMatrix &MO_0p, RealMatrix &MO_0m, 
  RealMatrix &S_00_p0, RealMatrix &S_00_m0, RealMatrix &S_00_0p, 
  RealMatrix &S_00_0m){

  std::vector<Eigen::VectorXd> NACME(2,Eigen::VectorXd(nFreq));

  for(auto i = 0; i < NACME.size(); i++) NACME[i].setZero();

  if(renorm) {
    T_p0 *= std::sqrt(0.5);
    T_m0 *= std::sqrt(0.5);
    T_0p *= std::sqrt(0.5);
    T_0m *= std::sqrt(0.5);
  }

  RealMatrix SWAPPED_p0(MO_00);
  RealMatrix SWAPPED_m0(MO_00);
  RealMatrix SWAPPED_0p(MO_00);
  RealMatrix SWAPPED_0m(MO_00);

  RealMatrix Prod_00_0p(MO_00);
  RealMatrix Prod_00_0m(MO_00);
  RealMatrix Prod_00_p0(MO_00);
  RealMatrix Prod_00_m0(MO_00);

  SWAPPED_p0.setZero();
  SWAPPED_m0.setZero();
  SWAPPED_0p.setZero();
  SWAPPED_0m.setZero();

  Prod_00_0p.setZero();
  Prod_00_0m.setZero();
  Prod_00_p0.setZero();
  Prod_00_m0.setZero();
  for(auto iSt = 0; iSt < nFreq; iSt++){
    for(auto i = 0, ia = 0; i < nocc; i++)
    for(auto a = 0; a < nvir; a++, ia++){
      if(
        (std::abs(T_p0(ia,iSt)) < 1e-8) &&
        (std::abs(T_m0(ia,iSt)) < 1e-8) &&
        (std::abs(T_0p(ia,iSt)) < 1e-8) &&
        (std::abs(T_0m(ia,iSt)) < 1e-8)
      ) continue;

      SWAPPED_p0.setZero();
      SWAPPED_m0.setZero();
      SWAPPED_0p.setZero();
      SWAPPED_0m.setZero();

      Prod_00_0p.setZero();
      Prod_00_0m.setZero();
      Prod_00_p0.setZero();
      Prod_00_m0.setZero();
      
      SWAPPED_p0 = MO_p0;
      SWAPPED_m0 = MO_m0;
      SWAPPED_0p = MO_0p;
      SWAPPED_0m = MO_0m;

      SWAPPED_p0.col(i).swap(SWAPPED_p0.col(nocc+a)); 
      SWAPPED_m0.col(i).swap(SWAPPED_m0.col(nocc+a)); 
      SWAPPED_0p.col(i).swap(SWAPPED_0p.col(nocc+a)); 
      SWAPPED_0m.col(i).swap(SWAPPED_0m.col(nocc+a)); 
      
      Prod_00_0p = MO_00.transpose() * S_00_0p * SWAPPED_0p;
      Prod_00_0m = MO_00.transpose() * S_00_0m * SWAPPED_0m;
      Prod_00_p0 = MO_00.transpose() * S_00_p0 * SWAPPED_p0;
      Prod_00_m0 = MO_00.transpose() * S_00_m0 * SWAPPED_m0;

      
      double OvLp_IASwap_00_0p =
        Prod_00_0p.block(0,0,nocc,nocc).determinant();
      double OvLp_IASwap_00_0m =
        Prod_00_0m.block(0,0,nocc,nocc).determinant();
      double OvLp_IASwap_00_p0 =
        Prod_00_p0.block(0,0,nocc,nocc).determinant();
      double OvLp_IASwap_00_m0 =
        Prod_00_m0.block(0,0,nocc,nocc).determinant();

      NACME[0](iSt) += (
        T_p0(ia,iSt)*OvLp_IASwap_00_p0 - 
        T_m0(ia,iSt)*OvLp_IASwap_00_m0
      ) / (2*step1);
      NACME[1](iSt) += (
        T_0p(ia,iSt)*OvLp_IASwap_00_0p - 
        T_0m(ia,iSt)*OvLp_IASwap_00_0m
      ) / (2*step2);
    } // loop ove ia
  } // loop over states
  NACME[0] *= 2;
  NACME[1] *= 2;

  return NACME;
};

std::vector<RealMatrix> ES_ES_NACME_CIS(int nFreq,bool renorm,
  int nocc, int nvir, double step1, double step2, RealMatrix &T_00, 
  RealMatrix &T_p0, RealMatrix &T_m0, RealMatrix &T_0p, 
  RealMatrix &T_0m, RealMatrix &MO_00, RealMatrix &MO_p0, 
  RealMatrix &MO_m0, RealMatrix &MO_0p, RealMatrix &MO_0m, 
  RealMatrix &S_00_p0, RealMatrix &S_00_m0, RealMatrix &S_00_0p, 
  RealMatrix &S_00_0m){

  bool doValidation = false;

  std::vector<RealMatrix> NACME(2,RealMatrix(nFreq,nFreq));

  for(auto i = 0; i < NACME.size(); i++) NACME[i].setZero();

  if(renorm) {
    T_p0 *= std::sqrt(0.5);
    T_m0 *= std::sqrt(0.5);
    T_0p *= std::sqrt(0.5);
    T_0m *= std::sqrt(0.5);
  }

  int nThreads = omp_get_max_threads();

  std::vector<RealMatrix> SWAPPED_IA_00(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> SWAPPED_JB_p0(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> SWAPPED_JB_m0(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> SWAPPED_JB_0p(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> SWAPPED_JB_0m(nThreads,RealMatrix(MO_00));

  std::vector<RealMatrix> Prod_00_0p(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_0m(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_p0(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_m0(nThreads,RealMatrix(MO_00));

// VALIDATION BEGIN **************
  RealMatrix T1Prod_00_0p(MO_00);
  RealMatrix T1Prod_00_0m(MO_00);
  RealMatrix T1Prod_00_p0(MO_00);
  RealMatrix T1Prod_00_m0(MO_00);
  RealMatrix T2Prod_00_0p(MO_00);
  RealMatrix T2Prod_00_0m(MO_00);
  RealMatrix T2Prod_00_p0(MO_00);
  RealMatrix T2Prod_00_m0(MO_00);
// VALIDATION END **************



  if(doValidation) {
// VALIDATION BEGIN **************
    T1Prod_00_0p.setZero();
    T1Prod_00_0m.setZero();
    T1Prod_00_p0.setZero();
    T1Prod_00_m0.setZero();
    T2Prod_00_0p.setZero();
    T2Prod_00_0m.setZero();
    T2Prod_00_p0.setZero();
    T2Prod_00_m0.setZero();
// VALIDATION END **************
  }

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < nFreq; iSt++)
  for(auto jSt = 0; jSt < nFreq; jSt++){
    cout << "Working on iSt = " << iSt <<" jSt = " << jSt << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto i = 0, ia = 0, iajb = 0; i < nocc; i++)
      for(auto a = 0; a < nvir; a++, ia++){
      //if(T_00(ia,iSt) < 1e-8) continue;   
     
        SWAPPED_IA_00[thread_id].setZero();
        SWAPPED_IA_00[thread_id] = MO_00;
     
        SWAPPED_IA_00[thread_id].col(i).swap(SWAPPED_IA_00[thread_id].col(nocc+a)); 
     
        for(auto j = 0, jb = 0; j < nocc; j++)
        for(auto b = 0; b < nvir; b++, jb++, iajb++){
          if(iajb % nThreads != thread_id) continue;
          if(
            (std::abs(T_p0(jb,jSt)) < 1e-8) &&
            (std::abs(T_m0(jb,jSt)) < 1e-8) &&
            (std::abs(T_0p(jb,jSt)) < 1e-8) &&
            (std::abs(T_0m(jb,jSt)) < 1e-8)
          ) continue;
       
          SWAPPED_JB_p0[thread_id].setZero();
          SWAPPED_JB_m0[thread_id].setZero();
          SWAPPED_JB_0p[thread_id].setZero();
          SWAPPED_JB_0m[thread_id].setZero();
       
          Prod_00_0p[thread_id].setZero();
          Prod_00_0m[thread_id].setZero();
          Prod_00_p0[thread_id].setZero();
          Prod_00_m0[thread_id].setZero();
          
          SWAPPED_JB_p0[thread_id] = MO_p0;
          SWAPPED_JB_m0[thread_id] = MO_m0;
          SWAPPED_JB_0p[thread_id] = MO_0p;
          SWAPPED_JB_0m[thread_id] = MO_0m;
       
          SWAPPED_JB_p0[thread_id].col(j).swap(SWAPPED_JB_p0[thread_id].col(nocc+b)); 
          SWAPPED_JB_m0[thread_id].col(j).swap(SWAPPED_JB_m0[thread_id].col(nocc+b)); 
          SWAPPED_JB_0p[thread_id].col(j).swap(SWAPPED_JB_0p[thread_id].col(nocc+b)); 
          SWAPPED_JB_0m[thread_id].col(j).swap(SWAPPED_JB_0m[thread_id].col(nocc+b)); 
          
          Prod_00_0p[thread_id] = 
            SWAPPED_IA_00[thread_id].transpose() * S_00_0p * SWAPPED_JB_0p[thread_id];
          Prod_00_0m[thread_id] = 
            SWAPPED_IA_00[thread_id].transpose() * S_00_0m * SWAPPED_JB_0m[thread_id];
          Prod_00_p0[thread_id] = 
            SWAPPED_IA_00[thread_id].transpose() * S_00_p0 * SWAPPED_JB_p0[thread_id];
          Prod_00_m0[thread_id] = 
            SWAPPED_IA_00[thread_id].transpose() * S_00_m0 * SWAPPED_JB_m0[thread_id];
       
         if(doValidation) {
//  VVALIDATION BEGIN **************
            T1Prod_00_0p.setZero();
            T1Prod_00_0m.setZero();
            T1Prod_00_p0.setZero();
            T1Prod_00_m0.setZero();
           
            T1Prod_00_0p = 
              MO_00.transpose() * S_00_0p * SWAPPED_JB_0p[thread_id];
            T1Prod_00_0m = 
              MO_00.transpose() * S_00_0m * SWAPPED_JB_0m[thread_id];
            T1Prod_00_p0 = 
              MO_00.transpose() * S_00_p0 * SWAPPED_JB_p0[thread_id];
            T1Prod_00_m0 = 
              MO_00.transpose() * S_00_m0 * SWAPPED_JB_m0[thread_id];
           
            T2Prod_00_0p.setZero();
            T2Prod_00_0m.setZero();
            T2Prod_00_p0.setZero();
            T2Prod_00_m0.setZero();
           
            T2Prod_00_0p = 
              SWAPPED_IA_00[thread_id].transpose() * S_00_0p * MO_0p;
            T2Prod_00_0m = 
              SWAPPED_IA_00[thread_id].transpose() * S_00_0m * MO_0m;
            T2Prod_00_p0 = 
              SWAPPED_IA_00[thread_id].transpose() * S_00_p0 * MO_p0;
            T2Prod_00_m0 = 
              SWAPPED_IA_00[thread_id].transpose() * S_00_m0 * MO_m0;
          }
//  VVALIDATION END **************
          
          double OvLp_IASwap_00_0p =
            Prod_00_0p[thread_id].block(0,0,nocc,nocc).determinant();
          double OvLp_IASwap_00_0m =
            Prod_00_0m[thread_id].block(0,0,nocc,nocc).determinant();
          double OvLp_IASwap_00_p0 =
            Prod_00_p0[thread_id].block(0,0,nocc,nocc).determinant();
          double OvLp_IASwap_00_m0 =
            Prod_00_m0[thread_id].block(0,0,nocc,nocc).determinant();
     
//  VVALIDATION BEGIN **************
          if(doValidation) {
            double T1OvLp_IASwap_00_0p =
              T1Prod_00_0p.block(0,0,nocc,nocc).determinant();
            double T1OvLp_IASwap_00_0m =
              T1Prod_00_0m.block(0,0,nocc,nocc).determinant();
            double T1OvLp_IASwap_00_p0 =
              T1Prod_00_p0.block(0,0,nocc,nocc).determinant();
            double T1OvLp_IASwap_00_m0 =
              T1Prod_00_m0.block(0,0,nocc,nocc).determinant();
           
            double T2OvLp_IASwap_00_0p =
              T2Prod_00_0p.block(0,0,nocc,nocc).determinant();
            double T2OvLp_IASwap_00_0m =
              T2Prod_00_0m.block(0,0,nocc,nocc).determinant();
            double T2OvLp_IASwap_00_p0 =
              T2Prod_00_p0.block(0,0,nocc,nocc).determinant();
            double T2OvLp_IASwap_00_m0 =
              T2Prod_00_m0.block(0,0,nocc,nocc).determinant();
           
            cout << "******** " << iSt << "," << jSt << "," << ia << "," << jb << endl;
            cout << OvLp_IASwap_00_0p  << endl;
            cout << OvLp_IASwap_00_0m  << endl;
            cout << OvLp_IASwap_00_p0  << endl;
            cout << OvLp_IASwap_00_m0  << endl;
            cout << endl;
            cout << T1OvLp_IASwap_00_0p  << endl;
            cout << T1OvLp_IASwap_00_0m  << endl;
            cout << T1OvLp_IASwap_00_p0  << endl;
            cout << T1OvLp_IASwap_00_m0  << endl;
            cout << endl;
            cout << T2OvLp_IASwap_00_0p  << endl;
            cout << T2OvLp_IASwap_00_0m  << endl;
            cout << T2OvLp_IASwap_00_p0  << endl;
            cout << T2OvLp_IASwap_00_m0  << endl;
            cout << endl;
            cout << OvLp_IASwap_00_0p - T1OvLp_IASwap_00_0p*T2OvLp_IASwap_00_0p << endl;
            cout << OvLp_IASwap_00_0m - T1OvLp_IASwap_00_0m*T2OvLp_IASwap_00_0m << endl;
            cout << OvLp_IASwap_00_p0 - T1OvLp_IASwap_00_p0*T2OvLp_IASwap_00_p0 << endl;
            cout << OvLp_IASwap_00_m0 - T1OvLp_IASwap_00_m0*T2OvLp_IASwap_00_m0 << endl;
          }
//  VVALIDATION END **************
       
          #pragma omp critical
          {
            NACME[0](iSt,jSt) += T_00(ia,iSt)*(
              T_p0(jb,jSt)*OvLp_IASwap_00_p0 - 
              T_m0(jb,jSt)*OvLp_IASwap_00_m0
            ) / (2*step1);
            NACME[1](iSt,jSt) += T_00(ia,iSt)*(
              T_0p(jb,jSt)*OvLp_IASwap_00_0p - 
              T_0m(jb,jSt)*OvLp_IASwap_00_0m
            ) / (2*step2);
          }
        } // loop ove jb
      } // loop ove ia
    }
  } // loop over states
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "ES -> ES NACME took " << elapsed.count() << " s" << endl;
  NACME[0] *= 2;
  NACME[1] *= 2;

  return NACME;
};

std::vector<RealMatrix> ES_ES_NACME_PPTDA(int nFreq,bool renorm,
  int nocc, int nvir, double step1, double step2, RealMatrix &T_00, 
  RealMatrix &T_p0, RealMatrix &T_m0, RealMatrix &T_0p, 
  RealMatrix &T_0m, RealMatrix &MO_00, RealMatrix &MO_p0, 
  RealMatrix &MO_m0, RealMatrix &MO_0p, RealMatrix &MO_0m, 
  RealMatrix &S_00_p0, RealMatrix &S_00_m0, RealMatrix &S_00_0p, 
  RealMatrix &S_00_0m){

  std::vector<RealMatrix> NACME(2,RealMatrix(nFreq,nFreq));

  for(auto i = 0; i < NACME.size(); i++) NACME[i].setZero();

  int nThreads = omp_get_max_threads();

  std::vector<RealMatrix> SWAPPED_AB_00(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> SWAPPED_CD_p0(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> SWAPPED_CD_m0(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> SWAPPED_CD_0p(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> SWAPPED_CD_0m(nThreads,RealMatrix(MO_00));

  std::vector<RealMatrix> Prod_00_0p(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_0m(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_p0(nThreads,RealMatrix(MO_00));
  std::vector<RealMatrix> Prod_00_m0(nThreads,RealMatrix(MO_00));

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < nFreq; iSt++)
  for(auto jSt = 0; jSt < nFreq; jSt++){
    cout << "Working on iSt = " << iSt <<" jSt = " << jSt << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto a = 0, ab = 0, abcd = 0; a < nvir; a++      )
      for(auto b = 0        ; b <= a ;  b++, ab++)
      for(auto c = 0, cd = 0; c < nvir; c++      )
      for(auto d = 0        ; d <= c ;  d++, cd++, abcd++){
     
        if(abcd % nThreads != thread_id) continue; 
        if(
          (std::abs(T_p0(cd,jSt)) < 1e-8) &&
          (std::abs(T_0p(cd,jSt)) < 1e-8) &&
          (std::abs(T_m0(cd,jSt)) < 1e-8) &&
          (std::abs(T_0m(cd,jSt)) < 1e-8)
        ) continue;
     
     
        // AC OVerlap
        SWAPPED_AB_00[thread_id].setZero();
        SWAPPED_CD_p0[thread_id].setZero();
        SWAPPED_CD_m0[thread_id].setZero();
        SWAPPED_CD_0p[thread_id].setZero();
        SWAPPED_CD_0m[thread_id].setZero();
     
        Prod_00_0p[thread_id].setZero();
        Prod_00_0m[thread_id].setZero();
        Prod_00_p0[thread_id].setZero();
        Prod_00_m0[thread_id].setZero();
     
        SWAPPED_AB_00[thread_id] = MO_00;
        SWAPPED_CD_p0[thread_id] = MO_p0;
        SWAPPED_CD_m0[thread_id] = MO_m0;
        SWAPPED_CD_0p[thread_id] = MO_0p;
        SWAPPED_CD_0m[thread_id] = MO_0m;
     
        SWAPPED_AB_00[thread_id].col(nocc).swap(SWAPPED_AB_00[thread_id].col(nocc+a));
        SWAPPED_CD_p0[thread_id].col(nocc).swap(SWAPPED_CD_p0[thread_id].col(nocc+c));
        SWAPPED_CD_m0[thread_id].col(nocc).swap(SWAPPED_CD_m0[thread_id].col(nocc+c));
        SWAPPED_CD_0p[thread_id].col(nocc).swap(SWAPPED_CD_0p[thread_id].col(nocc+c));
        SWAPPED_CD_0m[thread_id].col(nocc).swap(SWAPPED_CD_0m[thread_id].col(nocc+c));
     
     
        Prod_00_0p[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_0p * SWAPPED_CD_0p[thread_id];
        Prod_00_0m[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_0m * SWAPPED_CD_0m[thread_id];
        Prod_00_p0[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_p0 * SWAPPED_CD_p0[thread_id];
        Prod_00_m0[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_m0 * SWAPPED_CD_m0[thread_id];
       
        
        double OvLp_AC_00_0p =
          Prod_00_0p[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_AC_00_0m =
          Prod_00_0m[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_AC_00_p0 =
          Prod_00_p0[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_AC_00_m0 =
          Prod_00_m0[thread_id].block(0,0,nocc+1,nocc+1).determinant();
     
     
        // BD Overlap
        SWAPPED_AB_00[thread_id].setZero();
        SWAPPED_CD_p0[thread_id].setZero();
        SWAPPED_CD_m0[thread_id].setZero();
        SWAPPED_CD_0p[thread_id].setZero();
        SWAPPED_CD_0m[thread_id].setZero();
     
        Prod_00_0p[thread_id].setZero();
        Prod_00_0m[thread_id].setZero();
        Prod_00_p0[thread_id].setZero();
        Prod_00_m0[thread_id].setZero();
     
        SWAPPED_AB_00[thread_id] = MO_00;
        SWAPPED_CD_p0[thread_id] = MO_p0;
        SWAPPED_CD_m0[thread_id] = MO_m0;
        SWAPPED_CD_0p[thread_id] = MO_0p;
        SWAPPED_CD_0m[thread_id] = MO_0m;
     
        SWAPPED_AB_00[thread_id].col(nocc).swap(SWAPPED_AB_00[thread_id].col(nocc+b));
        SWAPPED_CD_p0[thread_id].col(nocc).swap(SWAPPED_CD_p0[thread_id].col(nocc+d));
        SWAPPED_CD_m0[thread_id].col(nocc).swap(SWAPPED_CD_m0[thread_id].col(nocc+d));
        SWAPPED_CD_0p[thread_id].col(nocc).swap(SWAPPED_CD_0p[thread_id].col(nocc+d));
        SWAPPED_CD_0m[thread_id].col(nocc).swap(SWAPPED_CD_0m[thread_id].col(nocc+d));
     
     
        Prod_00_0p[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_0p * SWAPPED_CD_0p[thread_id];
        Prod_00_0m[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_0m * SWAPPED_CD_0m[thread_id];
        Prod_00_p0[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_p0 * SWAPPED_CD_p0[thread_id];
        Prod_00_m0[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_m0 * SWAPPED_CD_m0[thread_id];
       
        
        double OvLp_BD_00_0p =
          Prod_00_0p[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_BD_00_0m =
          Prod_00_0m[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_BD_00_p0 =
          Prod_00_p0[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_BD_00_m0 =
          Prod_00_m0[thread_id].block(0,0,nocc+1,nocc+1).determinant();
     
        // AD Overlap
        SWAPPED_AB_00[thread_id].setZero();
        SWAPPED_CD_p0[thread_id].setZero();
        SWAPPED_CD_m0[thread_id].setZero();
        SWAPPED_CD_0p[thread_id].setZero();
        SWAPPED_CD_0m[thread_id].setZero();
     
        Prod_00_0p[thread_id].setZero();
        Prod_00_0m[thread_id].setZero();
        Prod_00_p0[thread_id].setZero();
        Prod_00_m0[thread_id].setZero();
     
        SWAPPED_AB_00[thread_id] = MO_00;
        SWAPPED_CD_p0[thread_id] = MO_p0;
        SWAPPED_CD_m0[thread_id] = MO_m0;
        SWAPPED_CD_0p[thread_id] = MO_0p;
        SWAPPED_CD_0m[thread_id] = MO_0m;
     
        SWAPPED_AB_00[thread_id].col(nocc).swap(SWAPPED_AB_00[thread_id].col(nocc+a));
        SWAPPED_CD_p0[thread_id].col(nocc).swap(SWAPPED_CD_p0[thread_id].col(nocc+d));
        SWAPPED_CD_m0[thread_id].col(nocc).swap(SWAPPED_CD_m0[thread_id].col(nocc+d));
        SWAPPED_CD_0p[thread_id].col(nocc).swap(SWAPPED_CD_0p[thread_id].col(nocc+d));
        SWAPPED_CD_0m[thread_id].col(nocc).swap(SWAPPED_CD_0m[thread_id].col(nocc+d));
     
     
        Prod_00_0p[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_0p * SWAPPED_CD_0p[thread_id];
        Prod_00_0m[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_0m * SWAPPED_CD_0m[thread_id];
        Prod_00_p0[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_p0 * SWAPPED_CD_p0[thread_id];
        Prod_00_m0[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_m0 * SWAPPED_CD_m0[thread_id];
       
        
        double OvLp_AD_00_0p =
          Prod_00_0p[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_AD_00_0m =
          Prod_00_0m[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_AD_00_p0 =
          Prod_00_p0[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_AD_00_m0 =
          Prod_00_m0[thread_id].block(0,0,nocc+1,nocc+1).determinant();
     
     
        // BC Overlap
        SWAPPED_AB_00[thread_id].setZero();
        SWAPPED_CD_p0[thread_id].setZero();
        SWAPPED_CD_m0[thread_id].setZero();
        SWAPPED_CD_0p[thread_id].setZero();
        SWAPPED_CD_0m[thread_id].setZero();
     
        Prod_00_0p[thread_id].setZero();
        Prod_00_0m[thread_id].setZero();
        Prod_00_p0[thread_id].setZero();
        Prod_00_m0[thread_id].setZero();
     
        SWAPPED_AB_00[thread_id] = MO_00;
        SWAPPED_CD_p0[thread_id] = MO_p0;
        SWAPPED_CD_m0[thread_id] = MO_m0;
        SWAPPED_CD_0p[thread_id] = MO_0p;
        SWAPPED_CD_0m[thread_id] = MO_0m;
     
        SWAPPED_AB_00[thread_id].col(nocc).swap(SWAPPED_AB_00[thread_id].col(nocc+b));
        SWAPPED_CD_p0[thread_id].col(nocc).swap(SWAPPED_CD_p0[thread_id].col(nocc+c));
        SWAPPED_CD_m0[thread_id].col(nocc).swap(SWAPPED_CD_m0[thread_id].col(nocc+c));
        SWAPPED_CD_0p[thread_id].col(nocc).swap(SWAPPED_CD_0p[thread_id].col(nocc+c));
        SWAPPED_CD_0m[thread_id].col(nocc).swap(SWAPPED_CD_0m[thread_id].col(nocc+c));
     
     
        Prod_00_0p[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_0p * SWAPPED_CD_0p[thread_id];
        Prod_00_0m[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_0m * SWAPPED_CD_0m[thread_id];
        Prod_00_p0[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_p0 * SWAPPED_CD_p0[thread_id];
        Prod_00_m0[thread_id] = 
          SWAPPED_AB_00[thread_id].transpose() * S_00_m0 * SWAPPED_CD_m0[thread_id];
       
        
        double OvLp_BC_00_0p =
          Prod_00_0p[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_BC_00_0m =
          Prod_00_0m[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_BC_00_p0 =
          Prod_00_p0[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_BC_00_m0 =
          Prod_00_m0[thread_id].block(0,0,nocc+1,nocc+1).determinant();
     
     
        double OvLp_Swap_00_p0 = OvLp_AC_00_p0*OvLp_BD_00_p0 + OvLp_AD_00_p0*OvLp_BC_00_p0;
        double OvLp_Swap_00_m0 = OvLp_AC_00_m0*OvLp_BD_00_m0 + OvLp_AD_00_m0*OvLp_BC_00_m0;
        double OvLp_Swap_00_0p = OvLp_AC_00_0p*OvLp_BD_00_0p + OvLp_AD_00_0p*OvLp_BC_00_0p;
        double OvLp_Swap_00_0m = OvLp_AC_00_0m*OvLp_BD_00_0m + OvLp_AD_00_0m*OvLp_BC_00_0m;
     
        double fact = 1;
        if( a == b ) fact *= std::sqrt(0.5);
        if( c == d ) fact *= std::sqrt(0.5);
     
        #pragma omp critical
        {
          NACME[0](iSt,jSt) += fact*T_00(ab,iSt)*(
            T_p0(cd,jSt)*OvLp_Swap_00_p0 - 
            T_m0(cd,jSt)*OvLp_Swap_00_m0
          ) / (2*step1);
          NACME[1](iSt,jSt) += fact*T_00(ab,iSt)*(
            T_0p(cd,jSt)*OvLp_Swap_00_0p - 
            T_0m(cd,jSt)*OvLp_Swap_00_0m
          ) / (2*step2);
        }
      }
    }
  }
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "ES -> ES NACME took " << elapsed.count() << " s" << endl;


  return NACME;
};
