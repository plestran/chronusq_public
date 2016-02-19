#include <response.h>
#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;
std::string genFName(double,double,int,int,std::string);
Eigen::VectorXd checkPhase(double*,double*,int,int);
RealMatrix genSpx(BasisSet&,BasisSet&);
std::vector<Eigen::VectorXd> ES_GS_NACME(int,bool,int,int,double,
  double,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&);
std::vector<Eigen::VectorXd> GS_ES_NACME(int,bool,int,int,double,
  double,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,
  RealMatrix&,RealMatrix&);

template<typename Derived>
double diffNorm( const Derived &A, const Derived &B ){
  return (A-B).norm();
}
template<typename Derived>
double diffNormI( const Derived &A ){
  return (A - Derived::Identity(A.rows(),A.cols())).norm();
}
template<typename Derived>
double selfInner( const Derived &A){
  return A.dot(A);
}



void twoDScan(std::vector<double>& scanX, std::vector<double>& scanY){

  std::string PREFIX = "H20";
  std::string BINPREFIX = "BIN";
  std::string XYZPREFIX = "XYZ";

  initCQ();
  Controls controls;
  CQSetNumThreads(1);
  controls.iniControls();

  BasisSet globalBasis;
  std::string basisName = "cc-pVDZ";
  int charge = 0;
  int nFreq = 8;
  bool debug = true;

  FileIO fileio_00("test_00.inp",PREFIX+".out_00",PREFIX+".rst_00");
  FileIO fileio_p0("test_p0.inp",PREFIX+".out_p0",PREFIX+".rst_p0");
  FileIO fileio_m0("test_m0.inp",PREFIX+".out_m0",PREFIX+".rst_m0");
  FileIO fileio_0p("test_0p.inp",PREFIX+".out_0p",PREFIX+".rst_0p");
  FileIO fileio_0m("test_0m.inp",PREFIX+".out_0m",PREFIX+".rst_0m");

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


  for(auto IX = 1; IX < (scanX.size()-1); IX++)
  for(auto JY = 1; JY < (scanY.size()-1); JY++){
    // Output header
    cout << bannerTop << endl;
    cout << "Starting Numerical Differentiation for:" << endl;
    cout << "  IX = " << IX << "  JY = " << JY << endl;

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
    std::string fBaseName00 = genFName(X ,Y ,3,3,PREFIX);
    std::string fBaseNamep0 = genFName(XP,Y ,3,3,PREFIX);
    std::string fBaseNamem0 = genFName(XM,Y ,3,3,PREFIX);
    std::string fBaseName0p = genFName(X ,YP,3,3,PREFIX);
    std::string fBaseName0m = genFName(X ,YM,3,3,PREFIX);

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
   
    basis00.makeMaps(1,&geom_00);
    basisp0.makeMaps(1,&geom_p0);
    basism0.makeMaps(1,&geom_m0);
    basis0p.makeMaps(1,&geom_0p);
    basis0m.makeMaps(1,&geom_0m);

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

    moints_00.communicate(geom_00,basis00,fileio_00,controls,aoints_00, ss_00);
    moints_p0.communicate(geom_p0,basisp0,fileio_p0,controls,aoints_p0, ss_p0);
    moints_m0.communicate(geom_m0,basism0,fileio_m0,controls,aoints_m0, ss_m0);
    moints_0p.communicate(geom_0p,basis0p,fileio_0p,controls,aoints_0p, ss_0p);
    moints_0m.communicate(geom_0m,basis0m,fileio_0m,controls,aoints_0m, ss_0m);
    
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
    cout << "Performing SCF (X,Y)" << endl;
    ss_00.formGuess();
    ss_00.formFock();
    ss_00.computeEnergy();
    ss_00.SCF();
    ss_00.computeProperties();
    ss_00.printProperties();

    // SCF (X+DX,Y)
    cout << "Performing SCF (X+DX,Y)" << endl;
    ss_p0.formGuess();
    ss_p0.formFock();
    ss_p0.computeEnergy();
    ss_p0.SCF();
    ss_p0.computeProperties();
    ss_p0.printProperties();

    // SCF (X-DX,Y)
    cout << "Performing SCF (X-DX,Y)" << endl;
    ss_m0.formGuess();
    ss_m0.formFock();
    ss_m0.computeEnergy();
    ss_m0.SCF();
    ss_m0.computeProperties();
    ss_m0.printProperties();

    // SCF (X,Y+DY)
    cout << "Performing SCF (X,Y+DY)" << endl;
    ss_0p.formGuess();
    ss_0p.formFock();
    ss_0p.computeEnergy();
    ss_0p.SCF();
    ss_0p.computeProperties();
    ss_0p.printProperties();

    // SCF (X,Y-DY)
    cout << "Performing SCF (X,Y-DY)" << endl;
    ss_0m.formGuess();
    ss_0m.formFock();
    ss_0m.computeEnergy();
    ss_0m.SCF();
    ss_0m.computeProperties();
    ss_0m.printProperties();

    if(debug) {
      cout << endl;
      cout << "  Checking | C - C' | Before Phase Check:" << endl;
      
      cout << "  | C(X,Y) - C(X+DX,Y) | = " << diffNorm((*ss_00.moA()),(*ss_p0.moA())) << endl;  
      cout << "  | C(X,Y) - C(X-DX,Y) | = " << diffNorm((*ss_00.moA()),(*ss_m0.moA())) << endl;  
      cout << "  | C(X,Y) - C(X,Y+DY) | = " << diffNorm((*ss_00.moA()),(*ss_0p.moA())) << endl;  
      cout << "  | C(X,Y) - C(X,Y-DY) | = " << diffNorm((*ss_00.moA()),(*ss_0m.moA())) << endl;  
    }

    Eigen::VectorXd Trans_mo_0p;
    Eigen::VectorXd Trans_mo_0m;
    Eigen::VectorXd Trans_mo_p0;
    Eigen::VectorXd Trans_mo_m0;
    Trans_mo_p0 = checkPhase(ss_00.moA()->data(),ss_p0.moA()->data(),basis00.nBasis(),basis00.nBasis());
    Trans_mo_m0 = checkPhase(ss_00.moA()->data(),ss_m0.moA()->data(),basis00.nBasis(),basis00.nBasis());
    Trans_mo_0p = checkPhase(ss_00.moA()->data(),ss_0p.moA()->data(),basis00.nBasis(),basis00.nBasis());
    Trans_mo_0m = checkPhase(ss_00.moA()->data(),ss_0m.moA()->data(),basis00.nBasis(),basis00.nBasis());

    if(debug) {
      cout << endl;
      cout << "  Checking | C - C' | After Phase Check:" << endl;
      
      cout << "  | C(X,Y) - C(X+DX,Y) | = " << diffNorm((*ss_00.moA()),(*ss_p0.moA())) << endl;  
      cout << "  | C(X,Y) - C(X-DX,Y) | = " << diffNorm((*ss_00.moA()),(*ss_m0.moA())) << endl;  
      cout << "  | C(X,Y) - C(X,Y+DY) | = " << diffNorm((*ss_00.moA()),(*ss_0p.moA())) << endl;  
      cout << "  | C(X,Y) - C(X,Y-DY) | = " << diffNorm((*ss_00.moA()),(*ss_0m.moA())) << endl;  
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

    resp_00.setMeth(RESPONSE_TYPE::CIS);
    resp_p0.setMeth(RESPONSE_TYPE::CIS);
    resp_m0.setMeth(RESPONSE_TYPE::CIS);
    resp_0p.setMeth(RESPONSE_TYPE::CIS);
    resp_0m.setMeth(RESPONSE_TYPE::CIS);

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

    cout << "Performing Response (X,Y)" << endl;
    resp_00.doResponse();
    cout << "Performing Response (X+DX,Y)" << endl;
    resp_p0.doResponse();
    cout << "Performing Response (X-DX,Y)" << endl;
    resp_m0.doResponse();
    cout << "Performing Response (X,Y+DY)" << endl;
    resp_0p.doResponse();
    cout << "Performing Response (X,Y-DY)" << endl;
    resp_0m.doResponse();





    // Start Differentiation things
    //


    
    double scf_00 = ss_00.totalEnergy;
    double scf_p0 = ss_p0.totalEnergy;
    double scf_m0 = ss_m0.totalEnergy;
    double scf_0p = ss_0p.totalEnergy;
    double scf_0m = ss_0m.totalEnergy;

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

    Eigen::VectorXd freq_00 = resp_00.frequencies()[0].head(nFreq);
    Eigen::VectorXd freq_p0 = resp_p0.frequencies()[0].head(nFreq);
    Eigen::VectorXd freq_m0 = resp_m0.frequencies()[0].head(nFreq);
    Eigen::VectorXd freq_0p = resp_0p.frequencies()[0].head(nFreq);
    Eigen::VectorXd freq_0m = resp_0m.frequencies()[0].head(nFreq);

    cout << endl << "  Excitation Frequencies:" << endl;

    for(auto iSt = 0; iSt < nFreq; iSt++){
      cout << "  W(" << iSt << ") (X,Y):    " << std::setprecision(8) 
           << freq_00[iSt] << " Eh" << endl;
      cout << "  W(" << iSt << ") (X+DX,Y): " << std::setprecision(8) 
           << freq_p0[iSt] << " Eh" << endl;
      cout << "  W(" << iSt << ") (X-DX,Y): " << std::setprecision(8) 
           << freq_m0[iSt] << " Eh" << endl;
      cout << "  W(" << iSt << ") (X,Y+DY): " << std::setprecision(8) 
           << freq_0p[iSt] << " Eh" << endl;
      cout << "  W(" << iSt << ") (X,Y-DY): " << std::setprecision(8) 
           << freq_0m[iSt] << " Eh" << endl;
      cout << endl;
    }





    RealMatrix T_00 = resp_00.transDen()[0].block(0,0,resp_00.nMatDim()[0],nFreq);
    RealMatrix T_p0 = resp_p0.transDen()[0].block(0,0,resp_p0.nMatDim()[0],nFreq);
    RealMatrix T_m0 = resp_m0.transDen()[0].block(0,0,resp_m0.nMatDim()[0],nFreq);
    RealMatrix T_0p = resp_0p.transDen()[0].block(0,0,resp_0p.nMatDim()[0],nFreq);
    RealMatrix T_0m = resp_0m.transDen()[0].block(0,0,resp_0m.nMatDim()[0],nFreq);



    if(debug){
      cout << endl;
      cout << "  Checking | T - T' | Before Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " << diffNorm(T_00,T_p0) << endl;  
      cout << "  | T(X,Y) - T(X-DX,Y) | = " << diffNorm(T_00,T_m0) << endl;  
      cout << "  | T(X,Y) - T(X,Y+DY) | = " << diffNorm(T_00,T_0p) << endl;  
      cout << "  | T(X,Y) - T(X,Y-DY) | = " << diffNorm(T_00,T_0m) << endl;  
    }


    Eigen::VectorXd Trans_t_0p;
    Eigen::VectorXd Trans_t_0m;
    Eigen::VectorXd Trans_t_p0;
    Eigen::VectorXd Trans_t_m0;
    Trans_t_p0 = checkPhase(T_00.data(),T_p0.data(),resp_p0.nMatDim()[0],nFreq);
    Trans_t_m0 = checkPhase(T_00.data(),T_m0.data(),resp_m0.nMatDim()[0],nFreq);
    Trans_t_0p = checkPhase(T_00.data(),T_0p.data(),resp_0p.nMatDim()[0],nFreq);
    Trans_t_0m = checkPhase(T_00.data(),T_0m.data(),resp_0m.nMatDim()[0],nFreq);


    if(debug){
      cout << endl;
      cout << "  Checking | T - T' | After Phase Check:" << endl;
      
      cout << "  | T(X,Y) - T(X+DX,Y) | = " << diffNorm(T_00,T_p0) << endl;  
      cout << "  | T(X,Y) - T(X-DX,Y) | = " << diffNorm(T_00,T_m0) << endl;  
      cout << "  | T(X,Y) - T(X,Y+DY) | = " << diffNorm(T_00,T_0p) << endl;  
      cout << "  | T(X,Y) - T(X,Y-DY) | = " << diffNorm(T_00,T_0m) << endl;  
    }

    if(debug) {
      cout << endl << " Checking < T | T > = 1 " << endl;
      for(auto iSt = 0; iSt < nFreq; iSt++){
        cout << "   State " << iSt << ":" << endl;
        cout << "    < T(X,Y) | T(X,Y) >       = " << selfInner(T_00.col(iSt)) << endl;
        cout << "    < T(X+DX,Y) | T(X+DX,Y) > = " << selfInner(T_p0.col(iSt)) << endl;
        cout << "    < T(X-DX,Y) | T(X-DX,Y) > = " << selfInner(T_m0.col(iSt)) << endl;
        cout << "    < T(X,Y+DY) | T(X,Y+DY) > = " << selfInner(T_0p.col(iSt)) << endl;
        cout << "    < T(X,Y-DY) | T(X,Y-DY) > = " << selfInner(T_0m.col(iSt)) << endl;
      }
    }



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







    // Calculate Derivatives
      
    // SCF energy gradient
    double gsdx = (scf_p0 - scf_m0) / (2*(XP-X));
    double gsdy = (scf_0p - scf_0m) / (2*(YP-Y));
    double gsnormd = std::sqrt(gsdx*gsdx + gsdy*gsdy);


    cout << endl << "  GS Gradient = (" << gsdx << "," << gsdy <<
         ")" << endl;

    // Excitation frequency gradient
    Eigen::VectorXd freqDX = (freq_p0 - freq_m0)/(2*(XP-X));
    Eigen::VectorXd freqDY = (freq_0p - freq_0m)/(2*(YP-Y));

    Eigen::VectorXd freqNorm = freqDX.cwiseProduct(freqDX);
    freqNorm += freqDY.cwiseProduct(freqDY);

    for(auto iFreq = 0; iFreq < nFreq; iFreq++)
      freqNorm(iFreq) = std::sqrt(freqNorm(iFreq));

    cout << endl << "  ES Gradients:" << endl;
    for(auto iSt = 0; iSt < nFreq; iSt++){
      cout << "   W(" << iSt << ")' = (" << freqDX(iSt) 
           << "," << freqDY(iSt) << ")" << endl;
    }





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



    if(debug){
      // Quantify deviation from C**H * S * C = I
      cout << endl;
      cout << "  Checking | C(X)**H * S(X,X') * C(X') - I |:" << endl;
     
      cout << "  < (X,Y) | (X,Y) >       = " << diffNormI(SMO_00_00) << endl;
      cout << "  < (X+DX,Y) | (X+DX,Y) > = " << diffNormI(SMO_p0_p0) << endl;
      cout << "  < (X-DX,Y) | (X-DX,Y) > = " << diffNormI(SMO_m0_m0) << endl;
      cout << "  < (X,Y+DY) | (X,Y+DY) > = " << diffNormI(SMO_0p_0p) << endl;
      cout << "  < (X,Y-DY) | (X,Y-DY) > = " << diffNormI(SMO_0m_0m) << endl;
     
     
     
      cout << "  < (X,Y) | (X+DX,Y) >    = " << diffNormI(SMO_00_p0) << endl;
      cout << "  < (X,Y) | (X-DX,Y) >    = " << diffNormI(SMO_00_m0) << endl;
      cout << "  < (X,Y) | (X,Y+DY) >    = " << diffNormI(SMO_00_0p) << endl;
      cout << "  < (X,Y) | (X,Y-DY) >    = " << diffNormI(SMO_00_0m) << endl;
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
   auto NAC = ES_GS_NACME(nFreq,true,NOCC,NVIR,XP-X,YP-Y,
     T_00,(*ss_00.moA()),(*ss_p0.moA()),(*ss_m0.moA()),(*ss_0p.moA()),
     (*ss_0m.moA()),S_00_p0,S_00_m0,S_00_0p,S_00_0m);
   auto NAC2 = GS_ES_NACME(nFreq,true,NOCC,NVIR,XP-X,YP-Y,
     T_00,T_p0,T_m0,T_0p,T_0m,(*ss_00.moA()),(*ss_p0.moA()),
     (*ss_m0.moA()),(*ss_0p.moA()),(*ss_0m.moA()),S_00_p0,S_00_m0,
     S_00_0p,S_00_0m);
   prettyPrint(cout,NAC[0] + NAC2[0],"DIFF DX");
   prettyPrint(cout,NAC[1] + NAC2[1],"DIFF DY");

  }
  finalizeCQ();
}



// Naming convention PREFIX_xptx_ypty.xyz
std::string genFName(double x, double y, int px, int py,
  std::string PREFIX) {
    std::stringstream ss1, ss2;

    ss1 << std::fixed << std::setprecision(3) << x;
    ss2 << std::fixed << std::setprecision(3) << y;

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
std::vector<Eigen::VectorXd> ES_GS_NACME(int nFreq,bool renorm,
  int nocc, int nvir, double step1, double step2, RealMatrix &T_00, 
  RealMatrix &MO_00, RealMatrix &MO_p0, RealMatrix &MO_m0, RealMatrix &MO_0p, 
  RealMatrix &MO_0m, RealMatrix &S_00_p0, RealMatrix &S_00_m0,
  RealMatrix &S_00_0p, RealMatrix &S_00_0m){

  std::vector<Eigen::VectorXd> NACME(2,Eigen::VectorXd(nFreq));

  for(auto i = 0; i < NACME.size(); i++) NACME[i].setZero();

  if(renorm) T_00 *= std::sqrt(0.5);
//T_00.transposeInPlace();
//prettyPrint(cout,T_00,"T");

  RealMatrix SWAPPED_00(MO_00);
  RealMatrix Prod_00_0p(MO_00);
  RealMatrix Prod_00_0m(MO_00);
  RealMatrix Prod_00_p0(MO_00);
  RealMatrix Prod_00_m0(MO_00);
  SWAPPED_00.setZero();
  Prod_00_0p.setZero();
  Prod_00_0m.setZero();
  Prod_00_p0.setZero();
  Prod_00_m0.setZero();
  for(auto iSt = 0; iSt < nFreq; iSt++){
    for(auto i = 0, ia = 0; i < nocc; i++)
    for(auto a = 0; a < nvir; a++, ia++){
      if(std::abs(T_00(ia,iSt)) < 1e-8) continue;
      SWAPPED_00.setZero();
      Prod_00_0p.setZero();
      Prod_00_0m.setZero();
      Prod_00_p0.setZero();
      Prod_00_m0.setZero();
      
      SWAPPED_00 = MO_00;
      SWAPPED_00.col(i).swap(SWAPPED_00.col(nocc+a)); 
      
      Prod_00_0p = SWAPPED_00.transpose() * S_00_0p * MO_0p;
      Prod_00_0m = SWAPPED_00.transpose() * S_00_0m * MO_0m;
      Prod_00_p0 = SWAPPED_00.transpose() * S_00_p0 * MO_p0;
      Prod_00_m0 = SWAPPED_00.transpose() * S_00_m0 * MO_m0;

      double OvLp_IASwap_00_0p =
        Prod_00_0p.block(0,0,nocc,nocc).determinant();
      double OvLp_IASwap_00_0m =
        Prod_00_0m.block(0,0,nocc,nocc).determinant();
      double OvLp_IASwap_00_p0 =
        Prod_00_p0.block(0,0,nocc,nocc).determinant();
      double OvLp_IASwap_00_m0 =
        Prod_00_m0.block(0,0,nocc,nocc).determinant();

      NACME[0](iSt) += T_00(ia,iSt) *
        (OvLp_IASwap_00_p0 - OvLp_IASwap_00_m0) /
        (2*step1);
      NACME[1](iSt) += T_00(ia,iSt) *
        (OvLp_IASwap_00_0p - OvLp_IASwap_00_0m) /
        (2*step2);
    } // loop ove ia
  } // loop over states

  prettyPrint(cout,2*NACME[0]*0.52918,"ES->GS NACME: DX");
  prettyPrint(cout,2*NACME[1]*0.52918,"ES->GS NACME: DY");
  return NACME;
};

std::vector<Eigen::VectorXd> GS_ES_NACME(int nFreq,bool renorm,
  int nocc, int nvir, double step1, double step2, RealMatrix &T_00, 
  RealMatrix &T_p0, RealMatrix &T_m0, RealMatrix &T_0p, RealMatrix &T_0m,
  RealMatrix &MO_00, RealMatrix &MO_p0, RealMatrix &MO_m0, RealMatrix &MO_0p, 
  RealMatrix &MO_0m, RealMatrix &S_00_p0, RealMatrix &S_00_m0,
  RealMatrix &S_00_0p, RealMatrix &S_00_0m){

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

  prettyPrint(cout,2*NACME[0]*0.52918,"GS->ES NACME: DX");
  prettyPrint(cout,2*NACME[1]*0.52918,"GS->ES NACME: DY");
  return NACME;
};
