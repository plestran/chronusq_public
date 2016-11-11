#include <global.h>
#include <pythonapi.h>
#include <singleslater.h>
#include <mointegrals.h>
#include <response.h>
#include <realtime.h>
using namespace ChronusQ;

// Molecule Specification Presets
enum MOLECULE_PRESETS {
  WATER, Li, O2, SingO2
};

template<MOLECULE_PRESETS T> void loadPresets(Molecule&);

template <MOLECULE_PRESETS T> std::string moleculeName();
template <MOLECULE_PRESETS T> std::string moleculeGeom(){
  Molecule mol; loadPresets<T>(mol);
  std::stringstream tmp;

  for(auto iAtm = 0; iAtm < mol.nAtoms(); iAtm++){
     tmp << " " << elements[mol.index(iAtm)].symbol;
     tmp << " " << std::setw(15) << std::setprecision(10) 
                << (*mol.cart())(0,iAtm);
     tmp << " " << std::setw(15) << std::setprecision(10) 
                << (*mol.cart())(1,iAtm);
     tmp << " " << std::setw(15) << std::setprecision(10) 
                << (*mol.cart())(2,iAtm);
     tmp << endl;
  }
  return tmp.str();
};

template<> std::string moleculeName<WATER>() { return std::string("Water"); }; 
template<> std::string moleculeName<Li>() { return std::string("Li"); }; 
template<> std::string moleculeName<O2>() { return std::string("O2"); }; 
template<> std::string moleculeName<SingO2>() { 
  return std::string("Singlet O2"); 
}; 


struct SCFSettings {
  std::array<double,3> staticField;

//SCFSettings() {
//  staticField = {0.0,0.0,0.0};
//};
};

struct RTSettings {
  std::array<double,3> fieldAmp;
};

struct CQJob {
  CQMemManager *memManager;
  std::string  jobTyp;
  std::string  basisSet;
  std::string  ref;
  int          numThreads;
  GUESS        guess;

  SCFSettings scfSett;
  RTSettings  rtSett;
};

template <typename T> void runSCF(SingleSlater<T> &ss) {
  ss.formGuess();
  ss.SCF3();
  ss.computeProperties();
};
template <typename T> void runRT(RealTime<T> &rt) {
  rt.alloc();
  rt.doPropagation();
};

template <typename T> void writeSCFRecord(std::string name, H5::Group &gp, 
  SingleSlater<T> &ss) {

  hsize_t one(1), three(3), nine(9), tseven(27);
  H5::DataSpace One(1,&one);
  H5::DataSpace Three(1,&three);
  H5::DataSpace Nine(1,&nine);
  H5::DataSpace TSeven(1,&tseven);

  H5::DataSet TotalEnergy(gp.createDataSet(name + "/TotalEnergy",
    H5PredType<double>(),One));
  H5::DataSet ElecDipole(gp.createDataSet(name + "/ElecDipole",
    H5PredType<double>(),Three));
  H5::DataSet ElecQuadrupole(gp.createDataSet(name + "/ElecQuadrupole",
    H5PredType<double>(),Nine));
  H5::DataSet ElecTLessQuadrupole(gp.createDataSet(
    name + "/ElecTLessQuadrupole",H5PredType<double>(),Nine));
  H5::DataSet ElecOctupole(gp.createDataSet(name + "/ElecOctupole",
    H5PredType<double>(),TSeven));

  TotalEnergy.write(&ss.totalEnergy(),H5PredType<double>());
  ElecDipole.write(&ss.elecDipole()[0],H5PredType<double>());
  ElecQuadrupole.write(&ss.elecQuadpole()[0],H5PredType<double>());
  ElecTLessQuadrupole.write(&ss.elecTracelessQuadpole()[0],
    H5PredType<double>());
  ElecOctupole.write(&ss.elecOctpole()[0],H5PredType<double>());


  std::ifstream outFile("test.out",ios::binary);
  std::vector<char> buffer((std::istreambuf_iterator<char>(outFile)),
    (std::istreambuf_iterator<char>()));

  hsize_t fileLen(buffer.size());
  H5::DataSpace FileLen(1,&fileLen);

  H5::DataSet OutFile(gp.createDataSet(name + "/Output",
    H5::PredType::NATIVE_CHAR,FileLen));

  OutFile.write(&buffer[0],H5::PredType::NATIVE_CHAR);
}

template <typename T> void writeRTRecord(std::string name, H5::Group &gp, 
  RealTime<T> &rt) {

  typedef struct {
    double timeStep;
    double energy;
    std::array<double,4> dipole;
    std::array<double,4> appliedfield;
  } timept_t;

  typedef struct {
    double x;
    double y;
    double z;
    double t;
  } fourvec_t;

  H5::CompType FourVector(sizeof(fourvec_t));
  FourVector.insertMember("X",HOFFSET(fourvec_t,x),H5::PredType::NATIVE_DOUBLE);
  FourVector.insertMember("Y",HOFFSET(fourvec_t,y),H5::PredType::NATIVE_DOUBLE);
  FourVector.insertMember("Z",HOFFSET(fourvec_t,z),H5::PredType::NATIVE_DOUBLE);
  FourVector.insertMember("T",HOFFSET(fourvec_t,t),H5::PredType::NATIVE_DOUBLE);

  H5::CompType timePt(sizeof(timept_t));

  timePt.insertMember("Time",HOFFSET(timept_t,timeStep),
    H5::PredType::NATIVE_DOUBLE);
  timePt.insertMember("Energy",HOFFSET(timept_t,energy),
    H5::PredType::NATIVE_DOUBLE);
  timePt.insertMember("Dipole",HOFFSET(timept_t,dipole),FourVector);
  timePt.insertMember("Field",HOFFSET(timept_t,appliedfield),FourVector);

  int nPts = rt.propInfo.size();
  hsize_t NPts(nPts);
  H5::DataSpace PropSpace(1,&NPts);


  H5::DataSet PropInfo(gp.createDataSet(name + "/TimePropagation",timePt,PropSpace));


  std::vector<timept_t> tmp;
  for(auto pt : rt.propInfo) {
    tmp.push_back(timept_t{pt.timeStep,pt.energy,pt.dipole,pt.appliedfield});
  }

  PropInfo.write(&tmp[0],timePt);

  std::ifstream outFile("test.out",ios::binary);
  std::vector<char> buffer((std::istreambuf_iterator<char>(outFile)),
    (std::istreambuf_iterator<char>()));

  hsize_t fileLen(buffer.size());
  H5::DataSpace FileLen(1,&fileLen);

  H5::DataSet OutFile(gp.createDataSet(name + "/Output",
    H5::PredType::NATIVE_CHAR,FileLen));

  OutFile.write(&buffer[0],H5::PredType::NATIVE_CHAR);
}

template<typename T,MOLECULE_PRESETS M>
void writeInput(const std::string &baseName, const std::string &jbTyp,
  const std::string &basisSet, const std::string &ref, const SCFSettings scfSett,
  const RTSettings rtSett, int numThreads, const std::string &guess) {

  std::ofstream input(baseName + ".inp");

  std::array<double,3> null = {0,0,0};

  input << "#" << endl;
  input << "#  " << baseName << " - " << moleculeName<M>() << " " << ref;
  input << "/" << basisSet << " : "<< jbTyp;
  if(!jbTyp.compare("SCF") and scfSett.staticField != null) 
    input << " in a field";
  input << endl << "#  ";
  if(numThreads == 1) input << "SERIAL";
  else                input << "SMP";
  input << endl;
  input << "#" << endl;
  input << "#  Molecule Specification " << endl;

  input << "[Molecule]" << endl;
  input << "charge = " << 0 << endl;
  input << "mult = ";
  if(M == WATER) input << 1;
  if(M == O2) input << 3;
  if(M == Li) input << 2;
  input << endl;

  input << "geom: " << endl << moleculeGeom<M>() << endl;
  

  input << "# \n#  Job Specification\n#" << endl;
  input << "[QM]" << endl;
  input << "reference = ";
  if(std::is_same<double,T>::value) input << "Real";
  else if(std::is_same<dcomplex,T>::value) input << "Complex";
  input << " " << ref << endl;

  input << "job = " << jbTyp << endl;
  input << endl;


  input << "[BASIS]" << endl;
  input << "basis = " << basisSet << endl;
  input << endl;  


  input << "[SCF]" << endl;
  if(guess.compare("SAD")) {
    input << "guess = " << guess << endl;
  }
  if(!jbTyp.compare("SCF") and scfSett.staticField != null){ 
    input << "field = ";
    for(auto x : scfSett.staticField)
      input << std::setw(10) << std::setprecision(5) << x;
    input << endl;
    
  }
  input << endl;

  if(numThreads != 0){
     input << "[MISC]" << endl;
     input << "nsmp = " << numThreads << endl;
  }
  input << endl;

}

template <typename T, MOLECULE_PRESETS M> 
void runCQJob(H5::Group &res, std::string &fName, CQMemManager &memManager, 
  const std::string &jbTyp, const std::string &basisSet, 
  const std::string &ref, const SCFSettings scfSett, const RTSettings rtSett,
  int numThreads, const std::string &guess) {

  Molecule molecule;
  BasisSet basis;
  AOIntegrals aoints;
  SingleSlater<T> singleSlater;
  RealTime<T> rt;
  FileIO fileio("test.inp","test.out");

  fileio.iniH5Files();

  CQSetNumThreads(numThreads);

  loadPresets<M>(molecule);
  molecule.convBohr(); molecule.computeNucRep(); molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(ref);
  singleSlater.setGuess(guess);
  if(!jbTyp.compare("SCF"))
    singleSlater.setField(scfSett.staticField);


  basis.findBasisFile(basisSet);
  basis.communicate(fileio);
  basis.parseGlobal();
  basis.constructLocal(&molecule);
  basis.makeMaps(&molecule);
  basis.renormShells();

  aoints.communicate(molecule,basis,fileio,memManager);
  singleSlater.communicate(molecule,basis,aoints,fileio,memManager);

  aoints.initMeta();
  singleSlater.initMeta();
  aoints.alloc();
  singleSlater.alloc();

  runSCF(singleSlater);
  if(!jbTyp.compare("RT")){
    rt.communicate(singleSlater);
    runRT(rt);
  }
  
  fileio.out.close();


  if(!jbTyp.compare("SCF"))     writeSCFRecord(fName,res,singleSlater);
  else if(!jbTyp.compare("RT")) writeRTRecord(fName,res,rt);
}

inline bool file_exists(const std::string& name) {
  std::ifstream f(name.c_str());
  return f.good();
};

inline bool copy_file(const std::string& from_name, const std::string& to_name){
  std::ifstream src(from_name.c_str(),std::ios::binary);
  std::ofstream dst(to_name.c_str()  ,std::ios::binary);

  dst << src.rdbuf();
}

int main(int argc, char **argv){

  CQMemManager memManager;
  memManager.setTotalMem(256e6);
  memManager.allocMem();
  initCQ(argc,argv);
  
  std::string refFileName("chronusq-ref.bin");
  std::string cpyFileName("chronusq-ref-old.bin");

  if(file_exists(refFileName)) copy_file(refFileName,cpyFileName);

  H5::H5File RefFile(refFileName,H5F_ACC_TRUNC);
/*
  H5::Group  SCFTests(RefFile.createGroup("/SCF"));
  H5::Group  RTTests(RefFile.createGroup("/RT"));

  H5::Group  SCFSTOTests(SCFTests.createGroup("/SCF/STO-3G")); 
  H5::Group  SCF631GTests(SCFTests.createGroup("/SCF/6-31G")); 
  H5::Group  SCFccpVDZTests(SCFTests.createGroup("/SCF/cc-pVDZ")); 
*/

  
  // Which job types to test
  std::vector<std::string> jobs = {"SCF","RT"};

  // DFT Functionals to test
  std::vector<std::string> KS {
    "SLATER","B88","LSDA","SVWN5","BLYP","B3LYP","BHANDH"
  };

  // Generate all variants of DFT functionals
  std::vector<std::string> RKSL,UKSL,GKSL,X2CKSL;
  for(auto f : KS) {
    RKSL.emplace_back("R" + f);
    UKSL.emplace_back("U" + f);
    GKSL.emplace_back("G" + f);
    X2CKSL.emplace_back("X2C-" + f);
  }

  // Which references to test
  std::vector<std::string> rRefs = {"RHF"};
  std::vector<std::string> uRefs = {"UHF"};
  std::vector<std::string> gRefs = {"GHF", "X2C"};
//std::vector<std::string> gRefs = { "X2C"};

  // Add the DFT counterparts into the reference lists
  rRefs.insert(rRefs.end(),RKSL.begin(),RKSL.end());
  uRefs.insert(uRefs.end(),UKSL.begin(),UKSL.end());

/*
  // Concatinate all reference lists into one list
  std::vector<std::string> refs;
  auto it = refs.end();
  it = refs.insert(it,rRefs.begin(),rRefs.end());
  it = refs.insert(it,uRefs.begin(),uRefs.end());
  it = refs.insert(it,gRefs.begin(),gRefs.end());
*/


  // Which basis sets to test
  std::vector<std::string> bases = {"STO-3G","6-31G","cc-pVDZ"};
//std::vector<std::string> bases = {"6-31G"};

  // Field Types
  std::vector<std::string> fieldtps = 
    {"NOFIELD","WEAKXFIELD","WEAKYFIELD","WEAKZFIELD"};
//  {"NOFIELD"};

  // SCF Settings
  SCFSettings defaultSCF{std::array<double,3>({0.0,0.0,0.0})};
  SCFSettings weakXFieldSCF{std::array<double,3>({0.0006,0.0,0.0})};
  SCFSettings weakYFieldSCF{std::array<double,3>({0.0,0.0006,0.0})};
  SCFSettings weakZFieldSCF{std::array<double,3>({0.0,0.0,0.0006})};


  int testNum = 1;
//std::vector<H5::Group> groups;  
//for( auto jbTyp : jobs ){
//  std::string curJbTyp = "/" + jbTyp;
//  groups.emplace_back(RefFile.createGroup(curJbTyp));
//for( auto fieldtyp : fieldtps ) {
//  std::string curField = curJbTyp + "/" + fieldtyp;
//  groups.emplace_back(RefFile.createGroup(curField));
//for( auto ref : refs ) {
//  std::string curRef = curField + "/" + ref;
//  groups.emplace_back(RefFile.createGroup(curRef));
//for( auto basis : bases ) {
//  std::string curBasis = curRef + "/" + basis;
//  groups.emplace_back(RefFile.createGroup(curBasis));

//  std::string lastStr = curBasis;

//  SCFSettings scfSett;
//  if(!fieldtyp.compare("NOFIELD")) scfSett = defaultSCF;
//  else if(!fieldtyp.compare("WEAKXFIELD")) scfSett = weakXFieldSCF;
//  else if(!fieldtyp.compare("WEAKYFIELD")) scfSett = weakYFieldSCF;
//  else if(!fieldtyp.compare("WEAKZFIELD")) scfSett = weakZFieldSCF;

//  /** Water Test **/
//  std::string fName = lastStr + "/" + moleculeName<WATER>();
//  groups.emplace_back(RefFile.createGroup(fName));
//  
//  cout << "Running Job " << fName;
//  runCQJob<double,WATER>(groups.back(),fName,memManager,jbTyp,basis,ref,
//    scfSett,1, "CORE");

//  std::stringstream linkName;
//  linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
//  H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
//    H5P_DEFAULT,H5P_DEFAULT);
//  testNum++;
//  cout << " -> " << linkName.str() << endl;

//  // Tests that only make sense using unrestricted references
//  if(std::find(rRefs.begin(),rRefs.end(),ref) == rRefs.end()){
//    // O2 Doesnt like to converge in the presence of a field
//    if(!fieldtyp.compare("NOFIELD")) {
//      /** O2 Triplet Test **/
//      fName = lastStr + "/" + moleculeName<O2>();
//      groups.emplace_back(RefFile.createGroup(fName));
//     
//      cout << "Running Job " << fName;
//      runCQJob<double,O2>(groups.back(),fName,memManager,jbTyp,basis,ref,
//        scfSett,1,"CORE");
//      
//      linkName.str("");
//      linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
//      H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
//        H5P_DEFAULT,H5P_DEFAULT);
//      testNum++;
//      cout << " -> " << linkName.str() << endl;
//    }

//    /** Li Test **/
//    fName = lastStr + "/" + moleculeName<Li>();
//    groups.emplace_back(RefFile.createGroup(fName));

//    cout << "Running Job " << fName;
//    runCQJob<double,Li>(groups.back(),fName,memManager,jbTyp,basis,ref,
//      scfSett,1,"CORE");
//    
//    linkName.str("");
//    linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
//    H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
//      H5P_DEFAULT,H5P_DEFAULT);
//    testNum++;
//    cout << " -> " << linkName.str() << endl;
//  }
//}
//}
//}
//}
  
  // Define reference to test for real and complex
  std::vector<std::string> realRefs,complexRefs;
  // RHF, UHF, RKS, UKS for Real
  auto it = realRefs.end();
  it = realRefs.insert(it,rRefs.begin(),rRefs.end());
  it = realRefs.insert(it,uRefs.begin(),uRefs.end());

  // RHF, UHF, GHF, X2C, RKS, UKS for Complex
  it = complexRefs.end();
  it = complexRefs.insert(it,rRefs.begin(),rRefs.end());
  it = complexRefs.insert(it,uRefs.begin(),uRefs.end());
  it = complexRefs.insert(it,gRefs.begin(),gRefs.end());


  H5::Group realTests(RefFile.createGroup("/REAL"));
  H5::Group complexTests(RefFile.createGroup("/COMPLEX"));

  // Create Groups
  for( auto fld : {std::string("REAL"),std::string("COMPLEX")} )
  for( auto jbTyp : jobs ){
    std::string curJbTyp = "/" + fld + "/" + jbTyp;
    H5::Group jbGroup(RefFile.createGroup(curJbTyp));
  for( auto para : {std::string("SERIAL"),std::string("SMP")}) {
    std::string curPara = curJbTyp + "/" + para;
    H5::Group paraGroup(RefFile.createGroup(curPara));
  for( auto fieldtyp : fieldtps ) {
    std::string curField = curPara + "/" + fieldtyp;
    H5::Group fldGroup(RefFile.createGroup(curField));
    std::vector<std::string> curRefs;
    if(!fld.compare("REAL")) curRefs = realRefs;
    else curRefs = complexRefs;
  for( auto ref : curRefs ) {
    std::string curRef = curField + "/" + ref;
    H5::Group refGroup(RefFile.createGroup(curRef));
  for( auto basis : bases ) {
    std::string curBasis = curRef + "/" + basis;
    H5::Group basisGroup(RefFile.createGroup(curBasis));

/*
    if(!basis.compare("cc-pVDZ") and !ref.compare("X2C"))
      continue;

    std::string lastStr = curBasis;

    SCFSettings scfSett;
    if(!fieldtyp.compare("NOFIELD")) scfSett = defaultSCF;
    else if(!fieldtyp.compare("WEAKXFIELD")) scfSett = weakXFieldSCF;
    else if(!fieldtyp.compare("WEAKYFIELD")) scfSett = weakYFieldSCF;
    else if(!fieldtyp.compare("WEAKZFIELD")) scfSett = weakZFieldSCF;

    int nCores = 1;
    if(!para.compare("SMP")) nCores = 2;

    std::string fName = lastStr + "/" + moleculeName<WATER>();
    H5::Group waterGroup(RefFile.createGroup(fName));
    
    cout << "Running Job " << fName;
    if(!fld.compare("REAL"))
      runCQJob<double,WATER>(waterGroup,fName,memManager,"SCF",basis,ref,
        scfSett,nCores, "CORE");
    else 
      runCQJob<dcomplex,WATER>(waterGroup,fName,memManager,"SCF",basis,ref,
        scfSett,nCores, "CORE");

    std::stringstream linkName;
    linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
    H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
      H5P_DEFAULT,H5P_DEFAULT);
    testNum++;
    cout << " -> " << linkName.str() << endl;
*/
  }
  }
  }
  }
  }


  // Water Tests
  for( auto fld : {std::string("REAL"),std::string("COMPLEX")} )
  for( auto jbTyp : {std::string("SCF"),std::string("RT")} )
  for( auto para : {std::string("SERIAL"),std::string("SMP")}) 
  for( auto fieldtyp : fieldtps ) 
  for( auto ref : rRefs )
  for( auto basis : bases ) {

    int nCores = 1;
    if(!para.compare("SMP")) nCores = 2;

    std::string fName = "/" + fld + "/" + jbTyp + "/" + para + "/" + fieldtyp 
                        + "/" + ref + "/" + basis + "/" + moleculeName<WATER>();
    H5::Group waterGroup(RefFile.createGroup(fName));

    if(!jbTyp.compare("RT") and fieldtyp.compare("NOFIELD")) continue;
    if(!jbTyp.compare("RT") and 
       (ref.compare("RHF") and ref.compare("RLSDA") and ref.compare("RBLYP")
         and ref.compare("RB3LYP"))) continue;

    SCFSettings scfSett;
    if(!fieldtyp.compare("NOFIELD")) scfSett = defaultSCF;
    else if(!fieldtyp.compare("WEAKXFIELD")) scfSett = weakXFieldSCF;
    else if(!fieldtyp.compare("WEAKYFIELD")) scfSett = weakYFieldSCF;
    else if(!fieldtyp.compare("WEAKZFIELD")) scfSett = weakZFieldSCF;

    RTSettings rtSett;

    cout << "Running Job " << fName;
    if(!fld.compare("REAL"))
      runCQJob<double,WATER>(waterGroup,fName,memManager,jbTyp,basis,ref,
        scfSett,rtSett,nCores, "SAD");
    else 
      runCQJob<dcomplex,WATER>(waterGroup,fName,memManager,jbTyp,basis,ref,
        scfSett,rtSett,nCores, "SAD");


    std::stringstream linkName;
    linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
    H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
      H5P_DEFAULT,H5P_DEFAULT);
    testNum++;
    cout << " -> " << linkName.str() << endl;
    if(!fld.compare("REAL"))
      writeInput<double,WATER>(linkName.str(),jbTyp,basis,ref,scfSett,rtSett,nCores,"SAD");
    else 
      writeInput<dcomplex,WATER>(linkName.str(),jbTyp,basis,ref,scfSett,rtSett,nCores,"SAD");
  }

  // Li Tests
  for( auto fld : {std::string("REAL"),std::string("COMPLEX")} )
  for( auto jbTyp : {std::string("SCF"),std::string("RT")} )
  for( auto para : {std::string("SERIAL"),std::string("SMP")}) 
  for( auto fieldtyp : {std::string("NOFIELD"),std::string("WEAKZFIELD")}) 
  for( auto ref : uRefs )
  for( auto basis : bases ) {

    int nCores = 1;
    if(!para.compare("SMP")) nCores = 2;

    std::string fName = "/" + fld + "/" + jbTyp + "/" + para + "/" + fieldtyp 
                        + "/" + ref + "/" + basis + "/" + moleculeName<Li>();
    H5::Group waterGroup(RefFile.createGroup(fName));

    if(!jbTyp.compare("RT") and fieldtyp.compare("NOFIELD")) continue;
    if(!jbTyp.compare("RT") and 
       (ref.compare("UHF") and ref.compare("ULSDA") and ref.compare("UBLYP")
         and ref.compare("UB3LYP"))) continue;

    SCFSettings scfSett;
    if(!fieldtyp.compare("NOFIELD")) scfSett = defaultSCF;
    else if(!fieldtyp.compare("WEAKXFIELD")) scfSett = weakXFieldSCF;
    else if(!fieldtyp.compare("WEAKYFIELD")) scfSett = weakYFieldSCF;
    else if(!fieldtyp.compare("WEAKZFIELD")) scfSett = weakZFieldSCF;

    RTSettings rtSett;

    cout << "Running Job " << fName;
    if(!fld.compare("REAL"))
      runCQJob<double,Li>(waterGroup,fName,memManager,jbTyp,basis,ref,
        scfSett,rtSett,nCores, "CORE");
    else 
      runCQJob<dcomplex,Li>(waterGroup,fName,memManager,jbTyp,basis,ref,
        scfSett,rtSett,nCores, "CORE");

    std::stringstream linkName;
    linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
    H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
      H5P_DEFAULT,H5P_DEFAULT);
    testNum++;
    cout << " -> " << linkName.str() << endl;
    if(!fld.compare("REAL"))
      writeInput<double,Li>(linkName.str(),jbTyp,basis,ref,scfSett,rtSett,nCores,"CORE");
    else 
      writeInput<dcomplex,Li>(linkName.str(),jbTyp,basis,ref,scfSett,rtSett,nCores,"CORE");
  }

  // O2 Tests
  for( auto fld : {std::string("REAL"),std::string("COMPLEX")} )
  for( auto jbTyp : {std::string("SCF"),std::string("RT")} )
  for( auto para : {std::string("SERIAL"),std::string("SMP")}) 
  for( auto fieldtyp : {std::string("NOFIELD"),std::string("WEAKXFIELD"),std::string("WEAKZFIELD")}) 
  for( auto ref : uRefs )
  for( auto basis : bases ) {

    int nCores = 1;
    if(!para.compare("SMP")) nCores = 2;

    std::string fName = "/" + fld + "/" + jbTyp + "/" + para + "/" + fieldtyp 
                        + "/" + ref + "/" + basis + "/" + moleculeName<O2>();
    H5::Group waterGroup(RefFile.createGroup(fName));

    if(!jbTyp.compare("RT") and fieldtyp.compare("NOFIELD")) continue;
    if(!jbTyp.compare("RT") and 
       (ref.compare("UHF") and ref.compare("ULSDA") and ref.compare("UBLYP")
         and ref.compare("UB3LYP"))) continue;

    SCFSettings scfSett;
    if(!fieldtyp.compare("NOFIELD")) scfSett = defaultSCF;
    else if(!fieldtyp.compare("WEAKXFIELD")) scfSett = weakXFieldSCF;
    else if(!fieldtyp.compare("WEAKYFIELD")) scfSett = weakYFieldSCF;
    else if(!fieldtyp.compare("WEAKZFIELD")) scfSett = weakZFieldSCF;

    RTSettings rtSett;

    cout << "Running Job " << fName;
    std::stringstream linkName;
    try {
    if(!fld.compare("REAL"))
      runCQJob<double,O2>(waterGroup,fName,memManager,jbTyp,basis,ref,
        scfSett,rtSett,nCores, "SAD");
    else 
      runCQJob<dcomplex,O2>(waterGroup,fName,memManager,jbTyp,basis,ref,
        scfSett,rtSett,nCores, "SAD");

    linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
    H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
      H5P_DEFAULT,H5P_DEFAULT);
    testNum++;
    cout << " -> " << linkName.str() << endl;
    } catch(...) {
      cout << " FAILED! " << endl;
    }
    if(!fld.compare("REAL"))
      writeInput<double,O2>(linkName.str(),jbTyp,basis,ref,scfSett,rtSett,nCores,"SAD");
    else 
      writeInput<dcomplex,O2>(linkName.str(),jbTyp,basis,ref,scfSett,rtSett,nCores,"SAD");
  }

/*
  H5::DataSet tmp(RefFile.openDataSet("test0001/Output"));
  std::vector<char> buffer(tmp.getStorageSize());
  
  tmp.read(&buffer[0],H5::PredType::NATIVE_CHAR);

  for(auto X : buffer) cout << X ;
*/
  
  finalizeCQ();
  return 0;
}



















template<>
void loadPresets<WATER>(Molecule &mol) {
  mol.setNAtoms(3);
  mol.setCharge(0);
  mol.setNTotalE(10);
  mol.setMultip(1);
  mol.alloc();
  mol.setIndex(0,HashAtom("O",0));
  mol.setIndex(1,HashAtom("H",0));
  mol.setIndex(2,HashAtom("H",0));
  mol.setCart(0,0.000000000 ,-0.07579184359, 0.0);
  mol.setCart(1,0.866811829 ,0.6014357793  ,0.0);
  mol.setCart(2,-0.866811829, 0.6014357793 ,0.0);
};

template<>
void loadPresets<Li>(Molecule &mol){
  mol.setNAtoms(1);
  mol.setCharge(0);
  mol.setNTotalE(3);
  mol.setMultip(2);
  mol.alloc();
  mol.setIndex(0,HashAtom("Li",0));
  mol.setCart(0,0.000000000 ,0.00000000000, 0.0);
};


template<>
void loadPresets<O2>(Molecule &mol){
  mol.setNAtoms(2);
  mol.setCharge(0);
  mol.setNTotalE(16);
  mol.setMultip(3);
  mol.alloc();
  mol.setIndex(0,HashAtom("O",0));
  mol.setIndex(1,HashAtom("O",0));
  mol.setCart(0,0.000000000 ,0.00000000000, 0.608586  );
  mol.setCart(1,0.000000000 ,0.00000000000, -0.608586 );
};

template<>
void loadPresets<SingO2>(Molecule &mol){
  mol.setNAtoms(2);
  mol.setCharge(0);
  mol.setNTotalE(16);
  mol.setMultip(1);
  mol.alloc();
  mol.setIndex(0,HashAtom("O",0));
  mol.setIndex(1,HashAtom("O",0));
  mol.setCart(0,0.000000000 ,0.00000000000, 0.608586  );
  mol.setCart(1,0.000000000 ,0.00000000000, -0.608586 );
};
