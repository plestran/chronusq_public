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

template <typename T, MOLECULE_PRESETS M> 
void runCQJob(H5::Group &res, std::string &fName, CQMemManager &memManager, 
  const std::string &jbTyp, const std::string &basisSet, 
  const std::string &ref, const SCFSettings scfSett, int numThreads, const std::string &guess) {

  Molecule molecule;
  BasisSet basis;
  AOIntegrals aoints;
  SingleSlater<T> singleSlater;
  FileIO fileio("test.inp","test.out");

  fileio.iniH5Files();

  CQSetNumThreads(numThreads);

  loadPresets<M>(molecule);
  molecule.convBohr(); molecule.computeNucRep(); molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(ref);
  singleSlater.setGuess(guess);
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
  
  fileio.out.close();


  if(!jbTyp.compare("SCF")) writeSCFRecord(fName,res,singleSlater);
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
  std::vector<std::string> jobs = {"SCF"};

  // DFT Functionals to test
  std::vector<std::string> KS {
    "SLATER","B88","LSDA","SVWN5","BLYP","B3LYP","BHandH"
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
  std::vector<std::string> gRefs = {"GHF"};

  // Add the DFT counterparts into the reference lists
  rRefs.insert(rRefs.end(),RKSL.begin(),RKSL.end());
  uRefs.insert(uRefs.end(),UKSL.begin(),UKSL.end());

  // Concatinate all reference lists into one list
  std::vector<std::string> refs;
  auto it = refs.end();
  it = refs.insert(it,rRefs.begin(),rRefs.end());
  it = refs.insert(it,uRefs.begin(),uRefs.end());
  it = refs.insert(it,gRefs.begin(),gRefs.end());


  // Which basis sets to test
  std::vector<std::string> bases = {"STO-3G","6-31G","cc-pVDZ"};

  // Field Types
  std::vector<std::string> fieldtps = 
//    {"NOFIELD","WEAKXFIELD","WEAKYFIELD","WEAKZFIELD"};
    {"NOFIELD"};

  // SCF Settings
  SCFSettings defaultSCF{std::array<double,3>({0.0,0.0,0.0})};
  SCFSettings weakXFieldSCF{std::array<double,3>({0.0006,0.0,0.0})};
  SCFSettings weakYFieldSCF{std::array<double,3>({0.0,0.0006,0.0})};
  SCFSettings weakZFieldSCF{std::array<double,3>({0.0,0.0,0.0006})};


  int testNum = 1;
  std::vector<H5::Group> groups;  
  for( auto jbTyp : jobs ){
    std::string curJbTyp = "/" + jbTyp;
    groups.emplace_back(RefFile.createGroup(curJbTyp));
  for( auto fieldtyp : fieldtps ) {
    std::string curField = curJbTyp + "/" + fieldtyp;
    groups.emplace_back(RefFile.createGroup(curField));
  for( auto ref : refs ) {
    std::string curRef = curField + "/" + ref;
    groups.emplace_back(RefFile.createGroup(curRef));
  for( auto basis : bases ) {
    std::string curBasis = curRef + "/" + basis;
    groups.emplace_back(RefFile.createGroup(curBasis));

    std::string lastStr = curBasis;

    SCFSettings scfSett;
    if(!fieldtyp.compare("NOFIELD")) scfSett = defaultSCF;
    else if(!fieldtyp.compare("WEAKXFIELD")) scfSett = weakXFieldSCF;
    else if(!fieldtyp.compare("WEAKYFIELD")) scfSett = weakYFieldSCF;
    else if(!fieldtyp.compare("WEAKZFIELD")) scfSett = weakZFieldSCF;

    /** Water Test **/
    std::string fName = lastStr + "/" + moleculeName<WATER>();
    groups.emplace_back(RefFile.createGroup(fName));
    
    cout << "Running Job " << fName;
    runCQJob<double,WATER>(groups.back(),fName,memManager,jbTyp,basis,ref,
      scfSett,1, "CORE");

    std::stringstream linkName;
    linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
    H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
      H5P_DEFAULT,H5P_DEFAULT);
    testNum++;
    cout << " -> " << linkName.str() << endl;

    // Tests that only make sense using unrestricted references
    if(std::find(rRefs.begin(),rRefs.end(),ref) == rRefs.end()){
      // O2 Doesnt like to converge in the presence of a field
      if(!fieldtyp.compare("NOFIELD")) {
        /** O2 Triplet Test **/
        fName = lastStr + "/" + moleculeName<O2>();
        groups.emplace_back(RefFile.createGroup(fName));
       
        cout << "Running Job " << fName;
        runCQJob<double,O2>(groups.back(),fName,memManager,jbTyp,basis,ref,
          scfSett,1,"CORE");
        
        linkName.str("");
        linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
        H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
          H5P_DEFAULT,H5P_DEFAULT);
        testNum++;
        cout << " -> " << linkName.str() << endl;
      }

      /** Li Test **/
      fName = lastStr + "/" + moleculeName<Li>();
      groups.emplace_back(RefFile.createGroup(fName));

      cout << "Running Job " << fName;
      runCQJob<double,Li>(groups.back(),fName,memManager,jbTyp,basis,ref,
        scfSett,1,"CORE");
      
      linkName.str("");
      linkName << "test" <<std::setfill('0') << std::setw(4) << testNum;
      H5Lcreate_soft(fName.c_str(),RefFile.getId(),linkName.str().c_str(),
        H5P_DEFAULT,H5P_DEFAULT);
      testNum++;
      cout << " -> " << linkName.str() << endl;
    }
  }
  }
  }
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
