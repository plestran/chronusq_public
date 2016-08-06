#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;
//using AOIntegrals::INTEGRAL_ALGORITHM;

enum MOLECULE_PRESETS {
  WATER, Li, O2, SingO2
};

template<MOLECULE_PRESETS T>
void loadPresets(Molecule&);

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

struct SCFSettings{
  bool  useDefaultSettings;
  std::array<double,3> field;

};

struct RTSettings{
  bool useDefaultSettings;
};

struct CQJob {
  Molecule *mol;
  int intAlg;
  SCFSettings scfSettings;
  RTSettings  rtSettings;
  std::string field;
  std::string jobType;
  std::string reference;
  std::size_t nProcShared;
  std::string basisSet;

  CQJob(Molecule &mol_, std::string field_, std::string jobtype_,
    std::string reference_, int intAlg_,
    std::size_t nProcShared_, std::string basis_) :
    mol(&mol_),field(std::move(field_)),jobType(std::move(jobtype_)),
    reference(std::move(reference_)),nProcShared(nProcShared_),
    basisSet(std::move(basis_)), scfSettings{true}, rtSettings{true},
    intAlg(intAlg_){ }

  CQJob(Molecule &mol_, std::string field_, std::string jobtype_,
    std::string reference_, int intAlg_,
    std::size_t nProcShared_, std::string basis_, SCFSettings scfset_) :
    CQJob(mol_,field_,jobtype_,reference_,intAlg_,nProcShared_,basis_) 
    { 
      scfSettings.useDefaultSettings = scfset_.useDefaultSettings;
      scfSettings.field = scfset_.field;
    };

  void writeInput(std::ostream &os){
    os << "#"                         << endl; 
    os << "#  Molecule Specification" << endl;
    os << "#"                         << endl ;

    os << "[Molecule]" << endl;
    os << "charge = " << mol->charge() << endl;
    os << "mult = " << mol->multip() << endl;
    os << "geom:" << endl;
    for(auto iAtm = 0; iAtm < mol->nAtoms(); iAtm++){
      os << " " << std::setw(4) << atom[mol->index(iAtm)].symbol;
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        os << std::setw(15) << std::setprecision(8) << std::fixed
           << (*mol->cart())(iXYZ,iAtm);
      }
      os << endl;
    }
    os << endl;
    
    os << "#"                    << endl;
    os << "#  Job Specification" << endl;
    os << "#"                    << endl;
    os << "[QM]" << endl;
    os << "reference = ";
    if(!field.compare("D"))
      os << "Real";
    else
      os << "Complex";
    os << " " << reference << endl;
  
    os << "job = " << jobType << endl;
    os << "basis = " << basisSet << endl;
    os << "ints = ";
    if(intAlg == AOIntegrals::INCORE)
      os << "incore";
    else if (intAlg == AOIntegrals::DIRECT)
      os << "direct";

    os << endl;
    os << endl;
    os << "#"                     << endl;
    os << "#  Misc Specification" << endl;
    os << "#"                     << endl;
    os << "[Misc]" << endl;
    os << "nsmp = " << nProcShared << endl;


  };

  void runJob() {
    CQMemManager memManager;
    BasisSet     basis;
    AOIntegrals  aoints;
    FileIO       fileio("test.inp","test.out");

    memManager.setTotalMem(256e6);
    fileio.iniH5Files();
    fileio.iniStdGroups();

    basis.findBasisFile(basisSet);
    basis.communicate(fileio);
    basis.parseGlobal();
    basis.constructLocal(mol);
    basis.makeMaps(mol);
    basis.renormShells();

    aoints.communicate(*mol,basis,fileio,memManager);
    aoints.initMeta();
    aoints.integralAlgorithm = intAlg;
    aoints.alloc();
    if(!field.compare("D")){
      SingleSlater<double> ss;
      if(!reference.compare("RHF")) {
        ss.setRef(SingleSlater<double>::RHF);
        ss.isClosedShell = true;
      } else if(!reference.compare("UHF")) {
        ss.setRef(SingleSlater<double>::UHF);
        ss.isClosedShell = false;
      }
      ss.communicate(*mol,basis,aoints,fileio,memManager);
      ss.initMeta();
      ss.genMethString();
      ss.alloc();
      ss.formGuess();
      ss.SCF3();
      ss.computeProperties();
    } else {
      SingleSlater<dcomplex> ss;
      if(!reference.compare("RHF")) {
        ss.setRef(SingleSlater<dcomplex>::RHF);
        ss.isClosedShell = true;
      } else if(!reference.compare("UHF")) {
        ss.setRef(SingleSlater<dcomplex>::UHF);
        ss.isClosedShell = false;
      }
      ss.communicate(*mol,basis,aoints,fileio,memManager);
      ss.initMeta();
      ss.genMethString();
      ss.alloc();
      ss.formGuess();
      ss.SCF3();
      ss.computeProperties();
    }

  };


};

std::ostream& operator<< (std::ostream &os, const CQJob& job) {
  return os;
}

int main() {
  std::vector<CQJob> jobs;

  initCQ(0,NULL);
  CQSetNumThreads(1);

  Molecule water,li,o2,singo2;
  loadPresets<WATER>(water);
  loadPresets<Li>(li);
  loadPresets<O2>(o2);
  loadPresets<SingO2>(singo2);


  // Real RHF Serial INCORE
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"6-31G");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ");

  // Real RHF Serial DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ");

  // Real RHF SMP DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ");


  // Real UHF Serial INCORE
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::INCORE,1,"STO-3G");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::INCORE,1,"6-31G");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::INCORE,1,"cc-pVDZ");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::INCORE,1,"STO-3G");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::INCORE,1,"6-31G");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::INCORE,1,"cc-pVDZ");

  // Real UHF Serial DIRECT
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,1,"STO-3G");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,1,"6-31G");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,1,"cc-pVDZ");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,1,"STO-3G");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,1,"6-31G");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,1,"cc-pVDZ");

  // Real UHF SMP DIRECT
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,2,"STO-3G");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,2,"6-31G");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,2,"cc-pVDZ");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,2,"STO-3G");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,2,"6-31G");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,2,"cc-pVDZ");

  
  // Complex RHF Serial INCORE
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::INCORE,1,"6-31G");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ");

  // Complex RHF Serial DIRECT
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ");

  // Complex RHF SMP DIRECT
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ");

  for(auto JOB : jobs){
    JOB.writeInput(cout);
  };

  water.convBohr();
  water.computeNucRep();
  water.computeRij();
  water.computeI();
  li.convBohr();
  li.computeNucRep();
  li.computeRij();
  li.computeI();
  o2.convBohr();
  o2.computeNucRep();
  o2.computeRij();
  o2.computeI();
  singo2.convBohr();
  singo2.computeNucRep();
  singo2.computeRij();
  singo2.computeI();

  for(auto JOB : jobs){
    JOB.runJob();
  };

  finalizeCQ();

  return 0;
}
