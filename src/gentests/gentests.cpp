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

  SCFSettings(bool def=true, std::array<double,3> f = {{0,0,0}}):
  useDefaultSettings(def), field(std::move(f)){ };
  

};

struct RTSettings{
  bool useDefaultSettings;
  int  maxStep;
  int  orthoType;
  int  irstrt;
  int  iEnvlp;
  double freq;
  double sigma;

  std::array<double,3> eField;

/*
  RTSettings(bool def = true, int mxStp = 20, 
    int orth = RealTime<double>::Lowdin, int iRstrt = 51):
    useDefaultSettings(def), maxStep(mxStp), orthoType(orth),
    irstrt(iRstrt),iEnvlp(RealTime<double>::Constant){ };
*/

  RTSettings(bool def = true, int mxStp = 20, 
    int orth = RealTime<double>::Lowdin, int iRstrt = 51,
    int iEnv = RealTime<double>::Constant, double freq_ = 0,
    std::array<double,3> fld = {{0,0,0}}, double sig = 0):
    useDefaultSettings(def), maxStep(mxStp), orthoType(orth),
    irstrt(iRstrt),iEnvlp(iEnv),freq(freq_),eField(fld),
    sigma(sig){ };
};

struct SCFResults{
  double Energy;
  std::array<double,3> elecDipole;
  std::array<double,6> elecQuadpole;
  std::array<double,10> elecOctpole;
};

struct RTResults{
  double LastEnergy;
  std::array<double,4> LastDipole;
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
  std::vector<std::string> func;

  std::string tag;

  SCFResults scfResults;
  RTResults  rtResults;

  CQJob(Molecule &mol_, std::string field_, std::string jobtype_,
    std::string reference_, int intAlg_,
    std::size_t nProcShared_, std::string basis_, std::string tag_) :
    mol(&mol_),
    field(std::move(field_)),
    jobType(std::move(jobtype_)),
    reference(std::move(reference_)),
    nProcShared(nProcShared_),
    basisSet(std::move(basis_)), 
    intAlg(intAlg_),
    tag(std::move(tag_)){ }

  CQJob(Molecule &mol_, std::string field_, std::string jobtype_,
    std::string reference_, int intAlg_,
    std::size_t nProcShared_, std::string basis_, std::vector<std::string> func_,
    std::string tag_) :
    CQJob(mol_,field_,jobtype_,reference_,intAlg_,nProcShared_,basis_,tag_){
      func = func_;
    } 
    

  CQJob(Molecule &mol_, std::string field_, std::string jobtype_,
    std::string reference_, int intAlg_,
    std::size_t nProcShared_, std::string basis_,std::string tag_, 
    SCFSettings scfset_) :
    CQJob(mol_,field_,jobtype_,reference_,intAlg_,nProcShared_,basis_,tag_) 
    { 
      scfSettings.useDefaultSettings = scfset_.useDefaultSettings;
      scfSettings.field = scfset_.field;
    };

  CQJob(Molecule &mol_, std::string field_, std::string jobtype_,
    std::string reference_, int intAlg_,
    std::size_t nProcShared_, std::string basis_,std::string tag_, 
    RTSettings rtset_) :
    CQJob(mol_,field_,jobtype_,reference_,intAlg_,nProcShared_,basis_,tag_) 
    { 
      rtSettings = rtset_;
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

    auto tmp = reference;
    if(func.size() != 0)
      tmp = tmp.replace(tmp.end()-2,tmp.end(),"KS");

    os << " " << tmp << endl;
  
    os << "job = " << jobType << endl;
    os << "basis = " << basisSet << endl;
    os << "ints = ";
    if(intAlg == AOIntegrals::INCORE)
      os << "incore";
    else if (intAlg == AOIntegrals::DIRECT)
      os << "direct";

    os << endl;
    if(func.size() != 0) {
      for(auto F : func) {
        if(!F.compare("SLATER"))      os << "exchange = SLATER";
        else if(!F.compare("B88"))    os << "exchange = B88";
        else if(!F.compare("LYP"))    os << "corr = LYP";
        else if(!F.compare("VWN5"))   os << "corr = VWN5";
        else if(!F.compare("VWN3"))   os << "corr = VWN3";
        os << endl;
      }
    }



    os << endl;
    os << endl;
    if(!jobType.compare("RT")){
      os << "[RT]" << endl;
      os << "maxstep = " << rtSettings.maxStep << endl;
      os << "ortho = ";
      if(rtSettings.orthoType == RealTime<double>::Cholesky)
        os << "Cholesky";
      else if(rtSettings.orthoType == RealTime<double>::Lowdin)
        os << "lowdin";
      os << endl;
      os << "IRSTRT = " << rtSettings.irstrt << endl;
      if(!rtSettings.useDefaultSettings) {
        os << "envelope = ";
        if(rtSettings.iEnvlp == RealTime<double>::Constant)
          os << "PW";
        else if(rtSettings.iEnvlp == RealTime<double>::LinRamp)
          os << "LinRamp";
        else if(rtSettings.iEnvlp == RealTime<double>::Gaussian)
          os << "Gaussian";
       

        os << endl;
        os << "edfield = ";
        os << std::setw(10) << std::setprecision(3);
        for(auto I : rtSettings.eField) 
          os << std::setw(10) << I;
        os << endl;

        os << "frequency = " << rtSettings.freq << endl;

        if(rtSettings.iEnvlp == RealTime<double>::Gaussian){
          os << "sigma = " << rtSettings.sigma;
        }
        
      }
    
    }
    os << endl;
    os << endl;
    os << "#"                     << endl;
    os << "#  Misc Specification" << endl;
    os << "#"                     << endl;
    os << "[Misc]" << endl;
    os << "nsmp = " << nProcShared << endl;


  };

  template<typename T> 
  void collectResults(SingleSlater<T> &ss){
    scfResults.Energy     = ss.totalEnergy;
    scfResults.elecDipole = ss.elecDipole();
    for(auto iXYZ = 0, run = 0; iXYZ < 3; iXYZ++)
    for(auto jXYZ = iXYZ; jXYZ < 3; jXYZ++, run++  ){
      scfResults.elecQuadpole[run] = ss.elecQuadpole()[iXYZ][jXYZ];
    }
    for(auto iXYZ = 0, run = 0; iXYZ < 3; iXYZ++)
    for(auto jXYZ = iXYZ; jXYZ < 3; jXYZ++) 
    for(auto kXYZ = jXYZ; kXYZ < 3; kXYZ++, run++  ){
      scfResults.elecOctpole[run] = ss.elecOctpole()[iXYZ][jXYZ][kXYZ];
    }
  };

  template<typename T> 
  void collectResults(RealTime<T> &rt){
    rtResults.LastEnergy     = rt.lastEnergy();
    rtResults.LastDipole     = rt.lastDipole();
  };

  void runJob() {
    CQSetNumThreads(nProcShared);
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

      if(func.size() != 0) {
        ss.isHF = false;
        ss.isDFT = true;
        for(auto F : func) {
          if(!F.compare("SLATER"))      ss.addSlater();
          else if(!F.compare("B88"))    ss.addB88();
          else if(!F.compare("LYP"))    ss.addLYP();
	  else if(!F.compare("VWN5"))    ss.addVWN5();
	  else if(!F.compare("VWN3"))    ss.addVWN3();
        }
      }

      if(mol->nAtoms() == 1) ss.setGuess(SingleSlater<double>::CORE);
      ss.communicate(*mol,basis,aoints,fileio,memManager);
      ss.initMeta();
      ss.genMethString();
      ss.alloc();
      ss.formGuess();
      ss.SCF3();
      ss.computeProperties();
      collectResults(ss);
      if(!jobType.compare("RT")){
        RealTime<double> rt;
        rt.communicate(fileio,aoints,ss);
        if(!rtSettings.useDefaultSettings) {
          rt.setMaxSteps(rtSettings.maxStep);
          rt.setOrthoTyp(rtSettings.orthoType);
          rt.setIRstrt(rtSettings.irstrt);
          rt.setEnvelope(rtSettings.iEnvlp);
          rt.setFreq(rtSettings.freq);
          rt.setFieldAmp(rtSettings.eField);
          rt.setSigma(rtSettings.sigma);
        }
        rt.initMeta();
        rt.alloc();
        rt.iniDensity();
        rt.doPropagation();
        collectResults(rt);
      }
    } else {
      SingleSlater<dcomplex> ss;
      if(!reference.compare("RHF")) {
        ss.setRef(SingleSlater<dcomplex>::RHF);
        ss.isClosedShell = true;
      } else if(!reference.compare("UHF")) {
        ss.setRef(SingleSlater<dcomplex>::UHF);
        ss.isClosedShell = false;
      }
      if(mol->nAtoms() == 1) ss.setGuess(SingleSlater<dcomplex>::CORE);
      ss.communicate(*mol,basis,aoints,fileio,memManager);
      ss.initMeta();
      ss.genMethString();
      ss.alloc();
      ss.formGuess();
      ss.SCF3();
      ss.computeProperties();
      collectResults(ss);
      if(!jobType.compare("RT")){
        RealTime<dcomplex> rt;
        if(!rtSettings.useDefaultSettings) {
          rt.setMaxSteps(rtSettings.maxStep);
          rt.setOrthoTyp(rtSettings.orthoType);
          rt.setIRstrt(rtSettings.irstrt);
          rt.setEnvelope(rtSettings.iEnvlp);
          rt.setFreq(rtSettings.freq);
          rt.setFieldAmp(rtSettings.eField);
          rt.setSigma(rtSettings.sigma);
        }
        rt.communicate(fileio,aoints,ss);
        rt.initMeta();
        rt.alloc();
        rt.iniDensity();
        rt.doPropagation();
        collectResults(rt);
      }
    }

  };


};

std::ostream& operator<< (std::ostream &os, const CQJob& job) {
  if(!job.field.compare("D"))
    os << "Real";
  else
    os << "Complex";

  auto tmp = job.reference;
  if(job.func.size() != 0)
    tmp = tmp.replace(tmp.end()-2,tmp.end(),"KS");

  os << " " << tmp << "-" << job.jobType;
  os << " - " << job.basisSet;
  os << " - " << job.tag;
  os << " ";
  if(job.nProcShared == 1) os << "Serial";
  else                 os << "OpenMP";

  os << " - ";
  if(job.intAlg == AOIntegrals::INCORE)
    os << "In-Core";
  else if(job.intAlg == AOIntegrals::DIRECT)
    os << "Direct";

  os << " Integrals";
  return os;
}

std::ostream& operator<< (std::ostream &os, const SCFResults &res){
  os << std::setprecision(12) << std::fixed;
  os << res.Energy << "/";
  for(auto I : res.elecDipole)
    os << I << "/";
  for(auto I : res.elecQuadpole)
    os << I << "/";
  for(auto I : res.elecOctpole)
    os << I << "/";
  return os;
}

std::ostream& operator<< (std::ostream &os, const RTResults &res){
  os << std::setprecision(12) << std::fixed;
  os << res.LastEnergy << "/";
  for(auto I : res.LastDipole)
    os << I << "/";
  return os;
}

// MOLECULE
  Molecule water,li,o2,singo2;

void RSmRHF_SCF(std::vector<CQJob> &jobs) {
  // Small Molecule Real RHF SCF Serial INCORE
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"6-31G",
    "Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small");

  // Small Molecule Real RHF SCF Serial DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small");

  // Small Molecule Real RHF SMP DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small");
}

void RSmUHF_SCF(std::vector<CQJob> &jobs) {
  // Small Molecule Real UHF SCF Serial INCORE
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::INCORE,1,"6-31G",
    "Small");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::INCORE,1,"6-31G",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small");

  // Small Molecule Real UHF SCF Serial DIRECT
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small");

  // Small Molecule Real UHF SCF SMP DIRECT
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small");
  jobs.emplace_back(li,"D","SCF","UHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small");
  jobs.emplace_back(o2,"D","SCF","UHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small");
}

void CSmRHF_SCF(std::vector<CQJob> &jobs) {
  // Small Molecule Complex RHF SCF Serial INCORE
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::INCORE,1,"6-31G",
    "Small");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small");

  // Small Molecule Complex RHF SCF Serial DIRECT
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small");

  // Small Molecule Complex RHF SCF SMP DIRECT
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small");
  jobs.emplace_back(singo2,"C","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small");
}

void RSmRSLATER_SCF(std::vector<CQJob> &jobs) {
  // Small Molecule Real R-SLATER SCF Serial INCORE
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G",
    std::vector<std::string>({"SLATER"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"6-31G",
    std::vector<std::string>({"SLATER"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    std::vector<std::string>({"SLATER"}),"Small");

  // Small Molecule Real R-SLATER SCF Serial DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    std::vector<std::string>({"SLATER"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G",
    std::vector<std::string>({"SLATER"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    std::vector<std::string>({"SLATER"}),"Small");

  // Small Molecule Real R-SLATER SCF SMP DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    std::vector<std::string>({"SLATER"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G",
    std::vector<std::string>({"SLATER"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    std::vector<std::string>({"SLATER"}),"Small");
}

void RSmRLSDA_SCF(std::vector<CQJob> &jobs) {
  // Small Molecule Real R-LSDA SCF Serial INCORE
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"6-31G",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");

  // Small Molecule Real R-LSDA SCF Serial DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");

  // Small Molecule Real R-LSDA SCF SMP DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    std::vector<std::string>({"SLATER","VWN5"}),"Small");
}

void RSmRSVWN3_SCF(std::vector<CQJob> &jobs) {
  // Small Molecule Real R-SVWN3 SCF Serial INCORE
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"6-31G",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");

  // Small Molecule Real R-SVWN3 SCF Serial DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");

  // Small Molecule Real R-SVWN3 SCF SMP DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    std::vector<std::string>({"SLATER","VWN3"}),"Small");
}

void RSmRB88_SCF(std::vector<CQJob> &jobs) {
  // Small Molecule Real R-B88 SCF Serial INCORE
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G",
    std::vector<std::string>({"B88"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"6-31G",
    std::vector<std::string>({"B88"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    std::vector<std::string>({"B88"}),"Small");

  // Small Molecule Real R-B88 SCF Serial DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    std::vector<std::string>({"B88"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G",
    std::vector<std::string>({"B88"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    std::vector<std::string>({"B88"}),"Small");

  // Small Molecule Real R-B88 SCF SMP DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    std::vector<std::string>({"B88"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G",
    std::vector<std::string>({"B88"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    std::vector<std::string>({"B88"}),"Small");
}

void RSmRBLYP_SCF(std::vector<CQJob> &jobs) {
  // Small Molecule Real R-BLYP SCF Serial INCORE
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"STO-3G",
    std::vector<std::string>({"B88","LYP"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"6-31G",
    std::vector<std::string>({"B88","LYP"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    std::vector<std::string>({"B88","LYP"}),"Small");

  // Small Molecule Real R-BLYP SCF Serial DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    std::vector<std::string>({"B88","LYP"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"6-31G",
    std::vector<std::string>({"B88","LYP"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    std::vector<std::string>({"B88","LYP"}),"Small");

  // Small Molecule Real R-BLYP SCF SMP DIRECT
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    std::vector<std::string>({"B88","LYP"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"6-31G",
    std::vector<std::string>({"B88","LYP"}),"Small");
  jobs.emplace_back(water,"D","SCF","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    std::vector<std::string>({"B88","LYP"}),"Small");
}

void RSmRHF_RT_NOFIELD(std::vector<CQJob> &jobs) {
  // Small Molecule Real RHF RT (LOWDIN) Serial INCORE
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small");
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"6-31G",
    "Small");
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small");

  // Small Molecule Real RHF RT (LOWDIN) Serial DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small");
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small");
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small");

  // Small Molecule Real RHF RT (LOWDIN) SMP DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small");
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small");
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small");

  // Small Molecule Real RHF RT (CHOLESKY) Serial INCORE
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));

  // Small Molecule Real RHF RT (CHOLESKY) Serial DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));

  // Small Molecule Real RHF RT (CHOLESKY) SMP DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky));
}

void RSmRHF_RT_PW(std::vector<CQJob> &jobs) {
  // Small Molecule Real RHF RT (LOWDIN) Serial INCORE
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (LOWDIN) Serial DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (LOWDIN) SMP DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (CHOLESKY) Serial INCORE
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (CHOLESKY) Serial DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (CHOLESKY) SMP DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::Constant,9.6,
    std::array<double,3>({0.01,0,0})));
}

void RSmRHF_RT_LIN(std::vector<CQJob> &jobs) {
  // Small Molecule Real RHF RT (LOWDIN) Serial INCORE
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (LOWDIN) Serial DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (LOWDIN) SMP DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Lowdin,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (CHOLESKY) Serial INCORE
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::INCORE,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (CHOLESKY) Serial DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,1,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));

  // Small Molecule Real RHF RT (CHOLESKY) SMP DIRECT
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"STO-3G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"6-31G",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
  jobs.emplace_back(water,"D","RT","RHF",AOIntegrals::DIRECT,2,"cc-pVDZ",
    "Small",RTSettings(false,20,RealTime<double>::Cholesky,51,
     RealTime<double>::LinRamp,9.6,
    std::array<double,3>({0.01,0,0})));
}


int main() {
  std::vector<CQJob> jobs;

  initCQ(0,NULL);

  loadPresets<WATER>(water);
  loadPresets<Li>(li);
  loadPresets<O2>(o2);
  loadPresets<SingO2>(singo2);


  RSmRHF_SCF(jobs);
  CSmRHF_SCF(jobs);
  RSmUHF_SCF(jobs);
  RSmRSLATER_SCF(jobs);
  RSmRLSDA_SCF(jobs);
  RSmRSVWN3_SCF(jobs);
  RSmRB88_SCF(jobs);
  RSmRBLYP_SCF(jobs);
  RSmRHF_RT_NOFIELD(jobs);
  RSmRHF_RT_PW(jobs);
  RSmRHF_RT_LIN(jobs);


  std::vector<std::string> fnames;
  std::ofstream index("test.index");
  for(auto iJob = 0; iJob < jobs.size(); iJob++){
    std::stringstream fname;
    cout << "Generating Job ";
    fname << "test" << std::setw(4) << std::setfill('0') << iJob+1;
    cout << fname.str() << " ";
    cout << "\"";
    cout << jobs[iJob];
    cout << "\"" << endl;

    fnames.emplace_back(fname.str() + ".inp");
    std::ofstream file(fnames.back());
    jobs[iJob].writeInput(file);

    index << fnames.back() << " - " << jobs[iJob] << endl;
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

  std::ofstream refs("chronus-ref.val");
  for(auto iJob = 0; iJob < jobs.size(); iJob++){
    cout << "Running Job ";
    cout << fnames[iJob] << endl;

    jobs[iJob].runJob();
//  cout << JOB.scfResults << endl;
    if(!jobs[iJob].jobType.compare("SCF"))
      refs << fnames[iJob] << "/" << jobs[iJob].scfResults  
           << jobs[iJob].jobType << endl;
    else if(!jobs[iJob].jobType.compare("RT"))
      refs << fnames[iJob] << "/" << jobs[iJob].rtResults  
           << jobs[iJob].jobType << endl;
  };

  finalizeCQ();

  return 0;
}
