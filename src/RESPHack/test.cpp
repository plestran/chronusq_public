#include <response.h>
#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;

int main(int argc, char*argv[]){
  auto inputXYZName = argv[1];
  Molecule molecule;
  BasisSet basis;
  Controls controls;
  AOIntegrals aoints;
  MOIntegrals<double> moints;
  SingleSlater<double> singleSlater;
  Response<double> resp;
  FileIO fileio("test.inp",std::string(basename(inputXYZName)) + ".out");


  initCQ(argc,argv);
  controls.iniControls();
  fileio.iniH5Files();
  fileio.iniStdGroups();
  CQSetNumThreads(1);

  std::ifstream inputXYZ(inputXYZName);
//cout << inputXYZName << endl;
  std::string inLine;
  int nAtoms(0);
  getline(inputXYZ,inLine);
  nAtoms = std::atoi(inLine.c_str());


  getline(inputXYZ,inLine);

//cout << nAtoms << endl;

  molecule.setNAtoms(nAtoms);
  molecule.setCharge(3);
  molecule.setNTotalE(14);
  molecule.setMultip(1);
  molecule.alloc(fileio.out); 

  for(auto iAtm = 0; iAtm < nAtoms; iAtm++){
    getline(inputXYZ,inLine);

    std::istringstream iss(inLine);
    std::vector<std::string> tokens{
      std::istream_iterator<std::string>{iss},
      std::istream_iterator<std::string>{ }
    };
    molecule.setIndex(iAtm,HashAtom(tokens[0],0));
    std::vector<double> xyz;
    xyz.push_back(std::atof(tokens[1].c_str()));
    xyz.push_back(std::atof(tokens[2].c_str()));
    xyz.push_back(std::atof(tokens[3].c_str()));

    molecule.setCart(iAtm,xyz[0],xyz[1],xyz[2]);
  }

  molecule.convBohr();
  molecule.computeNucRep();
  molecule.computeRij();
  molecule.computeI();

  singleSlater.setRef(SingleSlater<double>::RHF);
  singleSlater.isClosedShell = true;

  basis.findBasisFile("cc-pVTZ");
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
  resp.setMeth(RESPONSE_TYPE::PPTDA);
  resp.doSA();
  resp.setNSek(8);
  resp.doFull();
  resp.doResponse();

//cout << resp.frequencies()[0]*phys.eVPerHartree << endl;

  VectorXd freq = resp.frequencies()[0]*phys.eVPerHartree;
  RealMatrix den = resp.transDen()[0];

  H5::H5File outfile(std::string(basename(inputXYZName)) + ".bin",H5F_ACC_TRUNC);

  hsize_t dimFreq[] = {8,1};
  H5::DataSpace space(2,dimFreq,NULL);

  H5::DataSet dataSet = outfile.createDataSet("Freq",H5::PredType::NATIVE_DOUBLE,space); 
  dataSet.write(freq.data(),H5::PredType::NATIVE_DOUBLE,space);

  hsize_t dimDen[] = {8,den.rows()};
  H5::DataSpace spaceDen(2,dimDen,NULL);
//cout << dimDen[0] << " " << dimDen[1] << endl;

  H5::DataSet dataSetDen = outfile.createDataSet("Den",H5::PredType::NATIVE_DOUBLE,spaceDen);
  dataSetDen.write(den.data(),H5::PredType::NATIVE_DOUBLE,spaceDen);

  hsize_t nbsq[] = {basis.nBasis(),basis.nBasis()};
  H5::DataSpace spaceNBSQ(2,nbsq,NULL);

  H5::DataSet dataMO = outfile.createDataSet("MO",H5::PredType::NATIVE_DOUBLE,spaceNBSQ);
  dataMO.write(singleSlater.moA()->data(),H5::PredType::NATIVE_DOUBLE,spaceNBSQ);

  hsize_t one = {1};
  H5::DataSpace spaceOne(1,&one,NULL);

  H5::DataSet dataEne = outfile.createDataSet("SCF",H5::PredType::NATIVE_DOUBLE,spaceOne);
  dataEne.write(&singleSlater.totalEnergy,H5::PredType::NATIVE_DOUBLE,spaceOne);

// cout << singleSlater.totalEnergy << endl;
  fileio.out << "SUCESS" << endl;
//prettyPrint(cout,den,"Den");

/*
  int nOA = 8;
  int nVA = basis.nBasis() - nOA;

  for(auto i = 0, ia = 0; i < nOA; i++)
  for(auto a = 0, A = nOA; a < nVA; a++, A++, ia++){
    cout << "(" << i+1 << "," << nOA + a+1 << ") " << den(ia,0) << endl;
  }
*/
  finalizeCQ(); 

  return 0;
};
