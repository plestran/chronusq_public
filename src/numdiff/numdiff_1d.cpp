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
std::string genFName(double,int,std::string);
Eigen::VectorXd checkPhase(double*,double*,int,int);
RealMatrix genSpx(BasisSet&,BasisSet&);

Eigen::VectorXd GS_ES_NACME_CIS(int,bool,int,int,double,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&, 
  RealMatrix&,RealMatrix&);
Eigen::VectorXd ES_GS_NACME_CIS(int,bool,int,int,double,RealMatrix&,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&);
RealMatrix ES_ES_NACME_CIS(int,bool,int,int,double,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&, 
  RealMatrix&,RealMatrix&);
RealMatrix ES_ES_NACME_PPTDA(int,bool,int,int,double,
  RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&,RealMatrix&, 
  RealMatrix&,RealMatrix&);
Derivatives computeNAC2pt1D(Molecule&,Molecule&,Molecule&,FileIO&,FileIO&, 
  FileIO&,Controls&,std::string&,RESPONSE_TYPE,int,double);




void oneDScan(std::vector<double>& scanX, 
              std::string &basisName,
              RESPONSE_TYPE respType,
              int charge, int nFreq,
              int px,
              std::string &PREFIX){

//std::string PREFIX = "LiH";
  std::string BINPREFIX = "BIN";
  std::string XYZPREFIX = "XYZ";

  initCQ();
  Controls controls;
  CQSetNumThreads(3);
  controls.iniControls();

//std::string basisName = "cc-pVDZ";
//RESPONSE_TYPE respType = RESPONSE_TYPE::CIS;
//int charge = 0;
//int nFreq = 4;
//int px = 3;
//int py = 6;
  bool debug = true;

  std::string GS_ENERGY_FILE_NAME = PREFIX + "_gsEnergy.dat";
  std::string ES_ENERGY_FILE_NAME = PREFIX + "_esEnergy.dat";
  std::string GS_GRADIENT_DX_FILE_NAME = PREFIX + "_gsGrad_IX_1.dat";
  std::string ES_GRADIENT_DX_FILE_NAME = PREFIX + "_esGrad_IX_1.dat";
  std::string ES_GS_NACME_DX_FILE_NAME = PREFIX + "_esgsNAC_IX_1.dat";
  std::string ES_ES_NACME_DX_FILE_NAME = PREFIX + "_esesNAC_IX_1.dat";

  std::ofstream GS_ENERGY_FILE      (GS_ENERGY_FILE_NAME      ); 
  std::ofstream ES_ENERGY_FILE      (ES_ENERGY_FILE_NAME      );
  std::ofstream GS_GRADIENT_DX_FILE (GS_GRADIENT_DX_FILE_NAME );
  std::ofstream ES_GRADIENT_DX_FILE (ES_GRADIENT_DX_FILE_NAME );
  std::ofstream ES_GS_NACME_DX_FILE (ES_GS_NACME_DX_FILE_NAME );
  std::ofstream ES_ES_NACME_DX_FILE (ES_ES_NACME_DX_FILE_NAME );

  GS_ENERGY_FILE << std::setprecision(8) << std::fixed;
  ES_ENERGY_FILE << std::setprecision(8) << std::fixed;
  GS_GRADIENT_DX_FILE  << std::setprecision(8)    << std::fixed;
  ES_GRADIENT_DX_FILE  << std::setprecision(8)    << std::fixed;
  ES_GS_NACME_DX_FILE  << std::setprecision(8)    << std::fixed;
  ES_ES_NACME_DX_FILE  << std::setprecision(8)    << std::fixed;


  for(auto IX = 1; IX < (scanX.size()-1); IX++){
    // Output header
    cout << bannerTop << endl;
    cout << "Starting Numerical Differentiation for:" << endl;
    cout << "  IX = " << IX << endl;
    FileIO fileio_0(PREFIX+"_"+std::to_string(IX)+"_0.inp",PREFIX+"_"+std::to_string(IX)+".out_0",
      PREFIX+"_"+std::to_string(IX)+".rst_0");
    FileIO fileio_p(PREFIX+"_"+std::to_string(IX)+"_p.inp",PREFIX+"_"+std::to_string(IX)+".out_p",
      PREFIX+"_"+std::to_string(IX)+".rst_p");
    FileIO fileio_m(PREFIX+"_"+std::to_string(IX)+"_m.inp",PREFIX+"_"+std::to_string(IX)+".out_m",
      PREFIX+"_"+std::to_string(IX)+".rst_m");

    fileio_0.iniH5Files();
    fileio_p.iniH5Files();
    fileio_m.iniH5Files();

    fileio_0.iniStdGroups();
    fileio_p.iniStdGroups();
    fileio_m.iniStdGroups();

    // Get coordinate information
    double X = scanX[IX];


    double XP = scanX[IX+1];
    double XM = scanX[IX-1];


    // Output coordinate information
    cout << endl;
    cout << "  X      = " << X << endl;
    cout << "  X + DX = " << XP << endl;
    cout << "  X - DX = " << XM << endl;


    // Generate file names baed on convention
    std::string fBaseName0 = genFName(X ,px,PREFIX);
    std::string fBaseNamep = genFName(XP,px,PREFIX);
    std::string fBaseNamem = genFName(XM,px,PREFIX);

    std::string fNameGeom_0 = XYZPREFIX+"/"+fBaseName0 + ".xyz";
    std::string fNameGeom_p = XYZPREFIX+"/"+fBaseNamep + ".xyz";
    std::string fNameGeom_m = XYZPREFIX+"/"+fBaseNamem + ".xyz";

    // Output file names
    cout << endl;
    cout << "  Base File Name X:    " << fBaseName0 << endl;
    cout << "  Base File Name X+DX: " << fBaseNamep << endl;
    cout << "  Base File Name X-DX: " << fBaseNamem << endl;

    cout << endl;
    cout << "  XYZ File Name X:    " << fNameGeom_0 << endl;
    cout << "  XYZ File Name X+DX: " << fNameGeom_p << endl;
    cout << "  XYZ File Name X-DX: " << fNameGeom_m << endl;

    

    // Create CQ Objects
    Molecule geom_0;
    Molecule geom_p;
    Molecule geom_m;
/*
    BasisSet basis0;
    BasisSet basisp;
    BasisSet basism;

    AOIntegrals aoints_0;
    AOIntegrals aoints_p;
    AOIntegrals aoints_m;

    SingleSlater<double> ss_0;
    SingleSlater<double> ss_p;
    SingleSlater<double> ss_m;

    MOIntegrals<double> moints_0;
    MOIntegrals<double> moints_p;
    MOIntegrals<double> moints_m;

    Response<double> resp_0;
    Response<double> resp_p;
    Response<double> resp_m;
*/


    // Set up Geometries
    geom_0.setMultip(1);
    geom_p.setMultip(1);
    geom_m.setMultip(1);

    std::ifstream XYZ0File(fNameGeom_0);
    std::ifstream XYZpFile(fNameGeom_p);
    std::ifstream XYZmFile(fNameGeom_m);

    if(!XYZ0File.good()) CErr("Cannot find XYZFile0");
    if(!XYZpFile.good()) CErr("Cannot find XYZFilep");
    if(!XYZmFile.good()) CErr("Cannot find XYZFilem");

    // Get number of atoms
    int nAtoms(0);
    std::string inLine_0;
    std::string inLine_p;
    std::string inLine_m;

    getline(XYZ0File,inLine_0);
    getline(XYZpFile,inLine_p);
    getline(XYZmFile,inLine_m);

    nAtoms = std::atoi(inLine_0.c_str());
    geom_0.setNAtoms(nAtoms);
    geom_p.setNAtoms(nAtoms);
    geom_m.setNAtoms(nAtoms);

    // Allocate Molecule objects
    geom_0.alloc(fileio_0.out);
    geom_p.alloc(fileio_p.out);
    geom_m.alloc(fileio_m.out);

    // Skip line
    getline(XYZ0File,inLine_0);
    getline(XYZpFile,inLine_p);
    getline(XYZmFile,inLine_m);

    int nElec = -charge;
    for(auto iAtm = 0; iAtm < nAtoms; iAtm++){
      getline(XYZ0File,inLine_0);
      getline(XYZpFile,inLine_p);
      getline(XYZmFile,inLine_m);

      std::istringstream iss_0(inLine_0);
      std::istringstream iss_p(inLine_p);
      std::istringstream iss_m(inLine_m);

      std::vector<std::string> tokens_0{
        std::istream_iterator<std::string>{iss_0},
        std::istream_iterator<std::string>{ }
      };
      std::vector<std::string> tokens_p{
        std::istream_iterator<std::string>{iss_p},
        std::istream_iterator<std::string>{ }
      };
      std::vector<std::string> tokens_m{
        std::istream_iterator<std::string>{iss_m},
        std::istream_iterator<std::string>{ }
      };

      auto INDX = HashAtom(tokens_0[0],0);
      nElec += atom[INDX].atomicNumber;

      geom_0.setIndex(iAtm,INDX);
      geom_p.setIndex(iAtm,INDX);
      geom_m.setIndex(iAtm,INDX);


      geom_0.setCart(iAtm,
        std::atof(tokens_0[1].c_str()),
        std::atof(tokens_0[2].c_str()),
        std::atof(tokens_0[3].c_str())
      );
      geom_p.setCart(iAtm,
        std::atof(tokens_p[1].c_str()),
        std::atof(tokens_p[2].c_str()),
        std::atof(tokens_p[3].c_str())
      );
      geom_m.setCart(iAtm,
        std::atof(tokens_m[1].c_str()),
        std::atof(tokens_m[2].c_str()),
        std::atof(tokens_m[3].c_str())
      );
    }

    geom_0.setCharge(charge);
    geom_p.setCharge(charge);
    geom_m.setCharge(charge);

    geom_0.setNTotalE(nElec);
    geom_p.setNTotalE(nElec);
    geom_m.setNTotalE(nElec);

    geom_0.convBohr();
    geom_p.convBohr();
    geom_m.convBohr();

    geom_0.computeNucRep();
    geom_p.computeNucRep();
    geom_m.computeNucRep();

    geom_0.computeRij();
    geom_p.computeRij();
    geom_m.computeRij();

    geom_0.computeI();
    geom_p.computeI();
    geom_m.computeI();

    if(debug){
      cout << "GEOMETRY X" << endl;
      geom_0.printInfo();
      cout << "GEOMETRY X+DX" << endl;
      geom_p.printInfo();
      cout << "GEOMETRY X-DX" << endl;
      geom_m.printInfo();
    }

    Derivatives derv = computeNAC2pt1D(geom_0,geom_p,geom_m,fileio_0,fileio_p,
      fileio_m,controls,basisName,respType,nFreq,XP-X);


/*
    double          scf_0  = derv.GS_ENERGY;
    Eigen::VectorXd freq_0 = derv.ES_ENERGY;
    double          gsdx   = derv.GS_GRAD;
    GS_ENERGY_FILE << X << "," << scf_0 << endl;

    ES_ENERGY_FILE << X << ",";
    for(auto iSt = 0; iSt < nFreq; iSt++){
      ES_ENERGY_FILE << freq_0(iSt);
      if(iSt != (nFreq-1)) ES_ENERGY_FILE << ",";
    }
    ES_ENERGY_FILE << endl;

    GS_GRADIENT_DX_FILE << X << "," << gsdx << endl;

    ES_GRADIENT_DX_FILE << X << ",";
    for(auto iSt = 0; iSt < nFreq; iSt++){
      ES_GRADIENT_DX_FILE << freqDX(iSt);
      if(iSt != (nFreq-1)) 
        ES_GRADIENT_DX_FILE << ",";
    }
    ES_GRADIENT_DX_FILE << endl;

    ES_GS_NACME_DX_FILE << X << ",";
    for(auto iSt = 0; iSt < nFreq; iSt++){
      ES_GS_NACME_DX_FILE << NAC_ES_GS(iSt);
      if(iSt != (nFreq-1)) 
        ES_GS_NACME_DX_FILE << ",";
    }
    ES_GS_NACME_DX_FILE << endl;

    ES_ES_NACME_DX_FILE << X << ",";
    for(auto iSt = 0; iSt < nFreq; iSt++)
    for(auto jSt = 0; jSt < nFreq; jSt++){
      ES_ES_NACME_DX_FILE << NAC_ES_ES(iSt,jSt);
      if(iSt != (nFreq-1) || jSt != (nFreq-1)) 
        ES_ES_NACME_DX_FILE << ",";
    }
    ES_ES_NACME_DX_FILE << endl;
*/
  }
  finalizeCQ();
}



// Naming convention PREFIX_xptx_ypty.xyz
std::string genFName(double x, int px, std::string PREFIX) {
    std::stringstream ss1;

    ss1 << std::fixed << std::setprecision(px) << x;

    std::string str1 = ss1.str();

    auto ptPos = str1.find(".");
    str1.replace(str1.begin()+ptPos,str1.begin()+ptPos+1,"pt");

    std::string fName = PREFIX + "_" + str1;    

    return fName;
};


Eigen::VectorXd ES_GS_NACME_CIS(int nFreq,bool renorm, int nocc, int nvir, 
  double step, RealMatrix &T_0, RealMatrix &MO_0, RealMatrix &MO_p, 
  RealMatrix &MO_m, RealMatrix &S_0_p, RealMatrix &S_0_m){

  Eigen::VectorXd NACME(nFreq);

  NACME.setZero();

  if(renorm) T_0 *= std::sqrt(0.5);
//T_0.transposeInPlace();
//prettyPrint(cout,T_0,"T");

  int nThreads = omp_get_max_threads();
//std::vector<Eigen::VectorXd> TMPNACME(nThreads,Eigen::VectorXd(nFreq));
//for(auto iTh = 0; iTh < nThreads; iTh++) TMPNACME[iTh].setZero();


//RealMatrix SWAPPED_0(MO_0);
//RealMatrix Prod_0_p(MO_0);
//RealMatrix Prod_0_m(MO_0);
//SWAPPED_0.setZero();
//Prod_0_p.setZero();
//Prod_0_m.setZero();

  std::vector<RealMatrix> SWAPPED_0(nThreads,RealMatrix(MO_0));
  std::vector<RealMatrix> Prod_0_p(nThreads,RealMatrix(MO_0));
  std::vector<RealMatrix> Prod_0_m(nThreads,RealMatrix(MO_0));

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < nFreq; iSt++){
    cout << "Working on iSt = " << iSt << " ES->GS" << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto i = 0, ia = 0; i < nocc; i++)
      for(auto a = 0; a < nvir; a++, ia++){
        if(ia % nThreads != thread_id) continue;
        if(std::abs(T_0(ia,iSt)) < 1e-8) continue;
        SWAPPED_0[thread_id].setZero();
        Prod_0_p[thread_id].setZero();
        Prod_0_m[thread_id].setZero();
        
        SWAPPED_0[thread_id] = MO_0;
        SWAPPED_0[thread_id].col(i).swap(SWAPPED_0[thread_id].col(nocc+a)); 
        
        Prod_0_p[thread_id] = SWAPPED_0[thread_id].transpose() * S_0_p * MO_p;
        Prod_0_m[thread_id] = SWAPPED_0[thread_id].transpose() * S_0_m * MO_m;
     
        double OvLp_IASwap_0_p =
          Prod_0_p[thread_id].block(0,0,nocc,nocc).determinant();
        double OvLp_IASwap_0_m =
          Prod_0_m[thread_id].block(0,0,nocc,nocc).determinant();
     
      //TMPNACME[thread_id](iSt) += T_0(ia,iSt) *
      //  (OvLp_IASwap_0_p - OvLp_IASwap_0_m) /
      //  (2*step);
        #pragma omp critical
        {
          NACME(iSt) += T_0(ia,iSt) *
            (OvLp_IASwap_0_p - OvLp_IASwap_0_m) /
            (2*step);
        }
      } // loop ove ia
    }
  } // loop over states

  
//for(auto iTh = 0; iTh < nThreads; iTh++) NACME += TMPNACME[iTh];
  NACME = 2*NACME;
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "GS -> ES NACME took " << elapsed.count() << " s" << endl;
  return NACME;
};

Eigen::VectorXd GS_ES_NACME_CIS(int nFreq,bool renorm,
  int nocc, int nvir, double step, RealMatrix &T_0, RealMatrix &T_p, 
  RealMatrix &T_m, RealMatrix &MO_0, RealMatrix &MO_p, RealMatrix &MO_m, 
  RealMatrix &S_0_p, RealMatrix &S_0_m){

  Eigen::VectorXd NACME(nFreq);

  NACME.setZero();

  if(renorm) {
    T_p *= std::sqrt(0.5);
    T_m *= std::sqrt(0.5);
  }

  RealMatrix SWAPPED_p(MO_0);
  RealMatrix SWAPPED_m(MO_0);

  RealMatrix Prod_0_p(MO_0);
  RealMatrix Prod_0_m(MO_0);

  SWAPPED_p.setZero();
  SWAPPED_m.setZero();

  Prod_0_p.setZero();
  Prod_0_m.setZero();

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < nFreq; iSt++){
    cout << "Working on iSt = " << iSt << " GS->ES" << endl;
    for(auto i = 0, ia = 0; i < nocc; i++)
    for(auto a = 0; a < nvir; a++, ia++){
      if(
        (std::abs(T_p(ia,iSt)) < 1e-8) &&
        (std::abs(T_m(ia,iSt)) < 1e-8) 
      ) continue;

      SWAPPED_p.setZero();
      SWAPPED_m.setZero();

      Prod_0_p.setZero();
      Prod_0_m.setZero();
      
      SWAPPED_p = MO_p;
      SWAPPED_m = MO_m;

      SWAPPED_p.col(i).swap(SWAPPED_p.col(nocc+a)); 
      SWAPPED_m.col(i).swap(SWAPPED_m.col(nocc+a)); 
      
      Prod_0_p = MO_0.transpose() * S_0_p * SWAPPED_p;
      Prod_0_m = MO_0.transpose() * S_0_m * SWAPPED_m;

      
      double OvLp_IASwap_0_p =
        Prod_0_p.block(0,0,nocc,nocc).determinant();
      double OvLp_IASwap_0_m =
        Prod_0_m.block(0,0,nocc,nocc).determinant();

      NACME(iSt) += (
        T_p(ia,iSt)*OvLp_IASwap_0_p - 
        T_m(ia,iSt)*OvLp_IASwap_0_m
      ) / (2*step);
    } // loop ove ia
  } // loop over states
  NACME *= 2;
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "GS -> ES NACME took " << elapsed.count() << " s" << endl;

  return NACME;
};

RealMatrix ES_ES_NACME_CIS(int nFreq,bool renorm, int nocc, int nvir, 
  double step, RealMatrix &T_0, RealMatrix &T_p, RealMatrix &T_m, 
  RealMatrix &MO_0, RealMatrix &MO_p, RealMatrix &MO_m, 
  RealMatrix &S_0_p, RealMatrix &S_0_m){

  bool doValidation = false;

  RealMatrix NACME(nFreq,nFreq);

  NACME.setZero();

  if(renorm) {
    T_p *= std::sqrt(0.5);
    T_m *= std::sqrt(0.5);
  }

  int nThreads = omp_get_max_threads();


  std::vector<RealMatrix> SWAPPED_IA_0(nThreads,RealMatrix(MO_0));
  std::vector<RealMatrix> SWAPPED_JB_p(nThreads,RealMatrix(MO_0));
  std::vector<RealMatrix> SWAPPED_JB_m(nThreads,RealMatrix(MO_0));

  std::vector<RealMatrix> Prod_0_p(nThreads,RealMatrix(MO_0));
  std::vector<RealMatrix> Prod_0_m(nThreads,RealMatrix(MO_0));

// VALIDATION BEGIN **************
  RealMatrix T1Prod_0_p(MO_0);
  RealMatrix T1Prod_0_m(MO_0);
  RealMatrix T2Prod_0_p(MO_0);
  RealMatrix T2Prod_0_m(MO_0);
// VALIDATION END **************

//SWAPPED_IA_0.setZero();
//SWAPPED_JB_p.setZero();
//SWAPPED_JB_m.setZero();

//Prod_0_p.setZero();
//Prod_0_m.setZero();


  if(doValidation) {
// VALIDATION BEGIN **************
    T1Prod_0_p.setZero();
    T1Prod_0_m.setZero();
    T2Prod_0_p.setZero();
    T2Prod_0_m.setZero();
// VALIDATION END **************
  }

  auto NACMEStart = std::chrono::high_resolution_clock::now();
  for(auto iSt = 0; iSt < nFreq; iSt++)
  for(auto jSt = 0; jSt < nFreq; jSt++){
    cout << "Working on iSt = " << iSt <<" jSt = " << jSt << endl;
    #pragma omp parallel
    {
      int thread_id = omp_get_thread_num();
      for(auto i = 0, iajb = 0, ia = 0; i < nocc; i++)
      for(auto a = 0; a < nvir; a++, ia++){
        //if(T_0(ia,iSt) < 1e-8) continue;   
     
        SWAPPED_IA_0[thread_id].setZero();
        SWAPPED_IA_0[thread_id] = MO_0;
     
        SWAPPED_IA_0[thread_id].col(i).swap(SWAPPED_IA_0[thread_id].col(nocc+a)); 
     
        for(auto j = 0, jb = 0; j < nocc; j++)
        for(auto b = 0; b < nvir; b++, jb++, iajb++){
          if(iajb % nThreads != thread_id) continue;
          if(
            (std::abs(T_p(jb,jSt)) < 1e-8) &&
            (std::abs(T_m(jb,jSt)) < 1e-8) 
          ) continue;
       
          SWAPPED_JB_p[thread_id].setZero();
          SWAPPED_JB_m[thread_id].setZero();
       
          Prod_0_p[thread_id].setZero();
          Prod_0_m[thread_id].setZero();
          
          SWAPPED_JB_p[thread_id] = MO_p;
          SWAPPED_JB_m[thread_id] = MO_m;
       
          SWAPPED_JB_p[thread_id].col(j).swap(SWAPPED_JB_p[thread_id].col(nocc+b)); 
          SWAPPED_JB_m[thread_id].col(j).swap(SWAPPED_JB_m[thread_id].col(nocc+b)); 
          
          Prod_0_p[thread_id] = 
            SWAPPED_IA_0[thread_id].transpose() * S_0_p * SWAPPED_JB_p[thread_id];
          Prod_0_m[thread_id] = 
            SWAPPED_IA_0[thread_id].transpose() * S_0_m * SWAPPED_JB_m[thread_id];
       
         if(doValidation) {
//  VVALIDATION BEGIN **************
            T1Prod_0_p.setZero();
            T1Prod_0_m.setZero();
           
            T1Prod_0_p = 
              MO_0.transpose() * S_0_p * SWAPPED_JB_p[thread_id];
            T1Prod_0_m = 
              MO_0.transpose() * S_0_m * SWAPPED_JB_m[thread_id];
           
            T2Prod_0_p.setZero();
            T2Prod_0_m.setZero();
           
            T2Prod_0_p = 
              SWAPPED_IA_0[thread_id].transpose() * S_0_p * MO_p;
            T2Prod_0_m = 
              SWAPPED_IA_0[thread_id].transpose() * S_0_m * MO_m;
          }
//  VVALIDATION END **************
          
          double OvLp_IASwap_0_p =
            Prod_0_p[thread_id].block(0,0,nocc,nocc).determinant();
          double OvLp_IASwap_0_m =
            Prod_0_m[thread_id].block(0,0,nocc,nocc).determinant();
     
//  VVALIDATION BEGIN **************
          if(doValidation) {
            double T1OvLp_IASwap_0_p =
              T1Prod_0_p.block(0,0,nocc,nocc).determinant();
            double T1OvLp_IASwap_0_m =
              T1Prod_0_m.block(0,0,nocc,nocc).determinant();
           
            double T2OvLp_IASwap_0_p =
              T2Prod_0_p.block(0,0,nocc,nocc).determinant();
            double T2OvLp_IASwap_0_m =
              T2Prod_0_m.block(0,0,nocc,nocc).determinant();
           
            cout << "******** " << iSt << "," << jSt << "," << ia << "," << jb << endl;
            cout << OvLp_IASwap_0_p  << endl;
            cout << OvLp_IASwap_0_m  << endl;
            cout << endl;
            cout << T1OvLp_IASwap_0_p  << endl;
            cout << T1OvLp_IASwap_0_m  << endl;
            cout << endl;
            cout << T2OvLp_IASwap_0_p  << endl;
            cout << T2OvLp_IASwap_0_m  << endl;
            cout << endl;
            cout << OvLp_IASwap_0_p - T1OvLp_IASwap_0_p*T2OvLp_IASwap_0_p << endl;
            cout << OvLp_IASwap_0_m - T1OvLp_IASwap_0_m*T2OvLp_IASwap_0_m << endl;
          }
//  VVALIDATION END **************
       
          #pragma omp critical
          {
            NACME(iSt,jSt) += T_0(ia,iSt)*(
              T_p(jb,jSt)*OvLp_IASwap_0_p - 
              T_m(jb,jSt)*OvLp_IASwap_0_m
            ) / (2*step);
          }
        } // loop ove jb
      } // loop ove ia
    }
  } // loop over states
  NACME *= 2;
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "GS -> ES NACME took " << elapsed.count() << " s" << endl;

  return NACME;
};

RealMatrix ES_ES_NACME_PPTDA(int nFreq,bool renorm, int nocc, int nvir, 
  double step, RealMatrix &T_0, RealMatrix &T_p, RealMatrix &T_m, 
  RealMatrix &MO_0, RealMatrix &MO_p, RealMatrix &MO_m, 
  RealMatrix &S_0_p, RealMatrix &S_0_m){

  RealMatrix NACME(nFreq,nFreq);

  NACME.setZero();

  int nThreads = omp_get_max_threads();

  std::vector<RealMatrix> SWAPPED_AB_0(nThreads,RealMatrix(MO_0));
  std::vector<RealMatrix> SWAPPED_CD_p(nThreads,RealMatrix(MO_0));
  std::vector<RealMatrix> SWAPPED_CD_m(nThreads,RealMatrix(MO_0));

  std::vector<RealMatrix> Prod_0_p(nThreads,RealMatrix(MO_0));
  std::vector<RealMatrix> Prod_0_m(nThreads,RealMatrix(MO_0));

  
  cout << nvir << endl;
  cout << (nvir*(nvir+1)/2)*(nvir*(nvir+1)/2) << endl;
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
//      printf("abcd:%d\n",abcd);
        if(
          (std::abs(T_p(cd,jSt)) < 1e-8) &&
          (std::abs(T_m(cd,jSt)) < 1e-8)
        ) continue;
     
     
        // AC OVerlap
        SWAPPED_AB_0[thread_id].setZero();
        SWAPPED_CD_p[thread_id].setZero();
        SWAPPED_CD_m[thread_id].setZero();
     
        Prod_0_p[thread_id].setZero();
        Prod_0_m[thread_id].setZero();
     
        SWAPPED_AB_0[thread_id] = MO_0;
        SWAPPED_CD_p[thread_id] = MO_p;
        SWAPPED_CD_m[thread_id] = MO_m;
     
        SWAPPED_AB_0[thread_id].col(nocc).swap(SWAPPED_AB_0[thread_id].col(nocc+a));
        SWAPPED_CD_p[thread_id].col(nocc).swap(SWAPPED_CD_p[thread_id].col(nocc+c));
        SWAPPED_CD_m[thread_id].col(nocc).swap(SWAPPED_CD_m[thread_id].col(nocc+c));
     
     
        Prod_0_p[thread_id] = 
          SWAPPED_AB_0[thread_id].transpose() * S_0_p * SWAPPED_CD_p[thread_id];
        Prod_0_m[thread_id] = 
          SWAPPED_AB_0[thread_id].transpose() * S_0_m * SWAPPED_CD_m[thread_id];
       
        
        double OvLp_AC_0_p =
          Prod_0_p[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_AC_0_m =
          Prod_0_m[thread_id].block(0,0,nocc+1,nocc+1).determinant();
     
     
        // BD Overlap
        SWAPPED_AB_0[thread_id].setZero();
        SWAPPED_CD_p[thread_id].setZero();
        SWAPPED_CD_m[thread_id].setZero();
     
        Prod_0_p[thread_id].setZero();
        Prod_0_m[thread_id].setZero();
     
        SWAPPED_AB_0[thread_id] = MO_0;
        SWAPPED_CD_p[thread_id] = MO_p;
        SWAPPED_CD_m[thread_id] = MO_m;
     
        SWAPPED_AB_0[thread_id].col(nocc).swap(SWAPPED_AB_0[thread_id].col(nocc+b));
        SWAPPED_CD_p[thread_id].col(nocc).swap(SWAPPED_CD_p[thread_id].col(nocc+d));
        SWAPPED_CD_m[thread_id].col(nocc).swap(SWAPPED_CD_m[thread_id].col(nocc+d));
     
        Prod_0_p[thread_id] = 
          SWAPPED_AB_0[thread_id].transpose() * S_0_p * SWAPPED_CD_p[thread_id];
        Prod_0_m[thread_id] = 
          SWAPPED_AB_0[thread_id].transpose() * S_0_m * SWAPPED_CD_m[thread_id];
       
        
        double OvLp_BD_0_p =
          Prod_0_p[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_BD_0_m =
          Prod_0_m[thread_id].block(0,0,nocc+1,nocc+1).determinant();
     
        // AD Overlap
        SWAPPED_AB_0[thread_id].setZero();
        SWAPPED_CD_p[thread_id].setZero();
        SWAPPED_CD_m[thread_id].setZero();
     
        Prod_0_p[thread_id].setZero();
        Prod_0_m[thread_id].setZero();
     
        SWAPPED_AB_0[thread_id] = MO_0;
        SWAPPED_CD_p[thread_id] = MO_p;
        SWAPPED_CD_m[thread_id] = MO_m;
     
        SWAPPED_AB_0[thread_id].col(nocc).swap(SWAPPED_AB_0[thread_id].col(nocc+a));
        SWAPPED_CD_p[thread_id].col(nocc).swap(SWAPPED_CD_p[thread_id].col(nocc+d));
        SWAPPED_CD_m[thread_id].col(nocc).swap(SWAPPED_CD_m[thread_id].col(nocc+d));
     
     
        Prod_0_p[thread_id] = 
          SWAPPED_AB_0[thread_id].transpose() * S_0_p * SWAPPED_CD_p[thread_id];
        Prod_0_m[thread_id] = 
          SWAPPED_AB_0[thread_id].transpose() * S_0_m * SWAPPED_CD_m[thread_id];
       
        
        double OvLp_AD_0_p =
          Prod_0_p[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_AD_0_m =
          Prod_0_m[thread_id].block(0,0,nocc+1,nocc+1).determinant();
     
     
        // BC Overlap
        SWAPPED_AB_0[thread_id].setZero();
        SWAPPED_CD_p[thread_id].setZero();
        SWAPPED_CD_m[thread_id].setZero();
     
        Prod_0_p[thread_id].setZero();
        Prod_0_m[thread_id].setZero();
     
        SWAPPED_AB_0[thread_id] = MO_0;
        SWAPPED_CD_p[thread_id] = MO_p;
        SWAPPED_CD_m[thread_id] = MO_m;
     
        SWAPPED_AB_0[thread_id].col(nocc).swap(SWAPPED_AB_0[thread_id].col(nocc+b));
        SWAPPED_CD_p[thread_id].col(nocc).swap(SWAPPED_CD_p[thread_id].col(nocc+c));
        SWAPPED_CD_m[thread_id].col(nocc).swap(SWAPPED_CD_m[thread_id].col(nocc+c));
     
     
        Prod_0_p[thread_id] = 
          SWAPPED_AB_0[thread_id].transpose() * S_0_p * SWAPPED_CD_p[thread_id];
        Prod_0_m[thread_id] = 
          SWAPPED_AB_0[thread_id].transpose() * S_0_m * SWAPPED_CD_m[thread_id];
       
        
        double OvLp_BC_0_p =
          Prod_0_p[thread_id].block(0,0,nocc+1,nocc+1).determinant();
        double OvLp_BC_0_m =
          Prod_0_m[thread_id].block(0,0,nocc+1,nocc+1).determinant();
     
     
        double OvLp_Swap_0_p = OvLp_AC_0_p*OvLp_BD_0_p + OvLp_AD_0_p*OvLp_BC_0_p;
        double OvLp_Swap_0_m = OvLp_AC_0_m*OvLp_BD_0_m + OvLp_AD_0_m*OvLp_BC_0_m;
     
        double fact = 1;
        if( a == b ) fact *= std::sqrt(0.5);
        if( c == d ) fact *= std::sqrt(0.5);
     
        #pragma omp critical
        {
          NACME(iSt,jSt) += fact*T_0(ab,iSt)*(
            T_p(cd,jSt)*OvLp_Swap_0_p - 
            T_m(cd,jSt)*OvLp_Swap_0_m
          ) / (2*step);
        }
      }
    }
  }
  auto NACMEEnd = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = NACMEEnd - NACMEStart;
  cout << "ES -> ES NACME took " << elapsed.count() << " s" << endl;

  return NACME;
};

Derivatives computeNAC2pt1D(Molecule &geom_0, Molecule &geom_p,
  Molecule &geom_m, FileIO &fileio_0, FileIO &fileio_p, FileIO &fileio_m,
  Controls &controls, std::string &basisName, RESPONSE_TYPE respType,
  int nFreq, double step){

  bool debug = true;
 
  BasisSet basis0;
  BasisSet basisp;
  BasisSet basism;

  AOIntegrals aoints_0;
  AOIntegrals aoints_p;
  AOIntegrals aoints_m;

  SingleSlater<double> ss_0;
  SingleSlater<double> ss_p;
  SingleSlater<double> ss_m;

  MOIntegrals<double> moints_0;
  MOIntegrals<double> moints_p;
  MOIntegrals<double> moints_m;

  Response<double> resp_0;
  Response<double> resp_p;
  Response<double> resp_m;



  basis0.findBasisFile(basisName);
  basisp.findBasisFile(basisName);
  basism.findBasisFile(basisName);

  basis0.communicate(fileio_0);
  basisp.communicate(fileio_p);
  basism.communicate(fileio_m);

  basis0.parseGlobal();
  basisp.parseGlobal();
  basism.parseGlobal();

  basis0.constructLocal(&geom_0);
  basisp.constructLocal(&geom_p);
  basism.constructLocal(&geom_m);
  
  basis0.makeMaps(1,&geom_0);
  basisp.makeMaps(1,&geom_p);
  basism.makeMaps(1,&geom_m);


  RealMatrix S_0_p = genSpx(basis0,basisp);
  RealMatrix S_0_m = genSpx(basis0,basism);

  RealMatrix S_0_0,S_p_p,S_m_m;

  if(debug){
    S_0_0 = genSpx(basis0,basis0);
    S_p_p = genSpx(basisp,basisp);
    S_m_m = genSpx(basism,basism);
  }

  aoints_0.communicate(geom_0,basis0,fileio_0,controls);
  aoints_p.communicate(geom_p,basisp,fileio_p,controls);
  aoints_m.communicate(geom_m,basism,fileio_m,controls);

  ss_0.communicate(geom_0,basis0,aoints_0,fileio_0,controls);
  ss_p.communicate(geom_p,basisp,aoints_p,fileio_p,controls);
  ss_m.communicate(geom_m,basism,aoints_m,fileio_m,controls);

  moints_0.communicate(geom_0,basis0,fileio_0,controls,
    aoints_0, ss_0);
  moints_p.communicate(geom_p,basisp,fileio_p,controls,
    aoints_p, ss_p);
  moints_m.communicate(geom_m,basism,fileio_m,controls,
    aoints_m, ss_m);
  
  aoints_0.initMeta();
  aoints_p.initMeta();
  aoints_m.initMeta();

  aoints_0.integralAlgorithm = AOIntegrals::INCORE;
  aoints_p.integralAlgorithm = AOIntegrals::INCORE;
  aoints_m.integralAlgorithm = AOIntegrals::INCORE;

  ss_0.setRef(SingleSlater<double>::RHF);
  ss_p.setRef(SingleSlater<double>::RHF);
  ss_m.setRef(SingleSlater<double>::RHF);

  ss_0.isClosedShell = true;
  ss_p.isClosedShell = true;
  ss_m.isClosedShell = true;

  ss_0.initMeta();
  ss_p.initMeta();
  ss_m.initMeta();

  ss_0.genMethString();
  ss_p.genMethString();
  ss_m.genMethString();

  aoints_0.alloc();
  aoints_p.alloc();
  aoints_m.alloc();

  ss_0.alloc();
  ss_p.alloc();
  ss_m.alloc();

  // SCF (X,Y)
  cout << "  Performing SCF (X,Y)" << endl;
  ss_0.formGuess();
  ss_0.formFock();
  ss_0.computeEnergy();
  ss_0.SCF();
  ss_0.computeProperties();
  ss_0.printProperties();

  // SCF (X+DX,Y)
  cout << "  Performing SCF (X+DX,Y)" << endl;
  ss_p.formGuess();
  ss_p.formFock();
  ss_p.computeEnergy();
  ss_p.SCF();
  ss_p.computeProperties();
  ss_p.printProperties();

  // SCF (X-DX,Y)
  cout << "  Performing SCF (X-DX,Y)" << endl;
  ss_m.formGuess();
  ss_m.formFock();
  ss_m.computeEnergy();
  ss_m.SCF();
  ss_m.computeProperties();
  ss_m.printProperties();

  if(debug) {
    cout << endl;
    cout << "  Checking | C - C' | Before Phase Check:" << endl;
    
    cout << "  | C(X,Y) - C(X+DX,Y) | = " 
         << diffNorm((*ss_0.moA()),(*ss_p.moA())) 
         << endl;  
    cout << "  | C(X,Y) - C(X-DX,Y) | = " 
         << diffNorm((*ss_0.moA()),(*ss_m.moA())) 
         << endl;  
  }



  cout << endl << "  Performing Phase Check on MOs" << endl;

  RealMatrix SMO_0_0,SMO_p_p,SMO_m_m;

  if(debug){
  // < p (X,Y) | S[(X,Y),(X,Y)] | q(X,Y) >
   SMO_0_0 = ss_0.moA()->transpose() * S_0_0 * (*ss_0.moA());
  // < p (X+DX,Y) | S[(X+DX,Y),(X+DX,Y)] | q(X+DX,Y)(* >
   SMO_p_p = ss_p.moA()->transpose() * S_p_p * (*ss_p.moA());
  // < p (X-DX,Y) | S[(X-DX,Y),(X-DX,Y)] | q(X-DX,Y)(* >
   SMO_m_m = ss_m.moA()->transpose() * S_m_m * (*ss_m.moA());
  }

  // < p (X,Y) | S[(X,Y),(X+DX,Y)] | q(X+DX,Y) >
  RealMatrix SMO_0_p = 
    ss_0.moA()->transpose() * S_0_p * (*ss_p.moA());
  // < p (X,Y) | S[(X,Y),(X-DX,Y)] | q(X-DX,Y) >
  RealMatrix SMO_0_m = 
    ss_0.moA()->transpose() * S_0_m * (*ss_m.moA());



/*
  Eigen::VectorXd Trans_mo_p;
  Eigen::VectorXd Trans_mo_m;
  Trans_mo_p = checkPhase(ss_0.moA()->data(),
    ss_p.moA()->data(),basis0.nBasis(),basis0.nBasis());
  Trans_mo_m = checkPhase(ss_0.moA()->data(),
    ss_m.moA()->data(),basis0.nBasis(),basis0.nBasis());
*/

  RealMatrix O_0_p(SMO_0_p);
  RealMatrix O_0_m(SMO_0_m);

  prettyPrint(cout,*ss_p.epsA(),"EPS A(X+DX)");
  prettyPrint(cout,*ss_m.epsA(),"EPS A(X-DX)");
  prettyPrint(cout,*ss_0.epsA(),"EPS A(X)");

  cout << "  X -> X+DX MO Mapping" << endl;
  for(auto iMO = 0; iMO < SMO_0_p.cols(); iMO++){
    RealMatrix::Index maxRow;
    double maxVal = O_0_p.cwiseAbs().col(iMO).maxCoeff(&maxRow);
    std::stringstream mapStr, ovlpStr; 
    mapStr  << "    " << iMO << " -> " << maxRow;
    ovlpStr << "< " << iMO << "(X) | " << maxRow << "(X+DX) >"; 
            

    cout << std::left << std::setw(10) << mapStr.str();
    cout << "\t";
    cout << std::left << std::setw(25) << ovlpStr.str() << " = " << maxVal;
    cout << endl;
  }

  cout << endl;
  cout << "  X -> X-DX MO Mapping" << endl;
  for(auto iMO = 0; iMO < SMO_0_m.cols(); iMO++){
    RealMatrix::Index maxRow;
    O_0_m.cwiseAbs().col(iMO).maxCoeff(&maxRow);
    double maxVal = O_0_m(maxRow,iMO);
    std::stringstream mapStr, ovlpStr; 
    mapStr  << "    " << iMO << " -> " << maxRow;
    ovlpStr << "< " << iMO << "(X) | " << maxRow << "(X-DX) >"; 
            

    cout << std::left << std::setw(10) << mapStr.str();
    cout << "\t";
    cout << std::left << std::setw(25) << ovlpStr.str() << " = " << maxVal;
    cout << endl;
  }
  cout << endl;

  for(auto mu = 0; mu < SMO_0_p.rows(); mu++)
  for(auto nu = 0; nu < SMO_0_p.rows(); nu++){
    if(std::abs(O_0_p(mu,nu)) < 1e-2) O_0_p(mu,nu) = 0.0;
    else if(O_0_p(mu,nu) > 0.0)       O_0_p(mu,nu) = 1.0;
    else                                  O_0_p(mu,nu) = -1.0;
    if(std::abs(O_0_m(mu,nu)) < 1e-2) O_0_m(mu,nu) = 0.0;
    else if(O_0_m(mu,nu) > 0.0)       O_0_m(mu,nu) = 1.0;
    else                                  O_0_m(mu,nu) = -1.0;
  }

  RealMatrix TMP;
  TMP = (*ss_p.moA()) * O_0_p; 
  (*ss_p.moA()) = TMP;
  TMP = (*ss_m.moA()) * O_0_m; 
  (*ss_m.moA()) = TMP;
  

  if(debug) {
    cout << endl;
    cout << "  Checking | C - C' | After Phase Check:" << endl;
    
    cout << "  | C(X,Y) - C(X+DX,Y) | = " 
         << diffNorm((*ss_0.moA()),(*ss_p.moA())) 
         << endl;  
    cout << "  | C(X,Y) - C(X-DX,Y) | = " 
         << diffNorm((*ss_0.moA()),(*ss_m.moA())) 
         << endl;  
  }

  if(debug){
    // < p (X,Y) | S[(X,Y),(X+DX,Y)] | q(X+DX,Y) >
    SMO_0_p = 
      ss_0.moA()->transpose() * S_0_p * (*ss_p.moA());
    // < p (X,Y) | S[(X,Y),(X-DX,Y)] | q(X-DX,Y) >
    SMO_0_m = 
      ss_0.moA()->transpose() * S_0_m * (*ss_m.moA());
   
    // Quantify deviation from C**H * S * C = I
    cout << endl;
    cout << "  Checking | C(X)**H * S(X,X') * C(X') - I |:" 
         << endl;
   
    cout << "  < (X,Y) | (X,Y) >       = " 
         << diffNormI(SMO_0_0) 
         << endl;
    cout << "  < (X+DX,Y) | (X+DX,Y) > = " 
         << diffNormI(SMO_p_p) 
         << endl;
    cout << "  < (X-DX,Y) | (X-DX,Y) > = " 
         << diffNormI(SMO_m_m) 
         << endl;
   
   
    cout << "  < (X,Y) | (X+DX,Y) >    = " 
         << diffNormI(SMO_0_p) 
         << endl;
    cout << "  < (X,Y) | (X-DX,Y) >    = " 
         << diffNormI(SMO_0_m) 
         << endl;
  }

  moints_0.initMeta();
  moints_p.initMeta();
  moints_m.initMeta();

  resp_0.communicate(ss_0,moints_0,fileio_0);
  resp_p.communicate(ss_p,moints_p,fileio_p);
  resp_m.communicate(ss_m,moints_m,fileio_m);

  resp_0.setMeth(respType);
  resp_p.setMeth(respType);
  resp_m.setMeth(respType);

  resp_0.doSA();
  resp_p.doSA();
  resp_m.doSA();

  resp_0.setNSek(nFreq);
  resp_p.setNSek(nFreq);
  resp_m.setNSek(nFreq);

  resp_0.doFull();
  resp_p.doFull();
  resp_m.doFull();

  cout << endl;
  cout << "  Performing Response Calculations (";
  if(respType == RESPONSE_TYPE::CIS)       cout << "CIS" ;
  else if(respType == RESPONSE_TYPE::PPTDA) cout << "PPTDA" ;
  cout << ")" << endl;

  cout << "  Performing Response (X,Y)" << endl;
  resp_0.doResponse();
  cout << "  Performing Response (X+DX,Y)" << endl;
  resp_p.doResponse();
  cout << "  Performing Response (X-DX,Y)" << endl;
  resp_m.doResponse();





  // Start Differentiation things
  //


  
  double scf_0 = ss_0.totalEnergy;
  double scf_p = ss_p.totalEnergy;
  double scf_m = ss_m.totalEnergy;


  if(debug) {
    cout << endl;
    cout << "  SCF ENERGIES:" << endl;
    cout << "  SCF Energy (X,Y):    " << std::setprecision(8) 
         << scf_0 << " Eh" << endl;
    cout << "  SCF Energy (X+DX,Y): " << std::setprecision(8) 
         << scf_p << " Eh" << endl;
    cout << "  SCF Energy (X-DX,Y): " << std::setprecision(8) 
         << scf_m << " Eh" << endl;
  }
  

  Eigen::VectorXd freq_0,freq_p,freq_m; 

  if(this->respType_ == RESPONSE_TYPE::CIS){
    freq_0 = resp_0.template frequencies<SINGLETS>().head(nFreq);
    freq_p = resp_p.template frequencies<SINGLETS>().head(nFreq);
    freq_m = resp_m.template frequencies<SINGLETS>().head(nFreq);
  } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
    freq_0 = resp_0.template frequencies<A_PPTDA_SINGLETS>().head(nFreq);
    freq_p = resp_p.template frequencies<A_PPTDA_SINGLETS>().head(nFreq);
    freq_m = resp_m.template frequencies<A_PPTDA_SINGLETS>().head(nFreq);
  }

  if(debug) {
    cout << endl << "  Excitation template frequencies:" << endl;
   
    for(auto iSt = 0; iSt < nFreq; iSt++){
      cout << "  W(" << iSt << ") (X,Y):    " 
           << std::setprecision(8) << freq_0[iSt] << " Eh" 
           << endl;
      cout << "  W(" << iSt << ") (X+DX,Y): " 
           << std::setprecision(8) << freq_p[iSt] << " Eh" 
           << endl;
      cout << "  W(" << iSt << ") (X-DX,Y): " 
           << std::setprecision(8) << freq_m[iSt] << " Eh" 
           << endl;
      cout << endl;
    }
  }





  RealMatrix T_0,T_p,T_m;

  if(this->respType_ == RESPONSE_TYPE::CIS){
    T_0 = resp_0.template transDen<SINGLETS>().block(0,0,
     resp_0.template nMatDim<SINGLETS>(),nFreq);
    T_p = resp_p.template transDen<SINGLETS>().block(0,0,
     resp_p.template nMatDim<SINGLETS>(),nFreq);
    T_m = resp_m.template transDen<SINGLETS>().block(0,0,
     resp_m.template nMatDim<SINGLETS>(),nFreq);
  } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
    T_0 = resp_0.template transDen<A_PPTDA_SINGLETS>().block(0,0,
     resp_0.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
    T_p = resp_p.template transDen<A_PPTDA_SINGLETS>().block(0,0,
     resp_p.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
    T_m = resp_m.template transDen<A_PPTDA_SINGLETS>().block(0,0,
     resp_m.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
  }



  if(debug){
    cout << endl;
    cout << "  Checking | T - T' | Before Phase Check:" << endl;
    
    cout << "  | T(X,Y) - T(X+DX,Y) | = " << diffNorm(T_0,T_p) 
         << endl;  
    cout << "  | T(X,Y) - T(X-DX,Y) | = " << diffNorm(T_0,T_m) 
         << endl;  
  }


  Eigen::VectorXd Trans_t_p;
  Eigen::VectorXd Trans_t_m;
  if(this->respType_ == RESPONSE_TYPE::CIS){
    Trans_t_p = checkPhase(T_0.data(),T_p.data(),
      resp_p.template nMatDim<SINGLETS>(),nFreq);
    Trans_t_m = checkPhase(T_0.data(),T_m.data(),
      resp_m.template nMatDim<SINGLETS>(),nFreq);
  } else if(this->respType_ == RESPONSE_TYPE::PPTDA){
    Trans_t_p = checkPhase(T_0.data(),T_p.data(),
      resp_p.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
    Trans_t_m = checkPhase(T_0.data(),T_m.data(),
      resp_m.template nMatDim<A_PPTDA_SINGLETS>(),nFreq);
  }


  if(debug){
    cout << endl;
    cout << "  Checking | T - T' | After Phase Check:" << endl;
    
    cout << "  | T(X,Y) - T(X+DX,Y) | = " << diffNorm(T_0,T_p) 
         << endl;  
    cout << "  | T(X,Y) - T(X-DX,Y) | = " << diffNorm(T_0,T_m) 
         << endl;  
  }

  if(debug) {
    cout << endl << " Checking < T | T > = 1 " << endl;
    for(auto iSt = 0; iSt < nFreq; iSt++){
      cout << "   State " << iSt << ":" << endl;
      cout << "    < T(X,Y) | T(X,Y) >       = " 
           << selfInner(T_0.col(iSt)) 
           << endl;
      cout << "    < T(X+DX,Y) | T(X+DX,Y) > = " 
           << selfInner(T_p.col(iSt)) 
           << endl;
      cout << "    < T(X-DX,Y) | T(X-DX,Y) > = " 
           << selfInner(T_m.col(iSt)) 
           << endl;
    }
  }








  // Calculate Derivatives
    
  // SCF energy gradient
  double gsdx = (scf_p - scf_m) / (2*step);

  if(debug)
    cout << endl << "  GS Gradient = " << gsdx << endl;

  // Excitation frequency gradient
  Eigen::VectorXd freqDX = (freq_p - freq_m)/(2*step);

  if(debug){
    cout << endl << "  ES Gradients:" << endl;
    for(auto iSt = 0; iSt < nFreq; iSt++){
      cout << "   W(" << iSt << ")' = " << freqDX(iSt) << endl; 
    }
  }





  // Assemble ground to excited state couplings
  //
  
  int NOCC = ss_0.nOccA();
  int NVIR = ss_0.nVirA();

  if(debug){
    // Compute Overlaps of wavefunctions at different geometries
    double OvLp_0_0 = 
      SMO_0_0.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
    double OvLp_p_p = 
      SMO_p_p.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
    double OvLp_m_m = 
      SMO_m_m.block(0,0,NOCC,NOCC).eigenvalues().prod().real();

    double OvLp_0_p = 
      SMO_0_p.block(0,0,NOCC,NOCC).eigenvalues().prod().real();
    double OvLp_0_m = 
      SMO_0_m.block(0,0,NOCC,NOCC).eigenvalues().prod().real();


    cout << endl;
    cout << "  Checking Wavefunction Overlaps < X | X' >:" << endl;
    cout << "  < (X,Y) | (X,Y) >          = " << OvLp_0_0 << endl;
    cout << "  < (X+DX,Y) | (X+DX,Y) >    = " << OvLp_p_p << endl;
    cout << "  < (X-DX,Y) | (X-DX,Y) >    = " << OvLp_m_m << endl;

    cout << "  < (X,Y) | (X+DX,Y) >       = " << OvLp_0_p << endl;
    cout << "  < (X,Y) | (X-DX,Y) >       = " << OvLp_0_m << endl;
  }
  


  Eigen::VectorXd NAC_ES_GS;
  RealMatrix      NAC_ES_ES;

  if(respType == RESPONSE_TYPE::CIS){
    NAC_ES_GS = ES_GS_NACME_CIS(nFreq,true,NOCC,NVIR,step,
      T_0,(*ss_0.moA()),(*ss_p.moA()),(*ss_m.moA()),
      S_0_p,S_0_m);
    NAC_ES_ES = ES_ES_NACME_CIS(nFreq,true,NOCC,NVIR,step,
      T_0,T_p,T_m,(*ss_0.moA()),(*ss_p.moA()),
      (*ss_m.moA()),S_0_p,S_0_m);

  } else if(respType == RESPONSE_TYPE::PPTDA) {
    NAC_ES_ES = ES_ES_NACME_PPTDA(nFreq,true,NOCC,NVIR,step,
      T_0,T_p,T_m,(*ss_0.moA()),(*ss_p.moA()),
      (*ss_m.moA()),S_0_p,S_0_m);
  }
  

  if(debug){
    Eigen::VectorXd NAC_GS_ES;

    if(respType == RESPONSE_TYPE::CIS) {
      NAC_GS_ES = GS_ES_NACME_CIS(nFreq,false,NOCC,NVIR,step,
        T_0,T_p,T_m,(*ss_0.moA()),(*ss_p.moA()),
        (*ss_m.moA()),S_0_p,S_0_m);
     
      prettyPrint(cout,NAC_ES_GS*0.529177,"ES->GS NACME DX");
      prettyPrint(cout,NAC_GS_ES*0.529177,"GS->ES NACME DX");
      prettyPrint(cout,NAC_ES_GS + NAC_GS_ES,"DIFF DX");


    }
  }

  prettyPrint(cout,NAC_ES_ES*0.529177,"ES->ES DX");

  Derivatives derv;

  derv.GS_GRAD = gsdx;
  derv.ES_GRAD = freqDX; 
  derv.ES_GS_NACME = NAC_ES_GS;
  derv.ES_ES_NACME = NAC_ES_ES;

  return derv;
  
}
