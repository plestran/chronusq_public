/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
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
#include <classtools.h>
using ChronusQ::Molecule;
using ChronusQ::BasisSet;
using ChronusQ::Controls;

/************************************************/
/* read input files and initialize everything   */
/* memory allocations are done here too         */
/************************************************/
namespace ChronusQ {
void readInput(FileIO * fileio, Molecule * mol, BasisSet * basis, Controls * controls,
               BasisSet * dfBasis) {
  int i, j, n, readInt;
  std::string readString;
  fileio->in >> readString;
  while(!(fileio->in.eof())) {
    readString=stringupper(readString);
    if(!readString.compare("/*")) while(readString.compare("*/")) fileio->in >> readString;
    else if(!readString.compare("$CHARGE")) {
      fileio->in >> readInt;
      mol->setCharge(readInt);
    } else if(!readString.compare("$SPIN")) {
      fileio->in >> readInt;
      mol->setMultip(readInt);
    } else if(!readString.compare("$EXTRA")) {
      fileio->in>>readString;
      readString=stringupper(readString);
      if(!readString.compare("RESTART")) controls->restart=true;
    } else if(!readString.compare("$PRINT")) {
      fileio->in >> (controls->printLevel);
    } else if(!readString.compare("$METHOD")) {
      fileio->in >> readString;
      readString=stringupper(readString);
      if(!readString.compare("HF")) controls->HF=true;
      else if(!readString.compare("ROHF")){
        controls->HF     = true;
        controls->doCUHF = true;
      } else if(!readString.compare("GHF")){
        controls->HF = true;
        controls->doTCS = true;
      } else if(!readString.compare("DFT")){
	controls->HF=false;
	controls->DFT=true;
      };
    } else if(!readString.compare("$GEOM")) {
      fileio->in >> readString;
      readString=stringupper(readString);  
      fstream *geomRead;
      if(!readString.compare("READ")) {
        geomRead = &fileio->in;
      } else if(!readString.compare("FILE")) {
        fileio->in >> readString;
        geomRead = new fstream(readString,ios::in);
        if(geomRead->fail()) CErr("Unable to open "+std::string(readString),fileio->out); 
        else fileio->out << "Reading geometry from " << readString << endl;
      } else {
        CErr("Unrecognized GEOM option: " + std::string(readString),fileio->out);
      }
      mol->readMolecule(fileio,*geomRead);
      if(!(geomRead==&fileio->in)){
//      fileio->out << "Closing " << readString << endl;
        geomRead->close();
        delete geomRead;
      }
    } else if(!readString.compare("$BASIS")) {
      //basis->readBasisSet(fileio,mol);
      basis->basisSetRead(fileio,mol,controls);
    } else if(!readString.compare("$DFBASIS")) {
      fileio->in >> readString;
      readString = stringupper(readString);
      if(!readString.compare("ON")) controls->doDF = true;
      if(controls->doDF) {
        //basis->readBasisSet(fileio,mol);
        dfBasis->basisSetRead(fileio,mol,controls);
      }
//dbwys
    } else if(!readString.compare("$NSMP")) {
      fileio->in >> readInt;
      controls->readSMP(readInt);
    } else if(!readString.compare("$GUESS")) {
      fileio->in>>readString;
      readString=stringupper(readString);  
      if(!readString.compare("INPUT")) {
        controls->guess = 1;
      } else if(!readString.compare("GAUMATEL")) {
        controls->guess = 2;
        fileio->in>>controls->gauMatElName;
      } else if(!readString.compare("GAUFCHK")) {
        controls->guess = 3;
        fileio->in>>controls->gauFChkName;
      };
    } else if(!readString.compare("$SCF")) {
      fileio->in>>readString;
      readString=stringupper(readString);
      if(!readString.compare("OFF"))
        controls->optWaveFunction = false;
      else if(!readString.compare("ON"))
        controls->optWaveFunction = true;
      else if(!readString.compare("DIRECT"))
        controls->directTwoE = true;
      else if(!readString.compare("INCORE")){
        controls->directTwoE = false;
        controls->buildn4eri = true;
      } else if(!readString.compare("COMPLEX")){
        controls->doComplex = true;
      } else {
        CErr("Input SCF Option Not Recognized",fileio->out);
      }
    } else if(!readString.compare("$DEBUG")) {
      fileio->in>>readString;
      readString=stringupper(readString);
      controls->readDebug(readString);
    } else if(!readString.compare("$PSCF")) {
      controls->readPSCF(fileio->in,fileio->out);
    } else if(!readString.compare("$UNITTEST")) {
      controls->doUnit = true;
      fileio->in >> readString;
      readString = stringupper(readString);
      if(!readString.compare("SCF")) 
        controls->unitTest = Controls::UnitSCF;
      else if(!readString.compare("RESP")) 
        controls->unitTest = Controls::UnitResp;
      else if(!readString.compare("RT")) 
        controls->unitTest = Controls::UnitRT;
      else
        CErr("Unit Test Option: "+readString+" not recognized",fileio->out);
    } else if(!readString.compare("$FIELD")) {
      for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
        fileio->in >> controls->field_[iXYZ];
      } 
    } else if(!readString.compare("$RT")) {
      controls->doRealTime = true;
      fileio->in>>readString;
      readString=stringupper(readString);
      if(!readString.compare("MAXSTEP")) {
        fileio->in >> (controls->rtMaxSteps);
      } else if(!readString.compare("TIMESTEP")) {
        fileio->in >> (controls->rtTimeStep);
      } else if(!readString.compare("EDFIELD")) {
        for(auto iXYZ = 0; iXYZ < 3; iXYZ++){
          fileio->in >> controls->rtField_[iXYZ];
        } 
      } else if(!readString.compare("TIME_ON")) {
        fileio->in >> (controls->rtTOn_);
      } else if(!readString.compare("TIME_OFF")) {
        fileio->in >> (controls->rtTOff_);
      } else if(!readString.compare("FREQUENCY")) {
        fileio->in >> (controls->rtFreq_);
      } else if(!readString.compare("PHASE")) {
        fileio->in >> (controls->rtPhase_);
      } else if(!readString.compare("SIGMA")) {
        fileio->in >> (controls->rtSigma_);
      } else if(!readString.compare("ENVELOPE")) {
        fileio->in>>readString;
        readString=stringupper(readString);
        if(!readString.compare("PW")) {
          controls->rtEnvelope_ = 0; 
        } else if(!readString.compare("LINEAR_RAMP")) {
          controls->rtEnvelope_ = 1;
        } else if(!readString.compare("GAUSSIAN")) {
          controls->rtEnvelope_ = 2;
        } else if(!readString.compare("STEP")) {
          controls->rtEnvelope_ = 3;
        } else if(!readString.compare("SINE_SQUARE")) {
          CErr("Real Time Envelope Option: "+readString+" not yet implemented. \n",fileio->out); 
          controls->rtEnvelope_ = 4;
        } else { 
          CErr("Real Time Envelope Option: "+readString+" not recognized. \n"+
               "Try PW, LINEAR_RAMP, GAUSSIAN, STEP, or SINE_SQUARE",fileio->out); 
        }
      } else if(!readString.compare("ORTHO")) {
        fileio->in>>readString;
        readString=stringupper(readString);
        if(!readString.compare("LOWDIN"))
          controls->rtTypeOrtho = 0;
        else if(!readString.compare("CHOLESKY"))
          controls->rtTypeOrtho = 1;
        else if(!readString.compare("CANONICAL"))
          controls->rtTypeOrtho = 2;
        else 
          CErr("Real Time Orthogonalization Option: "+readString+" not recognized. \n"+
               "Try LOWDIN, CHOLESKY, or CANONICAL",fileio->out); 
      } else if(!readString.compare("INIDEN")) {
        fileio->in>>readString;
        readString=stringupper(readString);
        if(!readString.compare("SCF"))
          controls->rtInitDensity = 0;
        else if(!readString.compare("SWAP")) {
          controls->rtInitDensity = 1;
          fileio->in>>readString;
          readString=stringupper(readString);
          if(!readString.compare("ALPHA")) {
            fileio->in >> (controls->rtSwapMOA);
          } else if(!readString.compare("BETA")) {
            fileio->in >> (controls->rtSwapMOB);
          } else {
            CErr("Real Time Initial Density SWAP Option: "+readString+" not recognized. \n"+
                 "Try ALPHA or BETA",fileio->out); 
          } 
        } else if(!readString.compare("READ_AO"))
          controls->rtInitDensity = 2;
        else if(!readString.compare("READ_ORTHO"))
          controls->rtInitDensity = 3;
        else {
          CErr("Real Time Initial Density Option: "+readString+" not recognized. \n"+
                "Try SCF, SWAP, READ_AO, or READ_ORTHO",fileio->out); 
        }
      } else if(!readString.compare("UPROP")) {
        fileio->in>>readString;
        readString=stringupper(readString);
        if(!readString.compare("EIGEN"))
          controls->rtMethFormU = 0;
        else if(!readString.compare("TAYLOR"))
          controls->rtMethFormU = 1;
        else {
          CErr("Real Time U Matrix / Propagator Option: "+readString+" not recognized. \n"+
               "Try EIGEN or TAYLOR",fileio->out); 
        }
      } else if(!readString.compare("DEFAULT")) {
          continue;
      } else {
          CErr("Real Time Option: "+readString+" not recognized. \n"+
               "Valid options: \n"+
               "\t MAXSTEP:   Maximum number of time steps.       Default = 10 \n"+
               "\t TIMESTEP:  Size of time steps (au).            Default = 0.05 \n"+
               "\t EDFIELD:   Electric dipole Ex,Ey,Ez (au).      Default = 0.0 0.0 0.0 \n"+
               "\t TIME_ON:   Time (fs) field is applied.         Default = 0.0 \n"+ 
               "\t TIME_OFF:  Time (fs) field is removed.         Default = 10000.0 \n"+ 
               "\t FREQUENCY: Field frequency (eV).               Default = 0.0 \n"+ 
               "\t PHASE:     Field phase offset (rad).           Default = 0.0 \n"+ 
               "\t SIGMA:     For Gaussian envelope, this sets    Default = 0.0 \n"+ 
               "\t              the range of the frequency\n"+ 
               "\t              (FWHM, in eV)\n"+ 
               "\t ENVELOPE:  Type of field envelope function.    Default = PW \n"+
               "\t ORTHO:     Type of orthogonalization.          Default = LOWDIN \n"+
               "\t INIDEN:    Initial density for system.         Default = SCF \n"+
               "\t UPROP:     How the propagator is formed.       Default = EIGEN \n"+ 
               "\t DEFAULT:   RT-TDSCF with default settings.\n"+ 
               "\n Note: if you get 'Real Time Option: $RT not recognized',  try '$RT DEFAULT' instead of just '$RT'. \n",fileio->out); 
      
      }
    }
    fileio->in >> readString;
  };
};
/********************************/
/* trace of producto of two     */
/* symmetric RealMatrices       */
/********************************/
double traceSymm(RealMatrix* a, RealMatrix* b) {
  if(a->size()!=b->size()) CErr("Only able to trace matricies of the same size"); // FIXME this is depreciated
  double tmpVal = 0.0;
  int i;
  for(i=0;i<a->size();i++) tmpVal+=a->data()[i]*b->data()[i];

  return tmpVal;
 };

void printUnitInfo(Controls * controls, SingleSlater<double> * singleSlater, SDResponse<double> * sdResponse,
                   RealTime<double> * realTime){
  if(controls->unitTest == Controls::UnitSCF)
    cout << std::setprecision(10) << singleSlater->totalEnergy << "/"
         << std::setprecision(4)
         << (*singleSlater->dipole())(0)/phys.debye << "/"
         << (*singleSlater->dipole())(1)/phys.debye << "/"
         << (*singleSlater->dipole())(2)/phys.debye << "/"
         << (*singleSlater->quadpole())(0,0)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(1,1)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(2,2)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(0,1)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(0,2)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(1,2)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(0,0)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(1,1)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(2,2)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(0,1)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(0,2)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(1,2)*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,0,0)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(1,1,1)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(2,2,2)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,1,1)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,0,1)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,0,2)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,2,2)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(1,2,2)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(1,1,2)*phys.bohr*phys.bohr/phys.debye << "/"
    << (*singleSlater->octpole())(0,1,2)*phys.bohr*phys.bohr/phys.debye << endl;
  else if(controls->unitTest == Controls::UnitResp){
    for(auto iSt = 0; iSt < sdResponse->nSek(); iSt++){
      cout << (*sdResponse->omega())(iSt)*phys.eVPerHartree << "," << (*sdResponse->oscStrength())(0,iSt+1);
      if(iSt != (sdResponse->nSek() - 1)) cout << "/";
    }
    cout << endl;
  }
  else if(controls->unitTest == Controls::UnitRT){
    cout << "RT DOUBLE UNIT" << endl;
    }
}

void printUnitInfo(Controls * controls, SingleSlater<dcomplex> * singleSlater, SDResponse<double> * sdResponse, 
                   RealTime<dcomplex> * realTime){
  if(controls->unitTest == Controls::UnitSCF)
    cout << std::setprecision(10) << singleSlater->totalEnergy << "/"
         << std::setprecision(4)
         << (*singleSlater->dipole())(0)/phys.debye << "/"
         << (*singleSlater->dipole())(1)/phys.debye << "/"
         << (*singleSlater->dipole())(2)/phys.debye << "/"
         << (*singleSlater->quadpole())(0,0)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(1,1)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(2,2)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(0,1)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(0,2)*phys.bohr/phys.debye << "/"
         << (*singleSlater->quadpole())(1,2)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(0,0)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(1,1)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(2,2)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(0,1)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(0,2)*phys.bohr/phys.debye << "/"
        << (*singleSlater->tracelessQuadpole())(1,2)*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,0,0)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(1,1,1)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(2,2,2)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,1,1)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,0,1)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,0,2)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(0,2,2)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(1,2,2)*phys.bohr*phys.bohr/phys.debye << "/"
     << (*singleSlater->octpole())(1,1,2)*phys.bohr*phys.bohr/phys.debye << "/"
    << (*singleSlater->octpole())(0,1,2)*phys.bohr*phys.bohr/phys.debye << endl;
  else if(controls->unitTest == Controls::UnitResp){
    for(auto iSt = 0; iSt < sdResponse->nSek(); iSt++){
      cout << (*sdResponse->omega())(iSt)*phys.eVPerHartree << "," << (*sdResponse->oscStrength())(0,iSt+1);
      if(iSt != (sdResponse->nSek() - 1)) cout << "/";
    }
    cout << endl;
  }
  else if(controls->unitTest == Controls::UnitRT){
    cout <<  std::setprecision(4) << (realTime->maxTime())     << "/" 
         << std::setprecision(10) << (realTime->Energy())      << "/"
         <<  std::setprecision(6) << (realTime->EDx())         << "/"
                                  << (realTime->EDy())         << "/"
                                  << (realTime->EDz())         << "/"
                                  << (realTime->EDtot())       << "/"
         << endl;
    }
}
} // namespace ChronusQ
