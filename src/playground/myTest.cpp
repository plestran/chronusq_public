#include <global.h>
#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;


int main() {
    CQMemManager memManager;
    Molecule molecule;
    BasisSet basis;
    //Controls controls;
    AOIntegrals aoints;
    MOIntegrals<dcomplex> moints;
    SingleSlater<dcomplex> singleSlater;
    FileIO fileio("test.inp","test.out");

    memManager.setTotalMem(256e6);
    initCQ();
    //controls.iniControls(); //artifact from Xiaosong's code
    fileio.iniH5Files();  //Creates restart file and scratch file  
    fileio.iniStdGroups(); //Setting up and partitioning files
    CQSetNumThreads(1); //Sets up open MP threads

/*
// Li3
    molecule.setCharge(0);
    molecule.setNTotalE(9);
    molecule.setMultip(6);
    molecule.setNAtoms(3);
    molecule.alloc(); //allocates all memory for the class

    molecule.setIndex(0,HashAtom("Li",0));
    molecule.setIndex(1,HashAtom("Li",0));
    molecule.setIndex(2,HashAtom("Li",0));

    molecule.setCart(0,-1.05,0.0,0.0); //In Angstroms!
    molecule.setCart(1,1.05,0.0,0.0);
    molecule.setCart(2,0.0,1.81865,0.0);
*/

// WATER (H20 +)
    molecule.setCharge(0);
    molecule.setNTotalE(10);
    molecule.setMultip(1);
    molecule.setNAtoms(3);
    molecule.alloc(); //allocates all memory for the class

    molecule.setIndex(0,HashAtom("O",0));
    molecule.setIndex(1,HashAtom("H",0));
    molecule.setIndex(2,HashAtom("H",0));

    molecule.setCart(0,-0.464,0.177,0.0); //In Angstroms!
    molecule.setCart(1,-0.464,1.137,0.0);
    molecule.setCart(2,0.441,-0.143,0.0);
//

    molecule.convBohr();
    molecule.computeNucRep();
//    fileio.out << molecule.energyNuclei() << endl;

    molecule.computeRij();
    molecule.computeI();

    singleSlater.setRef(SingleSlater<dcomplex>::TCS); //TCS == GHF?
    singleSlater.setSCFEneTol(1e-12);
    singleSlater.setNTCS(2);
    singleSlater.isClosedShell = false;
    singleSlater.doDIIS = false;

    basis.findBasisFile("STO3G");
    basis.communicate(fileio);  // This function passes fileio reference 
    basis.parseGlobal(); //Reads entire basis set file into memory
    basis.constructLocal(&molecule);
    basis.makeMaps(&molecule); 
    basis.renormShells();

    aoints.communicate(molecule,basis,fileio,memManager);
    singleSlater.communicate(molecule,basis,aoints,fileio,memManager);
    moints.communicate(molecule,basis,fileio,aoints,singleSlater);

    aoints.setPrintLevel(2);
    aoints.initMeta();
    aoints.integralAlgorithm = AOIntegrals::INCORE;
    aoints.doX2C = true;
    aoints.useFiniteWidthNuclei = true;
    
    singleSlater.initMeta();
    singleSlater.genMethString();

//    fileio.out << "Allocate memory for aoints and singleSlater" << endl;
    aoints.alloc();
    singleSlater.alloc();

//    cout << "Form Guess" << endl;
    singleSlater.formGuess();
//    cout << "Form Fock" << endl;
    singleSlater.formFock();
//    fileio.out << "Compute Energy" << endl;
    singleSlater.computeEnergy();
//    fileio.out << "Do SCF" << endl;
    singleSlater.SCF2();
//    fileio.out << "Compute Properties" << endl;
    singleSlater.computeProperties();
    singleSlater.printProperties();

    finalizeCQ();

};
