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

// H Atom
    molecule.setCharge(0);
    molecule.setNTotalE(1);
    molecule.setMultip(2);
    molecule.setNAtoms(1);
    molecule.alloc();

    molecule.setIndex(0,HashAtom("H",0));
    molecule.setCart(0,0.0,0.0,0.0);


/*
// H2
    molecule.setCharge(0);
    molecule.setNTotalE(2);
    molecule.setMultip(1);
    molecule.setNAtoms(2);
    molecule.alloc(); //allocates all memory for the class

    molecule.setIndex(0,HashAtom("H",0));
    molecule.setIndex(1,HashAtom("H",0));

    molecule.setCart(0,-0.5,0.0,0.0); //In Angstroms!
    molecule.setCart(1,0.5,0.0,0.0);
*/

/*
// O2
    molecule.setCharge(0);
    molecule.setNTotalE(16);
    molecule.setMultip(1);
    molecule.setNAtoms(2);
    molecule.alloc(); //allocates all memory for the class

    molecule.setIndex(0,HashAtom("O",0));
    molecule.setIndex(1,HashAtom("O",0));

    molecule.setCart(0,-0.74,0.0,0.0); //In Angstroms!
    molecule.setCart(1,0.74,0.0,0.0);
*/

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

/*
// WATER (H2O)
    molecule.setCharge(0);
    molecule.setNTotalE(10);
    molecule.setMultip(1);
    molecule.setNAtoms(3);
    molecule.alloc();

    molecule.setIndex(0,HashAtom("O",0));
    molecule.setIndex(1,HashAtom("H",0));
    molecule.setIndex(2,HashAtom("H",0));

    molecule.setCart(0,0.0,0.110843,0.0);
    molecule.setCart(1,0.783809,-0.443452,0.0);
    molecule.setCart(2,-0.783809,-0.443452,0.0);
*/

/*
// Two WATER (H20)
    molecule.setCharge(0);
    molecule.setNTotalE(20);
    molecule.setMultip(1);
    molecule.setNAtoms(6);
    molecule.alloc(); //allocates all memory for the class

    molecule.setIndex(0,HashAtom("O",0));
    molecule.setIndex(1,HashAtom("H",0));
    molecule.setIndex(2,HashAtom("H",0));
    molecule.setIndex(3,HashAtom("O",0));
    molecule.setIndex(4,HashAtom("H",0));
    molecule.setIndex(5,HashAtom("H",0));

    molecule.setCart(0,0.0,0.110843,0.0);
    molecule.setCart(1,0.783809,-0.443452,0.0);
    molecule.setCart(2,-0.783809,-0.443452,0.0);
    molecule.setCart(3,0.0,0.110843,10.0);
    molecule.setCart(4,0.783809,-0.443452,10.0);
    molecule.setCart(5,-0.783809,-0.443452,10.0);
//    molecule.setCart(0,-0.464,0.177,0.0); //In Angstroms!
//    molecule.setCart(1,-0.464,1.137,0.0);
//    molecule.setCart(2,0.441,-0.143,0.0);

*/

/*
// Methanol
    molecule.setCharge(0);
    molecule.setNTotalE(18);
    molecule.setMultip(1);
    molecule.setNAtoms(6);
    molecule.alloc(); //allocates all memory for the class

    molecule.setIndex(0,HashAtom("C",0));
    molecule.setIndex(1,HashAtom("H",0));
    molecule.setIndex(2,HashAtom("O",0));
    molecule.setIndex(3,HashAtom("H",0));
    molecule.setIndex(4,HashAtom("H",0));
    molecule.setIndex(5,HashAtom("H",0));

    molecule.setCart(0,-1.013487,1.725956,1.257405);
    molecule.setCart(1,-0.069872,1.679306,0.755096);
    molecule.setCart(2,-1.488002,3.074931,1.256099);
    molecule.setCart(3,-1.168237,3.527873,2.039804);
    molecule.setCart(4,-1.718167,1.099499,0.751562);
    molecule.setCart(5,-0.897365,1.389691,2.266534);

*/

    molecule.convBohr();
    molecule.computeNucRep();
//    fileio.out << molecule.energyNuclei() << endl;

    molecule.computeRij();
    molecule.computeI();

    singleSlater.setRef(SingleSlater<dcomplex>::TCS); //TCS == GHF?
    singleSlater.setGuess(SingleSlater<dcomplex>::CORE);
    singleSlater.setSCFEneTol(1e-12);
    singleSlater.setNTCS(2);
    singleSlater.isClosedShell = false;
    singleSlater.doDIIS = false;

    basis.findBasisFile("CC-PVDZ");
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
    singleSlater.setPrintLevel(1);

//    fileio.out << "Allocate memory for aoints and singleSlater" << endl;
    aoints.alloc();
    singleSlater.alloc();

    singleSlater.formGuess();
    singleSlater.formFock();
    singleSlater.computeEnergy();
    singleSlater.SCF2();
    singleSlater.computeProperties();
    singleSlater.printProperties();

    finalizeCQ();

};
