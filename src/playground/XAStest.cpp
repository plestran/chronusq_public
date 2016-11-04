#include <global.h>
#include <workers.h>
#include <pythonapi.h>

using namespace ChronusQ;


int main(int argc,char **argv) {
    CQMemManager memManager;
    Molecule molecule;
    BasisSet basis;
    AOIntegrals aoints;
    SingleSlater<dcomplex> singleSlater;
    FileIO fileio("XAStest.inp","XAStest.out","MnAcAc.bin");

    memManager.setTotalMem(256e7);
    initCQ(argc,argv);
    CQSetNumThreads(1); //Sets up open MP threads


// Molecule Specification: Mn(acac)2
    molecule.setCharge(0);
    molecule.setNTotalE(87);
    molecule.setMultip(2);
    molecule.setNAtoms(15);
    molecule.alloc();

    molecule.setIndex(0,HashAtom("C",0));
    molecule.setIndex(1,HashAtom("O",0));
    molecule.setIndex(2,HashAtom("H",0));
    molecule.setIndex(3,HashAtom("C",0));
    molecule.setIndex(4,HashAtom("O",0));
    molecule.setIndex(5,HashAtom("H",0));
    molecule.setIndex(6,HashAtom("H",0));
    molecule.setIndex(7,HashAtom("C",0));
    molecule.setIndex(8,HashAtom("O",0));
    molecule.setIndex(9,HashAtom("H",0));
    molecule.setIndex(10,HashAtom("C",0));
    molecule.setIndex(11,HashAtom("O",0));
    molecule.setIndex(12,HashAtom("H",0));
    molecule.setIndex(13,HashAtom("H",0));
    molecule.setIndex(14,HashAtom("Mn",0));

    molecule.setCart(0,  0.000000,  1.230231, 1.676654);
    molecule.setCart(1,  0.000000,  0.836540, 0.491685);
    molecule.setCart(2,  0.000000,  2.268652, 1.888047);
    molecule.setCart(3, -0.000000, -1.230231, 1.676654);
    molecule.setCart(4, -0.000000, -0.836540, 0.491685);
    molecule.setCart(5, -0.000000, -2.268652, 1.888047);
    molecule.setCart(6,  0.000000,  0.000000, 2.514988);
    molecule.setCart(7,  1.230231, -0.000000,-1.676654);
    molecule.setCart(8,  0.836540, -0.000000,-0.491685);
    molecule.setCart(9,  2.268652, -0.000000,-1.888047);
    molecule.setCart(10,-1.230231,  0.000000,-1.676654);
    molecule.setCart(11,-0.836540,  0.000000,-0.491685);
    molecule.setCart(12,-2.268652,  0.000000,-1.888047);
    molecule.setCart(13, 0.000000,  0.000000,-2.514988);
    molecule.setCart(14, 0.000000,  0.000000, 0.000000);
// --------------------------------------

    molecule.convBohr();
    molecule.computeNucRep();
    molecule.computeRij();
    molecule.computeI();

    singleSlater.setRef("X2C");
    singleSlater.isClosedShell = false;
    singleSlater.doDIIS = true;
    singleSlater.doDamp = true;
    singleSlater.dampParam = 0.2;
    singleSlater.setGuess(READ);
    singleSlater.isDFT = false;
    singleSlater.isHF = true;

    fileio.doRestart = true;
    fileio.iniH5Files();
    
    basis.findBasisFile("6-31G");
    basis.communicate(fileio);  // This function passes fileio reference 
    basis.parseGlobal(); //Reads entire basis set file into memory
    basis.constructLocal(&molecule);
    basis.makeMaps(&molecule); 
    basis.renormShells();

    aoints.communicate(molecule,basis,fileio,memManager);
    singleSlater.communicate(molecule,basis,aoints,fileio,memManager);

    aoints.setPrintLevel(1);
    aoints.initMeta();
    aoints.integralAlgorithm = AOIntegrals::INCORE;
    
    singleSlater.initMeta();
    singleSlater.setPrintLevel(1);

    aoints.alloc();
    singleSlater.alloc();

    singleSlater.formGuess();
    singleSlater.computeProperties();
    singleSlater.printProperties();
    singleSlater.SCF3();
    singleSlater.computeProperties();
    singleSlater.printProperties();

    finalizeCQ();
    return 0;
};
