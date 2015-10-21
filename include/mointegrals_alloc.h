template<typename T>
void MOIntegrals<T>::iniMOIntegrals(Molecule * molecule, BasisSet * basisSet, 
  FileIO * fileio, Controls * controls, AOIntegrals * aointegrals, 
  SingleSlater<T> * singleSlater){

  this->communicate(*molecule,*basisSet,*fileio,*controls,*aointegrals,
    *singleSlater);
  this->initMeta();

}

