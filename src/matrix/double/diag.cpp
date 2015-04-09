
using ChronusQ::Matrix;
char UPLO = 'L';
int LWORK;
double *CPY,*WORK
// Check if it's already been diagonalized
// Check if it's square

this->cleanEigen();
this->allocEigen();
if(this->symm_=='G') CPY=this->allocScratch();
LWORK = this->eigWORK(1,UPLO);
WORK = new (nothrow) double *[LWORK];
if(WORK==NULL) throw 3105;

if(this->symm_=='G'){
  dgeev_(&this->JOBVL_,this->JOBVR_,&this->rows_,CPY,&this->rows_,
    this->eigenvalue_re_,this->eigenvalue_im,this->eigenvector_l_;&this->rows_,
    this->eigenvector_r_,&this->rows_,WORK,&LWORK,&INFO);}
}
else if(this->symm_'S'){
  dsyev_(&this->JOBVR_,&UPLO,&this->rows_,this->eigenvector_,&this->rows_,
    this->eigenvalue_,WORK,&LWORK,&INFO);
};

delete[] CPY; delete[] WORK;
this->haveEigen_=1;
}
