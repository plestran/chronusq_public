// Integer Length Parameters
int LenScr        ; 
int LenSigma      ; 
int LenRho        ;
int LenXTSigma    ;
int LenXTRho      ;
int LenSuper      ;
int LenU          ;
int LenRes        ;
int LenTVec       ;
int LEN_LAPACK_SCR; 
int LWORK         ;
// Templated pointers for scratch
T * SCR        ; 
T * SigmaRMem  ; 
T * SigmaLMem  ; 
T * XTSigmaRMem; 
T * XTSigmaLMem; 
T * RhoRMem    ; 
T * RhoLMem    ; 
T * XTRhoRMem  ; 
T * XTRhoLMem  ; 
T * ASuperMem  ;
T * SSuperMem  ;
T * URMem      ; 
T * ULMem      ; 
T * ResRMem    ; 
T * ResLMem    ; 
T * TVecRMem   ;         
T * TVecLMem   ;         
T * LAPACK_SCR ;
/* Old member definitions of Eigen::Maps (not used)
TCMMap SigmaR   ;
TCMMap NewSR    ;
TCMMap XTSigmaR ;
TCMMap RhoR     ;
TCMMap NewRhoR  ;
TCMMap XTRhoR   ;
TCMMap UR       ;
TCMMap ResR     ;
TCMMap TrialVecR;
TCMMap NewVecR  ;

TCMMap SigmaL   ;
TCMMap NewSL    ;
TCMMap XTSigmaL ;
TCMMap RhoL     ;
TCMMap NewRhoL  ;
TCMMap XTRhoL   ;
TCMMap UL       ;
TCMMap ResL     ;
TCMMap TrialVecL;
TCMMap NewVecL  ;

TCMMap ASuper;
TCMMap SSuper;

TCMMap QR; 
TCMMap QL;
TCMMap RR;
TCMMap RL;

TVecMap ER;
TVecMap EI;
TCMMap  VR;
TCMMap  VL;
*/
/**DETERMINE LENGTH OF SCRATCH SPACE
 * DETERMINE LENGTH OF SCRATCH SPACE
 *
 * (1) Linear transform of A onto right / gerade trial vectors 
 *
 * (2) Reduced dimension of A onto the right / gerade trial vector subspace
 *     (also stores the reduced dimension eigenvectors after diagonalization
 *      via DSYEV)
 *
 * (3) Approximate right / gerade eigenvectors 
 *
 * (4) Residuals of the right / gerade eigenvectors 
 *
 * (5) Right / gerade trial vectors 
 *
 * (6) Temp storage for a single precondictioned eigenvector 
 *
 * (7) Linear transform of the metric S onto the right / gerade trial vectors
 *
 * (8) Reduced dimension of S into the right / gerade trial vector subspace
 *
 * (9) Linear transform of A onto left / ungerade trial vectors 
 *
 * (10) Linear transform of the metric S onto the left /un gerade trial vectors
 *
 * (11) Reduced dimension of A onto the left / ungerade trial vector subspace
 *
 * (12) Reduced dimension of S into the left / ungerade trial vector subspace
 *
 * (13) Approximate left / ungerade eigenvectors 
 *
 * (14) Residuals of the left / ungerade eigenvectors 
 *
 * (15) Left / ungerade trial vectors 
 *
 * (16) Supermatrix of reduced dimension A onto the two subspaces
 *
 * (17) Supermatrix of reduced dimension S onto the two subspaces
 *
 * (18) Local copy of the real part of the eigenvalues (reused for Tau storage 
 *      for QR)
 *
 * (19) Space for the paired / imaginary part of the eigenvalues
 *
 * (20) Length of LAPACK workspace (used in all LAPACK Calls)
 *
 * (21) Total double precision words required for LAPACK
 */
inline void initScrLen(){
  this->LenSigma   = this->N_ * this->maxSubSpace_;
  this->LenRho     = this->N_ * this->maxSubSpace_;
  this->LenXTSigma = this->maxSubSpace_ * this->maxSubSpace_;
  this->LenXTRho   = this->maxSubSpace_ * this->maxSubSpace_;
  this->LenSuper   = 4 * this->maxSubSpace_ * this->maxSubSpace_;
  this->LenU       = this->N_ * this->maxSubSpace_;
  this->LenRes     = this->N_ * this->maxSubSpace_;
  this->LenTVec    = this->N_ * this->maxSubSpace_;
  this->LWORK      = 0;
  this->LEN_LAPACK_SCR = 0;

  // Init LenScr
  this->LenScr     = 0;

  this->LenScr += this->LenSigma;   // 1
  this->LenScr += this->LenXTSigma; // 2
  this->LenScr += this->LenU;       // 3
  this->LenScr += this->LenRes;     // 4 
  this->LenScr += this->LenTVec;    // 5


  if(!this->isHermetian_ || this->symmetrizedTrial_){
    this->LenScr += this->LenRho;     // 7
    this->LenScr += this->LenXTRho;   // 8 
    this->LenScr += this->LenSigma;   // 9
    this->LenScr += this->LenRho;     // 10
    this->LenScr += this->LenXTSigma; // 11
    this->LenScr += this->LenXTRho;   // 12
    this->LenScr += this->LenU;       // 13
    this->LenScr += this->LenRes;     // 14
    this->LenScr += this->LenTVec;    // 15
    this->LenScr += this->LenSuper;   // 16
    this->LenScr += this->LenSuper;   // 17
  }

  this->LWORK          = 6*this->N_;
  this->LEN_LAPACK_SCR += this->maxSubSpace_;   // 18
  if(!this->isHermetian_ || this->symmetrizedTrial_)
    this->LEN_LAPACK_SCR += this->maxSubSpace_; // 19
  this->LEN_LAPACK_SCR += this->LWORK;          // 20
  this->LenScr += this->LEN_LAPACK_SCR;         // 21
}
inline void initScrPtr(){
  this->SCR         = NULL; 
  this->SigmaRMem   = NULL; 
  this->SigmaLMem   = NULL; 
  this->XTSigmaRMem = NULL; 
  this->XTSigmaLMem = NULL; 
  this->RhoRMem     = NULL; 
  this->RhoLMem     = NULL; 
  this->XTRhoRMem   = NULL; 
  this->XTRhoLMem   = NULL; 
  this->ASuperMem   = NULL;
  this->SSuperMem   = NULL;
  this->URMem       = NULL; 
  this->ULMem       = NULL; 
  this->ResRMem     = NULL; 
  this->ResLMem     = NULL; 
  this->TVecRMem    = NULL;         
  this->TVecLMem    = NULL;         
  this->LAPACK_SCR  = NULL;
}
inline void allocScr(){
  // Initialize various parameters
//this->initScrLen();
//this->initScrPtr();

  // Allocate scratch space
  this->SCR = new double [this->LenScr]; 

  // Partition scratch space
  this->SigmaRMem     = this->SCR;
  this->XTSigmaRMem   = this->SigmaRMem   + this->LenSigma;
  this->URMem         = this->XTSigmaRMem + this->LenXTSigma;
  this->ResRMem       = this->URMem       + this->LenU; 
  this->TVecRMem      = this->ResRMem     + this->LenRes;
  this->LAPACK_SCR    = this->TVecRMem    + this->LenTVec;
  if(!this->isHermetian_ || this->symmetrizedTrial_){
    this->RhoRMem       = this->LAPACK_SCR  + this->LEN_LAPACK_SCR;
    this->XTRhoRMem     = this->RhoRMem     + this->LenRho;
    this->SigmaLMem     = this->XTRhoRMem   + this->LenXTRho;
    this->XTSigmaLMem   = this->SigmaLMem   + this->LenSigma;
    this->RhoLMem       = this->XTSigmaLMem + this->LenXTSigma;
    this->XTRhoLMem     = this->RhoLMem     + this->LenRho;
    this->ULMem         = this->XTRhoLMem   + this->LenXTRho;
    this->ResLMem       = this->ULMem       + this->LenU;
    this->TVecLMem      = this->ResLMem     + this->LenRes;
    this->ASuperMem     = this->TVecLMem    + this->LenTVec;
    this->SSuperMem     = this->ASuperMem   + this->LenSuper;
  }

/*
  // Initialize Eigen Maps (not sure if this is nessacary...)
  new (&SigmaR)    TCMMap(SigmaRMem,  0,0);
  new (&NewSR)     TCMMap(SigmaRMem,  0,0);
  new (&XTSigmaR)  TCMMap(XTSigmaRMem,0,0);
  new (&RhoR)      TCMMap(RhoRMem,    0,0);
  new (&NewRhoR)   TCMMap(RhoRMem,    0,0);
  new (&XTRhoR)    TCMMap(XTRhoRMem,  0,0);
  new (&UR)        TCMMap(URMem,      0,0);
  new (&ResR)      TCMMap(ResRMem,    0,0);
  new (&TrialVecR) TCMMap(TVecRMem,   0,0);
  new (&NewVecR)   TCMMap(TVecRMem,   0,0);

  new (&SigmaL)    TCMMap(SigmaLMem,  0,0);
  new (&NewSL)     TCMMap(SigmaLMem,  0,0);
  new (&XTSigmaL)  TCMMap(XTSigmaLMem,0,0);
  new (&RhoL)      TCMMap(RhoLMem,    0,0);
  new (&NewRhoL)   TCMMap(RhoLMem,    0,0);
  new (&XTRhoL)    TCMMap(XTRhoLMem,  0,0);
  new (&UL)        TCMMap(ULMem,      0,0);
  new (&ResL)      TCMMap(ResLMem,    0,0);
  new (&TrialVecL) TCMMap(TVecLMem,   0,0);
  new (&NewVecL)   TCMMap(TVecLMem,   0,0);

  new (&ASuper) TCMMap(ASuperMem,0,0);
  new (&SSuper) TCMMap(SSuperMem,0,0);

  new (&QR) TCMMap(TVecRMem,0,0); 
  new (&QL) TCMMap(TVecLMem,0,0);
  new (&RR) TCMMap(ResRMem ,0,0);
  new (&RL) TCMMap(ResLMem ,0,0);


  new (&ER) TVecMap(LAPACK_SCR,0);
  new (&EI) TVecMap(LAPACK_SCR,0);
  new (&VR) TCMMap( LAPACK_SCR,0,0);
  new (&VL) TCMMap( LAPACK_SCR,0,0);
*/
}
void cleanupScr(){
  delete [] this->SCR;
}
