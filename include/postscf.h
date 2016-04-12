class PostSCF {
  protected:
  // Single Dimensions
  int       nOA_; ///< Number of occupied orbitals (Alpha)
  int       nVA_; ///< Number of virtual orbitals (Alpha)
  int       nOB_; ///< Number of occupied orbitals (Beta)
  int       nVB_; ///< Number of virtual orbitals (Beta)
  int       nO_;  ///< Total number of occupied orbitals NOA + NOB
  int       nV_;  ///< Total number of virtual orbitals NVA + NVB

  // Quadratic Dimenstions
    
  // Occupied-Occupied
  int       nOAOA_;      ///< NOA * NOA
  int       nOBOB_;      ///< NOB * NOB
  int       nOAOB_;      ///< NOA * NOB
  int       nOO_;        ///< NO * NO
  int       nOAOA_SLT_;  ///< NOA * (NOA - 1) / 2
  int       nOBOB_SLT_;  ///< NOB * (NOB - 1) / 2
  int       nOO_SLT_;    ///< NO * (NO - 1) / 2
  int       nOAOA_LT_;   ///< NOA * (NOA + 1) / 2
  int       nOO_LT_;     ///< NO * (NO + 1) / 2

  // Occupied-Virtual
  int       nOAVA_;      ///< NOA * NVA
  int       nOBVB_;      ///< NOB * NVB
  int       nOAVB_;      ///< NOA * NVB
  int       nOBVA_;      ///< NOB * NVA
  int       nOV_;        ///< NO * NV

  // Virtual-Virtual
  int       nVAVA_;      ///< NVA * NVA
  int       nVBVB_;      ///< NVB * NVB
  int       nVAVB_;      ///< NVA * NVB
  int       nVV_;        ///< NV * NV
  int       nVAVA_SLT_;  ///< NVA * (NVA - 1) / 2
  int       nVBVB_SLT_;  ///< NVB * (NVB - 1) / 2
  int       nVV_SLT_;    ///< NV * (NV - 1) / 2
  int       nVAVA_LT_;   ///< NVA * (NVA + 1) / 2
  int       nVV_LT_;     ///< NV * (NV + 1) / 2

  public:
    PostSCF() {
      // Zero out PSCF dimensions to be built later
      this->nOA_        = 0;
      this->nVA_        = 0;
      this->nOB_        = 0;
      this->nVB_        = 0;
      this->nOAVA_      = 0;      
      this->nOBVB_      = 0;      
      this->nOAVB_      = 0;      
      this->nOBVA_      = 0;      
      this->nVAVA_SLT_  = 0;  
      this->nVBVB_SLT_  = 0;  
      this->nVAVA_LT_   = 0;   
      this->nVAVA_      = 0;      
      this->nVBVB_      = 0;      
      this->nOAOA_SLT_  = 0;  
      this->nOBOB_SLT_  = 0;  
      this->nOAOA_LT_   = 0;   
      this->nOAOA_      = 0;      
      this->nOBOB_      = 0;      
      this->nVAVB_      = 0;      
      this->nOAOB_      = 0;      
      this->nO_         = 0;         
      this->nV_         = 0;         
      this->nOV_        = 0;        
      this->nVV_SLT_    = 0;    
      this->nVV_LT_     = 0;     
      this->nVV_        = 0;        
      this->nOO_SLT_    = 0;    
      this->nOO_LT_     = 0;     
      this->nOO_        = 0;        
    };

    virtual void initPSCFDims() = 0;
    int nOA()         { return this->nOA_;} 
    int nVA()         { return this->nVA_;} 
    int nOB()         { return this->nOB_;} 
    int nVB()         { return this->nVB_;} 
    int nO()          { return this->nO_ ;} 
    int nV()          { return this->nV_ ;} 
                                                 
    int nOAOA()       { return this->nOAOA_     ;} 
    int nOBOB()       { return this->nOBOB_     ;} 
    int nOAOB()       { return this->nOAOB_     ;} 
    int nOO()         { return this->nOO_       ;} 
    int nOAOA_SLT()   { return this->nOAOA_SLT_ ;} 
    int nOBOB_SLT()   { return this->nOBOB_SLT_ ;} 
    int nOO_SLT()     { return this->nOO_SLT_   ;} 
    int nOAOA_LT()    { return this->nOAOA_LT_  ;} 
    int nOO_LT()      { return this->nOO_LT_    ;} 
                                                 
    int nOAVA()       { return this->nOAVA_;}      
    int nOBVB()       { return this->nOBVB_;}      
    int nOAVB()       { return this->nOAVB_;}      
    int nOBVA()       { return this->nOBVA_;}      
    int nOV()         { return this->nOV_  ;}      
                                                 
    int nVAVA()       { return this->nVAVA_     ;} 
    int nVBVB()       { return this->nVBVB_     ;} 
    int nVAVB()       { return this->nVAVB_     ;} 
    int nVV()         { return this->nVV_       ;} 
    int nVAVA_SLT()   { return this->nVAVA_SLT_ ;} 
    int nVBVB_SLT()   { return this->nVBVB_SLT_ ;} 
    int nVV_SLT()     { return this->nVV_SLT_   ;} 
    int nVAVA_LT()    { return this->nVAVA_LT_  ;} 
    int nVV_LT()      { return this->nVV_LT_    ;} 
};
