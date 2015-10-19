import libpythonapi as chronusQ
orthoMap = { 
  'LOWDIN'   :chronusQ.RealTime_ORTHO.Lowdin    ,
  'CHOLESKY' :chronusQ.RealTime_ORTHO.Cholesky  ,
  'CANONICAL':chronusQ.RealTime_ORTHO.Canonical 
}

formUMap = { 
  'EIGENDECOMP':chronusQ.RealTime_FORM_U.EigenDecomp ,
  'TAYLOR'     :chronusQ.RealTime_FORM_U.Taylor      
}

envMap   = {       
  'PW'      :chronusQ.RealTime_ENVELOPE.Constant ,
  'LINRAMP' :chronusQ.RealTime_ENVELOPE.LinRamp  ,
  'GAUSSIAN':chronusQ.RealTime_ENVELOPE.Gaussian ,
  'STEP'    :chronusQ.RealTime_ENVELOPE.Step     ,
  'SINSQ'   :chronusQ.RealTime_ENVELOPE.SinSq    
}




