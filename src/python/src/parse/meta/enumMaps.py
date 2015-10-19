import libpythonapi as chronusQ
orthoMap = { 
  'lowdin'   :chronusQ.RealTime_ORTHO.Lowdin    ,
  'cholesky' :chronusQ.RealTime_ORTHO.Cholesky  ,
  'canonical':chronusQ.RealTime_ORTHO.Canonical 
}

formUMap = { 
  'eigendecomp':chronusQ.RealTime_FORM_U.EigenDecomp ,
  'taylor'     :chronusQ.RealTime_FORM_U.Taylor      
}

envMap   = {       
  'pw'      :chronusQ.RealTime_ENVELOPE.Constant ,
  'linramp' :chronusQ.RealTime_ENVELOPE.LinRamp  ,
  'gaussian':chronusQ.RealTime_ENVELOPE.Gaussian ,
  'step'    :chronusQ.RealTime_ENVELOPE.Step     ,
  'sinsq'   :chronusQ.RealTime_ENVELOPE.SinSq    
}




