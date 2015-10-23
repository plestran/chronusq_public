import libpythonapi as chronusQ
#
# Map from input paramater for RT.ortho
# to the RT::ORTHO enum
#
orthoMap = { 
  'LOWDIN'   :chronusQ.RealTime_ORTHO.Lowdin    ,
  'CHOLESKY' :chronusQ.RealTime_ORTHO.Cholesky  ,
  'CANONICAL':chronusQ.RealTime_ORTHO.Canonical 
}

#
# Map from input paramater for RT.uprop
# to the RT::FORM_U enum
#
formUMap = { 
  'EIGENDECOMP':chronusQ.RealTime_FORM_U.EigenDecomp ,
  'TAYLOR'     :chronusQ.RealTime_FORM_U.Taylor      
}

#
# Map from input paramater for RT.envelope
# to the RT::ENVELOPE enum
#
envMap   = {       
  'PW'      :chronusQ.RealTime_ENVELOPE.Constant ,
  'LINRAMP' :chronusQ.RealTime_ENVELOPE.LinRamp  ,
  'GAUSSIAN':chronusQ.RealTime_ENVELOPE.Gaussian ,
  'STEP'    :chronusQ.RealTime_ENVELOPE.Step     ,
  'SINSQ'   :chronusQ.RealTime_ENVELOPE.SinSq    
}

#
# Map from input paramater for PSCF.METHOD
# to the SDResponse::METHOD enum
#
sdrMethodMap = {
  "INVALID":chronusQ.SDResponse_METHOD.INVALID,  
  "CIS"    :chronusQ.SDResponse_METHOD.CIS    ,
  "RPA"    :chronusQ.SDResponse_METHOD.RPA    ,
  "PPRPA"  :chronusQ.SDResponse_METHOD.PPRPA  ,
  "PPATDA" :chronusQ.SDResponse_METHOD.PPATDA ,
  "PPCTDA" :chronusQ.SDResponse_METHOD.PPCTDA ,
  "STAB"   :chronusQ.SDResponse_METHOD.STAB     
}

aointAlg = {
  "DIRECT":chronusQ.AOIntegrals_INTEGRAL_ALGORITHM.DIRECT,
  "INCORE":chronusQ.AOIntegrals_INTEGRAL_ALGORITHM.INCORE,
  "DENFIT":chronusQ.AOIntegrals_INTEGRAL_ALGORITHM.DENFIT
}
