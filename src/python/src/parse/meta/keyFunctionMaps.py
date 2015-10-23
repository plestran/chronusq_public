def genkeyFuncMap(workers):
  keyFuncMap = {}
  
  keyFuncMap['RT'] = {
    'MAXSTEP'  :workers['CQRealTime'].setMaxSteps ,
    'TIMESTEP' :workers['CQRealTime'].setStepSize ,
    'EDFIELD'  :workers['CQRealTime'].setFieldAmp ,
    'TIME_ON'  :workers['CQRealTime'].setTOn      ,
    'TIME_OFF' :workers['CQRealTime'].setTOff     ,
    'FREQUENCY':workers['CQRealTime'].setFreq     ,
    'PHASE'    :workers['CQRealTime'].setPhase    ,
    'SIGMA'    :workers['CQRealTime'].setSigma    ,
    'ENVELOPE' :workers['CQRealTime'].setEnvelope ,
    'ORTHO'    :workers['CQRealTime'].setOrthoTyp ,
    'INIDEN'   :workers['CQRealTime'].setInitDen  ,
    'UPROP'    :workers['CQRealTime'].setFormU    
  }

