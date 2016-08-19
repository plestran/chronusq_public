#
# The Chronus Quantum (ChronusQ) software package is high-performace 
# computational chemistry software with a strong emphasis on explicitly 
# time-dependent and post-SCF quantum mechanical methods.
# 
# Copyright (C) 2014-2016 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
# 
#
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
  'SINSQ'   :chronusQ.RealTime_ENVELOPE.SinSq    ,
  'ELLIPTIC':chronusQ.RealTime_ENVELOPE.Elliptic    
}

ellPolMap = {
  'LXY' :chronusQ.RealTime_ELL_POL.LXY,
  'LXZ' :chronusQ.RealTime_ELL_POL.LXZ,
  'LYZ' :chronusQ.RealTime_ELL_POL.LYZ,
  'RXY' :chronusQ.RealTime_ELL_POL.RXY,
  'RXZ' :chronusQ.RealTime_ELL_POL.RXZ,
  'RYZ' :chronusQ.RealTime_ELL_POL.RYZ,
}

#
# Map from input paramater for PSCF.METHOD
# to the SDResponse::METHOD enum
#
#sdrMethodMap = {
#  "INVALID":chronusQ.SDResponse_METHOD.INVALID,  
#  "CIS"    :chronusQ.SDResponse_METHOD.CIS    ,
#  "RPA"    :chronusQ.SDResponse_METHOD.RPA    ,
#  "PPRPA"  :chronusQ.SDResponse_METHOD.PPRPA  ,
#  "PPATDA" :chronusQ.SDResponse_METHOD.PPATDA ,
#  "PPCTDA" :chronusQ.SDResponse_METHOD.PPCTDA ,
#  "STAB"   :chronusQ.SDResponse_METHOD.STAB     
#}

aointAlg = {
  "DIRECT":chronusQ.AOIntegrals_INTEGRAL_ALGORITHM.DIRECT,
  "INCORE":chronusQ.AOIntegrals_INTEGRAL_ALGORITHM.INCORE,
  "DENFIT":chronusQ.AOIntegrals_INTEGRAL_ALGORITHM.DENFIT
}

guessMap = {
  "SAD"   :chronusQ.Guess.SAD    ,
  "CORE"  :chronusQ.Guess.CORE   ,
  "READ"  :chronusQ.Guess.READ   , 
  "RANDOM":chronusQ.Guess.RANDOM  
}

#exchMap = {
#  "NOEXCH":chronusQ.EXCH.NOEXCH,
#  "EXACT" :chronusQ.EXCH.EXACT ,
#  "SLATER":chronusQ.EXCH.SLATER,
#  "B88"   :chronusQ.EXCH.B88
#}
#
#corrMap = {
#  "NOCORR":chronusQ.CORR.NOCORR,
#  "VWN3"  :chronusQ.CORR.VWN3,
#  "VWN5"  :chronusQ.CORR.VWN5,
#  "LYP"   :chronusQ.CORR.LYP
#}

kernelMap = {
  "NODFT"      :chronusQ.DFT.NODFT,
  "USERDEFINED":chronusQ.DFT.USERDEFINED,
  "LSDA"       :chronusQ.DFT.LSDA
}

gridMap = {
  "EULERMACL" : chronusQ.GRID_TYPE.EULERMAC,
  "GAUSSCHEB" : chronusQ.GRID_TYPE.GAUSSCHEBFST
}

dftWeightScheme = {
 "BECKE"  : chronusQ.ATOMIC_PARTITION.BECKE,
 "FRISCH" : chronusQ.ATOMIC_PARTITION.FRISCH
}
