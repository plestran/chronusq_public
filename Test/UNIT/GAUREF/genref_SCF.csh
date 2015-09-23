#!/bin/tcsh

@ i = 1
foreach x ( h2o_sto3g.log h2o_631G.log h2o_ccpVDZ.log benzene_sto3g.log benzene_6-31G.log benzene_cc-pVDZ.log Li_sto3g.log Li_6-31G.log Li_cc-pVDZ.log O2_sto3g.log O2_6-31G.log O2_cc-pVDZ.log Li_sto3g_ROHF.log Li_6-31G_ROHF.log Li_cc-pVDZ_ROHF.log O2_sto3g_ROHF.log O2_6-31G_ROHF.log O2_cc-pVDZ_ROHF.log MnO2_3-21G.log MnO2_6-31G.log MnO2_cc-pVDZ.log MnO2_3-21G_ROHF.log MnO2_6-31G_ROHF.log MnO2_cc-pVDZ_ROHF.log)
	echo 'test'$i
	python grabinfo_SCF.py $x
	@ i ++
end

