#!/bin/tcsh

@ i = 25
foreach x ( h2o_sto3g_CIS.log h2o_631G_CIS.log h2o_ccpVDZ_CIS.log benzene_sto3g_CIS.log benzene_6-31G_CIS.log benzene_cc-pVDZ_CIS.log Li_sto3g_CIS.log Li_6-31G_CIS.log Li_cc-pVDZ_CIS.log O2_sto3g_CIS.log O2_6-31G_CIS.log O2_cc-pVDZ_CIS.log MnO2_3-21G_CIS.log MnO2_6-31G_CIS.log MnO2_cc-pVDZ_CIS.log)
	echo 'test'$i
	python grabinfo_RESP.py $x
	@ i ++
end

foreach x ( h2o_sto3g_TD.log h2o_631G_TD.log h2o_ccpVDZ_TD.log benzene_sto3g_TD.log benzene_6-31G_TD.log benzene_cc-pVDZ_TD.log Li_sto3g_TD.log Li_6-31G_TD.log Li_cc-pVDZ_TD.log O2_sto3g_TD.log O2_6-31G_TD.log O2_cc-pVDZ_TD.log MnO2_3-21G_TD.log MnO2_6-31G_TD.log MnO2_cc-pVDZ_TD.log)
	echo 'test'$i
	python grabinfo_RESP.py $x
	@ i ++
end
