###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Nov 19, 2014
# FileName: TestSSSLSSA.R
###############################################################################


debug(sssl.ssa.test)
sssl.ssa.test(exptl.design = "RCB",
	respvar = c("Yield", "TN"),
	geno="Genotype",
	block="Block",
	env="Env",
	data= sssl.data.balanced.3env
)
undebug(sssl.ssa.test)
