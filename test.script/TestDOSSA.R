###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Nov 20, 2014
# FileName: TestDOSSA.R
###############################################################################


#--- using create demo data ---#
sssl.data.balanced.3env
x <- read.pheno.data(sssl.data.balanced.3env, 
		pop.type = "SSSL", 
		exptl.design= "RCB",
		resp.var = c("Yield", "TN"),
		geno = "Genotype",
		block = "Block",
		env = "Env"
)
y <- restrict.pheno.data(x)
z <- doSSA(y)