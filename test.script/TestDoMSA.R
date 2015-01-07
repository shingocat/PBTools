###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Nov 24, 2014
# FileName: TestDoMSA.R
###############################################################################

sssl.data <- read.pheno.data(sssl.data.balanced.3env, 
		pop.type = "SSSL", 
		exptl.design= "RCB",
		resp.var = c("Yield", "TN"),
		geno = "Genotype",
		block = "Block",
		env = "Env"
)
sssl.data.restricted <- restrict.pheno.data(sssl.data)

py.data <- read.pheno.data(phenodata="e://data//sampledata//bigenes.csv", pop.type="PL", gene.num = 2, exptl.design="RCB", resp.var=c("Heading", "PH"), geno="Genotype", block = "Block", env="Env")
py.data.restricted <- restrict.pheno.data(py.data)


debug(doMSA)
sssl.data.restricted <- doSSA(sssl.data.restricted)
y <- doMSA(sssl.data.restricted, is.EnvFixed = T)
x <- doMSA(py.data.restricted, is.EnvFixed = F)
undebug(doMSA)