###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Jan 3, 2015
# FileName: TestDoMMA.R
###############################################################################


#--- Using xiao long TQ-IL population and snp data---#
geno.tq <- read.geno.data("e://data//xiaolong//TQ_IL_Geno.csv", dp.code = 2, rp.code = 0, ht.code= 1, na.code = -1, BCn = 2, Fn = 3);
geno.tq.restricted <- restrict.geno.data(geno.tq);
pheno.tq <- read.pheno.data("e://data//xiaolong//TQ_IL_Pheno.csv", geno = "L-Code", type="MEAN", pop.type="IL", resp.var = c("PH", "HD", "PN", "PL"), na.code = "." );

#--- for object PhenotypicData ---#
doMMA(pheno.tq, geno.tq.restricted,method="RIDGE_REGRESSION")
doMMA(pheno.tq, geno.tq.restricted,method="LASSO")
doMMA(pheno.tq, geno.tq.restricted,method="ELASTIC_NET")

#--- for object SingleEnvAnalysis ---#
pheno <- read.pheno.data("e://data//Small01.csv", geno = "Gname", type="RAW", pop.type="IL", block="Rep", resp.var = c("DFL", "PH", "TN"), exptl.design = "RCB", env= "Season", na.code = "." );
pheno.se <- doSEA(pheno)
markers <- paste("M", 1:5, sep = "");
set.seed(1);
genotypicData.sim <- sample(c(0,1,2), size = 1000, prob =c(0.8,0.05,0.15),replace = TRUE )
genotypicData.sim <- matrix(genotypicData.sim, ncol = 5, nrow = 200);
gene.code <- paste("Gen",1:200, sep = "");
genotypicData.sim <- cbind(gene.code, genotypicData.sim);
colnames(genotypicData.sim) <- c("Gname", markers)
genotypicData.sim <- as.data.frame(genotypicData.sim)
genodata.sim <- read.geno.data(genotypicData.sim, dp.code = 2, rp.code = 0, ht.code = 1, na.code = NA, BCn = 2, Fn = 3)

doMMA(pheno.se, genodata.sim);
doMMA(pheno.se, genodata.sim, method = "RIDGE_REGRESSION")
doMMA(pheno.se, genodata.sim, method = "ELASTIC_NET")


#--- for object MultiEnvAnalysis ---#
pheno.me <- doMEA(pheno)
doMMA(pheno.me, genodata.sim);
doMMA(pheno.me, genodata.sim, method = "RIDGE_REGRESSION")
doMMA(pheno.me, genodata.sim, method = "ELASTIC_NET")





