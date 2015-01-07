###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Dec 11, 2014
# FileName: TestDoMMAonMeans.R
###############################################################################

set.seed(1);
var <- var = rnorm(50, mean = 80, sd = 3)
genotypicdata <- matrix(sample(c(0,1,2), size = 100 * 50, replace = TRUE, prob=getExpectGenotypicFreq()), ncol = 100);
#--- reading data from xiao meng cold ILs C population---#
genotypicData <- read.geno.data("e://data//sampledata//sma//Col-.txt")

#--- Testing getElasticNetAlpha function---#
debug(getElasticNetAlpha);
alpha <- getElasticNetAlpha(genotypicdata, var)
undebug(getElasticNetAlpha);
#--- Testing doMMAonMeans function ---#
debug(doMMAonMeans);
x <- doMMAonMeans(var, genotypicdata)


#--- Using xiao long TQ-IL population and snp data---#
geno.tq <- read.geno.data("e://data//xiaolong//TQ_IL_Geno.csv", dp.code = 2, rp.code = 0, ht.code= 1, na.code = -1, BCn = 2, Fn = 3);
geno.tq.restricted <- restrict.geno.data(geno.tq);
pheno.tq <- read.csv("e://data//xiaolong//TQ_IL_Pheno.csv", header = T, encoding = "UTF-8", na.strings = ".");

debug(doMMAonMeans);
doMMAonMeans(pheno.tq, geno.tq.restricted, geno = "L.Code", resp.var = c("PH","HD"), method = "LASSO")
#Error in { : 
#			task 1 failed - "FUNLM not reporting p-values; possible overfit model
#			See help pages for more information"
doMMAonMeans(pheno.tq, geno.tq.restricted, geno = "L.Code", resp.var = c("PH","HD"), method = "RIDGE_REGRESSION")
#Error in { : 
#			task 1 failed - "FUNLM not reporting p-values; possible overfit model
#			See help pages for more information"
doMMAonMeans(pheno.tq, geno.tq.restricted, geno = "L.Code", resp.var = c("PH","HD"), method = "ELASTIC_NET")
#Error in elnet(x, is.sparse, ix, jx, y, weights, offset, type.gaussian,  : 
#				NA/NaN/Inf in foreign function call (arg 5)
undebug(doMMAonMeans);

#--- using xiao meng C population IL and SSR---#
geno.c <- read.geno.data("e://data//sampledata/sma//Geno_Cold_PopC.csv", dp.code = "B", rp.code = "A", ht.code ="H", na.code = "C");
geno.c.restricted <- restrict.geno.data(geno.c)
pheno.c <- read.csv("e://data//sampledata//sma//Pheno_Cold_PopC.csv")
doMMAonMeans(pheno.c, geno.c.restricted, geno = "Lines", resp.var =c("SF_08", "SF_09_POOL"), method = "LASSO" )


genetic.value <- geno.c.restricted$restricted$data
genetic.value <- genetic.value[, -1]; #removing the first column lines info.
genetic.value <- as.matrix(genetic.value)
genetic.value <- apply(genetic.value, 2, as.character);
genetic.value <- apply(genetic.value, 2, as.numeric);
genetic.value.subset <- genetic.value[ ,c(1:7,9:11)]
trait.value <- pheno.c[,3]
hdglm(trait.value ~ genetic.value.subset, alpha=1)

hdlm(trait.value ~ genetic.value.subset, alpha=0.0001)
glmnet(genetic.value.subset, trait.value)

#--- do imputation in GenotypicData object
#--- replacing missing value to major gene
set.seed(1)
x <- sample(c(0,1,2,NA), replace = TRUE, size = 100, prob=c(0.8,0.05,0.14,0.01));
x <- factor(x);
n <- length(x);
y <- table(x, useNA ="always");
prob <- y / n;
major.gene <- names(which(max(prob) == prob))
x <- replace.by.major.gene(x);

