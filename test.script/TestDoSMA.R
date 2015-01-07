###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Dec 5, 2014
# FileName: TestDoSMA.R
###############################################################################

#--- Meng Cold ILs data ---#
read.geno.data(source="e://data//sampledata/SMA/Geno_Cold_PopC.csv", dp.code = "B", rp.code = "A", ht.code="H", na.code= NA, BCn=2, Fn=3);
geno.cold.c <- read.csv("e://data//sampledata/SMA/Geno_Cold_PopC.csv", header = T)
geno.cold.c <- read.geno.data("e://data//sampledata/SMA/Geno_Cold_PopC.csv", dp.code="B", rp.code = "A", ht.code = "H", na.code = NA, BCn = 2, Fn = 3);
geno.cold.c <- restrict.geno.data(geno.cold.c)
pheno.cold.c <- read.pheno.data("e://data//sampledata/SMA/Pheno_Cold_PopC.csv",type="MEAN", pop.type = "IL", geno="Lines", resp.var = c("SF_08", "SF_09_POOL"), na.code = NA)
pheno.cold.c <- restrict.pheno.data(pheno.cold.c)
debug(doSMA)
x <- doSMA(pheno.cold.c, geno.cold.c)
undebug(doSMA)

#--- Testing Multiple Environmental phenotypic data
pheno.cold.c <- read.csv("e://data//sampledata/SMA/Pheno_Cold_PopC.csv", header = T)
pheno.cold.c.e1 <- cbind(pheno.cold.c, Env = "E1");
pheno.cold.c.e2 <- cbind(pheno.cold.c, Env = "E2");
pheno.cold.c.e3 <- cbind(pheno.cold.c, Env = "E3");
pheno.cold.c.e <- rbind(pheno.cold.c.e1, pheno.cold.c.e2, pheno.cold.c.e3);
pheno.cold.c.e <- read.pheno.data(pheno.cold.c.e,type="MEAN", pop.type = "IL", env="Env", geno="Lines", resp.var = c("SF_08", "SF_09_POOL"), na.code = NA)
pheno.cold.c.e <- restrict.pheno.data(pheno.cold.c.e)
doSMA(pheno.cold.c.e, geno.cold.c, include.env = TRUE, is.EnvFixed = FALSE);
doSMA(pheno.cold.c.e, geno.cold.c, include.env = TRUE, is.EnvFixed = TRUE);
doSMA(pheno.cold.c.e, geno.cold.c, include.env = FALSE)
#--- Testing doSMA.MEANS function---#
debug(doSMA.MEANS)
x <- doSMA.MEANS(pheno.cold.c, geno.cold.c, "Lines", resp.var = c("SF_08", "SF_09_POOL"))
y <- doSMA.MEANS(pheno.cold.c.e, geno.cold.c, "Lines", resp.var = c("SF_08", "SF_09_POOL"), env = "Env", is.EnvFixed = TRUE)
z <- doSMA.MEANS(pheno.cold.c.e, geno.cold.c, "Lines", resp.var = c("SF_08", "SF_09_POOL"), env = "Env", is.EnvFixed = FALSE)

undebug(doSMA.MEANS)

#--- 40 lines ---#
lines <- paste("L", 1:40, sep = "");
#--- RCB experimental design---#
blocks <- paste("B", 1:3, sep = "");
#--- three environments---#
envs <- paste("E", 1:3, sep = "");
#--- two traits---#
set.seed(10);
hd.env1 <- rnorm(120, mean = 100, sd = 5);
set.seed(10);
hd.env2 <- rnorm(120, mean = 90, sd = 5);
set.seed(10);
hd.env3 <- rnorm(120, mean = 80, sd = 5);
set.seed(10)
ph.env1 <- sample(8:10, 120, replace = T);
set


#--- 100 markers ---#
markers <- paste("M", 1:5, sep = "");
set.seed(1);
genotypicData.sim <- sample(c(0,1,2), size = 1000, prob =c(0.8,0.05,0.15),replace = TRUE )
genotypicData.sim <- matrix(genotypicData.sim, ncol = 5, nrow = 200);
gene.code <- paste("Gen",1:200, sep = "");
genotypicData.sim <- cbind(gene.code, genotypicData.sim);
colnames(genotypicData.sim) <- c("Gname", markers)
genotypicData.sim <- as.data.frame(genotypicData.sim)
#--- testing doSMA on SingleEnvAnalysis ---#
genodata.sim <- read.geno.data(genotypicData.sim, dp.code = 2, rp.code = 0, ht.code = 1, na.code = NA, BCn = 2, Fn = 3)
doSMA(y,genodata.sim)

aF <- function(
		a,
		b = c("1", "2"),
		c,
		d = NULL
		)
		{
			if(missing(a))
				stop("\tError: a could not be null!\n");
			if(missing(b))
				stop("\tError: b could not be null!\n");
			if(missing(c))
				stop("\tError: c could not be null!\n");
			b <- match.arg(b);
			UseMethod("aF");
		}
aF.a <- function(
		a,
		b = c("1", "2"),
		c,
		d = NULL
)
{
	print("class a\n");
	print(paste("a= ", a, ", b=", b, ", c=", c, ",d=",d));
}
aF.b <- function(
		a,
		b = c("1", "2"),
		c,
		d = NULL
)
{
	print("class b\n");
	print(paste("a= ", a, ", b=", b, ", c=", c, ",d=",d));
}