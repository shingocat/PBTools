#-------------------------------------------------
# This function computes several values first before calling ammi.analysis function
#
# ARGUMENTS:
# ENV - vector containing values of environment factor
# GEN - vector containing values of genotype factor
# REP - vector containing values of rep factor
# Y - vector containing the values of the response variable
# SIGMA2 - vector containing values of residual variance
# number - TRUE if levels of genotype and environment factors are numeric
# biplotPC12 - TRUE if biplot (PC1 vs PC2) will be created
# biplotPC13 - TRUE if biplot (PC1 vs PC3) will be created
# biplotPC23 - TRUE if biplot (PC2 vs PC3) will be created
# ammi1 - TRUE if ammi1 biplot will be created
# adaptMap - TRUE if adaptation map will be created
# respVar - name of response variable
#
# Author: Alaine A. Gulles, Nellwyn L. Sales
#-------------------------------------------------

ammi.analysis2 <- function (ENV, GEN, REP, Y, SIGMA2, number = TRUE, biplotPC12 = FALSE, biplotPC13 = FALSE, biplotPC23 = FALSE, ammi1 = FALSE, adaptMap = FALSE, respVar) UseMethod("ammi.analysis2")
  
ammi.analysis2.default <- function (ENV, GEN, REP, Y, SIGMA2, number = TRUE, biplotPC12 = FALSE, biplotPC13 = FALSE, biplotPC23 = FALSE, ammi1 = FALSE, adaptMap = FALSE, respVar) {

  numrepEnv <- tapply(REP, ENV, mean)
  sigma2Env <- tapply(SIGMA2, ENV, mean)
  MSE <- sum(numrepEnv*sigma2Env)/sum(numrepEnv)
  NREP <- 1/mean(1/numrepEnv)

  ammi12Out <- ammi.analysis(ENV, GEN, NREP, Y, MSE, number, biplotPC12 = biplotPC12, biplotPC13 = biplotPC13, biplotPC23 = biplotPC23, ammi1 = ammi1, adaptMap = adaptMap, yVar = respVar)

}

