#-------------------------------------------------
# This function performs gge analysis for multi-env two-stage. It computes several values first before calling gge.analysis function.
#
# PARAMETERS: 
# ENV - vector containing values of environment factor
# GEN - vector containing values of genotype factor
# REP - vector containing values of rep factor
# Y - vector containing the values of the response variable
# SIGMA2 - vector containing values of residual variance
# number - TRUE if levels of genotype and environment factors are numeric
# graphSym - TRUE if GGE biplot in symmetric view will be created
# graphEnv - TRUE if GGE biplot in environment view will be created
# graphGeno - TRUE if GGE biplot in genotype view will be created
# respVar - name of response variable
# f - 0.5 if symmetric view; 0 if environment view; 1 if genotype view
#
# Author: Alaine Gulles and Nellwyn Sales
#-------------------------------------------------

gge.analysis2 <- function (ENV, GEN, REP, Y, SIGMA2, number = TRUE, graphSym=FALSE, graphEnv=FALSE, graphGeno=FALSE, respVar, f=0.5) UseMethod("gge.analysis2")

gge.analysis2.default <- function (ENV, GEN, REP, Y, SIGMA2, number = TRUE, graphSym=FALSE, graphEnv=FALSE, graphGeno=FALSE, respVar,f=0.5) {
  
  numrepEnv <- tapply(REP, ENV, mean)
  sigma2Env <- tapply(SIGMA2, ENV, mean)
  MSE <- sum(numrepEnv*sigma2Env)/sum(numrepEnv)
  NREP <- 1/mean(1/numrepEnv)
  
  gge2Out <- gge.analysis(ENV, GEN, NREP, Y, MSE, number, graphSym=graphSym, graphEnv=graphEnv, graphGeno=graphGeno, yVar = respVar, f)
  
}