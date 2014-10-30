###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Oct 22, 2014
# FileName: hdlm.packages.R
###############################################################################


library(hdlm)
set.seed(1)
x <- matrix(rnorm(100*40), ncol = 100)
y <- x[,1] + x[,2] * 0.5 + rnorm(40, sd=0.1)
out <- hdlm(y ~ x, , bootstrap = 2, pval.method="fdr")

set.seed(42)
x <- matrix(rnorm(10*100),ncol=10)
mu <- exp(x[,1] + x[,2]*0.5) / (1 + exp(x[,1] + x[,2]*0.5))

y <- rbinom(100,1,prob=mu)

out <- hdglm(y ~ x, family='binomial')
summary(out)