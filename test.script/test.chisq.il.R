###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Sep 25, 2014
# FileName: test.chisq.il.R
###############################################################################


str(summary(geno.data.bc5))
aggregate(geno.data.bc5[1,],by=list((geno.data.bc5[1,])), FUN = length);
do.call(what=aggregate,args = list(x=geno.data.bc5[1,], by= list(geno.data.bc5[1,]), FUN = length))
outcome = chisq.test(c(96,4), p = getExpectGenotypicFreq(5,0));
colnames(geno.data.bc5)
M1 <- geno.data.bc5[,1];
M1 <- table(M1);
outcome.m1 <- chisq.test(M1,p = getExpectGenotypicFreq(5,0) )

M1.bcf.data <- geno.data.bc5f2[,1];
M1.bcf.data <- factor(M1.bcf.data);
M1.bcf <- table(M1.bcf.data);

outcome.m1.bcf <- chisq.test(M1.bcf, p = getExpectGenotypicFreq(5,2,2));
adjust.p <- getExpectGenotypicFreq(5,2,2)[c(3,1)]/sum(getExpectGenotypicFreq(5,2,2)[c(3,1)]);
test.adjust.p <- chisq.test(M1.bcf[c(1,3)], p = adjust.p, simulate.p.value = T, B = 1000);
test.adjust.p.sim <- chisq.test(M1.bcf[c(1,3)], p = adjust.p);

chisq.test(M1.bcf.data)


debug(chisqTestOnIL)
result <- chisqTestOnIL(geno.data.bc5f2, 5,2);
result2 <- chisqTestOnIL(geno.data.30.il, 5,2)
result3 <- chisqTestOnIL(geno.data.40.il,3,1)
write.csv(geno.data.30.il, file = "e://sampledata//geno.data.30.il.csv")

############################################
# Testing single marker analysis
############################################
library(lme4)
library(lsmeans);
library(phia)
demo.data.sm$line2 <- rep(c(1,2,3,4,5), each=4, times=2 )
demo.data.sm$line2 <- factor(demo.data.sm$line2);
#mod.lm <- lm(var ~  M3 , data = demo.data.sm)
mod.ran <- lmer(var ~  M3  + (1|lines) + (1|block), data = demo.data.sm)
#mod.fix <- lmer(var ~  M3  + lines + (1|block), data = demo.data.sm)
#mod.fix2 <- lmer(var ~  M3/lines + (1|block), data = demo.data.sm)
#mod.fix3 <- lmer(var ~  M3/line2 + (1|block), data = demo.data.sm)
#aov.fix3 <- aov(var ~  M3/line2 + block, data = demo.data.sm)


lsmeans(mod.ran, "M3");
testInteractions(mod.ran)


authors <- data.frame(
		surname = I(c("Tukey", "Venables", "Tierney", "Ripley", "McNeil")),
		nationality = c("US", "Australia", "US", "UK", "Australia"),
		deceased = c("yes", rep("no", 4)))
books <- data.frame(
		name = I(c("Tukey", "Venables", "Tierney",
						"Ripley", "Ripley", "McNeil", "R Core")),
		title = c("Exploratory Data Analysis",
				"Modern Applied Statistics ...",
				"LISP-STAT",
				"Spatial Statistics", "Stochastic Simulation",
				"Interactive Data Analysis",
				"An Introduction to R"),
		other.author = c(NA, "Ripley", NA, NA, NA, NA,
				"Venables & Smith"))

(m1 <- merge(authors, books, by.x = "surname", by.y = "name"))
(m2 <- merge(books, authors, by.x = "name", by.y = "surname"))


## Meng random population
e.ran.cold <- read.table("e://data//meng//il//E-Random-Cold.txt", header = T, row.names = 1);
marker.names <- rownames(e.ran.cold);
e.ran.cold <- t(e.ran.cold)
e.ran.cold <- rbind(marker.names, e.ran.cold)
## meng selected population
e.sele.cold <- read.table("e://data/meng//il//Cold-ILs-Geno.txt", header = T, row.names = 1);
e.sele.cold <- t(e.sele.cold);

debug(refFreq)

apply(e.ran.cold, 2, refFreq, "B", "H", "A")

debug(chisqTestOnIL)
chisqTestOnIL(e.ran.cold, 2, 6)
chisqTestOnIL(e.sele.cold, 2, 6)
x <- chisqTestOnIL(e.sele.cold, 2, 6, ref.matrix = e.ran.cold)
chisqTestOnIL(e.sele.cold, 2, 6, inc.ht = F,  ref.matrix = e.ran.cold)

undebug(chisqTestOnIL)


