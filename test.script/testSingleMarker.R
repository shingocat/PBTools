###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Nov 28, 2014
# FileName: testSingleMarker.R
###############################################################################



lines <- rep(1:10, each= 6);
block <- rep(1:3, each=2, times = 10);
rep <- rep(1:2, each = 1, times = 30);
var <- rnorm(60, mean = 100, sd = 5);
pheno <- cbind(lines, block, rep, var);
pheno <- as.data.frame(pheno);
# 0 for recurrent, 1 for hybrid, 2 for donor
M1 <- rep(c(0,2), each = 5);
M2 <- rep(c(2,1,0), times = c(6,1,3));
M3 <- rep(c(2,0, NA), times = c(4,5,1));
geno <- cbind(lines = c(1:5, letters[1:5]), M1, M2, M3);
geno <- as.data.frame(geno);
demo.data.sm <- merge(geno, pheno, by = "lines", incomparables=NA)
#demo.data.sm <- cbind(lines, M1, M2, M3, block, rep, var = var)
demo.data.sm <- as.data.frame(demo.data.sm);
demo.data.sm$lines <- factor(demo.data.sm$lines);
demo.data.sm$block <- factor(demo.data.sm$block);
demo.data.sm$M1 <- factor(demo.data.sm$M1);
demo.data.sm$M2 <- factor(demo.data.sm$M2);
demo.data.sm$M3 <- factor(demo.data.sm$M3);
demo.data.sm$rep <- factor(demo.data.sm$rep);
demo.data.sm$env <- "E1";
demo.data.sm$env <- factor(demo.data.sm$env);


#---one trait, one environment, 10 genotype, four block and three marker
demo.data.sm

library("lme4")
library("phia")
#--- single environment---#
#---lines nested in marker, and marker as fixed effect---#
lmer.m1 <- lmer(var ~ 1 + M1 + (1|M1:lines) + (1 | block), data= demo.data.sm)
aov.m1 <- aov(var ~ 1 +  M1/lines + Error(block), data = demo.data.sm )
lmer.m1.update.without.m1 <- update(lmer.m1, as.formula(paste(".~.-","M1")), evaluate = T)
lmer.m1.update.without.int <- update(lmer.m1, .~.-1);
m1.model.comp <- anova(lmer.m1.update, lmer.m1)
lmer.m1.anova <- anova(lmer.m1)
#Analysis of Variance Table
#	Df Sum Sq Mean Sq F value
#M1  1  71.94   71.94  2.5708
testInteractions(lmer.m1, fixed = "M1")
#Chisq Test: 
#		P-value adjustment method: holm
#Value Df  Chisq Pr(>Chisq)
#0  99.248  1 8780.4  < 2.2e-16
#2 101.649  1 9210.5  < 2.2e-16
testInteractions(lmer.m1,test=c("Chisq") )
#Chisq Test: 
#		P-value adjustment method: holm
#Value Df  Chisq Pr(>Chisq)
#0-2 -2.4017  1 2.5708     0.1089
testInteractions(lmer.m1,test=c("F") )
#F Test: 
#		P-value adjustment method: holm
#			Value Df      F Pr(>F)
#0-2       -2.4017  1 2.5708 0.1263
#Residuals         18   
lmer.m2 <- lmer(var ~ 1 + M2 + (1|M2:lines) + (1 | block), data= demo.data.sm)
aov.m2 <- aov(var ~ 1 + M2/lines + Error(block), data = demo.data.sm)
lmer.m2.anova <- anova(lmer.m2)
#Analysis of Variance Table
#	Df Sum Sq Mean Sq F value
#M2  2 100.31  50.155  1.8358
testInteractions(lmer.m2, fixed = "M2")
#Chisq Test: 
#		P-value adjustment method: holm
#	 Value Df  Chisq Pr(>Chisq)
#0 102.522  1 4524.3  < 2.2e-16
#1  99.750  1 4283.0  < 2.2e-16
#2  99.074  1 4225.1  < 2.2e-16
testInteractions(lmer.m2, test=c("Chisq"))
#Chisq Test: 
#		P-value adjustment method: holm
#     Value Df  Chisq Pr(>Chisq)
#0-1 2.7721  1 2.1128     0.2921
#0-2 3.4482  1 3.2690     0.2118
#1-2 0.6760  1 0.1257     0.7230
testInteractions(lmer.m2, test = "F")
#F Test: 
#		P-value adjustment method: holm
#			Value     Df      F Pr(>F)
#0-1       2.7721  1.000 1.6575 0.4176
#0-2       3.4482  1.000 2.5645 0.3625
#1-2       0.6760  1.000 0.0986 0.7559
#Residuals        27.205 
lmer.m3 <- lmer(var ~ 1 + M3 + (1|M3:lines) + (1 | block), data= demo.data.sm)
lmer.m3.anova <- anova(lmer.m3)
testInteractions(lmer.m3)

#---Multi Environment ---#
var1 <- rnorm(60, mean = 80, sd = 5);
demo.data.sm1 <- cbind(lines, M1, M2, M3, block, rep, var = var1)
demo.data.sm1 <- as.data.frame(demo.data.sm1);
demo.data.sm1$lines <- factor(demo.data.sm1$lines);
demo.data.sm1$block <- factor(demo.data.sm1$block);
demo.data.sm1$M1 <- factor(demo.data.sm1$M1);
demo.data.sm1$M2 <- factor(demo.data.sm1$M2);
demo.data.sm1$M3 <- factor(demo.data.sm1$M3);
demo.data.sm1$rep <- factor(demo.data.sm1$rep)
demo.data.sm1$env <- "E2";
demo.data.sm1$env <- factor(demo.data.sm1$env);
sm.data <- rbind(demo.data.sm, demo.data.sm1)


#--- include env factor and env is setted to be fixed---#
lmer.m1.env.fixed <- lmer(var ~ 1 + M1 + env + (1|M1:lines) + (1 | block), data= sm.data)
anova(lmer.m1.env.fixed)
#Analysis of Variance Table
#	 Df  Sum Sq Mean Sq F value
#M1   1     8.9     8.9   0.331
#env  1 13067.6 13067.6 485.405
testInteractions(lmer.m1.env.fixed, pairwise = "M1" , test = "F")
#F Test: 
#		P-value adjustment method: holm
#			  Value Df     F Pr(>F)
#0-2       -0.57727  1 0.331 0.5722
#Residuals          18  
testInteractions(lmer.m1.env.fixed, pairwise = "M1" , test = "Chisq")
#Chisq Test: 
#		P-value adjustment method: holm
#		Value Df Chisq Pr(>Chisq)
#0-2 -0.57727  1 0.331     0.5651
testInteractions(lmer.m1.env.fixed, fixed="env", across="M1")
#Chisq Test: 
#		P-value adjustment method: holm
#      Value Df Chisq Pr(>Chisq)
#E1 -0.57727  1 0.331          1
#E2 -0.57727  1 0.331          1
lmer.m2.env.fixed <- lmer(var ~ 1 + M2 + env + (1|M2:lines) + (1 | block), data = sm.data)
anova(lmer.m2.env.fixed)
#Analysis of Variance Table
#	 Df  Sum Sq Mean Sq  F value
#M2   2    70.7    35.3   1.3239
#env  1 13067.6 13067.6 489.7055
testInteractions(lmer.m2.env.fixed, pairwise = "M2", test="F")
#F Test: 
#		P-value adjustment method: holm
#			Value     Df      F Pr(>F)
#0-1       0.81671  1.000 0.3009 0.8782
#0-2       1.99573  1.000 1.7966 0.5920
#1-2       1.17902  1.000 0.6270 0.8782
#Residuals         17.415  
testInteractions(lmer.m2.env.fixed, pairwise = "M2", test="Chisq")
#Chisq Test: 
#		P-value adjustment method: holm
#		Value Df  Chisq Pr(>Chisq)
#0-1 0.81671  1 0.4386     0.6781
#0-2 1.99573  1 2.6190     0.3168
#1-2 1.17902  1 0.9140     0.6781
testInteractions(lmer.m2.env.fixed, fixed="env", across="M2")
#Chisq Test: 
#		P-value adjustment method: holm
#	   M21   M22 Df  Chisq Pr(>Chisq)
#E1 1.9957 1.179  2 2.6477     0.5322
#E2 1.9957 1.179  2 2.6477     0.5322
#
#--- include env factor and env is setted to be random---#
lmer.m1.env.random <- lmer(var ~ 1 + M1 + (1 | env) + (1|M1:env) + (1|M1:lines) + (1 | block), data = sm.data)
anova(lmer.m1.env.random)
#Analysis of Variance Table
#	Df Sum Sq Mean Sq F value
#M1  1 8.9109  8.9109   0.331
testInteractions(lmer.m1.env.random, pairwise = "M1", test="Chisq")
#Chisq Test: 
#		P-value adjustment method: holm
#       Value Df Chisq Pr(>Chisq)
#0-2 -0.57727  1 0.331     0.5651
# testInteractions(lmer.m1.env.random, fixed="M1", across="env", test="F") # Not working
#Error in parse(text = x) : <text>:2:0: unexpected end of input
#1: ~0 + 
#		^
#		In addition: Warning message:
#		In testInteractions(lmer.m1.env.random, fixed = "M1", across = "env",  :
#						Some factors with specified contrasts are not in the model and will be ignored.
lmer.m2.env.random <- lmer(var ~ 1 + M2 + (1 | env) + (1|M2:lines) + (1 | block), data = sm.data)
anova(lmer.m2.env.random)
#Analysis of Variance Table
#   Df Sum Sq Mean Sq F value
#M2  2 70.654  35.327  1.3239
testInteractions(lmer.m2.env.random, pairwise = "M2", test="Chisq")
#Chisq Test: 
#		P-value adjustment method: holm
#      Value Df  Chisq Pr(>Chisq)
#0-1 0.81671  1 0.4386     0.6781
#0-2 1.99573  1 2.6190     0.3168
#1-2 1.17902  1 0.9140     0.6781
testInteractions(lmer.m2.env.random, fixed = "env", across="M2", test="Chisq")
#Chisq Test: 
#		P-value adjustment method: holm
#		 M21   M22 Df  Chisq Pr(>Chisq)
#Mean 1.9957 1.179  2 2.6477     0.2661
