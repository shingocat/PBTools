###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Sep 8, 2014
# FileName: TestScript.R
###############################################################################

## --- prepare demo data of nine genotype of RCB design; --- #
nine.genotype <- c("AABB", "AABb",  "AAbb",	"AaBB",	"AaBb",	"Aabb",	"aaBB",	"aaBb",	"aabb");
three.block <- c("Block1", "Block2", "Block3");
three.rep <- c("Rep1", "Rep2", "Rep3");
five.env <- c("E1", "E2", "E3", "E4", "E5");
data <- data.frame("Env" = rep(five.env, each = 81),"Genotype" = rep(nine.genotype, each=9), "Block" =rep(three.block, each= 3), "Rep" = rep(three.rep));
genotype.means <- c(120,115,110,105,100,95,90,85,85);
env.effect <- c(3,2,-1,2,4);
yield <- c() ;
for(i in 1:5){
	for(j in 1:9) {
		yield <- c(yield,rnorm(9,genotype.means[j] + env.effect[i]));
	}
}
data$Yield <- yield;
tiller.number <- sample(9:15, 9, replace=T);
TN <- c() ;for(i in 1:5){
	for(j in 1:9) { 
		TN <- c(TN,rnorm(9,tiller.number[j] + env.effect[i]));
	}
}
data$TillerNum <- TN;

# --- testing source script of pyramided line multi site analysis --- #
# --- by far, it only considered only on trait first. --- #
library(PBTools);
result <- list();
geno="Genotype";
env="Env";
block="Block";
rep="Rep";
respvar <- c("Yield", "TillerNum");
result$respvar <- respvar[1];
# --- create temp.data without missing observations --- #
temp.data <- subset(data, subset = (is.na(data[,match(respvar[1], names(data))]) == FALSE));
temp.data[,match(geno, names(temp.data))] <- factor(trimStrings(temp.data[,match(geno, names(temp.data))]));
temp.data[,match(env, names(temp.data))] <- factor(trimStrings(temp.data[,match(env,names(temp.data))]));

# --- get levles of genotype and environment --- #
levelsGeno <- levels(temp.data[,match(geno, names(temp.data))]);
levelsEnv <- levels(temp.data[,match(env, names(temp.data))]);

# --- if max length of the character of levelsGeno or levelsEnv greater than 4, recode the levels --- #
if(max(nchar(levelsGeno)) > 4 || max(nchar(levelsEnv)) > 4)
{
	# --- recode genotype and environment levels --- #
	newCodingGeno <- data.frame(Genotype=levelsGeno, Code=paste("G", seq(1:length(levelsGeno)),sep=""));
	newCodingEnv <- data.frame(Envrironment = levelsEnv, Code=paste("E",seq(1:length(levelsEnv)),sep=""));
	
	result$newCodingGeno <- newCodingGeno;
	result$newCodingEnv <- newCodingEnv;
	recodedLevels <- TRUE;
	result$recodedLevels <- recodedLevels;
	
	# --- attach the new labels to temp.data --- #
	temp.data$CodedGeno <- newCodingGeno$Code[match(temp.data[,geno],newCodingGeno$Genotype)];
	temp.data$CodedEnv <- newCodingEnv$Code[match(temp.data[,env], newCodingEnv$Environment)];
	
} else
{
	temp.data$CodedGeno <- temp.data[ ,match(geno, names(temp.data))];
	temp.data$CodedEnv <- temp.data[ ,match(env, names(temp.data))];
	
	recodedLevels <- FALSE;
	result$recodedLevels <- recodedLevels;
}

# --- compute number of observations read, used and response rate --- #
obsread <- nrow(data);
obsused <- nrow(temp.data);
result$obsread <- obsread;
result$obsused <- obsused;
responseRate <- (obsused/obsread);
result$responseRate <- responseRate;

# --- compute summary statistics per environment --- #
sumStat.Env <- DescriptiveStatistics(temp.data, respvar[1], env, c("min", "mean", "max", "var","sd"));
sumStat.Env <- sumStat.Env[,c(2:ncol(sumStat.Env))];

is.envRandom = FALSE;

if(is.envRandom) 
{
	env.stmt <- paste("(1|", env, ") + (1|",geno, ":", env, ")", sep="") 
}else
{
	env.stmt <- paste(env," + ", geno, ":", env, sep ="");
}

exptl.design = "RCB";
# --- if design is Latinized Row-Column, check if the data follow case1 or case3 labeling --- #
if(exptl.design == "LatinRowCol")
{
	lengthPerCross <- tapply(temp.data[ , respvar[i]], temp.data[ ,c(row, column)], length);
	if(all(lengthPerCross <= nlevels(temp.data[,env]), na.rm=TRUE))
	{
		if(nlevels(temp.data[,row]) > nlevels(temp.data[,column]))
		{
			longerRow <- TRUE;
		} else
		{
			longerRow <- FALSE;
		}
	} else
	{
		stop("The levels of the row/column variable should be continuous across replicates.");
	}	
}

if(exptl.design == "RCB" || exptl.design =="AugRCB")
{
	myformula1 <- paste(respvar[1], " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , block , ":", env,")", sep="");
} else if(exptl.design == "AugLS")
{
	myformula1 <- paste(respvar[1], " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , row , ":", env, ") + (1|", column, ":", env, ")", sep="");
} else if(exptl.design == "Alpha")
{
	myformula1 <- paste(respvar[1], " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , rep , ":", env,") + (1|", rep, ":", block,":",env,")", sep = "" );
} else if(exptl.design == "RowCol")
{
	myformula1 <- paste(respvar[1], " ~ 1 + ", geno, " + " , env.stmt , " + (1|", rep , ":", env,") + (1|", rep,":", row,":", env,") + (1|", rep,":", column,":", env,")", sep = "");
} else if(exptl.design == "LatinAlpha")
{
	myformula1 <- paste(respvar[1], " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , rep , ":", env, ") + (1|", block, ":", env,") + (1|", rep, ":", block, ":", env,")", sep = "");
} else if(exptl.design == "LatinRowCol")
{
	if(longerRow)
	{
		myformula1 <- paste(respvar[1], " ~ 1 + ", geno, " + ", env.stmt , " + (1|" , rep, ":", env, ") + (1|", column, ":", env, ") + (1|", rep, ":", column, ":", env, ") + (1|", row, ":", env,")", sep = "");
	} else
	{
		myformula1 <- paste(respvar[1], " ~ 1 + ", geno, " + ", env.stmt , " + (1|" , rep, ":", env, ") + (1|", row, ":", env, ") + (1|", rep, ":" + row, ":", env,") + (1|", column, ":" , env, ")", sep = "");
	}
}


model <- lmer(formula(myformula1), data = temp.data);
result$formula1 <- myformula1;
result$model <- model;

# --- VARIANCE COMPONENTS --- #
varcomp <- NULL;
for(j in (1:length(VarCorr(model))))
{
	varcomp <- rbind(varcomp, data.frame(Groups = names(VarCorr(model))[j], Variance = VarCorr(model)[[j]][1], Std.Dev. = attr(VarCorr(model)[[j]], "stddev")[[1]]));
}
varcomp <- rbind(varcomp, data.frame(Groups = "Residual", Variance = attr(VarCorr(model), "sc")^2, Std.Dev. = attr(VarCorr(model), "sc")));
result$varcomp.table <- varcomp;

# --- TEST FOR SIGNIFICANCE OF GENOTYPIC EFFECT USING LRT && ANOVA TABLE IF GENO IS FIXED --- #

# --- myformula2 is full model minus the genotype term --- #
myformula2 <- sub(paste(" + ", geno, sep = ""), "", myformula1, fixed = TRUE);
model1 <- lmer(formula(myformula1), data = temp.data, REML = T);
model2 <- lmer(formula(myformula2), data = temp.data, REML = T);
result$formula2 <- myformula2;
# --- This is for genotype as Fixed --- #
anova.table1 <- anova(model2, model1);
rownames(anova.table1)<-c("Model2", "Model1");
attr(anova.table1, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF GENOTYPIC EFFECT USING LIKELIHOOD RATIO TEST:\n", sep = "");
attr(anova.table1, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
attr(anova.table1, "heading")[3] <- paste("Formula for Model2: ", myformula2, sep = "");
attr(anova.table1, "heading")[4] <- paste("", sep = "");
result$testsig.Geno <- anova.table1;

# --- TEST FOR SIGNIFICANCE OF ENVIRONMENT EFFECT USING LRT IF ENV IS FIXED--- #

# --- myformula3 is full model minus the environment term --- #
myformula3 <- gsub(paste(" + ", env, sep=""), "", myformula1, fixed = TRUE);
model3 <- lmer(formula(myformula3), data = temp.data, REML = TRUE);
result$formula3 <- myformula3;

anova.table2 <- anova(model3, model1);
rownames(anova.table2) <- c("Model3", "Model1");
attr(anova.table2, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF ENVIRONMENTAL EFFECT USING LIKELIHOOD RATIO TEST:\n", sep ="");
attr(anova.table2, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
attr(anova.table2, "heading")[3] <- paste("Formula for Model3: ", myformula3, sep = "");
attr(anova.table2, "heading")[4] <- paste("", sep="");
result$testsig.Env <- anova.table2;

# --- TEST OF SIGNIFICANCE OF GENOTYPE X ENVIRONMENT EFFECT USING LRT --- #

# --- myformula4 is full model minus the geno by environment interaction term --- #
myformula4 <- gsub(paste(" + ", geno, ":", env, sep = ""), "", myformula1, fixed = TRUE);
model4 <- lmer(formula(myformula4), data = temp.data, REML = TRUE);
result$formula4 <- myformula4;

anova.table3 <- anova(model4, model1);
rownames(anova.table3) <- c("Model4", "Model1");
attr(anova.table3, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF GENOTYPE BY ENVIRONMENT EFFECT USING LIKELIHOOD RATIO TEST:\n", sep="");
attr(anova.table3, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
attr(anova.table3, "heading")[3] <- paste("Formula for Model4: ", myformula4, sep = "");
attr(anova.table3, "heading")[4] <- paste("", sep = "");
result$testsig.GenoEnv <- anova.table3;


# --- PREDICTED MEANS/LSMEANS OF GENOTYPE AND SED STATISTICS--- #
# --- myformula5 is full model but without intercept --- #
myformula5 <- gsub("~ 1", "~ 0", myformula1, fixed = TRUE);
model.noint <- lmer(formula(myformula5), data = temp.data);

# the next line of code was deleted for R Version 3.0.2 by AAGulles 07.28.2014
# sumStat.Geno <- data.frame(summary(model.noint)@coef)[,1:2]
sumStat.Geno <- data.frame(summary(model.noint)$coefficients)[,1:2];
rownames(sumStat.Geno) <- gsub(geno,"", rownames(sumStat.Geno));
sumStat.Geno <- cbind(rownames(sumStat.Geno), sumStat.Geno);
colnames(sumStat.Geno) <- c(geno, "LSMean", "StdErrMean");

sumStat.Geno <- lsmeans(model1, geno);
sumStat.Geno <- summary(sumStat.Geno)[,1:3];
colnames(sumStat.Geno) <- c(geno, "LSMean", "StdErrMean");

library(lsmeans);
model.lmer <- lmer(Yield ~ 1 + Genotype + Env + Genotype:Env + (1|Block:Env), data = temp.data);
sumStat.Geno.lsmeans <- lsmeans(model, geno, cov.reduce = FALSE);
sumStat.Env.lsmeans <- lsmeans(model, env, cov.reduce = FALSE);

sumStat.GXE.lsmeans <- lsmeans(model, c(geno,env));


# --- display standard error of the differences --- #
noEntries<-nlevels(temp.data[,match(geno, names(temp.data))]);
covs <- as.matrix(vcov(model.noint)[1:noEntries, 1:noEntries]);
vars <- diag(covs);
vdiff <- outer(vars, vars, "+") - 2 * covs;
sed <- sqrt(vdiff[upper.tri(vdiff)]);

# --- display SED Table --- #
minSed<-formatC(as.numeric(format(min(sed), scientific=FALSE)), format="f");
meanSed<-formatC(as.numeric(format(mean(sed), scientific=FALSE)), format="f");
maxSed<-formatC(as.numeric(format(max(sed), scientific=FALSE)), format="f");
sedCol<-rbind(minSed, meanSed, maxSed);
rowNames<-rbind("Minimum  ", "Average  ", "Maximum  ");
sedTable<-as.table(cbind(rowNames, sedCol));
rownames(sedTable)<-c("","","");
colnames(sedTable)<-c("","Estimate");
result$sedTable <- sedTable;

rownames(sumStat.Geno) <- NULL;

# --- ESTIMATES OF EFFECTS --- #

# -- GENOTYPE EFFECT --- #
model.fix.effect <- fixef(model);
geno.with.levels.names <- paste(geno, levelsGeno, sep = "");
intercept <- fixef(model)[[1]];
geno.effect <- model.fix.effect[match(geno.with.levels.names, names(model.fix.effect))];
geno.effect <- as.data.frame(geno.effect[!is.na(geno.effect)]);
geno.effect <- data.frame(gsub(geno, "", rownames(geno.effect)), geno.effect);
colnames(geno.effect) <- c(geno, "geno_effect");
rownames(geno.effect) <- NULL;
# -- ENVIRONMENT EFFECT--- #
env.with.levels.names <- paste(env, levelsEnv, sep="");
env.effect <- model.fix.effect[match(env.with.levels.names, names(model.fix.effect))];
env.effect <- as.data.frame(env.effect[!is.na(env.effect)]);
env.effect <- data.frame(gsub(env, "", rownames(env.effect)), env.effect);
colnames(env.effect) <- c(env, "env_effect");
rownames(env.effect) <- NULL;

# -- G X E EFFECT --- #
GXE.with.levels.names <-  as.vector(outer(geno.with.levels.names, env.with.levels.names, paste, sep=":"));
GXE.effect <- model.fix.effect[match(GXE.with.levels.names, names(model.fix.effect))];
GXE.effect <- as.data.frame(GXE.effect[!is.na(GXE.effect)]);
names <- t(as.data.frame(strsplit(rownames(GXE.effect), ":")));
GXE.effect <- data.frame(gsub(geno, "", names[,1]), gsub(env, "", names[,2]), GXE.effect+intercept);
colnames(GXE.effect) <- c(geno, env, "ge_effect");
rownames(GXE.effect) <- NULL;

# -- G X E MEANS --- #
sumStat.GenoEnv <- merge(GXE.effect, env.effect, by = env, all = TRUE);
sumStat.GenoEnv <- merge(sumStat.GenoEnv, geno.effect, by = geno, all = TRUE);
sumStat.GenoEnv <- data.frame(sumStat.GenoEnv[,match(geno, names(sumStat.GenoEnv))], sumStat.GenoEnv[,match(env,names(sumStat.GenoEnv))], rowSums(subset(sumStat.GenoEnv, select = c(ge_effect, env_effect, geno_effect)), na.rm = TRUE));
colnames(sumStat.GenoEnv) <- c(geno, env, paste(respvar[1], "means", sep = "_"));

# -- create G x E MEANS with coded levels -- #
sumStat.GenoEnvCode<-sumStat.GenoEnv
sumStat.GenoEnvCode$CodedGeno <- sumStat.GenoEnvCode[,match(geno, names(sumStat.GenoEnvCode))]
sumStat.GenoEnvCode$CodedEnv <- sumStat.GenoEnvCode[,match(env, names(sumStat.GenoEnvCode))]


# --- display G x E means in 2-way table --- #
#wide.GenoEnv<-ToWide(sumStat.GenoEnv, paste(respvar[i], "means", sep = "_"), env, geno)
wide.GenoEnv<-reshape(sumStat.GenoEnv, v.names=paste(respvar[1], "means", sep = "_"), timevar=env, idv=geno, direction="wide")
colnames(wide.GenoEnv)<-gsub(paste(respvar[1], "means.", sep = "_"), "", colnames(wide.GenoEnv))
rownames(wide.GenoEnv)<-1:nrow(wide.GenoEnv)


result$means.Geno    <- sumStat.Geno
result$means.Env     <- sumStat.Env
result$means.GenoEnv <- sumStat.GenoEnv
result$wide.GenoEnv <- wide.GenoEnv
result$means.GenoEnvCode <- sumStat.GenoEnvCode
result$residuals <- resid(model1)
result$fitted.values <- fitted(model1)
result$data <- temp.data

# --- if genotype is fixed, output MSE and harmonic mean for AMMI --- #
# --- compute harmonic mean per environment level --- #
envgenorep <- as.data.frame.table(tapply(temp.data[, respvar[1]], temp.data[,c(env, geno)], length))
envgenorep <- envgenorep[(is.na(envgenorep[,"Freq"]) == FALSE),]
envgenorep$reciprocal <- 1/envgenorep$Freq
envgenorep2 <- as.data.frame.table(tapply(envgenorep[, "reciprocal"], envgenorep[,env], mean))
envgenorep2$reciprocal2 <- 1/envgenorep2$Freq
envgenorep3 <- merge(envgenorep, envgenorep2, by.x=env, by.y="Var1")

numrepEnv <- tapply(envgenorep3[,"reciprocal2"], envgenorep3[,env], mean)
no.reps <- 1/mean(1/numrepEnv)

result$harmonicMean <-no.reps
result$MSE <- varcomp[varcomp[,1] == "Residual", "Variance"]


# --- consolidate means and residuals --- #

for (m in 1:length(respvar)) {
	# --- consolidate genotype x environment means --- #
	if (m==1) {  means.GenoEnv.all <- as.data.frame(result$means.GenoEnv)
	} else {
		newGenoEnv <- as.data.frame(result$means.GenoEnv)
		if (nrow(newGenoEnv) > 0) {
			if (nrow(means.GenoEnv.all) == 0) { means.GenoEnv.all <- newGenoEnv
			} else { means.GenoEnv.all <- merge(means.GenoEnv.all, newGenoEnv, by=c(geno, env), all=TRUE)  }
			
		}
	}
	
	# --- consolidate genotype means --- #
	if (m==1) {
		newGeno <- as.data.frame(result$means.Geno)
		if (nrow(newGeno) > 0) { colnames(newGeno) <- c(geno, paste(respvar[m], "_LSMean", sep=""), paste(respvar[m], "_StdErrMean", sep="") ) }
		means.Geno.all <- newGeno
	} else {
		newGeno <- as.data.frame(result$means.Geno)
		if (nrow(newGeno) > 0) {
			colnames(newGeno) <- c(geno, paste(respvar[m], "_LSMean", sep=""), paste(respvar[m], "_StdErrMean", sep="") )
			if (nrow(means.Geno.all) == 0) { means.Geno.all <- newGeno
			} else { means.Geno.all <- merge(means.Geno.all, newGeno, by=c(geno), all=TRUE) }
		}
	}
}

# generate status of means.GenoEnv.all
if (nrow(means.GenoEnv.all) == 0) { meansGenoEnvWarning<-"empty"
} else { meansGenoEnvWarning<-"not empty" }

# generate status of means.Geno.all
if (nrow(means.Geno.all) == 0) { meansGenoWarning<-"empty"
} else { meansGenoWarning<-"not empty" }

options(warn = prev.opt)
detach("package:lme4")
return(list(output = result, means.GenoEnv.all=means.GenoEnv.all, meansGenoEnvWarning=meansGenoEnvWarning, means.Geno.all=means.Geno.all, meansGenoWarning=meansGenoWarning))










result.pyrmaidedline.envrandom <- pyramidedline.msa.onestage.test("RCB", data = data, respvar="Yield", geno = "Genotype", block="Block", env = "Env", is.envRandom = TRUE);
result.pyrmaidedline.onestage <- GEOneStage.test("RCB", data = data, respvar="Yield", geno = "Genotype", row="Block", env = "Env", is.genoRandom = FALSE);
result.pyrmaidedline.envfix <- pyramidedline.msa.onestage.test("RCB", data = data, respvar="Yield", geno = "Genotype", block="Block", env = "Env", is.envRandom = FALSE);

result.pyrmaidedline.envrandom.original <- pyramidedline.msa.onestage.test("RCB", data = data, respvar="Yield", geno = "Genotype", block="Block", env = "Env", is.envRandom = TRUE);
result.pyrmaidedline.envfix.original <- pyramidedline.msa.onestage.test("RCB", data = data, respvar="Yield", geno = "Genotype", block="Block", env = "Env", is.envRandom = FALSE);



# -- Testing the phia package for interaction contrast --- ##
library(phia);
some(Boik);
summary(Boik);

mod.boik <- lm(edr ~ therapy * medication, data = Boik);

par(mfocol = c(1,2));
plot(mod.boik, 1:2); # plot diagnostics for the model
Anova(mod.boik)

## contingency tabel
boik.means <- interactionMeans(mod.boik)

plot(boik.means, multiple = TRUE);# plot a graph of interaction and adjusted means


# define contrast
cntrl.vs.T1 <- list(therapy = c(1, -1, 0)) 
cntrl.vs.T2 <- list(therapy = c(1, 0, -1))
T1.vs.T2 <- list(therpy = c(0,1, -1));
plcb.vs.doses <- list(medication = c(1, -1/3, -1/3, -1/3))

testInteractions(mod.boik, fixed ="therapy", across="medication") ; # testing simple main effects for interaction on fixed factor by across factor;

testInteractions(mod.boik, pairwise="therapy", across="medication", adjustment="none")
testInteractions(mod.boik, custom=c(cntrl.vs.T1, plcb.vs.doses), adjustment="none")
testInteractions(mod.boik, custom=c(cntrl.vs.T1, cntrl.vs.T2), across="medication", adjustment="none")

#  Tesing simple effect
testInteractions(mod.boik, fixed ="therapy", across="medication") ; # testing simple main effects for interaction on fixed factor by across factor;

# Analysis of residual effects

boik.mtable <- xtabs(boik.means$"adjusted mean" ~ therapy + medication, boik.means);
boik.mtable <- addmargins(boik.mtable, FUN = mean, quiet = TRUE);
print(boik.mtable, digits=4);
# substract the mean
boik.resid <- boik.mtable - boik.mtable[4,5];
# substract row means
boik.resid <- sweep(boik.resid, 1, boik.resid[,5])
# substract column means
boik.resid <- sweep(boik.resid, 2, boik.resid[4,])
print(boik.resid, digits = 4)


# different between levels 2 and 3 on therapy
# by create a vector of the length of the coefficients list (including the intercetp)
# using multcomp pacakge method glht;
k <- matrix (c(0, 1, -1,0,0,0,0,0,0,0,0,0),1);
library(multcomp)
t <- glht(mod.boik, linfct = k);
summary(t);


# Testing my qm20 population 
qm20 <- read.csv("e://data/mine//qm20.new.txt",header = T,sep=" ")
qm20$Hd3a.SIND1 <- factor(qm20$Hd3a.SIND1);
qm20$PSM117 <- factor(qm20$PSM117);

summary(qm20);
# DaysToHeading column is not any missing value;
qm20.subset <- subset(qm20, subset=!is.na(qm20$DaysToFlowering)); 

# Define the model two marker and env as factor;
# linear model;
mod.qm20 <- lm(DaysToFlowering ~ Hd3a.SIND1 * PSM117 * Location, data = qm20);
# print anova table
anova(mod.qm20);
#Analysis of Variance Table
#
#Response: DaysToFlowering
#Df Sum Sq Mean Sq   F value    Pr(>F)    
#Hd3a.SIND1                    2  23747   11873  956.0588 < 2.2e-16 ***
#		PSM117                        2     53      26    2.1192  0.120352    
#Location                      1  68875   68875 5545.8185 < 2.2e-16 ***
#		Hd3a.SIND1:PSM117             4    197      49    3.9709  0.003245 ** 
#		Hd3a.SIND1:Location           2    297     149   11.9654 6.752e-06 ***
#		PSM117:Location               2     30      15    1.1928  0.303562    
#Hd3a.SIND1:PSM117:Location    4     15       4    0.2983  0.879188    
#Residuals                  2379  29545      12                        
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# print summary table;
summary(mod.qm20);
#Call:
#		lm(formula = DaysToFlowering ~ Hd3a.SIND1 * PSM117 * Location, 
#				data = qm20)
#
#Residuals:
#		Min       1Q   Median       3Q      Max 
#-14.2840  -1.8406   0.1257   2.0435  23.9366 
#
#Coefficients:
#		Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                    70.18519    0.33911 206.971  < 2e-16 ***
#		Hd3a.SIND12                     3.06337    0.43218   7.088 1.78e-12 ***
#		Hd3a.SIND13                    11.77134    0.80929  14.545  < 2e-16 ***
#		PSM1172                        -0.12175    0.40166  -0.303   0.7618    
#PSM1173                        -0.82698    0.45571  -1.815   0.0697 .  
#LocationJX                     10.75529    0.51268  20.979  < 2e-16 ***
#		Hd3a.SIND12:PSM1172             0.71384    0.51167   1.395   0.1631    
#Hd3a.SIND13:PSM1172            -1.55082    0.92446  -1.678   0.0936 .  
#Hd3a.SIND12:PSM1173             0.53675    0.58664   0.915   0.3603    
#Hd3a.SIND13:PSM1173            -1.42223    1.02496  -1.388   0.1654    
#Hd3a.SIND12:LocationJX         -0.04385    0.67736  -0.065   0.9484    
#Hd3a.SIND13:LocationJX          2.12819    1.13999   1.867   0.0620 .  
#PSM1172:LocationJX              0.05559    0.61404   0.091   0.9279    
#PSM1173:LocationJX              0.90420    0.68222   1.325   0.1852    
#Hd3a.SIND12:PSM1172:LocationJX -0.37239    0.80863  -0.461   0.6452    
#Hd3a.SIND13:PSM1172:LocationJX -0.13969    1.35123  -0.103   0.9177    
#Hd3a.SIND12:PSM1173:LocationJX -0.95193    0.90768  -1.049   0.2944    
#Hd3a.SIND13:PSM1173:LocationJX -0.30269    1.51087  -0.200   0.8412    
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#Residual standard error: 3.524 on 2379 degrees of freedom
#Multiple R-squared:  0.7593,	Adjusted R-squared:  0.7576 
#F-statistic: 441.5 on 17 and 2379 DF,  p-value: < 2.2e-16
# adjust means across each factor;
qm20.means <- interactionMeans(mod.qm20);
#Hd3a.SIND1 PSM117 Location adjusted mean
#1           1      1       GZ      70.18519
#2           2      1       GZ      73.24855
#3           3      1       GZ      81.95652
#4           1      2       GZ      70.06343
#5           2      2       GZ      73.84065
#6           3      2       GZ      80.28395
#7           1      3       GZ      69.35821
#8           2      3       GZ      72.95833
#9           3      3       GZ      79.70732
#10          1      1       JX      80.94048
#11          2      1       JX      83.96000
#12          3      1       JX      94.84000
#13          1      2       JX      80.87432
#14          2      2       JX      84.23529
#15          3      2       JX      93.08333
#16          1      3       JX      81.01770
#17          2      3       JX      83.62205
#18          3      3       JX      93.19231
# plot adjust means;
plot(qm20.means);

# adjusted means by factor Hd3a.SIND1
interactionMeans(mod.qm20, factors="Hd3a.SIND1");
# outcomes
#	Hd3a.SIND1 adjusted mean
#	1          1      75.40655
#	2          2      78.64415
#	3          3      87.17724
interactionMeans(mod.qm20, factors="PSM117");
#	Hd3a.SIND1 adjusted mean
#	1          1      75.40655
#	2          2      78.64415
#	3          3      87.17724
interactionMeans(mod.qm20, factors="Location");
#Location adjusted mean
#1       GZ      74.62246
#2       JX      86.19616
# Testing simple effects by fixed factor and across factor;
testInteractions(mod.qm20, fixed = "Hd3a.SIND1", across="PSM117");
#F Test: 
#		P-value adjustment method: holm
#			PSM1171 PSM1172   Df Sum of Sq      F  Pr(>F)  
#1         0.37488 0.28092    2      18.0 0.7252 0.48434  
#2         0.31409 0.74778    2     121.0 4.8703 0.01989 *
#3         1.94845 0.23383    2     124.9 5.0265 0.01989 *
#		Residuals                 2379   29545.3                 
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

testInteractions(mod.qm20, fixed = "PSM117", across="Hd3a.SIND1");
#F Test: 
#		P-value adjustment method: holm
#			Hd3a.SIND11 Hd3a.SIND12   Df Sum of Sq      F    Pr(>F)    
#1             -12.835     -9.7940    2    6321.9 254.52 < 2.2e-16 ***
#2             -11.215     -7.6457    2   12248.3 493.12 < 2.2e-16 ***
#3             -11.262     -8.1596    2    6502.3 261.78 < 2.2e-16 ***
#		Residuals                         2379   29545.3                     
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
testInteractions(mod.qm20, fixed = c("Hd3a.SIND1", "PSM117"), across="Location");
#F Test: 
#		P-value adjustment method: holm
#					Value   Df Sum of Sq       F    Pr(>F)    
#		1 : 1     -10.755    1    5465.7  440.10 < 2.2e-16 ***
#		2 : 1     -10.711    1    7270.8  585.44 < 2.2e-16 ***
#		3 : 1     -12.883    1    1988.4  160.10 < 2.2e-16 ***
#		1 : 2     -10.811    1   12709.6 1023.38 < 2.2e-16 ***
#		2 : 2     -10.395    1   16594.4 1336.19 < 2.2e-16 ***
#		3 : 2     -12.799    1    4937.6  397.58 < 2.2e-16 ***
#		1 : 3     -11.659    1    8333.8  671.05 < 2.2e-16 ***
#		2 : 3     -10.664    1    8692.2  699.90 < 2.2e-16 ***
#		3 : 3     -13.485    1    2893.2  232.96 < 2.2e-16 ***
#		Residuals         2379   29545.3                      
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# Anyalsis of residual effects;
qm20.mtable <- xtabs(qm20.means$"adjusted mean" ~ Hd3a.SIND1 + PSM117 + Location, qm20.means);
qm20.mtable <- addmargins(qm20.mtable, FUN = mean, quiet = TRUE);

# Create a new analysis based on nine genotype;

mod9.qm20 <- lm(DaysToFlowering ~ Genotype * Location, data = qm20);
anova(mod9.qm20);
summary(mod9.qm20);
qm20.9means <- interactionMeans(mod9.qm20);
#Genotype Location adjusted mean
#1      aabb       GZ      70.18519
#2      aaBb       GZ      70.06343
#3      aaBB       GZ      69.35821
#4      Aabb       GZ      73.24855
#5      AaBb       GZ      73.84065
#6      AaBB       GZ      72.95833
#7      AAbb       GZ      81.95652
#8      AABb       GZ      80.28395
#9      AABB       GZ      79.70732
#10     aabb       JX      80.94048
#11     aaBb       JX      80.87432
#12     aaBB       JX      81.01770
#13     Aabb       JX      83.96000
#14     AaBb       JX      84.23529
#15     AaBB       JX      83.62205
#16     AAbb       JX      94.84000
#17     AABb       JX      93.08333
#18     AABB       JX      93.19231
interactionMeans(mod9.qm20, factors ="Genotype");
#Genotype adjusted mean
#1     aabb      75.56283
#2     aaBb      75.46887
#3     aaBB      75.18795
#4     Aabb      78.60428
#5     AaBb      79.03797
#6     AaBB      78.29019
#7     AAbb      88.39826
#8     AABb      86.68364
#9     AABB      86.44981
interactionMeans(mod9.qm20, factors = "Location");
#Location adjusted mean
#1       GZ      74.62246
#2       JX      86.19616
plot(qm20.9means);

# Testing simple effect
testInteractions(mod9.qm20, fixed="Genotype", across="Location");
#F Test: 
#		P-value adjustment method: holm
#					Value   Df Sum of Sq       F    Pr(>F)    
#		aabb      -10.755    1    5465.7  440.10 < 2.2e-16 ***
#		aaBb      -10.811    1   12709.6 1023.38 < 2.2e-16 ***
#		aaBB      -11.659    1    8333.8  671.05 < 2.2e-16 ***
#		Aabb      -10.711    1    7270.8  585.44 < 2.2e-16 ***
#		AaBb      -10.395    1   16594.4 1336.19 < 2.2e-16 ***
#		AaBB      -10.664    1    8692.2  699.90 < 2.2e-16 ***
#		AAbb      -12.883    1    1988.4  160.10 < 2.2e-16 ***
#		AABb      -12.799    1    4937.6  397.58 < 2.2e-16 ***
#		AABB      -13.485    1    2893.2  232.96 < 2.2e-16 ***
#		Residuals         2379   29545.3                      
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


# Analysis of residual effects
qm20.mtable <- xtabs(qm20.9means$"adjusted mean" ~ Genotype + Location, qm20.9means);
qm20.mtable <- addmargins(qm20.mtable, FUN = mean, quiet = TRUE);
#				Location
#Genotype       GZ       JX     mean
#		aabb 70.18519 80.94048 75.56283
#		aaBb 70.06343 80.87432 75.46887
#		aaBB 69.35821 81.01770 75.18795
#		Aabb 73.24855 83.96000 78.60428
#		AaBb 73.84065 84.23529 79.03797
#		AaBB 72.95833 83.62205 78.29019
#		AAbb 81.95652 94.84000 88.39826
#		AABb 80.28395 93.08333 86.68364
#		AABB 79.70732 93.19231 86.44981
#		mean 74.62246 86.19616 80.40931

# substract the mean
qm20.resid <- qm20.mtable - qm20.mtable[10,3];
# substract the row means
qm20.resid <- sweep(qm20.resid, 1, qm20.resid[,3]);
# substract the column means
qm20.resid <- sweep(qm20.resid, 2, qm20.resid[10,]);
#					Location
#Genotype          GZ          JX        mean
#			aabb  0.40920579 -0.40920579  0.00000000
#			aaBb  0.38140924 -0.38140924  0.00000000
#			aaBB -0.04289378  0.04289378  0.00000000
#			Aabb  0.43112875 -0.43112875  0.00000000
#			AaBb  0.58952756 -0.58952756  0.00000000
#			AaBB  0.45499434 -0.45499434  0.00000000
#			AAbb -0.65488783  0.65488783  0.00000000
#			AABb -0.61284006  0.61284006  0.00000000
#			AABB -0.95564401  0.95564401  0.00000000
#			mean  0.00000000  0.00000000  0.00000000

# Testing directly by testInteractions via the argument residuals

testInteractions(mod9.qm20, residual = c("Genotype", "Location"));
#F Test: 
#		P-value adjustment method: holm
#								Value   Df Sum of Sq       F   Pr(>F)   
#aabb (resid.) : GZ (resid.)  0.40921    1      34.0  2.7391 0.588319   
#aaBb (resid.) : GZ (resid.)  0.38141    1      56.1  4.5133 0.404823   
#aaBB (resid.) : GZ (resid.) -0.04289    1       0.5  0.0372 1.000000   
#Aabb (resid.) : GZ (resid.)  0.43113    1      48.0  3.8617 0.413486   
#AaBb (resid.) : GZ (resid.)  0.58953    1     167.6 13.4973 0.004395 **
#AaBB (resid.) : GZ (resid.)  0.45499    1      61.8  4.9743 0.361464   
#AAbb (resid.) : GZ (resid.) -0.65489    1      25.2  2.0267 0.618743   
#AABb (resid.) : GZ (resid.) -0.61284    1      51.7  4.1662 0.413486   
#AABB (resid.) : GZ (resid.) -0.95564    1      70.1  5.6442 0.281465   
#aabb (resid.) : JX (resid.) -0.40921    1      34.0  2.7391 0.588319   
#aaBb (resid.) : JX (resid.) -0.38141    1      56.1  4.5133 0.404823   
#aaBB (resid.) : JX (resid.)  0.04289    1       0.5  0.0372 1.000000   
#Aabb (resid.) : JX (resid.) -0.43113    1      48.0  3.8617 0.413486   
#AaBb (resid.) : JX (resid.) -0.58953    1     167.6 13.4973 0.004395 **
#AaBB (resid.) : JX (resid.) -0.45499    1      61.8  4.9743 0.361464   
#AAbb (resid.) : JX (resid.)  0.65489    1      25.2  2.0267 0.618743   
#AABb (resid.) : JX (resid.)  0.61284    1      51.7  4.1662 0.413486   
#AABB (resid.) : JX (resid.)  0.95564    1      70.1  5.6442 0.281465   
#Residuals                            2379   29545.3                    
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# plot interaction
matplot(t(qm20.resid[-10,-3]), type = "b", xaxt ="n", ylab ="Interaction residuals");
axis(1, at=1:2, labels = levels(qm20$Location))

# pairwise contrast on Genotype and across Location
testInteractions(mod9.qm20, pairwise = "Genotype", across = "Location");
#F Test: 
#		P-value adjustment method: holm
#			  Value   Df Sum of Sq       F Pr(>F)  
#aabb-aaBb  0.05559    1       0.1  0.0082 1.0000  
#aabb-aaBB  0.90420    1      21.8  1.7566 1.0000  
#aabb-Aabb -0.04385    1       0.1  0.0042 1.0000  
#aabb-AaBb -0.36064    1       4.7  0.3784 1.0000  
#aabb-AaBB -0.09158    1       0.2  0.0197 1.0000  
#aabb-AAbb  2.12819    1      43.3  3.4851 1.0000  
#aabb-AABb  2.04409    1      76.9  6.1910 0.3485  
#aabb-AABB  2.72970    1      88.7  7.1413 0.2165  
#aaBb-aaBB  0.84861    1      28.2  2.2732 1.0000  
#aaBb-Aabb -0.09944    1       0.4  0.0319 1.0000  
#aaBb-AaBb -0.41624    1      11.0  0.8882 1.0000  
#aaBb-AaBB -0.14717    1       1.0  0.0783 1.0000  
#aaBb-AAbb  2.07259    1      46.4  3.7323 1.0000  
#aaBb-AABb  1.98850    1      93.3  7.5136 0.1851  
#aaBb-AABB  2.67411    1      99.3  7.9918 0.1564  
#aaBB-Aabb -0.94805    1      28.0  2.2551 1.0000  
#aaBB-AaBb -1.26484    1      70.1  5.6442 0.4574  
#aaBB-AaBB -0.99578    1      33.7  2.7162 1.0000  
#aaBB-AAbb  1.22399    1      15.0  1.2088 1.0000  
#aaBB-AABb  1.13989    1      26.3  2.1140 1.0000  
#aaBB-AABB  1.82550    1      42.1  3.3896 1.0000  
#Aabb-AaBb -0.31680    1       4.5  0.3625 1.0000  
#Aabb-AaBB -0.04773    1       0.1  0.0064 1.0000  
#Aabb-AAbb  2.17203    1      47.5  3.8271 1.0000  
#Aabb-AABb  2.08794    1      89.0  7.1698 0.2165  
#Aabb-AABB  2.77355    1      97.8  7.8773 0.1564  
#AaBb-AaBB  0.26907    1       3.7  0.2975 1.0000  
#AaBb-AAbb  2.48883    1      68.8  5.5425 0.4660  
#AaBb-AABb  2.40474    1     145.7 11.7316 0.0225 *
#AaBb-AABB  3.09034    1     137.7 11.0864 0.0309 *
#AaBB-AAbb  2.21976    1      51.0  4.1088 1.0000  
#AaBB-AABb  2.13567    1      98.6  7.9388 0.1564  
#AaBB-AABB  2.82128    1     104.8  8.4404 0.1259  
#AAbb-AABb -0.08410    1       0.1  0.0049 1.0000  
#AAbb-AABB  0.60151    1       2.5  0.1991 1.0000  
#AABb-AABB  0.68561    1       4.9  0.3941 1.0000  
#Residuals          2379   29545.3                 
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


custom.contr <- contrastCoefficients(
		Genotype ~ (AABB + AAbb)/4 - (aaBB + aabb)/4,
		Genotype ~ (AABB)/4 - (AAbb)/4 + (aaBB)/4 - (aabb)/4,
		Genotype ~ -(AABB)/4 - (AAbb)/4 + (AaBB)/2 + (Aabb)/2 -(aaBB)/4 -(aabb)/4,
		Genotype ~ -(AABB)/4 +  (AABb)/2  - (AAbb)/4 -(aaBB)/4 + (aaBb)/2 -(aabb)/4,
		Genotype ~ (AABB)/4 - (AAbb)/4 - (aaBB)/4 + (aabb) /4,
		Genotype ~ -(AABB)/4 + (AABb)/2 - (AAbb)/4 + (aaBB)/4 - (aaBb)/2 + (aabb)/4,
		Genotype ~ -(AABB)/4 + (AAbb)/4 + (AaBB)/2 - (Aabb)/2 - (aaBB)/4 + (aabb)/4,
		Genotype ~ (AABB)/4 - (AABb)/2 + (AAbb)/4 - (AaBB)/2 + (AaBb) - (Aabb)/2 + (aaBB)/4 - (aaBb)/2 + (aabb)/4,
		Location ~ GZ - JX, data = qm20, normalize = TRUE
)
names(custom.contr$Genotype) <- c("a1","a2","d1","d2","a1a2","a1d2","d1a2","d1d2");
#$Genotype
#		a1   a2         d1         d2 a1a2       a1d2       d1a2       d1d2
#aabb -0.5 -0.5 -0.2886751 -0.2886751  0.5  0.2886751  0.2886751  0.1666667
#aaBb  0.0  0.0  0.0000000  0.5773503  0.0 -0.5773503  0.0000000 -0.3333333
#aaBB -0.5  0.5 -0.2886751 -0.2886751 -0.5  0.2886751 -0.2886751  0.1666667
#Aabb  0.0  0.0  0.5773503  0.0000000  0.0  0.0000000 -0.5773503 -0.3333333
#AaBb  0.0  0.0  0.0000000  0.0000000  0.0  0.0000000  0.0000000  0.6666667
#AaBB  0.0  0.0  0.5773503  0.0000000  0.0  0.0000000  0.5773503 -0.3333333
#AAbb  0.5 -0.5 -0.2886751 -0.2886751 -0.5 -0.2886751  0.2886751  0.1666667
#AABb  0.0  0.0  0.0000000  0.5773503  0.0  0.5773503  0.0000000 -0.3333333
#AABB  0.5  0.5 -0.2886751 -0.2886751  0.5 -0.2886751 -0.2886751  0.1666667
#
#$Location
#Location
#GZ  0.7071068
#JX -0.7071068


# Testing gene effect by env interaction
testInteractions(mod9.qm20, custom=custom.contr);
#F Test: 
#		P-value adjustment method: holm
#					Value   Df Sum of Sq      F  Pr(>F)  
#  a1 : Location -1.39784    1      85.0 6.8478 0.06252 .
#  a2 : Location -0.53235    1      12.3 0.9932 1.00000  
#  d1 : Location  1.23147    1     121.6 9.7931 0.01418 *
#  d2 : Location  0.31899    1       6.9 0.5566 1.00000  
#a1a2 : Location  0.10702    1       0.5 0.0401 1.00000  
#a1d2 : Location -0.00476    1       0.0 0.0001 1.00000  
#d1a2 : Location  0.32684    1       8.6 0.6898 1.00000  
#d1d2 : Location -0.04608    1       0.3 0.0215 1.00000  
#Residuals                2379   29545.3                 
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Testing gene effect
custom.contr.geno <- contrastCoefficients(
		Genotype ~ (AABB + AAbb)/4 - (aaBB + aabb)/4,
		Genotype ~ (AABB)/4 - (AAbb)/4 + (aaBB)/4 - (aabb)/4,
		Genotype ~ -(AABB)/4 - (AAbb)/4 + (AaBB)/2 + (Aabb)/2 -(aaBB)/4 -(aabb)/4,
		Genotype ~ -(AABB)/4 +  (AABb)/2  - (AAbb)/4 -(aaBB)/4 + (aaBb)/2 -(aabb)/4,
		Genotype ~ (AABB)/4 - (AAbb)/4 - (aaBB)/4 + (aabb) /4,
		Genotype ~ -(AABB)/4 + (AABb)/2 - (AAbb)/4 + (aaBB)/4 - (aaBb)/2 + (aabb)/4,
		Genotype ~ -(AABB)/4 + (AAbb)/4 + (AaBB)/2 - (Aabb)/2 - (aaBB)/4 + (aabb)/4,
		Genotype ~ (AABB)/4 - (AABb)/2 + (AAbb)/4 - (AaBB)/2 + (AaBb) - (Aabb)/2 + (aaBB)/4 - (aaBb)/2 + (aabb)/4,
		 data = qm20, normalize = F
)
names(custom.contr.geno$Genotype) <- c("a1","a2","d1","d2","a1a2","a1d2","d1a2","d1d2");
testInteractions(mod9.qm20, custom=custom.contr.geno);
#F Test: 
#		P-value adjustment method: holm
#			 Value   Df Sum of Sq         F  Pr(>F)    
#  a1      12.0486    1   12636.8 1017.5182 < 2e-16 ***
#  a2      -1.1617    1     117.5    9.4586 0.01275 *  
#  d1      -3.4092    1    1864.3  150.1122 < 2e-16 ***
#  d2      -0.3735    1      19.0    1.5260 0.23632    
#a1a2      -0.7868    1      53.9    4.3389 0.14943    
#a1d2      -0.4814    1      31.5    2.5356 0.23632    
#d1a2       0.4893    1      38.4    3.0927 0.23632    
#d1d2       0.6095    1      93.4    7.5172 0.03078 *  
#		Residuals         2379   29545.3                      
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


custom.contr.env <- contrastCoefficients(
		Location ~ GZ - JX,
		data = qm20, normalize = TRUE
)
testInteractions(mod9.qm20, custom=custom.contr.env);
#F Test: 
#		P-value adjustment method: holm
#			 Value   Df Sum of Sq      F    Pr(>F)    
#Location  -8.1838    1     41483 3340.2 < 2.2e-16 ***
#Residuals         2379     29545                     
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Testing directly of contrast table;
# define genotype contrast
genotypeContrast <- list(Genotype = as.data.frame(t(getDefaultGenesContrast())));
testInteractions(mod9.qm20,custom = genotypeContrast);

temp.data1 <- qm20[,c("Location", "Genotype", "DaysToFlowering")];
temp.data1$Env <- NA;
temp.data1$Block <- NA;
temp.data2 <- qm20[,c("Location", "Genotype", "DaysToFlowering")];
temp.data2$Env <- NA;
temp.data2$Block <- NA;

temp.data1[temp.data1[1] == "GZ",]$Env <- c("E1");
temp.data1[temp.data1[1] == "JX",]$Env <- c("E2");

temp.data2[temp.data2[1] == "GZ",]$Env <- c("E3");
temp.data2[temp.data2[1] == "JX",]$Env <- c("E4");

temp.data <-rbind(temp.data1, temp.data2);

temp.data$Env <- factor(temp.data$Env);

names(temp.data) <- c("Location", "Genotype", "HD", "Env");
# define linear model of lm method;
mod.temp.data <- lm(HD ~ Genotype * Env, data = temp.data);
anova(mod.temp.data);



custom.contr.5env <- contrastCoefficients(
		Genotype ~ (AABB + AAbb)/4 - (aaBB + aabb)/4,
		Genotype ~ (AABB)/4 - (AAbb)/4 + (aaBB)/4 - (aabb)/4,
		Genotype ~ -(AABB)/4 - (AAbb)/4 + (AaBB)/2 + (Aabb)/2 -(aaBB)/4 -(aabb)/4,
		Genotype ~ -(AABB)/4 +  (AABb)/2  - (AAbb)/4 -(aaBB)/4 + (aaBb)/2 -(aabb)/4,
		Genotype ~ (AABB)/4 - (AAbb)/4 - (aaBB)/4 + (aabb) /4,
		Genotype ~ -(AABB)/4 + (AABb)/2 - (AAbb)/4 + (aaBB)/4 - (aaBb)/2 + (aabb)/4,
		Genotype ~ -(AABB)/4 + (AAbb)/4 + (AaBB)/2 - (Aabb)/2 - (aaBB)/4 + (aabb)/4,
		Genotype ~ (AABB)/4 - (AABb)/2 + (AAbb)/4 - (AaBB)/2 + (AaBb) - (Aabb)/2 + (aaBB)/4 - (aaBb)/2 + (aabb)/4,
		Env ~ E1 - E2,
		Env ~ E1 - E3,
		Env ~ E1 - E4,
		Env ~ E1 - E5,
		Env ~ E2 - E3,
		Env ~ E2 - E4,
		Env ~ E2 - E5,
		Env ~ E3 - E4,
		Env ~ E3 - E5,
		Env ~ E4 - E5,
		data = temp.data, normalize = TRUE
)
names(custom.contr.5env$Genotype) <- c("a1","a2","d1","d2","a1a2","a1d2","d1a2","d1d2");
names(custom.contr.5env$Env) <- c("E1 - E2", "E1 - E3", "E1 - E4", "E1 - E5", "E2 - E3", "E2 - E4", "E2 - E5", "E3 - E4","E3 - E5", "E4 - E5");
#$Genotype
#		a1   a2         d1         d2 a1a2       a1d2       d1a2       d1d2
#aabb -0.5 -0.5 -0.2886751 -0.2886751  0.5  0.2886751  0.2886751  0.1666667
#aaBb  0.0  0.0  0.0000000  0.5773503  0.0 -0.5773503  0.0000000 -0.3333333
#aaBB -0.5  0.5 -0.2886751 -0.2886751 -0.5  0.2886751 -0.2886751  0.1666667
#Aabb  0.0  0.0  0.5773503  0.0000000  0.0  0.0000000 -0.5773503 -0.3333333
#AaBb  0.0  0.0  0.0000000  0.0000000  0.0  0.0000000  0.0000000  0.6666667
#AaBB  0.0  0.0  0.5773503  0.0000000  0.0  0.0000000  0.5773503 -0.3333333
#AAbb  0.5 -0.5 -0.2886751 -0.2886751 -0.5 -0.2886751  0.2886751  0.1666667
#AABb  0.0  0.0  0.0000000  0.5773503  0.0  0.5773503  0.0000000 -0.3333333
#AABB  0.5  0.5 -0.2886751 -0.2886751  0.5 -0.2886751 -0.2886751  0.1666667
#
#$Env
#	  E1 - E2    E1 - E3    E1 - E4    E2 - E3    E2 - E4    E3 - E4
#E1  0.7071068  0.7071068  0.7071068  0.0000000  0.0000000  0.0000000
#E2 -0.7071068  0.0000000  0.0000000  0.7071068  0.7071068  0.0000000
#E3  0.0000000 -0.7071068  0.0000000 -0.7071068  0.0000000  0.7071068
#E4  0.0000000  0.0000000 -0.7071068  0.0000000 -0.7071068 -0.7071068
testInteractions(mod.temp.data, custom = custom.contr.4env);

#F Test: 
#		P-value adjustment method: holm
#					Value   Df Sum of Sq      F  Pr(>F)  
#  a1 : E1 - E2 -1.39784    1        85 6.8478 0.39174  
#  a2 : E1 - E2 -0.53235    1        12 0.9932 1.00000  
#  d1 : E1 - E2  1.23147    1       122 9.7931 0.08459 .
#  d2 : E1 - E2  0.31899    1         7 0.5566 1.00000  
#a1a2 : E1 - E2  0.10702    1         0 0.0401 1.00000  
#a1d2 : E1 - E2 -0.00476    1         0 0.0001 1.00000  
#d1a2 : E1 - E2  0.32684    1         9 0.6898 1.00000  
#d1d2 : E1 - E2 -0.04608    1         0 0.0215 1.00000  
#  a1 : E1 - E3  0.00000    1         0 0.0000 1.00000  
#  a2 : E1 - E3  0.00000    1         0 0.0000 1.00000  
#  d1 : E1 - E3  0.00000    1         0 0.0000 1.00000  
#  d2 : E1 - E3  0.00000    1         0 0.0000 1.00000  
#a1a2 : E1 - E3  0.00000    1         0 0.0000 1.00000  
#a1d2 : E1 - E3  0.00000    1         0 0.0000 1.00000  
#d1a2 : E1 - E3  0.00000    1         0 0.0000 1.00000  
#d1d2 : E1 - E3  0.00000    1         0 0.0000 1.00000  
#  a1 : E1 - E4 -1.39784    1        85 6.8478 0.39174  
#  a2 : E1 - E4 -0.53235    1        12 0.9932 1.00000  
#  d1 : E1 - E4  1.23147    1       122 9.7931 0.08459 .
#  d2 : E1 - E4  0.31899    1         7 0.5566 1.00000  
#a1a2 : E1 - E4  0.10702    1         0 0.0401 1.00000  
#a1d2 : E1 - E4 -0.00476    1         0 0.0001 1.00000  
#d1a2 : E1 - E4  0.32684    1         9 0.6898 1.00000  
#d1d2 : E1 - E4 -0.04608    1         0 0.0215 1.00000  
#  a1 : E2 - E3  1.39784    1        85 6.8478 0.39174  
#  a2 : E2 - E3  0.53235    1        12 0.9932 1.00000  
#  d1 : E2 - E3 -1.23147    1       122 9.7931 0.08459 .
#  d2 : E2 - E3 -0.31899    1         7 0.5566 1.00000  
#a1a2 : E2 - E3 -0.10702    1         0 0.0401 1.00000  
#a1d2 : E2 - E3  0.00476    1         0 0.0001 1.00000  
#d1a2 : E2 - E3 -0.32684    1         9 0.6898 1.00000  
#d1d2 : E2 - E3  0.04608    1         0 0.0215 1.00000  
#  a1 : E2 - E4  0.00000    1         0 0.0000 1.00000  
#  a2 : E2 - E4  0.00000    1         0 0.0000 1.00000  
#  d1 : E2 - E4  0.00000    1         0 0.0000 1.00000  
#  d2 : E2 - E4  0.00000    1         0 0.0000 1.00000  
#a1a2 : E2 - E4  0.00000    1         0 0.0000 1.00000  
#a1d2 : E2 - E4  0.00000    1         0 0.0000 1.00000  
#d1a2 : E2 - E4  0.00000    1         0 0.0000 1.00000  
#d1d2 : E2 - E4  0.00000    1         0 0.0000 1.00000  
#  a1 : E3 - E4 -1.39784    1        85 6.8478 0.39174  
#  a2 : E3 - E4 -0.53235    1        12 0.9932 1.00000  
#  d1 : E3 - E4  1.23147    1       122 9.7931 0.08459 .
#  d2 : E3 - E4  0.31899    1         7 0.5566 1.00000  
#a1a2 : E3 - E4  0.10702    1         0 0.0401 1.00000  
#a1d2 : E3 - E4 -0.00476    1         0 0.0001 1.00000  
#d1a2 : E3 - E4  0.32684    1         9 0.6898 1.00000  
#d1d2 : E3 - E4 -0.04608    1         0 0.0215 1.00000  
#Residuals               4758     59091                 
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

inter.geno.across.env <- testInteractions(mod.temp.data, custom = genotypeContrast, across = "Env");
inter.geno.across.env.col.length <- lengths(inter.geno.across.env)
inter.geno.across.env.extract <-inter.geno.across.env[,c("Df", "Sum of Sq", "F", "Pr(>F)")];
#F Test: 
#		P-value adjustment method: holm
#				Env1        Env2     Env3   Df Sum of Sq      F   Pr(>F)   
#  a1      -0.98842 -4.2966e-14 -0.98842    3       170 4.5652 0.023652 * 
#  a2      -0.37643 -1.3067e-13 -0.37643    3        25 0.6621 1.000000   
#  d1       1.50823 -1.1868e-13  1.50823    3       243 6.5287 0.001685 **
#  d2       0.39068 -8.0269e-14  0.39068    3        14 0.3710 1.000000   
#a1a2       0.07567  2.4980e-14  0.07567    3         1 0.0268 1.000000   
#a1d2      -0.00583  4.5630e-14 -0.00583    3         0 0.0001 1.000000   
#d1a2       0.40029  8.8707e-14  0.40029    3        17 0.4599 1.000000   
#d1d2      -0.09775  8.0269e-14 -0.09775    3         1 0.0143 1.000000   
#Residuals                               4758     59091                   
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# testing lmer model
library(lme4)
data <- rbind(data, data,data,data);

# mult site analysis
# fixed env
mod.lmer.temp.data.fixenv <- lmer(Yield ~ 1 + Genotype + Env + Genotype : Env + (1|Block : Env),data = data);
# if user do not submit env contrast
testInteractions(mod.lmer.temp.data.fixenv, custom = genotypeContrast, test="F");
#F Test: 
#		P-value adjustment method: holm
#			 Value  Df         F    Pr(>F)
#  a1      13.6934   1 34763.219 < 2.2e-16
#  a2       3.7939   1  2668.563 < 2.2e-16
#  d1      -1.2861   1   102.219 < 2.2e-16
#  d2      -1.2471   1    96.119 < 2.2e-16
#a1a2       1.1779   1   257.205 < 2.2e-16
#a1d2       1.3769   1   117.162 < 2.2e-16
#d1a2       1.3157   1   106.976 < 2.2e-16
#d1d2       1.4574   1    43.751 1.398e-10
#Residuals         350      
testInteractions(mod.lmer.temp.data.fixenv, custom = genotypeContrast, across = "Env", test="F");
#F Test: 
#		P-value adjustment method: holm
#			   Env1     Env2     Env3     Env4  Df      F Pr(>F)
#  a1      -0.41747 -0.22925  0.00959 -0.14257   4 1.1649      1
#  a2       0.04768  0.03374 -0.01308 -0.04079   4 0.0472      1
#  d1      -0.07897  0.38432 -0.23128 -0.33241   4 0.9411      1
#  d2      -0.35785 -0.04851 -0.04045 -0.42187   4 0.4916      1
#a1a2      -0.02305 -0.30141 -0.11923 -0.02786   4 0.5738      1
#a1d2       0.06506  0.24508 -0.46187 -0.00594   4 0.8426      1
#d1a2       0.16778  0.25461  0.17852  0.33903   4 0.1949      1
#d1d2       0.57356 -0.38124  0.11325  0.74524   4 0.8463      1
#Residuals                                     350  
# if user do submit env contrast
contr.5env <- list(Env = custom.contr.5env$Env);
testInteractions(mod.lmer.temp.data.fixenv, custom = contr.5env);
#Chisq Test: 
#		P-value adjustment method: holm
#          Value Df     Chisq Pr(>Chisq)
#E1 - E2  0.6171  1   31.7736  5.197e-08
#E1 - E3  2.5322  1  534.9335  < 2.2e-16
#E1 - E4  0.4242  1   15.0152  0.0002133
#E1 - E5 -0.9689  1   78.3128  < 2.2e-16
#E2 - E3  1.9151  1  305.9637  < 2.2e-16
#E2 - E4 -0.1929  1    3.1042  0.0780910
#E2 - E5 -1.5860  1  209.8518  < 2.2e-16
#E3 - E4 -2.1079  1  370.7045  < 2.2e-16
#E3 - E5 -3.5010  1 1022.5976  < 2.2e-16
#E4 - E5 -1.3931  1  161.9101  < 2.2e-16
# both geno and env contrast
testInteractions(mod.lmer.temp.data.fixenv, custom = custom.contr.5env);
#Chisq Test: 
#		P-value adjustment method: holm
#					Value Df  Chisq Pr(>Chisq)
#  a1 : E1 - E2 -0.26618  1 0.6568          1
#  a2 : E1 - E2  0.01971  1 0.0036          1
#  d1 : E1 - E2 -0.37828  1 1.3264          1
#  d2 : E1 - E2 -0.25257  1 0.5913          1
#a1a2 : E1 - E2  0.39367  1 1.4366          1
#a1d2 : E1 - E2 -0.14699  1 0.2003          1
#d1a2 : E1 - E2 -0.07090  1 0.0466          1
#d1d2 : E1 - E2  0.45009  1 1.8779          1
#  a1 : E1 - E3 -0.60394  1 3.3811          1
#  a2 : E1 - E3  0.08592  1 0.0684          1
#  d1 : E1 - E3  0.12436  1 0.1434          1
#  d2 : E1 - E3 -0.25915  1 0.6225          1
#a1a2 : E1 - E3  0.13602  1 0.1715          1
#a1d2 : E1 - E3  0.43023  1 1.7158          1
#d1a2 : E1 - E3 -0.00878  1 0.0007          1
#d1d2 : E1 - E3  0.21699  1 0.4365          1
#  a1 : E1 - E4 -0.38876  1 1.4010          1
#  a2 : E1 - E4  0.12511  1 0.1451          1
#  d1 : E1 - E4  0.20693  1 0.3969          1
#  d2 : E1 - E4  0.05227  1 0.0253          1
#a1a2 : E1 - E4  0.00680  1 0.0004          1
#a1d2 : E1 - E4  0.05797  1 0.0311          1
#d1a2 : E1 - E4 -0.13983  1 0.1812          1
#d1d2 : E1 - E4 -0.08093  1 0.0607          1
#  a1 : E1 - E5 -0.59039  1 3.2310          1
#  a2 : E1 - E5  0.06742  1 0.0421          1
#  d1 : E1 - E5 -0.06448  1 0.0385          1
#  d2 : E1 - E5 -0.29218  1 0.7914          1
#a1a2 : E1 - E5 -0.03260  1 0.0098          1
#a1d2 : E1 - E5  0.05312  1 0.0262          1
#d1a2 : E1 - E5  0.13699  1 0.1740          1
#d1d2 : E1 - E5  0.27038  1 0.6777          1
#  a1 : E2 - E3 -0.33776  1 1.0575          1
#  a2 : E2 - E3  0.06621  1 0.0406          1
#  d1 : E2 - E3  0.50264  1 2.3420          1
#  d2 : E2 - E3 -0.00658  1 0.0004          1
#a1a2 : E2 - E3 -0.25764  1 0.6153          1
#a1d2 : E2 - E3  0.57722  1 3.0885          1
#d1a2 : E2 - E3  0.06212  1 0.0358          1
#d1d2 : E2 - E3 -0.23310  1 0.5037          1
#  a1 : E2 - E4 -0.12258  1 0.1393          1
#  a2 : E2 - E4  0.10541  1 0.1030          1
#  d1 : E2 - E4  0.58521  1 3.1746          1
#  d2 : E2 - E4  0.30485  1 0.8615          1
#a1a2 : E2 - E4 -0.38687  1 1.3874          1
#a1d2 : E2 - E4  0.20495  1 0.3894          1
#d1a2 : E2 - E4 -0.06893  1 0.0440          1
#d1d2 : E2 - E4 -0.53102  1 2.6139          1
#  a1 : E2 - E5 -0.32421  1 0.9743          1
#  a2 : E2 - E5  0.04772  1 0.0211          1
#  d1 : E2 - E5  0.31380  1 0.9128          1
#  d2 : E2 - E5 -0.03961  1 0.0145          1
#a1a2 : E2 - E5 -0.42626  1 1.6843          1
#a1d2 : E2 - E5  0.20011  1 0.3712          1
#d1a2 : E2 - E5  0.20789  1 0.4006          1
#d1d2 : E2 - E5 -0.17972  1 0.2994          1
#  a1 : E3 - E4  0.21518  1 0.4292          1
#  a2 : E3 - E4  0.03920  1 0.0142          1
#  d1 : E3 - E4  0.08257  1 0.0632          1
#  d2 : E3 - E4  0.31142  1 0.8990          1
#a1a2 : E3 - E4 -0.12923  1 0.1548          1
#a1d2 : E3 - E4 -0.37227  1 1.2846          1
#d1a2 : E3 - E4 -0.13105  1 0.1592          1
#d1d2 : E3 - E4 -0.29792  1 0.8228          1
#  a1 : E3 - E5  0.01356  1 0.0017          1
#  a2 : E3 - E5 -0.01849  1 0.0032          1
#  d1 : E3 - E5 -0.18884  1 0.3306          1
#  d2 : E3 - E5 -0.03303  1 0.0101          1
#a1a2 : E3 - E5 -0.16862  1 0.2636          1
#a1d2 : E3 - E5 -0.37711  1 1.3183          1
#d1a2 : E3 - E5  0.14576  1 0.1970          1
#d1d2 : E3 - E5  0.05339  1 0.0264          1
#  a1 : E4 - E5 -0.20162  1 0.3768          1
#  a2 : E4 - E5 -0.05769  1 0.0308          1
#  d1 : E4 - E5 -0.27141  1 0.6829          1
#  d2 : E4 - E5 -0.34445  1 1.0998          1
#a1a2 : E4 - E5 -0.03940  1 0.0144          1
#a1d2 : E4 - E5 -0.00485  1 0.0002          1
#d1a2 : E4 - E5  0.27682  1 0.7103          1
#d1d2 : E4 - E5  0.35131  1 1.1440          1


# random env
mod.lmer.temp.data.ranenv <- lmer(Yield ~ 1 + Genotype + (1|Env) + (1|Genotype : Env) + (1|Block : Env),data = data);
testInteractions(mod.lmer.temp.data.ranenv, custom = genotypeContrast, test="F");
testInteractions(mod.lmer.temp.data.ranenv, custom = contr.5env ); # Not working
#Chisq Test: 
#		P-value adjustment method: holm
#		Value Df Chisq Pr(>Chisq)
#Mean 102.52  1 16034  < 2.2e-16
#Warning message:
#		In testInteractions(mod.lmer.temp.data.ranenv, custom = contr.5env) :
#		Some factors with specified contrasts are not in the model and will be ignored.
testInteractions(mod.lmer.temp.data.ranenv, custom = custom.contr.5env)
#Chisq Test: 
#		P-value adjustment method: holm
#		Value Df     Chisq Pr(>Chisq)
#  a1 27.3868  1 36072.762  < 2.2e-16
#  a2  7.5879  1  2769.089  < 2.2e-16
#  d1 -1.4851  1   106.070  < 2.2e-16
#  d2 -1.4401  1    99.739  < 2.2e-16
#a1a2  2.3557  1   266.894  < 2.2e-16
#a1d2  1.5899  1   121.575  < 2.2e-16
#d1a2  1.5192  1   111.005  < 2.2e-16
#d1d2  0.9716  1    45.400  1.607e-11
#Warning message:
#		In testInteractions(mod.lmer.temp.data.ranenv, custom = custom.contr.5env) :
#		Some factors with specified contrasts are not in the model and will be ignored.
# single site analysis
# no env
mod.lmer.temp.data.noenv <- lmer(Yield ~ 1 + Genotype + (1|Block), data = data);
testInteractions(mod.lmer.temp.data.noenv, custom = genotypeContrast, test="F");
#F Test: 
#		P-value adjustment method: holm
#					Value  Df        F    Pr(>F)    
#		  a1      13.6934   1 9355.646 < 2.2e-16 ***
#		  a2       3.7939   1  718.177 < 2.2e-16 ***
#		  d1      -1.2861   1   27.510 7.671e-07 ***
#		  d2      -1.2471   1   25.868 1.133e-06 ***
#		a1a2       1.1779   1   69.220 8.739e-15 ***
#		a1d2       1.3769   1   31.531 1.856e-07 ***
#		d1a2       1.3157   1   28.790 5.518e-07 ***
#		d1d2       1.4574   1   11.775  0.000664 ***
#		Residuals         394                       
#---
#		Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

testInteractions(lm(Yield ~ 1 + (Env / Genotype) , data = data), custom = genotypeContrast, across = "Env");

testInteractions(lmer(HD ~ 1 + Genotype + Env + Genotype : Env  , data = temp.data), custom = genotypeContrast, across = "Env");
model.site1 <- pyramidedline.ssa.outputs$output[[1]]$site[[1]]$model;
model.site2 <- pyramidedline.ssa.outputs$output[[1]]$site[[2]]$model;

testInteractions(model.site1, custom= genotypeContrast);
testInteractions(model.site2, custom= genotypeContrast);
model.msa <- pyramidedline.msa.outputs$output[[1]]$model;
testInteractions(model.msa, custom= genotypeContrast,across ="Env");

# Testing SSSL DEMO
sssl.data1.ssa <- sssl.ssa.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  data = sssl.data1);
sssl.without.rp <- sssl.ssa.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  data = sssl.data.without.rp);
sssl.with.na <- sssl.ssa.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  data = sssl.data.with.na);

contrastAnalysisOnIL.ssa(sssl.data1.ssa, contrastOpt = "RecurrentParent", recurrentParent = "HJX", acrossEnv = T);
sssl.ssa.without.rp <- contrastAnalysisOnIL.ssa(sssl.without.rp, contrastOpt = "RecurrentParent", recurrentParent = "HJX", acrossEnv = T);
sssl.ssa.with.na <- contrastAnalysisOnIL.ssa(sssl.with.na, contrastOpt = "RecurrentParent", recurrentParent = "HJX", acrossEnv = T);

# Testing PY DEMO

py.data.ssa.outputs <- pyramidedline.ssa.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  data = py.data5);

py.data.ssa.outputs.balanced <- pyramidedline.ssa.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  data = py.data.balanced);
py.data.ssa.outputs.unbalanced <- pyramidedline.ssa.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env", data = py.data.unbalanced);
py.data.ssa.outputs.balanced.3.genes <- pyramidedline.ssa.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env", data = py.data.balanced.3.genes);
py.data5.ssa.outputs <- pyramidedline.ssa.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  data = py.data5);


py.data.ssa.outputs.balanced.contrast <- contrastAnalysisOnIL.ssa(py.data.ssa.outputs.balanced, contrastOpt = "Custom", genoContrast = list(getDefaultGenesContrast()));
py.data.ssa.outputs.balanced.contrast.default <- contrastAnalysisOnIL.ssa(py.data.ssa.outputs.balanced, contrastOpt = "Default");
py.data.ssa.outputs.unbalanced.contrast <- contrastAnalysisOnIL.ssa(py.data.ssa.outputs.unbalanced, contrastOpt = "Custom", genoContrast = list(getDefaultGenesContrast(),getDefaultGenesContrast(2),getDefaultGenesContrast(2), getDefaultGenesContrast()),acrossEnv = F);
py.data.ssa.outputs.balanced.3.genes.contrast <- contrastAnalysisOnIL.ssa(py.data.ssa.outputs.balanced.3.genes, contrastOpt = "Default");
contrastAnalysisOnIL.ssa(py.data5.ssa.outputs, contrastOpt ="Default");


# --- MULTI SITE DEMO --- # 
# --- TESTING SSSL --- #
sssl.with.name.without.values <- rbind(sssl.data1, sssl.data6, sssl.data3);
sssl.with.name.without.values$Genotype <- factor(sssl.with.name.without.values$Genotype);
sssl.with.name.without.values$Env <- factor(sssl.with.name.without.values$Env)
sssl.with.name.without.values$Block <- factor(sssl.with.name.without.values$Block);
sssl.data1.msa <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  is.envFixed = FALSE, data = rbind(sssl.data1,sssl.data2));
sssl.data2.msa <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  is.envFixed = FALSE, data = sssl.with.name.without.values);

sssl.wit.name.without.values.msa.envFixed <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  is.envFixed = TRUE, data = sssl.with.name.without.values);
sssl.wit.name.without.values.msa.envRandom <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  is.envFixed = FALSE, data = sssl.with.name.without.values);

sssl.data.without.rp$Genotype <- factor(sssl.data.without.rp$Genotype);
sssl.data.without.rp$Env <- factor(sssl.data.without.rp$Env);
sssl.data.without.rp$Block <- factor(sssl.data.without.rp$Block);
sssl.without.rp.msa.envFixed <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env", is.envFixed = TRUE, data = sssl.data.without.rp);
sssl.without.rp.msa.envRandom <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env", is.envFixed = FALSE, data = sssl.data.without.rp);

sssl.with.na.msa <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno ="Genotype", block ="Block", env = "Env",  data = sssl.data.with.na);

# --- TESTING PYRAMIDEDLINE --- #
pl.balanced.msa.envrandom <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = FALSE, data = py.data.balanced);
pl.unbalanced.msa.envfixed <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = TRUE, data = py.data.unbalanced);
pl.unbalanced.msa.envrandom <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = FALSE, data = py.data.unbalanced);
pl.3genes.msa.envfixed <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = TRUE, data = py.data.balanced.3.genes);
pl.3genes.msa.envenvrandom <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = FALSE, data = py.data.balanced.3.genes);

# --- TESTING SSSL OF MSA --- #
library(lme4);
library(phia);
sssl.data.balanced.3env$Genotype <- factor(sssl.data.balanced.3env$Genotype);
sssl.data.balanced.3env$Env <- factor(sssl.data.balanced.3env$Env);
sssl.data.balanced.3env$Block <- factor(sssl.data.balanced.3env$Block);


sssl.envfixed.mod.with.names.without.values <- lmer(Yield ~ 1 + Genotype + Env + Genotype : Env + (1|Block : Env), data = sssl.with.name.without.values)
sssl.envfixed.mod <- lmer(Yield ~ 1 + Genotype + Env + Genotype : Env + (1|Block : Env), data = sssl.data.balanced.3env)
sssl.data.missing.2$Genotype <- factor(sssl.data.missing.2$Genotype)
sssl.data.missing.2$Env <- factor(sssl.data.missing.2$Env)
sssl.data.missing.2$Block <- factor(sssl.data.missing.2$Block)

sssl.envfixed.mod.missing.2 <- lmer(Yield ~ 1 + Genotype + Env + Genotype : Env + (1|Block : Env), data = sssl.data.missing.2)

sssl.envrandom.mod <- lmer(Yield ~ 1 + Genotype + (1|Env) + (1|Genotype : Env) + (1|Block : Env), data = sssl.with.name.without.values)
interactionMeans(sssl.envfixed.mod);
interactionMeans(sssl.envrandom.mod);
interactionMeans(sssl.envfixed.mod.missing.2 )
sssl.envrandom.mod.testInteractions.pairwise <- testInteractions(mod = sssl.envrandom.mod, pairwise = "Genotype",test = "F");
sssl.envfixed.mod.testInteractions.fixed <- testInteractions(mod = sssl.envfixed.mod, fixed = "Genotype");

# --- testing env is fixed --- #
sssl.envfixed.mod.testInteractions.pairwise <- testInteractions(mod = sssl.envfixed.mod, pairwise = "Genotype");
sssl.envfixed.mod.testInteractions.fixed <- testInteractions(mod = sssl.envfixed.mod, fixed = "Genotype");
sssl.contrast.21 <- diag(1,20,20);
sssl.contrast.21 <- cbind(-1,sssl.contrast.21);
colnames(sssl.contrast.21) <- c("HJX",levels(interaction("S", 1:20, sep="")));
rownames(sssl.contrast.21) <- interaction(c(levels(interaction("S", 1:20, sep=""))), " - HJX", sep =""); 
# -- testing custom contrast --- #
sssl.envfixed.mod.testInteractions.custom <- testInteractions(mod = sssl.envfixed.mod, custom = list(Genotype = as.data.frame(t(sssl.contrast.21))));
sssl.envfixed.mod.testInteractions.custom.acrossenv <- testInteractions(mod = sssl.envfixed.mod, custom = list(Genotype = as.data.frame(t(sssl.contrast.21))), across = "Env");

sssl.envrandom.mod.testInteractions <- testInteractions(mod = sssl.envrandom.mod, pairwise = "Genotype");
sssl.envrandom.mod.testInteractions.custom <- testInteractions(mod = sssl.envrandom.mod, custom = list(Genotype = as.data.frame(t(sssl.contrast.21))));
sssl.envrandom.mod.testInteractions.custom.acrossenv <- testInteractions(mod = sssl.envrandom.mod, custom = list(Genotype = as.data.frame(t(sssl.contrast.21))), across = "Env");

# -- testing function contrastAnalysisOnIL.msa --- #
sssl.envfixed.balanced.3env.msa <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno = "Genotype", block = "Block", env = "Env", is.envFixed = TRUE, data = sssl.data.balanced.3env );
sssl.envrandom.with.name.without.value.msa <- sssl.msa.onestage.test("RCB", respvar = c("Yield", "TN"), geno = "Genotype", block = "Block", env = "Env", is.envFixed = FALSE, data = sssl.data.with.name.without.values);

sssl.envfixed.balanced.3env.msa.contrast <- contrastAnalysisOnIL.msa(sssl.envfixed.balanced.3env.msa, contrastOpt="RecurrentParent", recurrentParent ="HJX");
sssl.envrandom.with.name.without.value.msa.contrast <- contrastAnalysisOnIL.msa(sssl.envrandom.with.name.without.value.msa, contrastOpt="RecurrentParent", recurrentParent ="HJX");

# --- testing recurrent parent options --- #
envContrast.1 <- matrix(c(1,-1,0,1,0,-1,0,1,-1), nrow = 3,byrow=T)
colnames(envContrast.1) <- c("E1","E3","E5")
rownames(envContrast.1) <- c("E1 - E3", "E1 - E5", "E3 - E5");
is.null(rownames(envContrast.1)) 
envContrast.1
sssl.envfixed.balanced.3env.msa.contrast.with.envcontrast <- contrastAnalysisOnIL.msa(sssl.envfixed.balanced.3env.msa, contrastOpt="RecurrentParent", recurrentParent ="HJX",envContrast = envContrast.1);


# --- testing Custom options on IL --- #
sssl.envfixed.balanced.3env.msa.contrast.with.envcontrast <-
		contrastAnalysisOnIL.msa(sssl.envfixed.balanced.3env.msa, 
				contrastOpt="Custom", genoContrast =sssl.contrast.21,envContrast = envContrast.1);

# --- testing Custom options on pyramided line ---- #
pl.balanced.msa.envfixed <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = TRUE, data = py.data.balanced);
pl.balanced.msa.envrandom <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = FALSE, data = py.data.balanced);
# No working # pl.unbalanced.msa.envfixed <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = TRUE, data = py.data.unbalanced);
pl.unbalanced.msa.envrandom <- pyramidedline.msa.onestage.test("RCB", respvar = c("Yield","TN"), geno ="Genotype", block="Block", env = "Env", is.envFixed = FALSE, data = py.data.unbalanced);
pl.balanced.msa.envfixed.contrast.without.envContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envfixed, contrastOpt="Custom", genoContrast=getDefaultGenesContrast());
envContrast.2 <- matrix(c(1,-1,0,1,0,-1,0,1,-1), nrow = 3,byrow=T)
colnames(envContrast.2) <- c("E1","E2","E3")
rownames(envContrast.2) <- c("E1 - E2", "E1 - E3", "E2 - E3");
pl.balanced.msa.envfixed.contrast.with.envContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envfixed, contrastOpt="Custom", genoContrast=getDefaultGenesContrast(),envContrast=envContrast.2);
pl.balanced.msa.envrandom.contrast.without.envContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envrandom, contrastOpt="Custom", genoContrast=getDefaultGenesContrast());
pl.balanced.msa.envrandom.contrast.with.envContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envrandom, contrastOpt="Custom", genoContrast=getDefaultGenesContrast(),envContrast=envContrast.2);

# --- testing default options on contrastAnalysisOnIL.msa 
# --- by far, it is only working on pyramidedline multi-site analysis outcomes
pl.balanced.msa.envfixed.contrast.without.envContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envfixed, contrastOpt="Default");
pl.balanced.msa.envfixed.contrast.with.envContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envfixed, contrastOpt="Default",envContrast=envContrast.2);
pl.unbalanced.msa.envrandom.contrast.without.envContrast <- contrastAnalysisOnIL.msa(pl.unbalanced.msa.envrandom, contrastOpt="Default")
pl.unbalanced.msa.envrandom.contrast.with.envContrast <- contrastAnalysisOnIL.msa(pl.unbalanced.msa.envrandom, contrastOpt="Default",envContrast=envContrast.2)
# No working # sssl.envfixed.balanced.3env.msa.contrast.with.envcontrast <- contrastAnalysisOnIL.msa(sssl.envfixed.balanced.3env.msa, contrastOpt="Default",envContrast = envContrast.1);
pl.balanced.msa.envrandom.contrast.without.envContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envrandom, contrastOpt="Default");
pl.balanced.msa.envrandom.contrast.with.envContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envrandom, contrastOpt="Default",envContrast=envContrast.2);
pl.balanced.msa.envfixed.contrast.with.genoAndEnvContrast <- contrastAnalysisOnIL.msa(pl.balanced.msa.envfixed, contrastOpt="Default",genoContrast=getDefaultGenesContrast(2),envContrast=envContrast.2);











