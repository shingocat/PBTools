###############################################################################
# TODO: MULTIPLE SITE ANALYSIS (ONE-STAGE) for Pyramided Line
# This function performs multi-environment one stage test for SSSL
#
# ARGUMENT:
# exptl.design - RCB, AugRCB, AugLS, Alpha, RowCol, LatinAlpha, LatinRowCol
# respvar - a vector of strings; variable names of the response variables
# geno - a string; variable name of the treatment/genotype variable
# block - a string; variable name of the block variable
# row - a string; variable name of the row variable
# column - a string; variable name of the column variable; NULL, if design is RCB, Alpha, LatinAlpha
# rep - a string; variable name of the replication variable; NULL, if design is RCB
# env - a string; variable name of the environment variable
# is.envFixed - logical; indicating whether environment factor is random or not; default value TRUE
# data - a string; name of the dataframe
#
# Author: mqin
# Date: Sep 17, 2014
# FileName: sssl.msa.onestage.R
###############################################################################

sssl.msa.onestage.test <- function(
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),
		respvar, 
		geno, 
		block = NULL,
		row = NULL, 
		column = NULL, 
		rep = NULL, 
		env, 
		is.envFixed = TRUE,
		data
	) UseMethod("sssl.msa.onestage.test");
sssl.msa.onestage.test.default <- function(
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),
		respvar, 
		geno, 
		block = NULL,
		row = NULL, 
		column = NULL, 
		rep = NULL, 
		env, 
		is.envFixed = TRUE,
		data
)
{
	library(lme4);
	library(lsmeans);
	options(show.signif.stars = FALSE);
	prev.opt <- options()$warn;
	#options(warn = -1); # -- By far not supress warnings
	
	# -- checking if columns specified are in the data set, if yes, factor all the factors
	if(sum(is.na(match(respvar,names(data))), na.rm = TRUE) > 0)
	{
		stop("At least one responsive variable does not match a column in the data frame.");
	}
	if(is.na(match(geno, names(data))))
	{
		stop("Genotypic variable does not match a column in the data frame.");
	} else
	{
		data[,match(geno, names(data))] <- factor(trimStrings(data[,match(geno, names(data))]));
		if(nlevels(data[,match(geno, names(data))]) <= 1)
		{
			stop("\tGenotypic variable cannot be less than two levels!\n");
		}
	}
	if(is.na(match(env, names(data))))
	{
		stop("Environmental variable does not match a column in the data frame.");
	}else
	{
		data[,match(env, names(data))] <- factor(trimStrings(data[,match(env, names(data))]));
		if(nlevels(data[,match(env, names(data))]) <= 1)
		{
			stop("\tEnvironmental variable cannot be less than two levels!\n");
		}
	}
	
	if(exptl.design == "RCB" || exptl.design == "AugRCB")
	{
		if(is.na(match(block, names(data))))
		{
			stop("The variable of Block does not match a column in the data frame.");
		} else
		{
			data[,match(block,names(data))] <- factor(trimStrings(data[,match(block,names(data))]));
			if(nlevels(data[,match(block,names(data))]) <= 1)
			{
				stop("\tThe variable of Block cannot be less than two levels!\n");
			}
		}
	} else if(exptl.design == "AugLS")
	{
		if(is.na(match(row, names(data))))
		{
			stop("Row variable does not match a column in the data frame.");
		}else
		{
			data[,match(row, names(data))] <- factor(trimStrings(data[,match(row, names(data))]));
			if(nlevels(data[,match(row, names(data))]) <= 1)
			{
				stop("\tThe variable of Row cannot be less than two levels!\n");
			}
		}
		if(is.na(match(column, names(data))))
		{
			stop("Column variable does not match a column in the data frame.");
		} else
		{
			data[,match(column, names(data))] <- factor(trimStrings(data[,match(column, names(data))]));
			if(nlevels(data[,match(column, names(data))]) <= 1)
			{
				stop("\tThe variable of Column cannot be less than two levels!\n");
			}
		}
	} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
	{
		if(is.na(match(block, names(data))))
		{
			stop("Block variable does not match a column in the data frame.");
		} else
		{
			data[,match(block, names(data))] <- factor(trimStrings(data[,match(block, names(data))]));
			if(nlevels(data[,match(block, names(data))]) <= 1)
			{
				stop("\tThe variable of Block cannot be less than two levels!\n");
			}
		}
		if(is.na(match(rep, names(data))))
		{
			stop("Rep variable does not match a column in the data frame.");
		} else
		{
			data[,match(rep, names(data))] <- factor(trimStrings(data[,match(rep, names(data))]));
			if(nlevels(data[,match(rep, names(data))]) <= 1)
			{
				stop("\tThe variable of Rep cannot be less than two levels!\n");
			}
		}
	} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
	{
		if(is.na(match(row, names(data))))
		{
			stop("Row variable does not match a column in the data frame.");
		} else
		{
			data[,match(row, names(data))] <- factor(trimStrings(data[,match(row, names(data))]));
			if(nlevels(data[,match(row, names(data))]) <= 1)
			{
				stop("\tThe variable of Row cannot be less than two levels!\n");
			}
		}
		if(is.na(match(column, names(data))))
		{
			stop("Column variable does not match a column in the data frame.");
		} else
		{
			data[,match(column, names(data))] <- factor(trimStrings(data[,match(column, names(data))]));
			if(nlevels(data[,match(column, names(data))]) <= 1)
			{
				stop("\tThe variable of Column cannot be less than two levels!\n");
			}
		}
		if(is.na(match(rep, names(data))))
		{
			stop("Rep variable does not match a column in the data frame.");
		} else
		{
			data[,match(rep, names(data))] <- factor(trimStrings(data[,match(rep, names(data))]));
			if(nlevels(data[,match(rep, names(data))]) <= 1)
			{
				stop("\tThe variable of Rep cannot be less than two levels!\n");
			}
		}
	}
	
	# checking whether missing whole levels of this respvar
	for(i in (1 : length(respvar)))
	{
		if(countMissingLevelsByGroupFactor(respvar[i], groupfactor = c(geno, env), data = data) && is.envFixed)
		{
			stop(paste("\tError: There are missing one levels of combination between Genotype and Env on " 
							, respvar[i] , " when the env factor is setted to fixed.",
							" Please checking your data!\n", sep = ""));
		}
	}
	
	result <- list();
	for(i in (1 : length(respvar)))
	{
		result[[i]] <- list();
		result[[i]]$respvar <- respvar[i];
		
		
		# --- create temp.data without missing observations --- #
		temp.data <- subset(data, subset = (is.na(data[,match(respvar[i], names(data))]) == FALSE));
		temp.data[,match(geno, names(temp.data))] <- factor(trimStrings(temp.data[,match(geno, names(temp.data))]));
		temp.data[,match(env, names(temp.data))] <- factor(trimStrings(temp.data[,match(env,names(temp.data))]));
		
		# --- get levles of genotype and environment --- #
		levelsGeno <- levels(temp.data[,match(geno, names(temp.data))]);
		levelsEnv <- levels(temp.data[,match(env, names(temp.data))]);
		
		result[[i]]$nlevelsGeno <- length(levelsGeno);
		result[[i]]$nlevelsEnv <- length(levelsEnv);
		
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
		
		# --- compute number of observations read, used and response rate --- #
		obsread <- nrow(data);
		obsused <- nrow(temp.data);
		result[[i]]$obsread <- obsread;
		result[[i]]$obsused <- obsused;
		responseRate <- (obsused/obsread);
		result[[i]]$responseRate <- responseRate;
		
		if(responseRate < 0.80)
		{
			result[[i]]$manyNAWarning <- "Too many missing observation. Canot proceed with the analysis.";
			warning("Too many missing observation. Canot proceed with the analysis.");
			next;
		} else
		{
			# --- compute summary statistics per environment --- #
			sumStat.Env <- DescriptiveStatistics(temp.data, respvar[i], env, c("min", "mean", "max", "var","sd"));
			sumStat.Env <- sumStat.Env[,c(2:ncol(sumStat.Env))];
			
			# --- CONSTRUCT THE MODEL --- #
			if(is.envFixed) 
			{
				env.stmt <- paste(env," + ", geno, ":", env, sep ="");
			}else
			{
				env.stmt <- paste("(1|", env, ") + (1|",geno, ":", env, ")", sep="");
			}
			if(exptl.design == "RCB" || exptl.design =="AugRCB")
			{
				myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , block , ":", env,")", sep="");
			} else if(exptl.design == "AugLS")
			{
				myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , row , ":", env, ") + (1|", column, ":", env, ")", sep="");
			} else if(exptl.design == "Alpha")
			{
				myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , rep , ":", env,") + (1|", rep, ":", block,":",env,")", sep = "" );
			} else if(exptl.design == "RowCol")
			{
				myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + " , env.stmt , " + (1|", rep , ":", env,") + (1|", rep,":", row,":", env,") + (1|", rep,":", column,":", env,")", sep = "");
			} else if(exptl.design == "LatinAlpha")
			{
				myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , rep , ":", env, ") + (1|", block, ":", env,") + (1|", rep, ":", block, ":", env,")", sep = "");
			} else if(exptl.design == "LatinRowCol")
			{
				if(longerRow)
				{
					myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + ", env.stmt , " + (1|" , rep, ":", env, ") + (1|", column, ":", env, ") + (1|", rep, ":", column, ":", env, ") + (1|", row, ":", env,")", sep = "");
				} else
				{
					myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + ", env.stmt , " + (1|" , rep, ":", env, ") + (1|", row, ":", env, ") + (1|", rep, ":" + row, ":", env,") + (1|", column, ":" , env, ")", sep = "");
				}
			}
			
			model <- lmer(formula(myformula1), data = temp.data);
			result[[i]]$formula1 <- myformula1;
			result[[i]]$model <- model;
			result[[i]]$envfixed <- is.envFixed;
			
			# --- VARIANCE COMPONENTS --- #
			varcomp <- NULL;
			for(j in (1:length(VarCorr(model))))
			{
				varcomp <- rbind(varcomp, 
						data.frame(
								Groups = names(VarCorr(model))[j], 
								Variance = VarCorr(model)[[j]][1], 
								Std.Dev. = attr(VarCorr(model)[[j]], "stddev")[[1]]
						)
				);
			}
			varcomp <- rbind(varcomp, 
					data.frame(
							Groups = "Residual", 
							Variance = attr(VarCorr(model), "sc")^2, 
							Std.Dev. = attr(VarCorr(model), "sc")
					)
			);
			result[[i]]$varcomp.table <- varcomp;
			
			# --- TEST FOR SIGNIFICANCE OF GENOTYPIC EFFECT USING LRT && ANOVA TABLE IF GENO IS FIXED --- #
			
			# --- myformula2 is full model minus the genotype term --- #
			myformula2 <- sub(paste(" + ", geno, sep = ""), "", myformula1, fixed = TRUE); ## There is one question on here, whether it should be remove the Genotype:Env term when remove the Genotype term in the model
			model1 <- lmer(formula(myformula1), data = temp.data, REML = T);
			model2 <- lmer(formula(myformula2), data = temp.data, REML = T);
			result[[i]]$formula2 <- myformula2;
			
			# --- This is for genotype as Fixed --- #
			anova.table1 <- anova(model2, model1);
			rownames(anova.table1)<-c("Model2", "Model1");
			attr(anova.table1, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF GENOTYPIC EFFECT USING LIKELIHOOD RATIO TEST:\n", sep = "");
			attr(anova.table1, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
			attr(anova.table1, "heading")[3] <- paste("Formula for Model2: ", myformula2, sep = "");
			attr(anova.table1, "heading")[4] <- paste("", sep = "");
			result[[i]]$testsig.Geno <- anova.table1;
			
			
			# --- TEST FOR SIGNIFICANCE OF ENVIRONMENT EFFECT USING LRT IF ENV IS FIXED--- #
			
			# --- myformula3 is full model minus the environment term --- #
			if(is.envFixed)
			{
				myformula3 <- gsub(paste(" + ", env, sep=""), "", myformula1, fixed = TRUE);
				model3 <- lmer(formula(myformula3), data = temp.data, REML = TRUE);
				result[[i]]$formula3 <- myformula3;
				
				anova.table2 <- anova(model3, model1);
				rownames(anova.table2) <- c("Model3", "Model1");
				attr(anova.table2, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF ENVIRONMENTAL EFFECT USING LIKELIHOOD RATIO TEST:\n", sep ="");
				attr(anova.table2, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
				attr(anova.table2, "heading")[3] <- paste("Formula for Model3: ", myformula3, sep = "");
				attr(anova.table2, "heading")[4] <- paste("", sep="");
				result[[i]]$testsig.Env <- anova.table2;
				
			} else
			{
				myformula3 <- gsub(paste(" + (1|", env, ")", sep = ""), "", myformula1, fixed = TRUE);
				model3 <- lmer(formula(myformula3), data = temp.data, REML = TRUE);
				result[[i]]$formula3 <- myformula3;
				
				models.table2 <- modelComparisonTable(model1, model3);
				result[[i]]$testsig.Env <- models.table2;
			}
			
			# --- TEST OF SIGNIFICANCE OF GENOTYPE X ENVIRONMENT EFFECT USING LRT --- #
			
			# --- myformula4 is full model minus the geno by environment interaction term --- #
			if(is.envFixed)
			{
				myformula4 <- gsub(paste(" + ", geno, ":", env, sep = ""), "", myformula1, fixed = TRUE);
				model4 <- lmer(formula(myformula4), data = temp.data, REML = TRUE);
				result[[i]]$formula4 <- myformula4;
				
				anova.table3 <- anova(model4, model1);
				rownames(anova.table3) <- c("Model4", "Model1");
				attr(anova.table3, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF GENOTYPE BY ENVIRONMENT EFFECT USING LIKELIHOOD RATIO TEST:\n", sep="");
				attr(anova.table3, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
				attr(anova.table3, "heading")[3] <- paste("Formula for Model4: ", myformula4, sep = "");
				attr(anova.table3, "heading")[4] <- paste("", sep = "");
				result[[i]]$testsig.GenoEnv <- anova.table3;
				
			} else 
			{
				myformula4 <- gsub(paste(" + (1|", geno, ":", env, ")", sep = ""), "", myformula1, fixed = TRUE);
				model4 <- lmer(formula(myformula4), data = temp.data, REML = TRUE);
				
				models.table3<-modelComparisonTable(model1, model4);
				
				result[[i]]$formula4 <- myformula4;
				result[[i]]$testsig.GenoEnv <- models.table3;
			}
			
			# --- PREDICTED MEANS/LSMEANS OF GENOTYPE AND SED STATISTICS--- #
			# --- myformula5 is full model but without intercept --- #
			myformula5 <- gsub("~ 1", "~ 0", myformula1, fixed = TRUE);
			model.noint <- lmer(formula(myformula5), data = temp.data);
			# --- Modify by QIN MAO FOR USING LSMENAS PACKAGE TO GET LSMEANS --- #
			sumStat.Geno <- lsmeans(model, geno, weights = "cells");
			sumStat.Geno <- summary(sumStat.Geno)[,1:3];
			colnames(sumStat.Geno) <- c(geno, "LSMean", "StdErrMean");
			rownames(sumStat.Geno) <- NULL;
			
			# --- display standard error of the differences --- #
			# ??? how to caculated standard error of the differences using with intercept model ??? #
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
			result[[i]]$sedTable <- sedTable;
			
			# --- ESTIMATES OF EFFECTS --- #
			if(is.envFixed)
			{
				# --- ENVIRONMENT LSMEANS --- #
				sumStat.Env <- lsmeans(model1,env);
				sumStat.Env <- summary(sumStat.Env)[,1:3];
				colnames(sumStat.Env) <- c(env, "LSMean", "StdErrMean");
				rownames(sumStat.Env) <- NULL;
				
				# --- G X E MEANS --- #
				# IF ENVIRONMENT IS FIXED, GXE LSMEANS COULD BE RETRIEVED DIRECTLY BY LSMEANS;
				sumStat.GenoEnv <- lsmeans(model1, c(geno, env));
				sumStat.GenoEnv <- summary(sumStat.GenoEnv)[,1:4];
				colnames(sumStat.GenoEnv) <- c(geno, env, paste(respvar[i],"LSMean",sep ="_"), "StdErrMean");
				rownames(sumStat.GenoEnv) <- NULL;
			} else
			{
				# -- GENOTYPE EFFECT --- #
				intercept <- fixef(model)[[1]];
				geno.effect <- as.data.frame(fixef(model)[-1]);
				geno.effect <- data.frame(gsub(geno, "", rownames(geno.effect)), geno.effect);
				colnames(geno.effect) <- c(geno, "geno_effect")
				rownames(geno.effect) <- NULL
				
				# -- ENVIRONMENT EFFECT--- #
				env.effect <- eval(parse(text = paste("ranef(model)$", env, sep = "")));
				env.effect <- data.frame(gsub(env, "", rownames(env.effect)), env.effect);
				colnames(env.effect) <- c(env, "env_effect");
				rownames(env.effect) <- NULL;
				
				# -- G X E EFFECT --- #
				GXE.effect <- as.data.frame(eval(parse(text = paste("ranef(model)$'", geno, ":", env, "'",sep = ""))));
				names <- t(as.data.frame(strsplit(rownames(GXE.effect), ":")));
				GXE.effect <- data.frame(names[,1], names[,2], GXE.effect+intercept);
				colnames(GXE.effect) <- c(geno, env, "ge_effect");
				rownames(GXE.effect) <- NULL;
				
				# -- G X E MEANS --- #
				sumStat.GenoEnv <- merge(GXE.effect, env.effect, by = env, all = TRUE);
				sumStat.GenoEnv <- merge(sumStat.GenoEnv, geno.effect, by = geno, all = TRUE);
				sumStat.GenoEnv <- data.frame(sumStat.GenoEnv[,match(geno, names(sumStat.GenoEnv))], sumStat.GenoEnv[,match(env,names(sumStat.GenoEnv))], rowSums(subset(sumStat.GenoEnv, select = c(ge_effect, env_effect, geno_effect)), na.rm = TRUE));
				colnames(sumStat.GenoEnv) <- c(geno, env, paste(respvar[i], "LSMean", sep = "_"));
			}
			
			# -- create G x E MEANS with coded levels -- #
			sumStat.GenoEnvCode <- sumStat.GenoEnv;
			sumStat.GenoEnvCode$CodedGeno <- sumStat.GenoEnvCode[,match(geno, names(sumStat.GenoEnvCode))];
			sumStat.GenoEnvCode$CodedEnv <- sumStat.GenoEnvCode[,match(env, names(sumStat.GenoEnvCode))];
			
			# --- display G x E means in 2-way table --- #
			#wide.GenoEnv<-ToWide(sumStat.GenoEnv, paste(respvar[i], "means", sep = "_"), env, geno)
			wide.GenoEnv <- reshape(
					sumStat.GenoEnv[,1:3], 
					v.names=paste(respvar[i], "LSMean", sep = "_"), 
					timevar=env, 
					idv=geno, 
					direction="wide"
			);
			colnames(wide.GenoEnv) <- gsub(
					paste(respvar[i], "LSMean.", sep = "_"), 
					"", 
					colnames(wide.GenoEnv));
			rownames(wide.GenoEnv) <- 1:nrow(wide.GenoEnv);
			
			result[[i]]$means.Geno <- sumStat.Geno;
			result[[i]]$means.Env <- sumStat.Env;
			result[[i]]$means.GenoEnv <- sumStat.GenoEnv;
			result[[i]]$wide.GenoEnv <- wide.GenoEnv;
			result[[i]]$means.GenoEnvCode <- sumStat.GenoEnvCode;
			result[[i]]$residuals <- resid(model1);
			result[[i]]$fitted.values <- fitted(model1);
			result[[i]]$data <- temp.data;
			
			# --- if genotype is fixed, output MSE and harmonic mean for AMMI --- #
			# --- compute harmonic mean per environment level --- #
			envgenorep <- as.data.frame.table(tapply(temp.data[, respvar[i]], temp.data[,c(env, geno)], length));
			envgenorep <- envgenorep[(is.na(envgenorep[,"Freq"]) == FALSE),];
			envgenorep$reciprocal <- 1/envgenorep$Freq;
			envgenorep2 <- as.data.frame.table(tapply(envgenorep[, "reciprocal"], envgenorep[,env], mean));
			envgenorep2$reciprocal2 <- 1/envgenorep2$Freq;
			envgenorep3 <- merge(envgenorep, envgenorep2, by.x=env, by.y="Var1");
			
			numrepEnv <- tapply(envgenorep3[,"reciprocal2"], envgenorep3[,env], mean);
			no.reps <- 1/mean(1/numrepEnv);
			
			result[[i]]$harmonicMean <- no.reps;
			result[[i]]$MSE <- varcomp[varcomp[,1] == "Residual", "Variance"];
			
		} # end of statement if(responseRate < 0.80)
		
	} # end of statement for(i in (1 : length(respvar)))	
	
	# --- consolidate means and residuals --- #
	
	for (m in 1:length(respvar)) {
		# --- consolidate genotype x environment means --- #
		if (m==1) { 
			means.GenoEnv.all <- as.data.frame(result[[m]]$means.GenoEnv);
		} else {
			newGenoEnv <- as.data.frame(result[[m]]$means.GenoEnv);
			if (nrow(newGenoEnv) > 0) {
				if (nrow(means.GenoEnv.all) == 0) { 
					means.GenoEnv.all <- newGenoEnv;
				} else { 
					means.GenoEnv.all <- merge(means.GenoEnv.all, newGenoEnv, by=c(geno, env), all=TRUE); 
				}
				
			}
		}
		
		# --- consolidate genotype means --- #
		if (m==1) {
			newGeno <- as.data.frame(result[[m]]$means.Geno);
			if (nrow(newGeno) > 0) { 
				colnames(newGeno) <- c(
						geno, 
						paste(respvar[m], "_LSMean", sep=""), 
						paste(respvar[m], "_StdErrMean", sep="")
				); 
			}
			means.Geno.all <- newGeno;
		} else {
			newGeno <- as.data.frame(result[[m]]$means.Geno);
			if (nrow(newGeno) > 0) {
				colnames(newGeno) <- c(
						geno, 
						paste(respvar[m], "_LSMean", sep=""), 
						paste(respvar[m], "_StdErrMean", sep=""));
				if (nrow(means.Geno.all) == 0) { 
					means.Geno.all <- newGeno;
				} else { 
					means.Geno.all <- merge(means.Geno.all, newGeno, by=c(geno), all=TRUE);
				}
			}
		}
	}
	
	# generate status of means.GenoEnv.all
	if (nrow(means.GenoEnv.all) == 0) {
		meansGenoEnvWarning<-"empty";
	} else { 
		meansGenoEnvWarning<-"not empty"; 
	}
	
	# generate status of means.Geno.all
	if (nrow(means.Geno.all) == 0) {
		meansGenoWarning<-"empty";
	} else { 
		meansGenoWarning<-"not empty"; 
	}
	
	options(warn = prev.opt);
	detach("package:lme4");
	detach("package:lsmeans");
	return(list(
					output = result, 
					means.GenoEnv.all=means.GenoEnv.all, 
					meansGenoEnvWarning=meansGenoEnvWarning, 
					means.Geno.all=means.Geno.all, 
					meansGenoWarning=meansGenoWarning,
					analysisType = "MultiSiteAnalysis",
					populationType = "SSSL"
			)
	)
} # end of function






