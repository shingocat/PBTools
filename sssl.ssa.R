###############################################################################
# TODO: SSS Single Site Analysis, 

# --------------------------------------------------------

# --------------------------------------------------------
# ARGUMENTS:
# exptl.design = RCB, AugRCB, AugLS, Alpha, RowCol, LatinAlpha, LatinRowCol
# respvar - a string; variable name of the response variable
# geno - a string; variable name of the treatment/genotype variable
# block - a string; variable name of the blocking variable
# row - a string; variable name of the row variable
# column - a string; variable name of the column variable
#        - NULL, if design is RCB, Aug RCB, Alpha Lattice
# rep - a string; variable name of the replication variable
#       - NULL, if design is RCB, Aug RCB, Aug LS, 
# env - a string; variable name of the environment variable
# data - a string; name of the dataframe
# --------------------------------------------------------

# --------------------------------------------------------
# Author: mqin
# Date: Sep 4, 2014
# FileName: pyramidedline.ssa.R
###############################################################################


sssl.ssa.test <- function(
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),	
		respvar,
		geno,
		block,
		row = NULL,
		column = NULL,
		rep = NULL,
		env = NULL,
		data
)
{
	UseMethod("sssl.ssa.test");
};

sssl.ssa.test.default <- function(
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),
		respvar,
		geno,
		block = NULL,
		row = NULL,
		column = NULL,
		rep = NULL,
		env = NULL,
		data
)
{
	options(show.signif.stars=FALSE);
	library(lme4);
	
	# --- check if columns specified are in the data set --- #	
	if(sum(is.na(match(respvar, names(data))), na.rm=TRUE) > 0) # if there are at least one variable name not match data frame name, after sum, it will lager than zero
	{
		stop("At least one resonpse variable name does not match a column in the data frame.");
	}
	if(is.na(match(geno, names(data))))
	{
		stop("Genotype variable name does not match a column in the data frame.");
	}
	if(!(is.null(env)) && is.na(match(env, names(data))))
	{
		stop("Environment variable name does not match a column in the data frame.");
	}
	if(exptl.design == "RCB" || exptl.design == "AugRCB")
	{
		if(block != NULL && is.na(match(block, names(data))))
		{
			stop("Block variable names does not match a column in the data frame.");
		}
	} else if(exptl.design == "AugLS")
	{
		if(row != NULL && is.na(match(row, names(data))))
		{
			stop("Row variable names does not match a column in the data frame.");
		}
		if(column != NULL && is.na(match(column, names(data))))
		{
			stop("Column variable names does not match a column in the data frame.");
		}
	} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
	{
		if(block != NULL && is.na(match(block, names(data))))
		{
			stop("Block variable names does not match a column in the data frame.");
		}
		if(rep != NULL && is.na(match(rep, names(data))))
		{
			stop("Rep variable names does not match a column in the data frame.");
		}
	} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
	{
		if(row != NULL && is.na(match(row, names(data))))
		{
			stop("Row variable names does not match a column in the data frame.");
		}
		if(column != NULL && is.na(match(column, names(data))))
		{
			stop("Column variable names does not match a column in the data frame.");
		}
		if(rep != NULL && is.na(match(rep, names(data))))
		{
			stop("Rep variable names does not match a column in the data frame.");
		}
	}
	
	# --- if env column is not specified, create EnvLevel column --- #
	if (is.null(env)) {
		env = "EnvLevel"
		data <- cbind(data, EnvLevel=1);
	}
	
	# --- set environment to factor --- #
	data[,match(env, names(data))] <- factor(trimStrings(data[ ,match(env, names(data))]));
	
	
	# --- save all output to this list. --- #
	result <- list();
	for(i in (1 : length(respvar)))
	{
		result[[i]] <- list();
		result[[i]]$respvar <- respvar[i];
		for( j in (1 : nlevels(data[,match(env, names(data))])))
		{
			# --- create temp.data with one respvar only --- #
			temp.data <- data;
			
			result[[i]]$site[[j]] <- list();
			result[[i]]$site[[j]]$env <- levels(temp.data[ ,match(env, names(temp.data))])[j];
			
			# --- create temp.data with one environment level only --- #
			temp.data <- temp.data[temp.data[ ,match(env, names(temp.data))] == levels(temp.data[ ,match(env, names(temp.data))])[j], ];	
			
			# --- count number of observation read and used --- #
			obsread <- nrow(temp.data);
			result[[i]]$site[[j]]$obsread <- obsread;
			temp.data <- subset(temp.data, subset = (is.na(temp.data[, match(respvar[i], names(temp.data))]) == FALSE));
			obsused <- nrow(temp.data);
			result[[i]]$site[[j]]$obsused <- obsused;
			
			# --- define all factors --- #
			temp.data[ ,match(geno, names(temp.data))] <- factor(trimStrings(temp.data[ ,match(geno, names(temp.data))]));
			if(nlevels(temp.data[ ,match(geno, names(temp.data))]) <= 1)
			{
				stop("\tThe genotypic variable cannot be less than two levels!\n");
			}
			if(exptl.design == "RCB" || exptl.design == "AugRCB")
			{
				temp.data[ ,match(block, names(temp.data))] <- factor(trimStrings(temp.data[,match(block, names(temp.data))]));
				if(nlevels(temp.data[ ,match(block, names(temp.data))]) <= 1)
				{
					stop("\tThe variable of Block cannot be less than two levels!\n");
				}
			} else if(exptl.design == "AugLS")
			{
				temp.data[ , match(row, names(temp.data))] <- factor(trimStrings(temp.data[,match(row, names(temp.data))]));
				if(nlevels(temp.data[ , match(row, names(temp.data))]) <= 1)
				{
					stop("\tThe variable of Row cannot be less than two levels!\n");
				}
				temp.data[ , match(column, names(temp.data))] <- factor(trimStrings(temp.data[,match(column,names(temp.data))]));
				if(nlevels(temp.data[ , match(column, names(temp.data))]) <= 1)
				{
					stop("\tThe variable of Column cannot be less than two levels!\n");
				}
			} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
			{
				temp.data[ , macth(block, names(temp.data))] <- factor(trimStrings(temp.data[,match(block, names(temp.data))]));
				if(nlevels(temp.data[ , macth(block, names(temp.data))]) <= 1)
				{
					stop("\tThe variable of Block cannot be less than two levels!\n");
				}
				temp.data[ , match(rep, names(temp.data))] <- factor(trimStrings(temp.data[,match(rep, names(temp.data))]));
				if(nlevels(temp.data[ , match(rep, names(temp.data))]) <= 1)
				{
					stop("\tThe variable of Rep cannot be less than two levels!\n");
				}
			} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
			{
				temp.data[ , match(row, names(temp.data))] <- factor(trimStrings(temp.data[,match(row, names(temp.data))]));
				if(nlevels(temp.data[ , match(row, names(temp.data))]) <= 1)
				{
					stop("\tThe variable of Row cannot be less than two levels!\n");
				}
				temp.data[ , match(column, names(temp.data))] <- factor(trimStrings(temp.data[,match(column,names(temp.data))]));
				if(nlevels(temp.data[ , match(column, names(temp.data))]) <= 1)
				{
					stop("\tThe variable of Column cannot be less than two levels!\n");
				}
				temp.data[ , match(rep, names(temp.data))] <- factor(trimStrings(temp.data[,match(column,names(temp.data))]));
				if(nlevels(temp.data[ , match(rep, names(temp.data))]) <= 1)
				{
					stop("\tThe variable of Rep cannot be less than two levels!\n");
				}
			}
			
			# --- checking geno factor should have correct coded labels and levels --- #
			genoLevels <- levels(temp.data[,match(geno, names(temp.data))]);
			if(length(genoLevels) <= 1)
			{
				stop("\tError: The levels of geno argument is less than two!\n");
			}
			
			# --- compute harmonic mean --- #
			no.reps <- data.frame(n = tapply(eval(parse(text = paste("temp.data$", respvar[i], sep = ""))),
							eval(parse(text = paste("temp.data$", geno, sep = ""))), FUN = length));
			no.reps <- as.numeric(1/sapply(1/no.reps, mean));
			result[[i]]$site[[j]]$numreps <- no.reps;
			
			# --- if design is Latinized Row-Column, check if the data follow case1 or case3 labeling --- #
			if(exptl.design == "LatinRowCol")
			{
				lengthPerCross <- tapply(temp.data[,respvar[i]], temp.data[ ,c(row, column)], length);
				if(all(lengthPerCross <= 1, na.rm = TRUE))
				{
					if(nlevels(temp.data[ , row]) > nlevels(temp.data[ , column]))
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
			
			# --- compute response rate --- #
			responseRate <- (obsused / obsread);
			result[[i]]$site[[j]]$responseRate <- responseRate;
			
			if(responseRate < 0.80)
			{
				result[[i]]$site[[j]]$manyNAWarning <- "Too many missing observations. Cannot proceed with analysis.";
				warning(paste("\tWarning: Too many missing observations. Cannot proceed with analysis on ", result[[i]]$site[[j]]$env, " of ",
								result[[i]]$respvar, ".\n", sep = ""))
				next;
			} else
			{
				# --- check if items in checkList are in temp.data, if not adjust checkList and create warnings --- #
				
				# --- CONSTRUCT THE MODEL --- #
				
				# --- CONSTRUCT THE MODEL: GENOTYPE FIXED FACTOR --- #
				if(exptl.design == "RCB" || exptl.design == "AugRCB")
				{
					myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + (1|",block, ")", sep="");
				} else if(exptl.design == "AugLS")
				{
					myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + (1|", row , "), + (1|", column, ")", sep ="");
				} else if(exptl.design == "Alpha")
				{
					myformula1 <- paste(respvar[i], " ~ 1 + ", geno," + (1|", rep,"/", block,")", sep = "");
				} else if(exptl.design == "RowCol")
				{
					myformula1 <- paste(respvar[i], " ~ 1 + ", geno," + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,")", sep = "");
				} else if(exptl.design == "LatinAlpha")
				{
					myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + (1|", rep,") + (1|", block,") + (1|", rep, ":", block, ")", sep = "");
				} else if(exptl.design == "LatinRowCol")
				{
					if(longerRow)
					{
						myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + (1|", rep,") + (1|", column, ") + (1|", rep,":", column,") + (1|", row, ")", sep = ""); 
					} else
					{
						myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + (1|", rep,") + (1|", row, ") + (1|", rep,":", row,") + (1|", column, ")", sep = "") ;
					}	
				}
				
				# --- call lmer function using myformula1 to get variance components table--- #
				#model <- lmer(formula(myformula1), data = temp.data)
				model <- try(lmer(formula(myformula1), data = temp.data), silent=TRUE);
				
				if (!is.null(model) && class(model)=="try-error") {  
					msg <- trimStrings(strsplit(model, ":")[[1]]);
					msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "));
					msg <- gsub("\"", "", msg);
					
					result[[i]]$site[[j]]$lmerRun <- "ERROR";
					result[[i]]$site[[j]]$lmerError <- msg;
					warning(paste("\tWarning: There are error occured on lmer model computed. Cannot proceed with analysis on ", result[[i]]$site[[j]]$env, " of ",
									result[[i]]$respvar, ".\n", sep = ""));
					next;
				} 
				
				result[[i]]$site[[j]]$lmerRun <- "NO ERROR";
				result[[i]]$site[[j]]$formula1 <- myformula1;
				result[[i]]$site[[j]]$model <- model;
				
				# --- get variance components --- #
				varcomp <- NULL;
				for( k in (1:length(VarCorr(model))))
				{
					varcomp <- rbind(varcomp, data.frame(Groups = names(VarCorr(model))[k], Variance = VarCorr(model)[[k]][1], Std.Dev. = attr(VarCorr(model)[[k]], "stddev")[[1]]));
				}
				varcomp <- rbind(varcomp, data.frame(Groups = "Residual", Variance = attr(VarCorr(model),"sc")^2, Std.Dev. = attr(VarCorr(model), "sc")));
				attr(varcomp, "heading") <- "Variance Components for Random Effects\n";
				result[[i]]$site[[j]]$varcomp.table <- varcomp;
				
				# --- for saving variance and num of reps --- #
				result[[i]]$site[[j]]$varcompnRep <- as.data.frame(attr(VarCorr(model), "sc")^2);
				result[[i]]$site[[j]]$varcompnRep$numRep <- result[[i]]$site[[j]]$numreps;
				result[[i]]$site[[j]]$varcompnRep$env <- result[[i]]$site[[j]]$env[[1]];
				colnames(result[[i]]$site[[j]]$varcompnRep) <- c(paste(respvar[i], "sigma2", sep="_"), paste(respvar[i],"No.Rep", sep="_"), env);
				if(j == 1) 
				{
					result[[i]]$out.sigma2 <- result[[i]]$site[[j]]$varcompnRep;
				} else
				{
					result[[i]]$out.sigma2 <- rbind(result[[i]]$out.sigma2, result[[i]]$site[[j]]$varcompnRep);
				}
				
				# --- TEST SIGNIFICANCE OF GENOTYPIC EFFECT USING MAXIMUM LIKELIHOOD RATIO TEST --- #
				myformula2 <- gsub(paste(" + ", geno, sep = ""), "", myformula1, fixed = TRUE);
				model1 <- lmer(formula(myformula1), data = temp.data, REML = T);
				model2 <- lmer(formula(myformula2), data = temp.data, REML = T);
				
				result[[i]]$site[[j]]$formula2 <- myformula2;
				
				# --- COMPUTE GENOTYPE MEANS --- #
				myformula3 <- gsub("~ 1", "~ 0", myformula1);
				model3 <- lmer(formula(myformula3), data = temp.data);
				sumStat.table <- data.frame(summary(model3)$coefficients)[, 1:2]; # modified for R 3.0.2
				rownames(sumStat.table) <- gsub(geno, "", rownames(sumStat.table));
				sumStat.table <- cbind(rownames(sumStat.table), sumStat.table);
				rownames(sumStat.table) <- NULL;
				colnames(sumStat.table) <- c(geno, "LSMean", "StdErrMean");
				result[[i]]$site[[j]]$summary.statistic <- sumStat.table;
				
				# --- display standard erro of the differents --- #
				noEntries <- nlevels(temp.data[,match(geno, names(temp.data))]);
				covs <- as.matrix(vcov(model3)[1:noEntries, 1:noEntries]);
				vars <- diag(covs);
				vdiff <- outer(vars, vars, "+") - 2 * covs;
				sed <- sqrt(vdiff[upper.tri(vdiff)]);
				
				# --- display SED table --- #
				minSed <- formatC(as.numeric(format(min(sed), scientific=FALSE)), format="f");
				meanSed <- formatC(as.numeric(format(mean(sed), scientific=FALSE)), format="f");
				maxSed <- format(as.numeric(format(max(sed), scientific=FALSE)), format="f");
				sedCol <- rbind(minSed, meanSed, maxSed);
				rowNames <- rbind("Minimu ", "Average ", "Maximum ");
				sedTable <- as.table(cbind(rowNames, sedCol));
				rownames(sedTable) <- c("", "", "");
				colnames(sedTable) <- c("", "Estimate");
				result[[i]]$site[[j]]$sedTable <- sedTable;
				
				# --- For saving to file --- #
				result[[i]]$site[[j]]$sum.out <- sumStat.table;
				result[[i]]$site[[j]]$sum.out$Env <- result[[i]]$site[[j]]$env[[1]];
				colnames(result[[i]]$site[[j]]$sum.out) <- c(geno, paste(result[[i]]$respvar, "Mean" , sep="_"), paste(result[[i]]$respvar, "StdErrMean", sep="_"), env);
				
				if(j==1) 
				{
					result[[i]]$meansse.out <- result[[i]]$site[[j]]$sum.out;
					result[[i]]$means.out <- result[[i]]$site[[j]]$sum.out[-3];
				} else
				{
					result[[i]]$meansse.out <- rbind(result[[i]]$meansse.out, result[[i]]$site[[j]]$sum.out);
					result[[i]]$means.out <- rbind(result[[i]]$means.out, result[[i]]$site[[j]]$sum.out[-3]);
				}
				
				result[[i]]$site[[j]]$residuals <- resid(model1);
				result[[i]]$site[[j]]$fitted.values <- fitted(model1);
				result[[i]]$site[[j]]$data <- temp.data;
			} ## --- end of else stmt --- if (responseRate < 0.80) --- #
			
		} ## --- end stmt --- for (j in (1:nlevels(data[,match(env, names(data))])))
		
		if(exptl.design == "RCB" || exptl.design == "AugRCB")
		{
			byVariables <- c(env, geno, block);
		} else if(exptl.design == "AugLS")
		{
			byVaribles <- c(env, geno, row, column);
		} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
		{
			byVaribles <- c(env, geno, block, rep);
		} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
		{
			byVariables <- c(env, geno, row, column, rep);
		}
		
		# --- to consolidate means and variances for fixed --- #
		if(i == 1)
		{
			means.out.all <- result[[i]]$means.out;
			meansse.out.all <- result[[i]]$meansse.out;
			varrep.out.all <- result[[i]]$out.sigma2;
			meansWarning <- "not empty";
			meansseWarning <- "not empty";
			varrepWarning <- "not emapty";
		} else
		{
			meansOut2 <- result[[i]]$means.out;
			if(!is.null(meansOut2))
			{
				if(is.null(means.out.all))
				{
					means.out.all <- meansOut2;
					meansWarning <- "empty";
				} else 
				{
					means.out.all <- merge(means.out.all, meansOut2, by = c(env, geno), all= TRUE);
					meansWarning <- "not empty";
				}
			}
			
			meansseOut2 <- result[[i]]$meansse.out;
			if(!is.null(meansseOut2))
			{
				if(is.null(meansse.out.all))
				{
					meansse.out.all <- meansseOut2;
					meansseWarning <- "empty";
				} else
				{
					meansse.out.all <- merge(meansse.out.all, meansseOut2, by=c(env, geno), all = TRUE);
					meansseWarning <- "not empty";
				}
			}
			
			sigmaOut2 <- result[[i]]$out.sigma2;
			if(!is.null(sigmaOut2))
			{
				if(is.null(varrep.out.all))
				{
					varrep.out.all <- sigmaOut2;
					varrepWarning <- "empty";
				} else
				{
					varrep.out.all <- merge(varrep.out.all, sigmaOut2, by=c(env), all= TRUE);
					varrepWarning <- "not empty";
				}
			}
		}	
	} ## -- end stmt -- for (i in (1:length(respvar)))
	
	detach("package:lme4");
	return(list(output = result, 
					means = means.out.all, meansWarning = meansWarning,
					meansse = meansse.out.all, meansseWarning = meansseWarning,
					varrep = varrep.out.all, varrepWarning = varrepWarning,
					byVars = byVariables,
					analysisType = "SingleSiteAnalysis",
					populationType = "SSSL"
			));
}
