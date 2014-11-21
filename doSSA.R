###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Nov 14, 2014
# FileName: doSSA.R
###############################################################################


doSSA <- function(
		data
)
{
	UseMethod("doSSA");
}
doSSA.PhenotypicData <- function(
		data
)
{
	if(!inherits(data, "PhenotypicData"))
		stop("\tError: The argument of data must be of class PhenotypicData!\n");
	if(is.null(data$doRestricted))
	{
		warning(paste("\tWarning: It will use the default parameter to restrict phenotypic data;\n",
						"\tMissing rate should be not larger than 0.2 of each trait on each environment!\n",
						"\tAnd used the grand mean of each trait on each environment to instead of missing all obseravtion!\n", sep=""));
		data <- restrict.pheno.data(data);
	}
	
	options(show.signif.stars=FALSE);
	library(lme4);
	
	for(i in (1:length(data$traits)))
	{
		#--- saving all analysis output to this list.---#
		data$traits[[i]]$analysis <- list();
		#--- saving single site analysis output---#
		data$traits[[i]]$analysis$ssa <- list();
		data$traits[[i]]$analysis$ssa$sites <- list();
		respvar <- data$traits[[i]]$name;
		
		for(j in (1:length(data$traits[[i]]$sites)))
		{
			data$traits[[i]]$analysis$ssa$sites[[j]] <- list();
			env.name <- data$traits[[i]]$sites[[j]]$name;
			data$traits[[i]]$analysis$ssa$sites[[j]]$name <- env.name;
			temp.data <- data$traits[[i]]$sites[[j]]$data;
			#--- retrive all design factor name and experimental design---#
			exptl.design <- data$traits[[i]]$sites[[j]]$design$exptl.design;
			geno <- data$traits[[i]]$sites[[j]]$design$geno;
			env <- data$traits[[i]]$sites[[j]]$design$env;
			if(exptl.design == "RCB" || exptl.design == "AugRCB")
			{
				block <- data$traits[[i]]$sites[[j]]$design$block;
			} else if(exptl.design == "AugLS")
			{
				row <- data$traits[[i]]$sites[[j]]$design$row;
				column <- data$traits[[i]]$sites[[j]]$design$column;
			} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
			{
				block <- data$traits[[i]]$sites[[j]]$design$block;
				rep <- data$traits[[i]]$sites[[j]]$design$rep;
			} else if(exptl.design == "RowCol")
			{
				row <- data$traits[[i]]$sites[[j]]$design$row;
				column <- data$traits[[i]]$sites[[j]]$design$column;
				rep <- data$traits[[i]]$sites[[j]]$design$rep;
			} else if(exptl.design == "LatinRowCol")
			{
				row <- data$traits[[i]]$sites[[j]]$design$row;
				column <- data$traits[[i]]$sites[[j]]$design$column;
				rep <- data$traits[[i]]$sites[[j]]$design$rep;
				longerRow <- data$traits[[i]]$sites[[j]]$design$longerRow
			}
			
			
			
			#--- checking whether this site is restricted
			if(data$traits[[i]]$sites[[j]]$restricted$isTRUE)
			{
				warning(paste("\tWarning: Too many missing observations. Cannot proceed with analysis on ", env.name, " of ",
								respvar, ".\n", sep = ""))
				next;
			} else
			{
				#--- compute harmonic mean---#
				no.reps <- data.frame(n = tapply(eval(parse(text = paste("temp.data$", respvar, sep = ""))), eval(parse(text = paste("temp.data$", geno, sep = ""))), FUN = length));
				no.reps <- as.numeric(1/sapply(1/no.reps, mean));
				data$traits[[i]]$analysis$ssa$sites[[j]]$numreps <- no.reps;
				
				# --- CONSTRUCT THE MODEL: GENOTYPE FIXED FACTOR --- #
				if(exptl.design == "RCB" || exptl.design == "AugRCB")
				{
					myformula1 <- paste(respvar, " ~ 1 + ", geno, " + (1|",block, ")", sep="");
				} else if(exptl.design == "AugLS")
				{
					myformula1 <- paste(respvar, " ~ 1 + ", geno, " + (1|", row , "), + (1|", column, ")", sep ="");
				} else if(exptl.design == "Alpha")
				{
					myformula1 <- paste(respvar, " ~ 1 + ", geno," + (1|", rep,"/", block,")", sep = "");
				} else if(exptl.design == "RowCol")
				{
					myformula1 <- paste(respvar, " ~ 1 + ", geno," + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,")", sep = "");
				} else if(exptl.design == "LatinAlpha")
				{
					myformula1 <- paste(respvar, " ~ 1 + ", geno, " + (1|", rep,") + (1|", block,") + (1|", rep, ":", block, ")", sep = "");
				} else if(exptl.design == "LatinRowCol")
				{
					if(longerRow)
					{
						myformula1 <- paste(respvar, " ~ 1 + ", geno, " + (1|", rep,") + (1|", column, ") + (1|", rep,":", column,") + (1|", row, ")", sep = ""); 
					} else
					{
						myformula1 <- paste(respvar, " ~ 1 + ", geno, " + (1|", rep,") + (1|", row, ") + (1|", rep,":", row,") + (1|", column, ")", sep = "") ;
					}	
				}
				#--- call lmer function using myformula1 to get variance components table ---#
				#--- model <- lmer(formula(myformula1), data = temp.data) ---#
				model <- try(lmer(formula(myformula1), data = temp.data), silent=TRUE);
				
				if (!is.null(model) && class(model)=="try-error") {  
					msg <- trimStrings(strsplit(model, ":")[[1]]);
					msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "));
					msg <- gsub("\"", "", msg);
					
					data$traits[[i]]$analysis$ssa$sites[[j]]$lmerRun <- "ERROR";
					data$traits[[i]]$analysis$ssa$sites[[j]]$lmerError <- msg;
					warning(paste("\tWarning: There are error occured on lmer model computed. Cannot proceed with analysis on ", env.name, " of ",
									respvar, ".\n", sep = ""));
					next;
				} 
				
				data$traits[[i]]$analysis$ssa$sites[[j]]$lmerRun <- "NO ERROR";
				data$traits[[i]]$analysis$ssa$sites[[j]]$formula1 <- myformula1;
				data$traits[[i]]$analysis$ssa$sites[[j]]$model <- model;
				
				#--- get variance components ---#
				varcomp <- NULL;
				for(k in (1:length(VarCorr(model))))
				{
					varcomp <- rbind(varcomp, data.frame(
									Groups = names(VarCorr(model))[k], 
									Variance = VarCorr(model)[[k]][1],
									Std.Dev. = attr(VarCorr(model)[[k]],"stddev")[[1]]
					));
				}
				varcomp <- rbind(varcomp, data.frame(
								Groups = "Residual",
								Variance = attr(VarCorr(model),"sc")^2,
								Std.Dev. = attr(VarCorr(model), "sc")
								));
				attr(varcomp, "heading") <- "Variance Components for Random Effects\n";
				data$traits[[i]]$analysis$ssa$sites[[j]]$varcomp.table <- varcomp;
				
				#--- for saving variance and sum of reps---#
				data$traits[[i]]$analysis$ssa$sites[[j]]$varcompnRep <- as.data.frame(attr(VarCorr(model), "sc")^2);
				data$traits[[i]]$analysis$ssa$sites[[j]]$varcompnRep$No.Rep <- no.reps;
				data$traits[[i]]$analysis$ssa$sites[[j]]$varcompnRep$Env <- env.name;
				colnames(data$traits[[i]]$analysis$ssa$sites[[j]]$varcompnRep) <- c(paste(respvar, "sigma2", sep="_"),
						paste(respvar, "No.Rep", sep="_"),
						env);
				#--- test significance of genotype effect using maximum likelihood ratio test ---#
				myformula2 <- gsub(paste("+", geno, sep=""), "", myformula1, fixed = TRUE);
				model1 <- lmer(formula(myformula1), data = temp.data, REML = T);
				model2 <- lmer(formula(myformula2), data = temp.data, REML = T);
				
				data$traits[[i]]$analysis$ssa$sites[[j]]$formula2 <- myformula2;
				
				#--- conmpute genotype means---#
				myformula3 <- gsub("~ 1", "~ 0", myformula1);
				model3 <- lmer(formula(myformula3), data = temp.data);
				sumStat.table <- data.frame(summary(model3)$coefficients)[ , 1:2];
				rownames(sumStat.table) <- gsub(geno, "", rownames(sumStat.table));
				sumStat.table <- cbind(rownames(sumStat.table), sumStat.table);
				rownames(sumStat.table) <- NULL;
				colnames(sumStat.table) <- c(geno, "LSMean", "StdErrMean");
				data$traits[[i]]$analysis$ssa$sites[[j]]$summary.statistic <- sumStat.table;
				
				#--- display standard error of the differents ---#
				noEntries <- nlevels(temp.data[,geno]);
				covs <- as.matrix(vcov(model3)[1:noEntries, 1:noEntries]);
				vars <- diag(covs);
				vdiff <- outer(vars, vars, "+") - 2 * covs;
				sed <- sqrt(vdiff[upper.tri(vdiff)]);
				
				#--- display sed table ---#
				minSed <- formatC(as.numeric(format(min(sed), scientific = FALSE)), format = "f");
				meanSed <- formatC(as.numeric(format(mean(sed), scientific = FALSE)), format = "f");
				maxSed <- formatC(as.numeric(format(max(sed), scientific = FALSE)), format = "f");
				sedCol <- rbind(minSed, meanSed, maxSed);
				rowNames <- rbind("Minimu", "Average", "Maximum");
				sedTable <- as.table(cbind(rowNames, sedCol));
				rownames(sedTable) <- c("", "", "");
				colnames(sedTable) <- c("", "Estimate");
				data$traits[[i]]$analysis$ssa$sites[[j]]$sedTable <- sedTable;
				
				#--- for saving to file ---#
				data$traits[[i]]$analysis$ssa$sites[[j]]$sum.out <- sumStat.table;
				data$traits[[i]]$analysis$ssa$sites[[j]]$sum.out$Env <- env.name;
				colnames(data$traits[[i]]$analysis$ssa$sites[[j]]$sum.out) <- c(geno, paste(respvar, "Mean", sep = "_"), paste(respvar, "StdErrMean", sep = "_"), env);
				
				data$traits[[i]]$analysis$ssa$sites[[j]]$residuals <- resid(model1);
				data$traits[[i]]$analysis$ssa$sites[[j]]$fitted.values <- fitted(model1);
		
			} #--- end statement of if(data$traits[[i]]$sites[[j]]$restricted$isTRUE) ---#
		}#--- end statement of for(j in (1:length(data$triatis[[i]]$sites))) ---#
		
	}#--- end statement of for(i in (1:length(data$traits)))---#
	
	class(data) <- c("SingleSiteAnalysis", class(data));
	return(data);
}
