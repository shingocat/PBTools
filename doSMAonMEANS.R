###############################################################################
# TODO: Do single marker analysis on Means
#
# [Arguments]
#	phenotypicData - means of each genotypic lines;
#	genoetypicData - GenotypicData after do restrict.geno.data;
#	geno - the geno variable to indicate the genotypic lines;
#	env - the env variable to indicate the environment if it have, otherwise keep it be null;
#	resp.var - the respone variable in phenotypic data;
#	is.EnvFixed - TRUE for specified env factor to be fixed if env variable is not null;
#	include.ht - TRUE for specified including heterozygous on single marker analysis;
#	test - specified the testing method for each marker;
#	digits - specified the digits for the compute result;
# 
# Author: mqin
# Date: Dec 10, 2014
# FileName: doSMAonMEANS.R
###############################################################################


#--- The is method is to implement using means of each lines to process single marker analysis---#
doSMA.MEANS <- function(
		phenotypicData,
		genotypicData,
		geno,
		env,
		resp.var,
		is.EnvFixed = TRUE,
		include.ht = FALSE,
		test = c("F", "Chisq"),
		digits = 4
)
{
	#--- not checking any conditions by far---#
	genotypicData <- genotypicData$restricted$data;
	marker.number <- ncol(genotypicData) - 1;
	marker.names <- colnames(genotypicData)[-1];
	#---suppress warn message---#
	old.options <- options();
	options(warn = -1);
	#library("lme4");
	#library("phia");
	library("lsmeans")
	
	#--- reformat all genotypic data using dp.code instead of ht.code basing on include.ht ----#
	#--- Because it will recode to default coding sytem, it could be define such code directly!---#
	dp.code <- 2;
	rp.code <- 0;
	ht.code <- 1;
	na.code <- NA;
	if(include.ht)
	{
		genotypicData[ , -1][genotypicData[,-1] == ht.code] <- dp.code;
	}
	genotypicData <- apply(genotypicData, 2, factor);
	
	include.env <- TRUE;
	if(missing(env))
		include.env <- FALSE;
	data <- merge(phenotypicData, genotypicData, by = geno);
	data[, geno] <- as.factor(data[ , geno]);
	result <- list();
	result$trait.number <- length(resp.var);
	result$trait.names <- resp.var;
	result$include.env <- include.env;
	result$include.ht <- include.ht;
	result$is.EnvFixed <- is.EnvFixed;
	result$marker.number <- marker.number;
	result$marker.names <- marker.names;
	result$traits <- list();
	
	for(i in 1:length(resp.var))
	{
		trait.name <- resp.var[i];
		result$traits[[i]] <- list();
		result$traits[[i]]$name <- trait.name;
		result$traits[[i]]$markers <- list();
		means.out <- c();
		var.table <- c();
		p.table <- c();
		for(j in 1:marker.number)
		{		
			marker.name <- marker.names[j];
			myformula1 <- c();
			if(include.env)
			{
				if(is.EnvFixed)
				{
					env.stmt <- paste(env ," + ", marker.name, " : ", env, sep ="");
					myformula1 <- paste(trait.name, " ~ 1 + ", marker.name, " + ", env.stmt, sep = "");
					model <- try(lm(formula(myformula1), data = data), silent = TRUE);					
				} else
				{
					env.stmt <- paste("(1|", env ,") + (1|", marker.name, " : ", env, ") ",sep ="");
					myformula1 <- paste(trait.name, " ~ 1 + ", marker.name, " + ", env.stmt, sep = "");
					model <- try(lmer(formula(myformula1), data = data), silent = TRUE);
				}
			} else
			{
				myformula1 <- paste(trait.name, " ~ 1 + ", marker.name, sep = "");
				model <- try(lm(formula(myformula1), data = data), silent = TRUE);
			}
			if(!is.null(model) && class(model) == "try-error")
			{
				warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
				next;
			}
			result$traits[[i]]$markers[[j]] <- list();
			result$traits[[i]]$markers[[j]]$name <- marker.name;
			result$traits[[i]]$markers[[j]]$myformula1 <- myformula1;
			result$traits[[i]]$markers[[j]]$model <- model;
			if(include.env)
			{
				if(is.EnvFixed)
				{
					adjust.means <- lsmeans(model, specs = marker.name, by = env, weights = "cells" );
					adjust.means <- summary(adjust.means)[ , c(1:3)];
					result$traits[[i]]$markers[[j]]$adjust.means <- adjust.means;
					env.levels <- levels(adjust.means[ , env]);
					temp.means <- c();
					for(k in 1:length(env.levels))
					{	
						env.name <- env.levels[k];
						env.means <- adjust.means[adjust.means[ , env] == env.name, ];
						rp.means <- if(rp.code %in% env.means[ , marker.name]) signif(env.means[env.means[ , marker.name] == rp.code,3], digits) else NA;
						ht.means <- if(ht.code %in% env.means[ , marker.name]) signif(env.means[env.means[ , marker.name] == ht.code,3], digits) else NA;
						dp.means <- if(dp.code %in% env.means[ , marker.name]) signif(env.means[env.means[ , marker.name] == dp.code,3], digits) else NA;					
						temp <- c(marker.name, env.name, rp.means, ht.means, dp.means);
						temp.means <- rbind(temp.means, temp);
					}
					
					#--- a formatted variance table ---#	
					model.aov <- anova(model);
					model.aov.rownames <- rownames(model.aov[1]);
					marker.sq <- if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Sum Sq", digits) else NA;
					marker.p <-  if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Pr(>F)", digits) else NA;
					env.sq <- if(env %in% model.aov.rownames) signif(model.aov[model.aov.rownames == env, ]$"Sum Sq", digits) else NA;
					env.p <-  if(env %in% model.aov.rownames) signif(model.aov[model.aov.rownames == env, ]$"Pr(>F)", digits) else NA;
					marker.env.sq <- if(paste(marker.name, ":", env, sep="") %in% model.aov.rownames) signif(model.aov[model.aov.rownames == paste(marker.name, ":", env, sep=""), ]$"Sum Sq", digits) else NA;
					marker.env.p <- if(paste(marker.name, ":", env, sep="") %in% model.aov.rownames) signif(model.aov[model.aov.rownames == paste(marker.name, ":", env, sep=""), ]$"Pr(>F)", digits) else NA;
					residuals.sq <- if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Sum Sq", digits) else NA;
					residuals.p <-  if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Pr(>F)", digits) else NA; 
					
					temp.var.table <- rbind(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""), paste(env.sq, " (", env.p, ")", sep = ""), paste(marker.env.sq, " (", marker.env.p, ")", sep = ""), residuals.sq);
					temp.p.table <- c(marker.name, marker.p, marker.env.p);
					
				} else
				{
					#--- when env is random, contrust a means variance table ---#
					#--- marker effect---#
					intercept <- fixef(model)[[1]];
					marker.effect <- as.data.frame(fixef(model)[-1]);
					marker.effect <- data.frame(gsub(marker.name, "", rownames(marker.effect)), marker.effect);
					colnames(marker.effect) <- c(marker.name, "M_Effect");
					rownames(marker.effect) <- NULL;
					#--- env effect ---#
					env.effect <- eval(parse(text = paste("ranef(model)$'", env, "'", sep = "")));
					env.effect <- data.frame(gsub(env, "", rownames(env.effect)), env.effect);
					colnames(env.effect) <- c(env, "E_Effect");
					rownames(env.effect) <- NULL;
					#--- marker by env effect ---#
					mxe.effect <- eval(parse(text = paste("ranef(model)$'", marker.name, ":", env,"'", sep = "")));
					names <- t(as.data.frame(strsplit(rownames(mxe.effect), ":")));
					mxe.effect <- data.frame(names[,1], names[,2], mxe.effect + intercept);
					colnames(mxe.effect) <- c(marker.name, env, "MxE.Effect");
					rownames(mxe.effect) <- NULL;					
					#--- marker by env means---#
					mxe.means <- merge(mxe.effect, env.effect, by = env, all = TRUE);
					mxe.means <- merge(mxe.means, marker.effect, by = marker.name, all = TRUE);
					mxe.means <- data.frame(mxe.means[,match(marker.name, names(mxe.means))], mxe.means[,match(env, names(mxe.means))], rowSums(subset(mxe.means, select = c(MxE.Effect, E_Effect, M_Effect)), na.rm = TRUE));
					colnames(mxe.means) <- c(marker.name, env, "LSMean");
					adjust.means <- mxe.means;
					result$traits[[i]]$markers[[j]]$adjust.means <- adjust.means;
					env.levels <- levels(adjust.means[ , env]);
					temp.means <- c();
					for(k in 1:length(env.levels))
					{	
						env.name <- env.levels[k];
						env.means <- adjust.means[adjust.means[ , env] == env.name, ];
						rp.means <- if(rp.code %in% env.means[ , marker.name]) signif(env.means[env.means[ , marker.name] == rp.code,3], digits) else NA;
						ht.means <- if(ht.code %in% env.means[ , marker.name]) signif(env.means[env.means[ , marker.name] == ht.code,3], digits) else NA;
						dp.means <- if(dp.code %in% env.means[ , marker.name]) signif(env.means[env.means[ , marker.name] == dp.code,3], digits) else NA;					
						temp <- c(marker.name, env.name, rp.means, ht.means, dp.means);
						temp.means <- rbind(temp.means, temp);
					}
					#--- a formatted variance table ---#	
					model2 <- update(model, formula(paste(". ~ . - (1|", marker.name, ":", env, ")", sep = "")));
					model.aov <- anova(model);
					marker.sq <- signif(model.aov$"Sum Sq", digits);
					marker.p <- signif(model.aov$"Pr(>F)", digits);
					model.com.aov <- anova(model2, model);
					marker.env.p <- signif(model.com.aov$"Pr(>Chisq)"[2], digits);
					temp.var.table <- c(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""));
					temp.p.table <- c(marker.name, marker.p, marker.env.p);
				} #--- end of statement if(is.EnvFixed)---#
			} else
			{			
				adjust.means <- lsmeans(model, specs = marker.name, weights = "cells");
				adjust.means <- summary(adjust.means)[,c(1:3)];
				result$traits[[i]]$markers[[j]]$adjust.means <- adjust.means;
				rp.means <- if(rp.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == rp.code, 2], digits) else NA;
				ht.means <- if(ht.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == ht.code, 2], digits) else NA;
				dp.means <- if(dp.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == dp.code, 2], digits) else NA;
				temp.means <- c(marker.name, rp.means, ht.means, dp.means);
#				marker.test <- testInteractions(model, pairwise = marker.name, test = test);
#				result$traits[[i]]$markers[[j]]$marker.test <- marker.test;
#				#--- a formatted means table---#
#				rp.marker.means <- if(length(adjust.means[adjust.means[,1] == rp.code, 2])) signif(adjust.means[adjust.means[,1] == rp.code, 2], digits) else NA;
#				ht.marker.means <- if(length(adjust.means[adjust.means[,1] == ht.code, 2])) signif(adjust.means[adjust.means[,1] == ht.code, 2], digits) else NA;
#				dp.marker.means <- if(length(adjust.means[adjust.means[,1] == dp.code, 2])) signif(adjust.means[adjust.means[,1] == dp.code, 2], digits) else NA;
#				temp.means <- c(marker.name, rp.marker.means, ht.marker.means, dp.marker.means);
#				marker.test.rownames <- rownames(marker.test["Value"]);
#				test0.1.value <- if("0-1" %in% marker.test.rownames) signif(marker.test[marker.test.rownames == "0-1", ]$"Value", digits) else NA;
#				test0.1.p <- if("0-1" %in% marker.test.rownames) signif(marker.test[marker.test.rownames == "0-1", ]$"Pr(>F)", digits) else NA;
#				test0.2.value <- if("0-2" %in% marker.test.rownames) signif(marker.test[marker.test.rownames == "0-2", ]$"Value", digits) else NA;
#				test0.2.p <- if("0-2" %in% marker.test.rownames) signif(marker.test[marker.test.rownames == "0-2", ]$"Pr(>F)", digits) else NA;
#				test1.2.value <- if("1-2" %in% marker.test.rownames) signif(marker.test[marker.test.rownames == "1-2", ]$"Value", digits) else NA;
#				test1.2.p <- if("1-2" %in% marker.test.rownames) signif(marker.test[marker.test.rownames == "1-2", ]$"Pr(>F)", digits) else NA;
#				temp.means <- c(temp.means, test0.1.value, test0.1.p, test0.2.value, test0.2.p, test1.2.value, test1.2.p);
				
				#--- a formatted variance table ---#	
				model.aov <- anova(model);
				model.aov.rownames <- rownames(model.aov[1]);
				marker.sq <- if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Sum Sq", digits) else NA;
				marker.p <-  if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Pr(>F)", digits) else NA;
				residuals.sq <- if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Sum Sq", digits) else NA;
				residuals.p <-  if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Pr(>F)", digits) else NA; 
				
				temp.var.table <- rbind(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""), residuals.sq);
				temp.p.table <- c(marker.name, marker.p);
			}#--- end of statement if(include.env) ---#
	
			means.out <- rbind(means.out, temp.means);
			var.table <- cbind(var.table, temp.var.table);
			p.table <- rbind(p.table, temp.p.table);
		}
		rownames(means.out) <- NULL;
		colnames(var.table) <- NULL;
		rownames(p.table) <- NULL;
		if(include.env)
		{
			colnames(means.out) <- c("Marker","Env", "mm", "Mm", "MM");
			if(is.EnvFixed) 
			{
				rownames(var.table) <- c("Marker", "M SS(P.value)", "E SS (P.value)", "MxE SS (P.value)", "Residuals SS");
			} else
			{
				rownames(var.table) <- c("Marker", "M SS(P.value)");
			}
			colnames(p.table) <- c("Marker", "P.value(M)", "P.value(MxE)")
		} else
		{
			colnames(means.out) <- c("Marker", "mm", "Mm", "MM");
			rownames(var.table) <- c("Marker", "M SS (p.value)", "Residuals SS");
			colnames(p.table) <- c("Marker", "P.value(M)");
		}
	
		result$traits[[i]]$means.table <- means.out;
		result$traits[[i]]$var.table <- var.table;
		result$traits[[i]]$p.table <- p.table;
	} #--- end of statement of for(i in 1:length(resp.var)) ---#
	
	class(result) <- "SMAonMeans";
	#detach("package:lme4", unload = TRUE);
	#detach("package:phia", unload = TRUE);
	detach("package:lsmeans", unload = TRUE);
	options(old.options);
	return(result);
}

print.SMAonMeans <- function(
		data, 
		levels = c("SIG", "ALL"), 
		p.value = 0.05
)
{	
	levels <- match.arg(levels);
	if(identical(levels, "SIG"))
	{
		if(missing(p.value))
			p.value <- 0.05;
	}
	cat(rep("=", 50), "\n", sep ="");
	cat("Single Marker Analysis\n");
	cat(rep("=", 50), "\n", sep ="");
	cat("The number of trait(s): ", data$trait.number, ".\n", sep = "");
	cat("Trait Name:", data$trait.names, ".\n", sep = " ");
	cat("The number of marker(s): ", data$marker.number, ".\n", sep = "");
	cat("Marker:", data$marker.names, ".\n", sep = " ")
	cat(rep("=", 50), "\n", sep ="");
	for(i in 1:data$trait.number)
	{
		cat("Trait Name: ", data$traits[[i]]$name, ".\n", sep = "");
		p.table <- as.data.frame(data$traits[[i]]$p.table);
		means.table <- as.data.frame(data$traits[[i]]$means.table);
		if(identical(levels, "SIG"))
		{
			if(data$include.env)
			{
				p.table.m <- p.table[as.numeric(as.character(p.table[,2])) <= p.value, ];
				p.table.me <- p.table[as.numeric(as.character(p.table[,3])) <= p.value, ];
				p.table<- merge(p.table.m, p.table.me, by = "Marker", all = TRUE)[, 1:3];
				colnames(p.table) <- c("Marker", "P.value(M)", "P.value(MxE)");
				marker.names <- as.character(p.table[,1]);
				means.table <- means.table[which(as.character(means.table[,1]) %in% marker.names), ]; 
				cat("Statistics test on each marker:\n");
				print(as.data.frame(p.table), quote = FALSE,  row.names = FALSE);
				cat("\n");
				cat("Means of each marker across genotype:\n");
				print(as.data.frame(means.table), quote = FALSE,  row.names = FALSE);
			} else
			{
				p.table <- p.table[as.numeric(as.character(p.table[,2])) <= p.value, ];
				marker.names <- as.character(p.table[,1]);
				means.table <- means.table[which(as.character(means.table[,1]) %in% marker.names), ]; 
				cat("Statistics test on each marker:\n");
				print(as.data.frame(p.table), quote = FALSE,  row.names = FALSE);
				cat("\n");
				cat("Means of each marker across genotype:\n");
				print(as.data.frame(means.table), quote = FALSE,  row.names = FALSE);
			}
		} else
		{
			cat("Statistics test on each marker:\n");
			print(as.data.frame(p.table), quote = FALSE,  row.names = FALSE);
			cat("\n");
			cat("Means of each marker across genotype:\n");
			print(as.data.frame(means.table), quote = FALSE,  row.names = FALSE);
		}
		cat(rep("-", 50), "\n", sep = "");
		cat("\n")
	}
}
