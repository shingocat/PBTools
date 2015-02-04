###############################################################################
# TODO: compute single marker analysis
# 
# Argument:
#	phenotypicData - 
#	genotypicData -
#	geno - the indicator of lines/materials shared between phenotypic and genotypic Data;
# 	resp.var - the responsed variable on phenotypic data;
#	include.env -
#	is.EnvFixed -
#	include.ht - if it is true then the heterozygous will change to donor type else it will remove from data
#	test - 
# Author: mqin
# Date: Nov 28, 2014
# FileName: doSM.R
###############################################################################



doSMA <- function(
		phenotypicData,
		genotypicData,
		#	geno,
		#	env,
		resp.var,
		include.env = TRUE,
		is.EnvFixed = TRUE,
		include.ht = FALSE,
		#	test = c("F", "Chisq"),
		digits = 4,
		...
)
{
	if(missing(phenotypicData))
		stop("\tError: The phenotypicData argument could not be null!\n");
	if(missing(genotypicData))
		stop("\tError: The genotypicData argument could not be null!\n");
#	if(missing(geno))
#		stop("\tError: The genotypicData argument could not be null!\n");
#	if(!missing(env))
#	{
#		if(missing(is.EnvFixed))
#			stop("\tError: The is.EnvFixed argument could not be null when it specified env argument!\n");
#	}
	if(missing(resp.var))
		warning("\tWarnings: It will compute all respone variables of phenotypicData");
#	if(missing(test))
#		test = match.arg(test);
	if(missing(digits))
		digits = 4;
	
	UseMethod("doSMA");
}

#--- It is used to process single marker analysis on PhenotypicData---#
doSMA.PhenotypicData <- function(
		phenotypicData,
		genotypicData,
#		geno,
#		env,
		resp.var,
		include.env = TRUE,
		is.EnvFixed = TRUE,
		include.ht = FALSE,
#		test = c("F", "Chisq"),
		digits = 4,
		...
)
{
	#--- checking resp.var whether exists in the phenotypic data, ortherwise using all traits in the phenotypicData ---#
	if(!missing(resp.var))
	{
		trait.names <- phenotypicData$trait.names;
		temp.resp.var <- c();
		for(i in 1:length(resp.var))
		{
			if(resp.var[i] %in% trait.names)
			{
				temp.resp.var <- c(temp.resp.var, resp.var[i]);
			} else
			{
				warnings(paste("\tWarnings: It will be omitted the ", resp.var[i], " because of it is not exist!\n", sep = ""));
			}
		}
		if(length(temp.resp.var) == 0)
			stop("\tError: All the specified resp.var are not exists in the phenotypicdata!\n");
		resp.var <- temp.resp.var;
	} else
	{
		resp.var <- phenotypicData$trait.names;
	}
	
	if(phenotypicData$isMean)
	{
		if(!phenotypicData$isRestricted)
		{
			warnings("\tWarning: It will be used default argument to restrict phenodtypic data!\n");
			phenotypicData <- restrict.pheno.data(phenotypicData);
		}
		if(!genotypicData$isRestricted)
		{
			warnings("\tWarning: It will be used default argument to restrict genotypic data!\n");
			genotypicData <- restrict.geno.data(genotypicData);
		}
		
		#---getting all the requried informations---#
		genotypicData.restricted <- genotypicData$restricted$data;
		marker.number <- ncol(genotypicData.restricted) - 1;
		marker.names <- colnames(genotypicData.restricted)[-1];
		#---suppress warn message---#
		old.options <- options();
		options(warn = -1);
		#library("lsmeans");
		
		#--- reformat all genotypic data using dp.code instead of ht.code basing on include.ht ----#
		#--- Because it will recode to default coding sytem, it could be define such code directly!---#
		dp.code <- 2;
		rp.code <- 0;
		ht.code <- 1;
		na.code <- NA;
		if(include.ht)
		{
			genotypicData.restricted[,-1][genotypicData.restricted[,-1] == ht.code] <- dp.code;
		}
		genotypicData.restricted <- apply(genotypicData.restricted, 2, factor);
		#--- checking whether is single environment---#
		if(include.env && phenotypicData$isSingleEnv)
		{	
			include.env <- FALSE;
			warning("\tWarnings: The phenotypic data is single environmental data and The is.EnvFixed paramter will be omitted!\n");
		}
		
		for(i in 1:length(resp.var))
		{
			#--- checking whether the resp.var name is the trait name of phenotypicData---#
			for(j in 1:length(phenotypicData$traits))
			{
				if(identical(resp.var[i], phenotypicData$traits[[j]]$name))
				{
					trait.name <- resp.var[i];
					#--- checking the analysis structure of phenotypicData--#
					if(is.null(phenotypicData$traits[[j]]$analysis))
						phenotypicData$traits[[j]]$analysis <- list();
					if(is.null(phenotypicData$traits[[j]]$analysis$sma))
					{
						phenotypicData$traits[[j]]$analysis$sma <- list();
					} else
					{
						phenotypicData$traits[[j]]$analysis$sma <- NULL;
						phenotypicData$traits[[j]]$analysis$sma <- list();
					}
					phenotypicData$traits[[j]]$analysis$sma$include.env <- include.env;
					phenotypicData$traits[[j]]$analysis$sma$is.EnvFixed <- is.EnvFixed;
					phenotypicData$traits[[j]]$analysis$sma$include.ht <- include.ht;
					phenotypicData$traits[[j]]$analysis$sma$marker.number <- marker.number;
					
					means.out <- c();
					var.table <- c();
					p.table <- c();
					
					#---according to whether including env to determinate how to analysis single marker ---#
					if(include.env)
					{
						#--- combined all environmental data into one data frame---#
						temp.phenotypicData <- c();
						env <- c();
						for(k in 1:length(phenotypicData$traits[[j]]$envs))
						{
							temp.phenotypicData <- rbind(temp.phenotypicData, phenotypicData$traits[[j]]$envs[[k]]$data);
							if(length(env) == 0)
								env <- phenotypicData$traits[[j]]$envs[[k]]$design$env;
						}
						#--- merge phenotypicData and genotoypicData into one data frame ---#
						geno.name.phenotypicData <- phenotypicData$traits[[j]]$envs[[1]]$design$geno;
						geno.name.genotypicData <- genotypicData$geno.name;
						data <- merge(temp.phenotypicData, genotypicData.restricted, by.x = geno.name.phenotypicData, by.y = geno.name.genotypicData);
						data[,trait.name] <- as.numeric(as.character(data[,trait.name]));
						
						library("lme4");
						
						for(k in 1:marker.number)
						{
							marker.name <- marker.names[k];
							myformula <- c();
							
							if(is.EnvFixed)
							{
								myformula <- paste(trait.name, " ~ 1 + ", marker.name, " + ", env, " + ", marker.name, " : ", env, sep = "");
								model <- try(lm(formula(myformula), data = data), silent = TRUE);
								if(!is.null(model) && class(model) == "try-error")
								{
									warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
									next;
								}
								adjust.means <- lsmeans(model, specs = marker.name, by = env, weights = "cells");
								adjust.means <- summary(adjust.means)[ , c(1:3)];
								env.levels <- levels(adjust.means[,env]);
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
								myformula <- paste(trait.name, " ~ 1 + ", marker.name, " + ", "(1|", env, ") + (1|", marker.name, " : ", env, ") ", sep = "");
								model <- try(lmer(formula(myformula), data = data), silent = TRUE);
								if(!is.null(model) && class(model) == "try-error")
								{
									warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
									next;
								}
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
								#result$traits[[i]]$markers[[j]]$adjust.means <- adjust.means;
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
								marker.p <- signif(model.aov$"F value", digits);
								model.com.aov <- anova(model2, model);
								marker.env.p <- signif(model.com.aov[[8]][2], digits);
								temp.var.table <- c(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""));
								temp.p.table <- c(marker.name, marker.p, marker.env.p);
							} #--- end of if(is.EnvFixed) ---#
							
							means.out <- rbind(means.out, temp.means);
							var.table <- cbind(var.table, temp.var.table);
							p.table <- rbind(p.table, temp.p.table);
						} #--- end of for(k in 1:marker.number)---#
						
						rownames(means.out) <- NULL;
						colnames(var.table) <- NULL;
						rownames(p.table) <- NULL;
						
						colnames(means.out) <- c("Marker","Env", "mm", "Mm", "MM");
						if(is.EnvFixed) 
						{
							rownames(var.table) <- c("Marker", "M SS(P.value)", "E SS (P.value)", "MxE SS (P.value)", "Residuals SS");
						} else
						{
							rownames(var.table) <- c("Marker", "M SS(P.value)");
						}
						colnames(p.table) <- c("Marker", "P.value(M)", "P.value(MxE)")
						
						phenotypicData$traits[[j]]$analysis$sma$means.table <- means.out;
						phenotypicData$traits[[j]]$analysis$sma$var.table <- var.table;
						phenotypicData$traits[[j]]$analysis$sma$p.table <- p.table;
						
						detach("package:lme4", unload=TRUE);
					} else
					{
						#--- if not include environment, separated process single marker analysis ---#
						#--- and then merge all the outcomes into table---#
						for(k in 1:length(phenotypicData$traits[[j]]$envs))
						{
							env.name <- phenotypicData$traits[[j]]$envs[[k]]$name;
							geno.name.phenotypicData <- phenotypicData$traits[[j]]$envs[[k]]$design$geno;
							geno.name.genotypicData <- genotypicData$geno.name;
							data <- merge(phenotypicData$traits[[j]]$envs[[k]]$data, genotypicData.restricted, by.x = geno.name.phenotypicData, by.y = geno.name.genotypicData);
							data[,geno.name.phenotypicData] <- as.factor(data[,geno.name.phenotypicData]);
							data[,trait.name] <- as.numeric(as.character(data[,trait.name]));
							#---do each marker anlysis ---#
							for(n in 1:marker.number)
							{
								marker.name <- marker.names[n];
								myformula <- c();
								myformula <- paste(trait.name, " ~ 1 + ", marker.name, sep = "");
								model <- try(lm(formula(myformula), data = data), silent = TRUE);
								if(!is.null(model) && class(model) == "try-error")
								{
									warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
									next;
								}
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$name <- marker.name;
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$myformula <- myformula;
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$model <- model;
								
								adjust.means <- lsmeans(model, specs = marker.name, weights = "cells");
								adjust.means <- summary(adjust.means)[,c(1:3)];
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$adjust.means <- adjust.means;
								rp.means <- if(rp.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == rp.code, 2], digits) else NA;
								ht.means <- if(ht.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == ht.code, 2], digits) else NA;
								dp.means <- if(dp.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == dp.code, 2], digits) else NA;
								temp.means <- c(marker.name, rp.means, ht.means, dp.means, env.name);
								#--- a formatted variance table ---#	
								model.aov <- anova(model);
								model.aov.rownames <- rownames(model.aov[1]);
								marker.sq <- if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Sum Sq", digits) else NA;
								marker.p <-  if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Pr(>F)", digits) else NA;
								residuals.sq <- if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Sum Sq", digits) else NA;
								residuals.p <-  if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Pr(>F)", digits) else NA; 
								
								temp.var.table <- rbind(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""), residuals.sq, env.name);
								temp.p.table <- c(marker.name, marker.p, env.name);
								
								means.out <- rbind(means.out, temp.means);
								var.table <- cbind(var.table, temp.var.table);
								p.table <- rbind(p.table, temp.p.table);
								
							} #--- end stmt of for(n in 1:marker.number)---#
						} #--- end stmt of for(k in 1:length(phenotypicData$traits[[j]]$envs)) ---#
						rownames(means.out) <- NULL;
						colnames(var.table) <- NULL;
						rownames(p.table) <- NULL;
						
						colnames(means.out) <- c("Marker", "mm", "Mm", "MM", "Env");
						rownames(var.table) <- c("Marker", "M SS (p.value)", "Residuals SS", "Env");
						var.table <- t(var.table);
						colnames(p.table) <- c("Marker", "P.value(M)", "Env");
					} #--end stmt of if(include.env) ---#
					phenotypicData$traits[[j]]$analysis$sma$means.table <- means.out;
					phenotypicData$traits[[j]]$analysis$sma$var.table <- var.table;
					phenotypicData$traits[[j]]$analysis$sma$p.table <- p.table;
				} else
				{
					next;
				} #--end stmt of if(identical(resp.var[i], phenotypicData$traits[[j]]$name))---#
			} #--- end stmt of for( j in 1:length(phenotypicData$traits)) ---#
		} #--- end stmt of for(i in 1:length(resp.var)---#
		
		phenotypicData$genotypicData <- genotypicData;
		#detach("package:lsmeans", unload = TRUE);
		options(old.options);
		class(phenotypicData) <- c("SingleMarkerAnalysis", class(phenotypicData));
		return(phenotypicData);
	} else
	{
		stop("\tError: It could not process single marker analysis on raw phenotypic data!\n");
	} #--- end stmt of if(phenotypicData$isMean)---#
}

#--- It is used to process single marker analysis on SingleEnvAnalysis---#
doSMA.SingleEnvAnalysis <- function(
		phenotypicData,
		genotypicData,
#		geno,
#		env,
		resp.var,
		include.env = TRUE,
		is.EnvFixed = TRUE,
		include.ht = FALSE,
#		test = c("F", "Chisq"),
		digits = 4,
		...
)
{
	#--- checking resp.var whether exists in the phenotypic data, ortherwise using all traits in the phenotypicData ---#
	if(!missing(resp.var))
	{
		trait.names <- phenotypicData$trait.names;
		temp.resp.var <- c();
		for(i in 1:length(resp.var))
		{
			if(resp.var[i] %in% trait.names)
			{
				temp.resp.var <- c(temp.resp.var, resp.var[i]);
			} else
			{
				warnings(paste("\tWarnings: It will be omitted the ", resp.var[i], " because of it is not exist!\n", sep = ""));
			}
		}
		if(length(temp.resp.var) == 0)
			stop("\tError: All the specified resp.var are not exists in the phenotypicdata!\n");
		resp.var <- temp.resp.var;
	} else
	{
		resp.var <- phenotypicData$trait.names;
	}
	
	if(!genotypicData$isRestricted)
	{
		warnings("\tWarning: It will be used default argument to restrict genotypic data!\n");
		genotypicData <- restrict.geno.data(genotypicData);
	}
	
	#---getting all the requried informations---#
	genotypicData.restricted <- genotypicData$restricted$data;
	marker.number <- ncol(genotypicData.restricted) - 1;
	marker.names <- colnames(genotypicData.restricted)[-1];
	#---suppress warn message---#
	old.options <- options();
	options(warn = -1);
	#library("lsmeans");
	
	#--- reformat all genotypic data using dp.code instead of ht.code basing on include.ht ----#
	#--- Because it will recode to default coding sytem, it could be define such code directly!---#
	dp.code <- 2;
	rp.code <- 0;
	ht.code <- 1;
	na.code <- NA;
	if(include.ht)
	{
		genotypicData.restricted[,-1][genotypicData[,-1] == ht.code] <- dp.code;
	}
	
	genotypicData.restricted <- apply(genotypicData.restricted, 2, factor);
	#--- checking whether is single environment---#
	if(include.env && phenotypicData$isSingleEnv)
	{	
		include.env <- FALSE;
		warning("\tWarnings: The phenotypic data is single environmental data and The is.EnvFixed paramter will be omitted!\n");
	}
	
	for(i in 1:length(resp.var))
	{
		#--- checking whether the resp.var name is the trait name of phenotypicData---#
		for(j in 1:length(phenotypicData$traits))
		{
			if(identical(resp.var[i], phenotypicData$traits[[j]]$name))
			{
				
				trait.name <- resp.var[i];
				#--- checking whether this trait is done single environment analysis, if not, try next ---#
				if(is.null(phenotypicData$traits[[j]]$analysis$sea))
				{
					warning("\tThis trait", trait.name, "does not conduct single environment analysis yet! It will be omitted!\n");
					next;
				}
				#--- adding Mean postfix on trait.name, becuase of single environment analysis sum out of data ---#
				trait.name <- paste(trait.name, "_Mean", sep = "");
				#--- checking the analysis structure of phenotypicData--#
				if(is.null(phenotypicData$traits[[j]]$analysis))
					phenotypicData$traits[[j]]$analysis <- list();
				if(is.null(phenotypicData$traits[[j]]$analysis$sma))
				{
					phenotypicData$traits[[j]]$analysis$sma <- list();
				} else
				{
					phenotypicData$traits[[j]]$analysis$sma <- NULL;
					phenotypicData$traits[[j]]$analysis$sma <- list();
				}
				phenotypicData$traits[[j]]$analysis$sma$include.env <- include.env;
				phenotypicData$traits[[j]]$analysis$sma$is.EnvFixed <- is.EnvFixed;
				phenotypicData$traits[[j]]$analysis$sma$include.ht <- include.ht;
				phenotypicData$traits[[j]]$analysis$sma$marker.number <- marker.number;
				
				means.out <- c();
				var.table <- c();
				p.table <- c();
				
				#---according to whether including env to determinate how to analysis single marker ---#
				if(include.env)
				{
					#--- combined all single environmental analysis data into one data frame---#
					temp.phenotypicData <- c();
					env <- c();
					for(k in 1:length(phenotypicData$traits[[j]]$analysis$sea$envs))
					{
						temp.phenotypicData <- rbind(temp.phenotypicData, phenotypicData$traits[[j]]$analysis$sea$envs[[k]]$sum.out);
						if(length(env) == 0)
							env <- phenotypicData$traits[[j]]$envs[[k]]$design$env;
					}
					#--- merge phenotypicData and genotoypicData into one data frame ---#
					geno.name.phenotypicData <- phenotypicData$traits[[j]]$envs[[1]]$design$geno;
					geno.name.genotypicData <- genotypicData$geno.name;
					data <- merge(temp.phenotypicData, genotypicData.restricted, by.x = geno.name.phenotypicData, by.y = geno.name.genotypicData);
					data[,env] <- as.factor(data[,env]);
					data[,trait.name] <- as.numeric(as.character(data[,trait.name]));
					
				#	library("lme4");
					
					for(k in 1:marker.number)
					{
						marker.name <- marker.names[k];
						myformula <- c();
						
						if(is.EnvFixed)
						{
							myformula <- paste(trait.name, " ~ 1 + ", marker.name, " + ", env, " + ", marker.name, " : ", env, sep = "");
							model <- try(lm(formula(myformula), data = data), silent = TRUE);
							if(!is.null(model) && class(model) == "try-error")
							{
								warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
								next;
							}
							adjust.means <- lsmeans(model, specs = marker.name, by = env, weights = "cells");
							adjust.means <- summary(adjust.means)[ , c(1:3)];
							env.levels <- levels(adjust.means[,env]);
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
							myformula <- paste(trait.name, " ~ 1 + ", marker.name, " + ", "(1|", env, ") + (1|", marker.name, " : ", env, ") ", sep = "");
							model <- try(lmer(formula(myformula), data = data), silent = TRUE);
							if(!is.null(model) && class(model) == "try-error")
							{
								warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
								next;
							}
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
							#result$traits[[i]]$markers[[j]]$adjust.means <- adjust.means;
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
							marker.p <- signif(model.aov$"F value", digits);
							model.com.aov <- anova(model2, model);
							marker.env.p <- signif(model.com.aov[[8]][2], digits);
							temp.var.table <- c(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""));
							temp.p.table <- c(marker.name, marker.p, marker.env.p);
						} #--- end of if(is.EnvFixed) ---#
						
						means.out <- rbind(means.out, temp.means);
						var.table <- cbind(var.table, temp.var.table);
						p.table <- rbind(p.table, temp.p.table);
					} #--- end of for(k in 1:marker.number)---#
					
					rownames(means.out) <- NULL;
					colnames(var.table) <- NULL;
					rownames(p.table) <- NULL;
					
					colnames(means.out) <- c("Marker","Env", "mm", "Mm", "MM");
					if(is.EnvFixed) 
					{
						rownames(var.table) <- c("Marker", "M SS(P.value)", "E SS (P.value)", "MxE SS (P.value)", "Residuals SS");
					} else
					{
						rownames(var.table) <- c("Marker", "M SS(P.value)");
					}
					colnames(p.table) <- c("Marker", "P.value(M)", "P.value(MxE)")
					
					phenotypicData$traits[[j]]$analysis$sma$means.table <- means.out;
					phenotypicData$traits[[j]]$analysis$sma$var.table <- var.table;
					phenotypicData$traits[[j]]$analysis$sma$p.table <- p.table;
					
				#	detach("package:lme4", unload=TRUE);
				} else
				{
					#--- if not include environment, separated process single marker analysis ---#
					#--- and then merge all the outcomes into table---#
					for(k in 1:length(phenotypicData$traits[[j]]$analysis$sea$envs))
					{
						env.name <- phenotypicData$traits[[j]]$analysis$sea$envs[[k]]$name;
						geno.name.phenotypicData <- phenotypicData$traits[[j]]$envs[[k]]$design$geno;
						geno.name.genotypicData <- genotypicData$geno.name;
						data <- merge(phenotypicData$traits[[j]]$analysis$sea$envs[[k]]$sum.out, genotypicData.restricted, by.x = geno.name.phenotypicData, by.y = geno.name.genotypicData);
						data[,geno.name.phenotypicData] <- as.factor(data[,geno.name.phenotypicData]);
						data[,trait.name] <- as.numeric(as.character(data[,trait.name]));
						#---do each marker anlysis ---#
						for(n in 1:marker.number)
						{
							marker.name <- marker.names[n];
							myformula <- c();
							myformula <- paste(trait.name, " ~ 1 + ", marker.name, sep = "");
							model <- try(lm(formula(myformula), data = data), silent = TRUE);
							if(!is.null(model) && class(model) == "try-error")
							{
								warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
								next;
							}
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$name <- marker.name;
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$myformula <- myformula;
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$model <- model;
							
							adjust.means <- lsmeans(model, specs = marker.name, weights = "cells");
							adjust.means <- summary(adjust.means)[,c(1:3)];
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$adjust.means <- adjust.means;
							rp.means <- if(rp.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == rp.code, 2], digits) else NA;
							ht.means <- if(ht.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == ht.code, 2], digits) else NA;
							dp.means <- if(dp.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == dp.code, 2], digits) else NA;
							temp.means <- c(marker.name, rp.means, ht.means, dp.means, env.name);
							#--- a formatted variance table ---#	
							model.aov <- anova(model);
							model.aov.rownames <- rownames(model.aov[1]);
							marker.sq <- if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Sum Sq", digits) else NA;
							marker.p <-  if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Pr(>F)", digits) else NA;
							residuals.sq <- if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Sum Sq", digits) else NA;
							residuals.p <-  if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Pr(>F)", digits) else NA; 
							
							temp.var.table <- rbind(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""), residuals.sq, env.name);
							temp.p.table <- c(marker.name, marker.p, env.name);
							
							means.out <- rbind(means.out, temp.means);
							var.table <- cbind(var.table, temp.var.table);
							p.table <- rbind(p.table, temp.p.table);
							
						} #--- end stmt of for(n in 1:marker.number)---#
					} #--- end stmt of for(k in 1:length(phenotypicData$traits[[j]]$envs)) ---#
					rownames(means.out) <- NULL;
					colnames(var.table) <- NULL;
					rownames(p.table) <- NULL;
					
					colnames(means.out) <- c("Marker", "mm", "Mm", "MM", "Env");
					rownames(var.table) <- c("Marker", "M SS (p.value)", "Residuals SS", "Env");
					var.table <- t(var.table);
					colnames(p.table) <- c("Marker", "P.value(M)", "Env");
				} #--end stmt of if(include.env) ---#
				phenotypicData$traits[[j]]$analysis$sma$means.table <- means.out;
				phenotypicData$traits[[j]]$analysis$sma$var.table <- var.table;
				phenotypicData$traits[[j]]$analysis$sma$p.table <- p.table;
			} else
			{
				next;
			} #--end stmt of if(identical(resp.var[i], phenotypicData$traits[[j]]$name))---#
		} #--- end stmt of for( j in 1:length(phenotypicData$traits)) ---#
	} #--- end stmt of for(i in 1:length(resp.var)---#
	
	phenotypicData$genotypicData <- genotypicData;
	#detach("package:lsmeans", unload = TRUE);
	options(old.options);
	class(phenotypicData) <- c("SingleMarkerAnalysis", class(phenotypicData));
	return(phenotypicData);
	
}

#--- It is used to process single marker analysis on MultiEnvAnalysis---#
doSMA.MultiEnvAnalysis <- function(
		phenotypicData,
		genotypicData,
#		geno,
#		env,
		resp.var,
		include.env = TRUE,
		is.EnvFixed = TRUE,
		include.ht = FALSE,
#		test = c("F", "Chisq"),
		digits = 4,
		...
)
{
	#--- checking resp.var whether exists in the phenotypic data, ortherwise using all traits in the phenotypicData ---#
	if(!missing(resp.var))
	{
		trait.names <- phenotypicData$trait.names;
		temp.resp.var <- c();
		for(i in 1:length(resp.var))
		{
			if(resp.var[i] %in% trait.names)
			{
				temp.resp.var <- c(temp.resp.var, resp.var[i]);
			} else
			{
				warnings(paste("\tWarnings: It will be omitted the ", resp.var[i], " because of it is not exist!\n", sep = ""));
			}
		}
		if(length(temp.resp.var) == 0)
			stop("\tError: All the specified resp.var are not exists in the phenotypicdata!\n");
		resp.var <- temp.resp.var;
	} else
	{
		resp.var <- phenotypicData$trait.names;
	}
	
	if(!genotypicData$isRestricted)
	{
		warnings("\tWarning: It will be used default argument to restrict genotypic data!\n");
		genotypicData <- restrict.geno.data(genotypicData);
	}
	
	#---getting all the requried informations---#
	genotypicData.restricted <- genotypicData$restricted$data;
	marker.number <- ncol(genotypicData.restricted) - 1;
	marker.names <- colnames(genotypicData.restricted)[-1];
	#---suppress warn message---#
	old.options <- options();
	options(warn = -1);
#	library("lsmeans");
	
	#--- reformat all genotypic data using dp.code instead of ht.code basing on include.ht ----#
	#--- Because it will recode to default coding sytem, it could be define such code directly!---#
	dp.code <- 2;
	rp.code <- 0;
	ht.code <- 1;
	na.code <- NA;
	if(include.ht)
	{
		genotypicData.restricted[,-1][genotypicData[,-1] == ht.code] <- dp.code;
	}
	genotypicData.restricted <- apply(genotypicData.restricted, 2, factor);
	
	#--- checking whether is single environment---#
	if(include.env && phenotypicData$isSingleEnv)
	{	
		include.env <- FALSE;
		warning("\tWarnings: The phenotypic data is single environmental data and The is.EnvFixed paramter will be omitted!\n");
	}
	
	for(i in 1:length(resp.var))
	{
		#--- checking whether the resp.var name is the trait name of phenotypicData---#
		for(j in 1:length(phenotypicData$traits))
		{
			if(identical(resp.var[i], phenotypicData$traits[[j]]$name))
			{
				trait.name <- resp.var[i];
				#--- check whether this trait is conducted multi environment analysis---#
				if(is.null(phenotypicData$traits[[j]]$analysis$mea))
				{
					warning("\tWarning: This trait", trait.name, "does not conduct multi environment analysis! It will be omitted!\n");
					next;
				}
				trait.name <- paste(trait.name, "_LSMean", sep = ""	);
				#--- checking the analysis structure of phenotypicData--#
				if(is.null(phenotypicData$traits[[j]]$analysis))
					phenotypicData$traits[[j]]$analysis <- list();
				if(is.null(phenotypicData$traits[[j]]$analysis$sma))
				{
					phenotypicData$traits[[j]]$analysis$sma <- list();
				} else
				{
					phenotypicData$traits[[j]]$analysis$sma <- NULL;
					phenotypicData$traits[[j]]$analysis$sma <- list();
				}
				phenotypicData$traits[[j]]$analysis$sma$include.env <- include.env;
				phenotypicData$traits[[j]]$analysis$sma$is.EnvFixed <- is.EnvFixed;
				phenotypicData$traits[[j]]$analysis$sma$include.ht <- include.ht;
				phenotypicData$traits[[j]]$analysis$sma$marker.number <- marker.number;
				
				means.out <- c();
				var.table <- c();
				p.table <- c();
				
				#---according to whether including env to determinate how to analysis single marker ---#
				if(include.env)
				{
					#--- combined all environmental data into one data frame---#
					temp.phenotypicData <- c();
					env <- c();
#					for(k in 1:length(phenotypicData$traits[[j]]$envs))
#					{
#						temp.phenotypicData <- rbind(temp.phenotypicData, phenotypicData$traits[[j]]$analysis$mea$data);
#						if(length(env) == 0)
#							env <- phenotypicData$traits[[j]]$envs[[k]]$design$env;
#					}
					temp.phenotypicData <- phenotypicData$traits[[j]]$analysis$mea$means.GenoEnv;
					env <- phenotypicData$traits[[j]]$envs[[1]]$design$env;
					
					#--- merge phenotypicData and genotoypicData into one data frame ---#
					geno.name.phenotypicData <- phenotypicData$traits[[j]]$envs[[1]]$design$geno;
					geno.name.genotypicData <- genotypicData$geno.name;
					data <- merge(temp.phenotypicData, genotypicData.restricted, by.x = geno.name.phenotypicData, by.y = geno.name.genotypicData);
					data[,trait.name] <- as.numeric(as.character(data[,trait.name]));
					data[,env] <- as.factor(data[,env]);
					
#					library("lme4");
					
					for(k in 1:marker.number)
					{
						marker.name <- marker.names[k];
						myformula <- c();
						
						if(is.EnvFixed)
						{
							myformula <- paste(trait.name, " ~ 1 + ", marker.name, " + ", env, " + ", marker.name, " : ", env, sep = "");
							model <- try(lm(formula(myformula), data = data), silent = TRUE);
							if(!is.null(model) && class(model) == "try-error")
							{
								warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
								next;
							}
							adjust.means <- lsmeans(model, specs = marker.name, by = env, weights = "cells");
							adjust.means <- summary(adjust.means)[ , c(1:3)];
							env.levels <- levels(adjust.means[,env]);
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
							myformula <- paste(trait.name, " ~ 1 + ", marker.name, " + ", "(1|", env, ") + (1|", marker.name, " : ", env, ") ", sep = "");
							model <- try(lmer(formula(myformula), data = data), silent = TRUE);
							if(!is.null(model) && class(model) == "try-error")
							{
								warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
								next;
							}
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
							#result$traits[[i]]$markers[[j]]$adjust.means <- adjust.means;
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
							marker.p <- signif(model.aov$"F value", digits);
							model.com.aov <- anova(model2, model);
							marker.env.p <- signif(model.com.aov[[8]][2], digits);
							temp.var.table <- c(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""));
							temp.p.table <- c(marker.name, marker.p, marker.env.p);
						} #--- end of if(is.EnvFixed) ---#
						
						means.out <- rbind(means.out, temp.means);
						var.table <- cbind(var.table, temp.var.table);
						p.table <- rbind(p.table, temp.p.table);
					} #--- end of for(k in 1:marker.number)---#
					
					rownames(means.out) <- NULL;
					colnames(var.table) <- NULL;
					rownames(p.table) <- NULL;
					
					colnames(means.out) <- c("Marker","Env", "mm", "Mm", "MM");
					if(is.EnvFixed) 
					{
						rownames(var.table) <- c("Marker", "M SS(P.value)", "E SS (P.value)", "MxE SS (P.value)", "Residuals SS");
					} else
					{
						rownames(var.table) <- c("Marker", "M SS(P.value)");
					}
					colnames(p.table) <- c("Marker", "P.value(M)", "P.value(MxE)")
					
					phenotypicData$traits[[j]]$analysis$sma$means.table <- means.out;
					phenotypicData$traits[[j]]$analysis$sma$var.table <- var.table;
					phenotypicData$traits[[j]]$analysis$sma$p.table <- p.table;
					
#					detach("package:lme4", unload=TRUE);
				} else
				{
					#--- if not include environment, separated process single marker analysis ---#
					#--- and then merge all the outcomes into table---#
					data.total <- phenotypicData$traits[[j]]$analysis$mea$means.GenoEnv;
					env <- phenotypicData$traits[[j]]$envs[[1]]$design$env;
					data.total[, env] <- as.factor(data.total[, env]);
					nlevels.env <- nlevels(data.total[,env]);
					names.env <- as.character(levels(data.total[,env]));
					for(k in 1:nlevels.env)
					{
						env.name <- names.env[k];
						geno.name.phenotypicData <- phenotypicData$traits[[j]]$envs[[1]]$design$geno;
						geno.name.genotypicData <- genotypicData$geno.name;
						temp.data <- subset(data.total, subset = (data.total[,env] == env.name));
						data <- merge(temp.data, genotypicData.restricted, by.x = geno.name.phenotypicData, by.y = geno.name.genotypicData);
						data[,geno.name.phenotypicData] <- as.factor(data[,geno.name.phenotypicData]);
						data[,trait.name] <- as.numeric(as.character(data[,trait.name]));
						#---do each marker anlysis ---#
						for(n in 1:marker.number)
						{
							marker.name <- marker.names[n];
							myformula <- c();
							myformula <- paste(trait.name, " ~ 1 + ", marker.name, sep = "");
							model <- try(lm(formula(myformula), data = data), silent = TRUE);
							if(!is.null(model) && class(model) == "try-error")
							{
								warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
								next;
							}
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$name <- marker.name;
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$myformula <- myformula;
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$model <- model;
							
							adjust.means <- lsmeans(model, specs = marker.name, weights = "cells");
							adjust.means <- summary(adjust.means)[,c(1:3)];
#								phenotypicData$traits[[j]]$analysis$sma$markers[[n]]$adjust.means <- adjust.means;
							rp.means <- if(rp.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == rp.code, 2], digits) else NA;
							ht.means <- if(ht.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == ht.code, 2], digits) else NA;
							dp.means <- if(dp.code %in% adjust.means[ , marker.name]) signif(adjust.means[adjust.means[ , marker.name] == dp.code, 2], digits) else NA;
							temp.means <- c(marker.name, rp.means, ht.means, dp.means, env.name);
							#--- a formatted variance table ---#	
							model.aov <- anova(model);
							model.aov.rownames <- rownames(model.aov[1]);
							marker.sq <- if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Sum Sq", digits) else NA;
							marker.p <-  if(marker.name %in% model.aov.rownames) signif(model.aov[model.aov.rownames == marker.name, ]$"Pr(>F)", digits) else NA;
							residuals.sq <- if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Sum Sq", digits) else NA;
							residuals.p <-  if("Residuals" %in% model.aov.rownames) signif(model.aov[model.aov.rownames == "Residuals", ]$"Pr(>F)", digits) else NA; 
							
							temp.var.table <- rbind(marker.name, paste(marker.sq, " (", marker.p, ")", sep = ""), residuals.sq, env.name);
							temp.p.table <- c(marker.name, marker.p, env.name);
							
							means.out <- rbind(means.out, temp.means);
							var.table <- cbind(var.table, temp.var.table);
							p.table <- rbind(p.table, temp.p.table);
							
						} #--- end stmt of for(n in 1:marker.number)---#
					} #--- end stmt of for(k in 1:length(phenotypicData$traits[[j]]$envs)) ---#
					rownames(means.out) <- NULL;
					colnames(var.table) <- NULL;
					rownames(p.table) <- NULL;
					
					colnames(means.out) <- c("Marker", "mm", "Mm", "MM", "Env");
					rownames(var.table) <- c("Marker", "M SS (p.value)", "Residuals SS", "Env");
					var.table <- t(var.table);
					colnames(p.table) <- c("Marker", "P.value(M)", "Env");
				} #--end stmt of if(include.env) ---#
				phenotypicData$traits[[j]]$analysis$sma$means.table <- means.out;
				phenotypicData$traits[[j]]$analysis$sma$var.table <- var.table;
				phenotypicData$traits[[j]]$analysis$sma$p.table <- p.table;
			} else
			{
				next;
			} #--end stmt of if(identical(resp.var[i], phenotypicData$traits[[j]]$name))---#
		} #--- end stmt of for( j in 1:length(phenotypicData$traits)) ---#
	} #--- end stmt of for(i in 1:length(resp.var)---#
	
	phenotypicData$genotypicData <- genotypicData;
#	detach("package:lsmeans", unload = TRUE);
	options(old.options);
	class(phenotypicData) <- c("SingleMarkerAnalysis", class(phenotypicData));
	return(phenotypicData);
}

print.SingleMarkerAnalysis <- function(
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
	for(i in 1:length(data$traits))
	{
		if(is.null(data$traits[[i]]$analysis$sma))
			next;
		cat("Trait Name:", data$traits[[i]]$name, ".\n", sep = "");
		cat("Including Environment:", ifelse(data$traits[[i]]$analysis$sma$include.env, 
						"TRUE", "FALSE"),".\n", sep = "");
		cat("Environment as Fixed:", ifelse(data$traits[[i]]$analysis$sma$include.env, 
						ifelse(data$traits[[i]]$analysis$sma$isEnvFixed,"TRUE", "FALSE"),
						"NA"), ".\n", sep ="");
		cat("Marker Number:", data$traits[[i]]$analysis$sma$marker.number, ".\n", sep = "");
		cat("\n");
		p.table <- as.data.frame(data$traits[[i]]$analysis$sma$p.table);
		
		means.table <- as.data.frame(data$traits[[i]]$analysis$sma$means.table);
		if(identical(levels, "SIG"))
		{
			if(data$traits[[i]]$analysis$sma$include.env)
			{
#				p.table.m <- p.table[as.numeric(as.character(p.table[,2])) <= p.value, ];
#				p.table.me <- p.table[as.numeric(as.character(p.table[,3])) <= p.value, ];
#				p.table<- merge(p.table.m, p.table.me, by = "Marker", all = TRUE)[, 1:3];
				p.table[ , "P.value(M)"] <- as.numeric(as.character(p.table[ , "P.value(M)"]));
				p.table[ , "P.value(MxE)"] <- as.numeric(as.character(p.table[ , "P.value(MxE)"]));
				p.table <- subset(p.table, p.table$"P.value(M)" <= p.value | p.table$"P.value(MxE)" <= p.value);
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
				#p.table <- p.table[as.numeric(as.character(p.table[,2])) <= p.value, ];
				p.table[ , "P.value(M)"] <- as.numeric(as.character(p.table[ , "P.value(M)"]));
				p.table <- subset(p.table, p.table$"P.value(M)" <= p.value);
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
		#cat("\n")
	}
}



##--- This s3 method is for processed data after single site analysis and using their lsmeans of each lines---#
#doSMA.SingleEnvAnalysis(
#		phenotypicData,
#		genotypicData,
#		include.env = FALSE,
#		is.EnvFixed = TRUE,
#		include.ht = FALSE,
#		test = c("F", "Chisq")
#)
#{
#	if(missing(pheno))
#		stop("\tError: The pheno could not be null!\n");
#	if(missing(geno))
#		stop("\tError: The geno could not be null!\n");
#	if(!inherits(pheno, "SingleSiteAnalysis"))
#		stop("\tError: The pheno should be of class SingleSiteAnalysis!\n");
#	#--- It should be done before do single site analysis---#
#	if(is.null(pheno$doRestricted))
#	{
#		waring(paste("\tWarning: It will use default parameters to restrict phenotypic data!\n",
#						"\tBecause it do not do restrict.pheno.data function!\n", sep=""));
#		pheno <- restrict.pheno.data(pheno);
#	}
#	if(!inherits(geno, "GenotypicData"))
#		stop("\tError: The pheno should be of class GenotypicData!\n");
#	if(is.null(geno$restricted))
#	{
#		warning(paste("\tWarning: It will use default parameters to restrict genotypic data!\n", 
#						"\tBecause it do not do restrict.geno.data function!\n", sep= ""));
#		geno <- restrict.geno.data(geno);
#	}
#	#--- checking the marker.num after restricted should be larger than 1---#
#	#--- if it is not, stop execution ---#
#	#--- because the first column is lines info and the rest is marker info---#
#	#--- retrieved the restricted genotypic data frame ---#
#	genotypicData <- geno$restricted$data;
#	marker.number <- ncol(genotypicData) - 1;
#	marker.names <- colnames(genotypicData)[-1];
#	
#	if(marker.number < 1)
#	{
#		stop("\tWarning: There are no more markers after restricted genotypic data!\n");
#	}
#	if(!is.logical(include.env))
#		stop("\tError: The include.env argument should be of type logical!\n");
#	if(!is.logical(is.EnvFixed))
#		stop("\tError: The is.EnvFixed argument should be of type logical!\n");
#	test <- match.arg(test);
#	
#	trait.number <- length(pheno$traits);
#	if(trait.number < 1)
#		stop("\tError: There is no trait on phenotypic data\n");
#	library("lme4");
#	library("phia");
#	
#	#--- reformat all genotypic data using dp.code instead of ht.code basing on include.ht ----#
#	#--- Because it will recode to default coding sytem, it could be define such code directly!---#
#	dp.code <- 2;
#	rp.code <- 0;
#	ht.code <- 1;
#	na.code <- NA;
#	if(include.ht)
#	{
#		genotypicData[ , -1][genotypic[,-1] == ht.code] <- dp.code;
#	}
#	genotypicData <- apply(genotypicData, 2, factor);
#	
#	for(i in 1:trait.number)
#	{
#		
#		trait.name <- pheno$traits[[i]]$name;
#		site.number <- length(pheno$traits[[i]]$sites);
#		site.number.check <- length(pheno$traits[[i]]$analysis$ssa$sites);
#		geno.phenotypicData <- pheno$traits[[i]]$sites[[1]]$design$geno;
#		env <- pheno$traits[[i]]$sites[[1]]$design$env;
#		geno.genotypicData <- names(genotypicData)[1];	
#		#--- if the single marker analysis element of pheno is not exist, then created it ---#
#		if(is.null(pheno$traits[[i]]$analysis$sma))
#			pheno$traits[[i]]$analysis$sma <- list();
#		pheno$traits[[i]]$analysis$sma$site.used <- site.number.check;
#		pheno$traits[[i]]$analysis$sma$site.read <- site.number;
#		pheno$traits[[i]]$analysis$sma$marerks.used <- marker.number;
#		pheno$traits[[i]]$analysis$sma$markers.names <- marker.names;
#		pheno$traits[[i]]$analysis$sma$markers <- list();
#		pheno$traits[[i]]$analysis$sma$include.env <- include.env;
#		pheno$traits[[i]]$analysis$sma$include.ht <- include.ht;
#		pheno$traits[[i]]$analysis$sma$is.EnvFixed <- is.EnvFixed;
#		
#		
#		if(is.null(site.number.check) || site.number.check < 1)
#		{
#			stop(paste("\tError: There are no any single site analysis outcomes on ", trait.name, "!\n", sep = ""));
#			next;
#		}
#		
#		if(include.env && site.number.check > 1)
#		{
#			#--- according site number and parameter include.env to do single marker analysis with/without env---#
#			#--- combind all not restricted data of each site--- #
#			phenotypicData <- c();		
#			for(j in 1:site.number.check)
#			{	
#				#--- sum.out[ , c(1,2,4)] stands for Genotype, triat_mean, env---#
#				phenotypicData <- rbind(phenotypicData, pheno$traits[[i]]$analysis$ssa$sites[[j]]$sum.out[ , c(1,2,4)]);
#			}			
#			#--- meger marker data with phenotypic data one by one ---#
#			data <- merge(phenotypicdata, genotypicData, by.x = geno.phenotypicData, by.y = geno.genotypicData);
#			
#			for(k in 1:marker.number)
#			{	
#				marker.name <- marker.names[k];
#				env.stmt <- paste(env ," + ", marker.name, " : ", env, sep ="");
#				if(!is.EnvFixed)
#				{
#					env.stmt <- paste("(1|", env ,") + (1|", marker.name, " : ", env, ") ",sep ="");
#				}
#				#myformula1 <- paste(trait.name, " ~ 1 + ", marker.name, " + (1|", marker.name, " : ", geno, ") + ", env.stmt, sep = "");
#				myformula1 <- paste(trait.name, " ~ 1 + ", marker.name, " + ", env.stmt, sep = "");
#				
#				model <- try(lmer(formula(myformula1), data = data), silent = TRUE);
#				if(!is.null(model) && class(model == "try-error"))
#				{
#					warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
#					next;
#				}
#				
#				pheno$traits[[i]]$analysis$sma$markers[[k]] <- list();
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$name <- marker.name;
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites <- list();
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$model <- model;
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$formula <- myformula1;
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$adjust.means <- interactionMeans(model, factor= marker.name);
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$marker.test <- testInteractions(model, pairwise = marker.name, test = test);
#				if(is.EnvFixed)
#				{
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$marker.across.env <- testInteractions(model, fixed = env, across = marker.name, test = test);
#				} else
#				{
#					model2 <- update(mode, formula(paste(". ~ . - ", env.stmt, sep = "")));
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$model2 <- model2
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$marker.across.env <- anova(model, model2);					
#				}
#				
#			} #--end of statment for(k in 1:marker.number) ---#
#		} else
#		{
#			if(include.env && site.number.check == 1)
#				warning("\tWarning: The include.env argument would be omitted because of only one site available!\n");
#			for(j in 1:site.number.check)
#			{
#				phenotypicData <- pheno$traits[[i]]$analysis$ssa$sites[[j]]$sum.out[ , c(1,2,4)];
#				
#				#--- meger marker data with phenotypic data one by one ---#
#				data <- merge(phenotypicdata, genotypicData, by.x = geno.phenotypicData, by.y = geno.genotypicData);
#				for(k in 1:marker.number)
#				{	
#					marker.name <- marker.names[k];
#					#myformula1 <- paste(trait.name, " ~ 1 + ", marker.name, " + (1|", marker.name, " : ", geno, ")", sep = "");
#					myformula1 <- paste(trait.name, " ~ 1 + ", marker.name, sep = "");
#					
#					model <- try(lmer(formula(myformula1), data = data), silent = TRUE);
#					if(!is.null(model) && class(model == "try-error"))
#					{
#						warning(paste("\tWarning: There are error occured on this marker, i.e. ", marker.name, ".\n", sep = ""));
#						next;
#					}
#					
#					pheno$traits[[i]]$analysis$sma$markers[[k]] <- list();
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$name <- marker.name;
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites <- list();
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$model <- model;
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$formula <- myformula1;
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$adjust.means <- interactionMeans(model, factor= marker.name);
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$marker.test <- testInteractions(model, pairwise = marker.name, test = test);
#					if(is.EnvFixed)
#					{
#						pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$marker.across.env <- testInteractions(model, fixed = env, across = marker.name, test = test);
#					} else
#					{
#						model2 <- update(mode, formula(paste(". ~ . - ", env.stmt, sep = "")));
#						pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$model2 <- model2
#						pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$marker.across.env <- anova(model, model2);					
#					}
#					
#				} #--end of statment for(k in 1:marker.number) ---#
#				
#			}
#		} #--- end of statement if(include.env && site.number.check > 1) ---#
#		
#	}#--- end of statement for( i in 1:trait.number) ---#
#}
#
##--- This s3 method is for raw data after read the phenotypic data and using their raw observations---#
#doSMA.PhenotypicData(
#		phenotypicData,
#		genotypicData,
#		include.env = FALSE,
#		is.EnvFixed = TRUE,
#		include.ht = FALSE,
#		test = c("F", "Chisq")
#)
#{
#	if(missing(phenotypicData))
#		stop("\tError: The phenotypicData could not be null!\n");
#	if(missing(genotypicData))
#		stop("\tError: The genotypicData could not be null!\n");
#	if(!inherits(phenotypicData, "PhneotypicData"))
#		stop("\tError: The pheno should be of class PhenotypicData!\n");
#	if(phenotypicData$doRestricted)
#	{
#		waring(paste("\tWarning: It will use default parameters to restrict phenotypic data!\n",
#						"\tBecause it do not have restrict.pheno.data function!\n", sep=""));
#		phenotypicData <- restrict.pheno.data(phenotypicData);
#	}
#	if(!inherits(geno, "GenotypicData"))
#		stop("\tError: The pheno should be of class GeneticData!\n");
#	if(is.null(geno$restricted))
#	{
#		warning(paste("\tWarning: It will use default parameters to restrict genotypic data!\n", 
#						"\tBecause it do not do restrict.geno.data function!\n", sep= ""));
#		geno <- restrict.geno.data(geno);
#	}
#	#--- checking the marker.num after restricted should be larger than 1---#
#	#--- if it is not, stop execution ---#
#	#--- because the first column is lines info and the rest is marker info---#
#	#--- retrieved the restricted genotypic data frame ---#
#	genotypicData <- geno$restricted$data;
#	marker.number <- ncol(genotypicData) - 1;
#	marker.names <- colnames(genotypicData)[-1];
#	
#	if(marker.number < 1)
#	{
#		stop("\tWarning: There are no more markers after restricted genotypic data!\n");
#	}
#	if(!is.logical(include.env))
#		stop("\tError: The include.env argument should be of type logical!\n");
#	if(!is.logical(is.EnvFixed))
#		stop("\tError: The is.EnvFixed argument should be of type logical!\n");
#	test <- match.arg(test);
#	
#	trait.number <- length(pheno$traits);
#	if(trait.number < 1)
#		stop("\tError: There is no trait on phenotypic data\n");
#	library("lme4");
#	library("phia");
#	
#	#--- reformat all genotypic data using dp.code instead of ht.code basing on include.ht ----#
#	#--- Because it will recode to default coding sytem, it could be define such code directly!---#
#	dp.code <- 2;
#	rp.code <- 0;
#	ht.code <- 1;
#	na.code <- NA;
#	if(include.ht)
#	{
#		genotypicData[ , -1][genotypic[,-1] == ht.code] <- dp.code;
#	}
#	genotypicData <- apply(genotypicData, 2, factor);
#	
#	for(i in 1:trait.number)
#	{
#		#--- checking whether the site number is large than 1, if not, do next trait---#
#		site.number <- length(pheno$traits[[i]]$sites);	
#		if(site.number < 1)
#		{
#			warning(paste("\tWarning: There is no available site on ", trait.name, " after restricted!\n", sep = ""));
#			next;
#		}
#		#--- get the available site number after restricted---#
#		site.number.check <- site.number
#		for(j in 1:site.number)
#		{
#			if(pheno$traits[[i]]$sites[[j]]$restricted$isTRUE)
#				site.number.check <- site.number.check -1;
#		}
#		#--- rechecking whether the site number is large than 1 after restricted ---#
#		if(site.number.check < 1)
#		{
#			warning(paste("\tWarning: There is no available site on ", trait.name,  " after restricted!\n", sep = ""));
#			next;
#		}
#		
#		#--- if the analysis element of pheno is not exist, created it---#
#		if(is.null(pheno$traits[[i]]$analysis))
#			pheno$traits[[i]]$analysis <- list();
#		#--- if the single marker analysis element of pheno is not exist, then created it ---#
#		if(is.null(pheno$traits[[i]]$analysis$sma))
#			pheno$traits[[i]]$analysis$sma <- list();
#		trait.name <- pheno$traits$name;
#		geno.phenotypicData <- pheno$traits[[i]]$sites[[1]]$design$geno;
#		exptl.design <- pheno$traits[[i]]$sites[[1]]$design$exptl.design;
#		env <- pheno$traits[[i]]$sites[[1]]$design$env;
#		geno.genotypicData <- names(genotypicData)[1];
#		pheno$traits[[i]]$analysis$sma$site.used <- site.number.check;
#		pheno$traits[[i]]$analysis$sma$site.read <- site.number;
#		pheno$traits[[i]]$analysis$sma$marerks.used <- marker.number;
#		pheno$traits[[i]]$analysis$sma$markers.names <- marker.names;
#		pheno$traits[[i]]$analysis$sma$markers <- list();
#		pheno$traits[[i]]$analysis$sma$include.env <- include.env;
#		pheno$traits[[i]]$analysis$sma$include.ht <- include.ht;
#		pheno$traits[[i]]$analysis$sma$is.EnvFixed <- is.EnvFixed;
#		
#		if(include.env && (site.number > 1))
#		{	
#			#--- according site number and parameter include.env to do single marker analysis with/without env---#
#			#--- combind all not restricted data of each site--- #
#			phenotypicData <- c();	
#			for(j in 1:site.number)
#			{
#				if(!pheno$traits[[i]]$sites[[j]]$restricted$isTRUE)
#					phenotypicData <- rbind(phenotypicData, pheno$traits[[i]]$sites[[j]]$data);
#			}
#			
#			#--- meger marker data with phenotypic data one by one ---#
#			data <- merge(phenotypicdata, genotypicData, by.x = geno.phenotypicData, by.y = geno.genotypicData);
#			
#			for(k in 1:marker.number)
#			{	
#				marker.name <- marker.names[k];
#				pheno$traits[[i]]$analysis$sma$markers[[k]] <- list();
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$name <- marker.name;
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites <- list();
#				env.stmt <- paste(env ," + ", marker.name, " : ", env, sep ="");
#				if(!is.EnvFixed)
#				{
#					env.stmt <- paste("(1|", env ,") + (1|", marker.name, " : ", env, ") ",sep ="");
#				}
#				if(exptl.design %in% c("RCB", "AugRCB"))
#				{
#					block <- pheno$traits[[i]]$sites[[1]]$design$block;
#					myformula1 <- paste(trait.name,  " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno.phenotypicData, ") + ", env.stmt, " + (1|", block, ")", sep = "");
#				} else if(exptl.design == "AugLS")
#				{
#					row <- pheno$traits[[i]]$sites[[1]]$design$row;
#					column <- pheno$traits[[i]]$sites[[1]]$design$column;
#					myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ") + ", env.stmt, " + (1|", row , "), + (1|", column, ")", sep ="");
#				} else if(exptl.design == "Alpha")
#				{
#					rep <- pheno$traits[[i]]$sites[[1]]$design$rep;
#					block <- pheno$traits[[i]]$sites[[1]]$design$block;
#					myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ") + ", env.stmt, " + (1|", row , "), + (1|", column, ")", sep ="");					
#				} else if(exptl.design == "RowCol")
#				{
#					rep <- pheno$traits[[i]]$sites[[1]]$design$rep;
#					row <- pheno$traits[[i]]$sites[[1]]$design$row;
#					column <- pheno$traits[[i]]$sites[[1]]$design$column;
#					myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ") + ", env.stmt, " + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,")", sep = "");
#				} else if(exptl.design == "LatinAlpha")
#				{
#					rep <- pheno$traits[[i]]$sites[[1]]$design$rep;
#					block <- pheno$traits[[i]]$sites[[1]]$design$block;
#					myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ") + ", env.stmt, " + (1|", rep,") + (1|", block,") + (1|", rep, ":", block, ")", sep = "");
#					
#				} else if(exptl.design == "LatinRowCol")
#				{	
#					row <- pheno$traits[[i]]$sites[[1]]$design$row;
#					column <- pheno$traits[[i]]$sites[[1]]$design$column;
#					rep <- pheno$traits[[i]]$sites[[1]]$design$rep;
#					longerRow <- pheno$traits[[i]]$sites[[1]]$design$longerRow;
#					if(longerRow)
#					{
#						myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ") + ", env.stmt, " + (1|", rep,") + (1|", column, ") + (1|", rep,":", column,") + (1|", row, ")", sep = ""); 
#					} else
#					{
#						myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ") + ", env.stmt, " + (1|", rep,") + (1|", row, ") + (1|", rep,":", row,") + (1|", column, ")", sep = "") ;
#					}			
#				}
#				
#				model <- try(lmer(formula(myformula1), data = data), silent = TRUE);
#				if(!is.null(model) && class(model == "try-error"))
#				{
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$error <- TRUE;
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$message <- paste("There are error occured on this marker, i.e. ", marker.name, ".\n", sep = "");
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$model <- model;
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$formula <- myformula1;
#					next;
#				}
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$error <- FALSE;
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$message <- paste("There are nothing error on this marker, i.e. ", marker.name, sep = "");
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$model <- model;
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$formula <- myformula1;
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$adjust.means <- interactionMeans(model, factor= marker.name);
#				pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$marker.test <- testInteractions(model, pairwise = marker.name, test = test);
#				if(is.EnvFixed)
#				{
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$marker.across.env <- testInteractions(model, fixed = env, across = marker.name, test = test);
#				} else
#				{
#					model2 <- update(mode, formula(paste(". ~ . - ", env.stmt, sep = "")));
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$model2 <- model2
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[1]]$marker.across.env <- anova(model, model2);					
#				}
#				
#			} #--end of statment for(k in 1:marker.number) ---#
#		} else
#		{
#			if(include.env && (site.number == 1))
#				warning("\tWarning: There is only one available site, omitted the include.env paramter!\n");
#			for(j in 1:site.number)
#			{
#				if(pheno$traits[[i]]$sites[[j]]$restricted$isTRUE)
#					next;
#				phenotypicData <- pheno$traits[[i]]$sites[[j]]$data;
#				
#				#--- meger marker data with phenotypic data one by one ---#
#				data <- merge(phenotypicdata, genotypicData, by.x = geno.phenotypicData, by.y = geno.genotypicData);
#				for(k in 1:marker.number)
#				{
#					if(exptl.design %in% c("RCB", "AugRCB"))
#					{
#						block <- pheno$traits[[i]]$sites[[1]]$design$block;
#						myformula1 <- paste(trait.name,  " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ")", " + (1|", block, ")", sep = "");
#					} else if(exptl.design == "AugLS")
#					{
#						row <- pheno$traits[[i]]$sites[[1]]$design$row;
#						column <- pheno$traits[[i]]$sites[[1]]$design$column;
#						myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ")", " + (1|", row , "), + (1|", column, ")", sep ="");
#					} else if(exptl.design == "Alpha")
#					{
#						rep <- pheno$traits[[i]]$sites[[1]]$design$rep;
#						block <- pheno$traits[[i]]$sites[[1]]$design$block;
#						myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ")", " + (1|", row , "), + (1|", column, ")", sep ="");					
#					} else if(exptl.design == "RowCol")
#					{
#						rep <- pheno$traits[[i]]$sites[[1]]$design$rep;
#						row <- pheno$traits[[i]]$sites[[1]]$design$row;
#						column <- pheno$traits[[i]]$sites[[1]]$design$column;
#						myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ")", " + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,")", sep = "");
#					} else if(exptl.design == "LatinAlpha")
#					{
#						rep <- pheno$traits[[i]]$sites[[1]]$design$rep;
#						block <- pheno$traits[[i]]$sites[[1]]$design$block;
#						myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ")", " + (1|", rep,") + (1|", block,") + (1|", rep, ":", block, ")", sep = "");
#						
#					} else if(exptl.design == "LatinRowCol")
#					{	
#						row <- pheno$traits[[i]]$sites[[1]]$design$row;
#						column <- pheno$traits[[i]]$sites[[1]]$design$column;
#						rep <- pheno$traits[[i]]$sites[[1]]$design$rep;
#						longerRow <- pheno$traits[[i]]$sites[[1]]$design$longerRow;
#						if(longerRow)
#						{
#							myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ")", " + (1|", rep,") + (1|", column, ") + (1|", rep,":", column,") + (1|", row, ")", sep = ""); 
#						} else
#						{
#							myformula1 <- paste(respvar, " ~ 1 + ", marker.name, " + (1|", marker.name, ":", geno, ")", " + (1|", rep,") + (1|", row, ") + (1|", rep,":", row,") + (1|", column, ")", sep = "") ;
#						}			
#					}
#					model <- try(lmer(formula(myformula1), data = data), silent = TRUE);
#					if(!is.null(model) && class(model == "try-error"))
#					{
#						pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$error <- TRUE;
#						pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$message <- paste("There are error occured on this marker, i.e. ", marker.name, ".\n", sep = "");
#						pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$model <- model;
#						pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$formula <- myformula1;
#						next;
#					}
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$error <- FALSE;
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$message <- paste("There are nothing error on this marker, i.e. ", marker.name, sep = "");
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$model <- model;
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$formula <- myformula1;
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$adjust.means <- interactionMeans(model, factor= marker.name);
#					pheno$traits[[i]]$analysis$sma$markers[[k]]$sites[[j]]$marker.test <- testInteractions(model, pairwise = marker.name, test = test);
#				} #--- end of statement for(k in 1:marker.number) ---#
#			}#--- end of statement for(j in 1:site.number)---#
#			
#		} #--- end of statement if(include.env && (site.number > 1)) ---#
#	}#--- end statement of for(i in 1:trait.number)---#
#}
#
