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
		is.EnvFixed = TRUE,
		include.ht = FALSE,
		test = c("F", "Chisq"),
		digits = 4,
		...
)
{
	if(missing(PhenotypicData))
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
		warings("\tWarnings: It could compute all respone variables of phenotypicData");
	if(missing(test))
		test = match.arg(test);
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
		is.EnvFixed = TRUE,
		include.ht = FALSE,
		test = c("F", "Chisq"),
		digits = 4,
		...
)
{
	#--- checking resp.var whether exists in the phenotypi data ---#
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
				warings(paste("\tWarnings: It will be omitted the ", resp.var[i], " because of it is not exist!\n", sep = ""));
			}
		}
		if(length(temp.resp.var) == 0)
			stop("\tError: All the specified resp.var are not exists in the phenotypicdata!\n");
	}
	if(phenotypicData$isMean)
	{
		if(!phenotypicData$isRestricted)
		{
			warnings("\tWarning: It will be used default argument to restrict phenodtypic data!\n");
			phenotypicData <- restict.pheno.data(phenotypicData);
		}
		if(!genotypicData$isRestricted)
		{
			warnings("\tWarning: It will be used default argument to restrict genotypic data!\n");
			genotypicData <- restrict.geno.data(genotypicData);
		}
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
		is.EnvFixed = TRUE,
		include.ht = FALSE,
		test = c("F", "Chisq"),
		digits = 4,
		...
)
{
	
}

#--- It is used to process single marker analysis on MultiEnvAnalysis---#
doSMA.MultiEnvAnalysis <- function(
		phenotypicData,
		genotypicData,
#		geno,
#		env,
		resp.var,
		is.EnvFixed = TRUE,
		include.ht = FALSE,
		test = c("F", "Chisq"),
		digits = 4,
		...
)
{
	
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
