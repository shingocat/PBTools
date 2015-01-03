###############################################################################
# TODO: compute multi marker analysis on PhenotypicData, SingleEnvAnalysis and 
#		MultiEnvAnalysis object. And if it is PhenotypicData object, it will 
#		only support MEAN type data of PhenotypicData object right now. The RAW
#		type of PhenotypicData should be conducte Single or Multi environment
#		analysis first.
# [Arguments]
#	phenotypicData - The phenotypicData is one of PhenothypicData, SingleEnvAnalysis and MultiEnvAnalysis object
#	genotypicData - The genotypicData is the GenotypicData object;
#	resp.var - respone variable of PhenotypicData you want to conduct multiple marker analysis
#	method - the regression method, there are three types of method  "LASSO", "RIDGE_REGRESSION", "ELASTIC_NET"
#	siglevel - default 0.05, the significant level
#	bootstrap - default 10;
#	pval.method - the method for computing p value, four method now "median", "fdr", "holm", "QA"
#	family - the family of response variable, three types family "gaussian", "binomial", "possion"
#	include.ht - default TRUE, include heterozygous of GenotypicData;
#	nfolds - default 3, separated data into how-much fold
#	step - default 0.1, the step used to check ELASTIC_NET which lambda is good;
#	max.try - default 10, the maximum of try to do high-dimensional regression;
# 
# Author: mqin
# Date: Nov 28, 2014
# FileName: doMM.R
###############################################################################


doMMA <- function(
		phenotypicData,
		genotypicData,
		resp.var,
		method = c("LASSO", "RIDGE_REGRESSION", "ELASTIC_NET"),
		siglevel = 0.05,
		bootstrap = 10,
		pval.method = c("median","fdr", "holm", "QA"),
		family = c("gaussian", "binomial", "possion"),
		include.ht = TRUE,
		nfolds = 3,
		step = 0.1,
		max.try = 10,
		...
)
{
	if(missing(phenotypicData))
		stop("\tThe phenotypicData argument could not be null!\n");
	if(missing(genotypicData))
		stop("\tThe genotypicData argument could not be null!\n");
	if(missing(siglevel))
		siglevel <- 0.05;
	if(missing(bootstrap))
		bootstrp <- 10;
	if(missing(include.ht))
		include.ht <- TRUE;
	if(missing(nfolds))
		nfolds <- 3;
	if(missing(step))
		step <-  0.1;
	if(missing(max.try))
		max.try <- 10;
	
	UseMethod("doMMA");
}

#--- the S3 method for PhenotypicData Object---#
doMMA.PhenotypicData <- function(
		phenotypicData,
		genotypicData,
		resp.var,
		method = c("LASSO", "RIDGE_REGRESSION", "ELASTIC_NET"),
		siglevel = 0.05,
		bootstrap = 10,
		pval.method = c("median","fdr", "holm", "QA"),
		family = c("gaussian", "binomial", "possion"),
		include.ht = TRUE,
		nfolds = 3,
		step = 0.1,
		max.try = 10,
		...
)
{
	method <- match.arg(method);
	pval.method <- match.arg(pval.method);
	family <- match.arg(family);
	#--- only supported PhenotypicData is MEAN type, if not, it should do the Single or Multi Environment Analysis first---#
	if(phenotypicData$isMean)
	{
		if(!phenotypicData$isRestricted)
		{
			warning("\tWarning: It will be used default parameters to conduct restricting phenotypic data!\n");
			phenotypicData <- restrict.pheno.data(phenotypicData);
		}
		if(!genotypicData$isRestricted)
		{
			warning("\tWarning: It will be used default parameters to conduct restricting genotypic data!\n");
			genotypicData <- restrict.geno.data(genotypicData);
		}
		
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
					warings(paste("\tWarnings: It will be omitted the ", resp.var[i], " because of it is not exist!\n", sep = ""));
				}
			}
			if(length(temp.resp.var) == 0)
				stop("\tError: All the specified resp.var are not exists in the phenotypicdata!\n");
			resp.var <- temp.resp.var;
		} else
		{
			resp.var <- phenotypicData$trait.names;
		}
		
		genotypicData.restricted <- genotypicData$restricted$data;
		marker.number <- ncol(genotypicData.restricted) - 1;
		marker.names <- colnames(genotypicData.restricted)[-1];
		
		#--- suppress warning message---#
		old.options <- options();
		options(warn = 0);
		library(hdlm);
		
		#--- reformat all genotypic data using dp.code instead of ht.code basing on include.ht ----#
		#--- Because it will recode to default coding sytem, it could be define such code directly!---#
		dp.code <- 2;
		rp.code <- 0;
		ht.code <- 1;
		na.code <- NA;
		if(include.ht)
		{
			genotypicData.restricted[ , -1][genotypicData.restricted[,-1] == ht.code] <- dp.code;
		}
		genotypicData.restricted <- apply(genotypicData.restricted, 2, factor);
		
		for(i in 1:length(resp.var))
		{
			#--- checking whether the resp.var name is the trait name of phenotypicData---#
			for(j in 1:length(phenotypicData$traits))
			{
				if(!identical(resp.var[i], phenotypicData$traits[[j]]$name))
					next;
				#--- checking the analysis structure of phenotypicData ---#
				if(is.null(phenotypicData$traits[[j]]$analysis))
					phenotypicData$traits[[j]]$analysis <- list();
				if(is.null(phenotypicData$traits[[j]]$analysis$mma))
				{
					phenotypicData$traits[[j]]$analysis$mma <- list();
				} else
				{
					phenotypicData$traits[[j]]$analysis$mma <- NULL;
					phenotypicData$traits[[j]]$analysis$mma <- list();
				}
				
				trait.name <- resp.var[i];
				phenotypicData$traits[[j]]$analysis$mma$method <- ifelse(is.null(method) | is.na(method), "NA", method);
				phenotypicData$traits[[j]]$analysis$mma$siglevel <- ifelse(is.null(siglevel) | is.na(siglevel), "NA", siglevel);
				phenotypicData$traits[[j]]$analysis$mma$bootstrap <- ifelse(is.null(bootstrap) | is.na(bootstrap), "NA", bootstrap);
				phenotypicData$traits[[j]]$analysis$mma$pval.method <- ifelse(is.null(pval.method) | is.na(pval.method), "NA", pval.method);
				phenotypicData$traits[[j]]$analysis$mma$family <- ifelse(is.null(family) | is.na(family), "NA", family);
				phenotypicData$traits[[j]]$analysis$mma$include.ht <- ifelse(is.null(include.ht) | is.na(include.ht), "NA", family);
				phenotypicData$traits[[j]]$analysis$mma$nfolds <- ifelse(is.null(nfolds) | is.na(nfolds), "NA", nfolds);
				phenotypicData$traits[[j]]$analysis$mma$step <- ifelse(is.null(step) | is.na(step), "NA", step);
				phenotypicData$traits[[j]]$analysis$mma$max.try <- ifelse(is.null(max.try) | is.na(max.try), "NA", step);
				phenotypicData$traits[[j]]$analysis$mma$envs <- list();
				
				#--- separated env on multi marker analysis ---#
				for(k in 1:length(phenotypicData$traits[[j]]$envs))
				{
					env.name <- phenotypicData$traits[[j]]$envs[[k]]$name;
					#--- checking this env data of this trait whether is restricted, if it is true, try next env--#
					if(phenotypicData$traits[[j]]$envs[[k]]$restricted$isTRUE)
					{
						warning("\tWarning: The phenotypic data on ", env.name, " of ", trait.name, " is restricted! It will conduct next environment!\n");
						next;
					}
					phenotypicData$traits[[j]]$analysis$mma$envs[[k]] <- list();
					phenotypicData$traits[[j]]$analysis$mma$envs[[k]]$name <- env.name;
					#--- reformating genotypic data to numeric type---#
#					genotypicData.restricted <- apply(genotypicData.restricted,2,trimStrings);
#					genotypicData.restricted[,-1] <- apply(genotypicData.restricted[,-1], 2, as.character);
#					genotypicData.restricted[,-1] <- apply(genotypicData.restricted[,-1], 2, as.numeric);
					#--- merge phenotypicData and genotoypicData into one data frame ---#
					geno.name.phenotypicData <- phenotypicData$traits[[j]]$envs[[k]]$design$geno;
					geno.name.genotypicData <- genotypicData$geno.name;
					data <- merge(phenotypicData$traits[[j]]$envs[[k]]$data, genotypicData.restricted, by.x = geno.name.phenotypicData, by.y = geno.name.genotypicData);
					data[,trait.name] <- as.numeric(as.character(data[,trait.name]));
					#--- remove missing value on trait and corresponding records on genotypic data ---#
					na.index <- which(is.na(data[,trait.name]));
					no.ind <- nrow(data);
					if(length(na.index) > 0)
						data <- data[-na.indx,];
					trait.value <- data[,trait.name];
					trait.value <- as.numeric(as.character(data[,trait.name]));
					genetic.value <- data[,marker.names];
					genetic.value <- apply(genetic.value, 2, trimStrings);
					genetic.value <- apply(genetic.value, 2, as.character);
					genetic.value <- apply(genetic.value, 2, as.numeric);
					if(method == "LASSO")
					{
						#---Try ten times whether hdlm could fit the model, if not, show error message!---#
						n <- 0;
						while(n <= max.try)
						{
							out <- try(hdglm(trait.value ~ genetic.value, M = no.ind, bootstrap = bootstrap, pval.method = pval.method, alpha = 1, siglevel = siglevel, family = family), silent = TRUE);
							if(!identical(class(out),"try-error"))
								break;
							n <- n + 1;
						}
					} else if(method == "RIDGE_REGRESSION")
					{ 
						#---Try ten times whether hdlm could fit the model, if not, show error message!---#
						n <- 0;
						while(n <= max.try)
						{
							out <- try(hdglm(trait.value ~ genetic.value, M = no.ind, bootstrap = bootstrap, pval.method = pval.method, alpha = 0.000001, siglevel = siglevel, family = family), silent = TRUE);
							if(!identical(class(out), "try-error"))
								break;
							n <- n + 1;
						}
					} else if(method == "ELASTIC_NET")
					{
						elastic.net.value <- getElasticNetAlpha(genetic.value, trait.value, nfolds = nfolds, step = step);
						
						#---Try ten times whether hdlm could fit the model, if not, show error message!---#
						n <- 0;
						while(n <= max.try)
						{
							out <- try(hdglm(trait.value ~ genetic.value, M = no.ind, bootstrap = bootstrap, pval.method = pval.method, alpha = elastic.net.value, siglevel = siglevel, family = family), silent = TRUE);
							if(!identical(class(out), "try-error"))
								break;
							n <- n + 1;
						}
						
					}
					
					phenotypicData$traits[[j]]$analysis$mma$envs[[k]]$outcomes <- out;
				} #--- end stmt of for(k in 1:length(phenotypicData$traits[[j]]$envs)) ---#
			} #--- end stmt of for(i in 1:length(phenotpyicData$traits)) ---#
		} #--- end stmt of for(i in 1:length(resp.var)) ---#
		
		phenotypicData$genotypicData <- genotypicData;
		options(old.options);
		detach("package:hdlm", unload = TRUE);
		class(phenotypicData) <- c("MultiMarkerAnalysis", class(phenotypicData));
		return(phenotypicData);
		
	} else
	{
		stop("\tError: It could not conduct multi-marker analysis on RAW type of PhenotypicData Object. Please do Single or Multi Environment Analysis at first!\n");
	} #--- end of if(phenotypicData$isMean)
}

#--- the S3 method for SingleEnvAnalysis object---#
doMMA.SingleEnvAnalysis <- function(
		phenotypicData,
		genotypicData,
		resp.var,
		method = c("LASSO", "RIDGE_REGRESSION", "ELASTIC_NET"),
		siglevel = 0.05,
		bootstrap = 10,
		pval.method = c("median","fdr", "holm", "QA"),
		family = c("gaussian", "binomial", "possion"),
		include.ht = TRUE,
		nfolds = 3,
		step = 0.1,
		max.try = 10,
		...
)
{
	
}

#--- the s3 method for MultiEnvAnalysis object ---#

doMMA.MultiEnvAnalysis <- function(
		phenotypicData,
		genotypicData,
		resp.var,
		method = c("LASSO", "RIDGE_REGRESSION", "ELASTIC_NET"),
		siglevel = 0.05,
		bootstrap = 10,
		pval.method = c("median","fdr", "holm", "QA"),
		family = c("gaussian", "binomial", "possion"),
		include.ht = TRUE,
		nfolds = 3,
		step = 0.1,
		max.try = 10,
		...
)
{
	
}

getElasticNetAlpha <- function
		(
		x, 
		y, 
		nfolds = 3,
		step = 0.1,
		...
)
{
	#---suppress warning messages ---#
	old.options <- options();
	options(warn = -1);
	if(missing(x))
		stop("\tError: The argument of x could not be null!\n");
	if(missing(y))
		stop("\tError: The argument of y could not be null!\n");
	if(missing(nfolds))
		nfolds <- 3;
	if(!is.numeric(nfolds))
		stop("\tError: The argument of nfolds should be of type integer!\n");
	if(missing(step))
		step <- 0.1; # set step by default value 0.1.
	if(!is.numeric(step))
		stop("\tError: The argument of step should be of numeric value!\n");
	if(step < 0 || step > 1)
		stop("\tError: The argument of step should be not bigger than 1 and samller than 0!\n");
	
	set.seed(1);
	ind.number <- length(y);
	foldids <- sample(1:nfolds, size = ind.number, replace = TRUE);
	alphas <- seq(from = 0, to = 1, by = step);
	#library(glmnet);
	cv.lambdas <- c();
	for(i in 1:length(alphas))
	{
		cv.lambdas <- c(cv.lambdas, cv.glmnet(x, y, foldid = foldids, alpha = alphas[i])$lambda.min);
	}
	alpha <- alphas[which(cv.lambdas == min(cv.lambdas))][1]; # get the first one of smallest lambda
	if(is.null(alpha) || length(alpha) == 0)
		stop("\tError: It could not be found promised alpha for Elastic Net!\n");
	#detach("package:glmnet", unload = TRUE);
	options(old.options);
	return(alpha);
}

print.MultiMarkerAnalysis <- function(
		data, 
		p.value = 0.05, 
		level = 1)
{
	if(!inherits(data, "MultiMarkerAnalysis"))
		stop("\tError: The argument of data should be of class MultiMarkerAnalysis!\n");
	if(!is.numeric(p.value))
		stop("\tError: The argument of p.value should be of value between 0 and 1.\n");
	if(p.value > 1 || p.value < 0)
		stop("\tError: The argument of p.value should be of value between 0 and 1.\n");
	if(!level %in% c(1,2,3))
		stop("\tError: The argument of level should be one of value 1, 2 and 3.\t");
	cat("Multiple Marker Analysis (p.value < ", p.value, "):\n", sep = "");
	
	library(hdlm);
	for(i in 1:length(data$traits))
	{
		if(is.null(data$traits[[i]]$analysis$mma))
			next;
		cat("The phenotype is ", data$traits[[i]]$name, sep = "");
		cat("\n")
		cat(rep("-", times = 50));
		cat("\n");
		for(j in 1:length(data$traits[[i]]$analysis$mma$envs))
		{
			cat("The environment is ", data$traits[[i]]$analysis$mma$envs[[j]]$name, sep = "");
			cat("\n");
			outcomes <- data$traits[[i]]$analysis$mma$envs[[j]]$outcomes;
			if(identical(class(outcomes), "try-error"))
			{
				cat("It could not fit the model with multiple markers on this trait.\n");
			}else
			{
				out.summary <- summary(outcomes, level = level);
				coef <- out.summary$coefficients;
				lower.bound <- out.summary$lower.bound;
				upper.bound <- out.summary$upper.bound;
				PValue <- out.summary$p.value;
				temp <- data.frame(Marker = names(coef), Estimate = coef, CI.Lower = lower.bound, CI.Upper = upper.bound, P.value = PValue);
				#--- remove intercept in the temp data frame ---#
				temp <- temp[-1,];
				print(temp[which(temp$P.value <= p.value),], row.names = FALSE);
			}
			cat("\n");
		}
	}
	detach("package:hdlm", unload = TRUE);
}