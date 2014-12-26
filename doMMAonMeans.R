###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Dec 10, 2014
# FileName: doMMSAonMeans.R
###############################################################################


doMMAonMeans <- function(
		phenotypicData,
		genotypicData,
		geno,
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
	if(missing(siglevel))
		siglevel <- 0.05;
	if(missing(bootstrap))
		bootstrp <- 10;
	
	#---not checking any conditions by far ---#
	genotypicData <- genotypicData$restricted$data;
	marker.number <- ncol(genotypicData) - 1;
	marker.names <- colnames(genotypicData)[-1];
	no.ind <- nrow(phenotypicData);
	#---The maximum number of trying hdlm function when the model is overfit! Default is 10.---#
	max.try <- max.try;
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
		genotypicData[ , -1][genotypicData[,-1] == ht.code] <- dp.code;
	}
	genotypicData <- apply(genotypicData, 2, factor);
	#---not merge genotypic and phenotypic data---#
	#---only checking the genotypic and phenotypic data have same lines informations---#
	#data <- merge(phenotypicData, genotypicData, by = geno);
	if(length(genotypicData[,geno]) != length(phenotypicData[ ,geno]))
		stop("\tError: The phenotypic data and genotypic data have different lines number!\n");
	if(!all(as.character(genotypicData[,geno]) %in% as.character(phenotypicData[,geno])) ||
			!all(as.character(phenotypicData[,geno]) %in% as.character(genotypicData[,geno])))
		stop("\tError: The phenotypic data and genotypic data have different lines names!\n");	
	result <- list();
	result$trait.number <- length(resp.var);
	result$trait.name <- resp.var;
	result$marker.number <- marker.number;
	result$marker.names <- marker.names;
	result$traits <- list();
	
	for(i in 1:length(resp.var))
	{
		genetic.value <- genotypicData[ , -1];
		genetic.value <- apply(genetic.value, 2, trimStrings);
		genetic.value <- apply(genetic.value, 2, as.character);
		genetic.value <- apply(genetic.value, 2, as.numeric);
		
		trait.name <- resp.var[i];
		result$traits[[i]] <- list();
		result$traits[[i]]$name <- resp.var[i];
		trait.value <- phenotypicData[,trait.name];
		#--- remove missing value on trait.value and corresponding records on genetic value---#
		na.index <- which(is.na(trait.value));
		if(length(na.index) > 0)
		{
			trait.value <- trait.value[-na.index];
			genetic.value <- genetic.value[-na.index,];
		}
		if(method == "LASSO")
		{
			#---Try ten times whether hdlm could fit the model, if not, show error message!---#
			j <- 0;
			while(j <= max.try)
			{
				out <- try(hdglm(trait.value ~ genetic.value, M = no.ind, bootstrap = bootstrap, pval.method = pval.method, alpha = 1, siglevel = siglevel, family = family), silent = TRUE);
				if(!identical(class(out),"try-error"))
					break;
				j <- j + 1;
			}
		} else if(method == "RIDGE_REGRESSION")
		{ 
			#---Try ten times whether hdlm could fit the model, if not, show error message!---#
			j <- 0;
			while(j <= max.try)
			{
				out <- try(hdglm(trait.value ~ genetic.value, M = no.ind, bootstrap = bootstrap, pval.method = pval.method, alpha = 0.000001, siglevel = siglevel, family = family), silent = TRUE);
				if(!identical(class(out), "try-error"))
					break;
				j <- j + 1;
			}
		} else if(method == "ELASTIC_NET")
		{
			elastic.net.value <- getElasticNetAlpha(genetic.value, trait.value, nfolds = nfolds, step = step);
			
			#---Try ten times whether hdlm could fit the model, if not, show error message!---#
			j <- 0;
			while(j <= max.try)
			{
				out <- try(hdglm(trait.value ~ genetic.value, M = no.ind, bootstrap = bootstrap, pval.method = pval.method, alpha = elastic.net.value, siglevel = siglevel, family = family), silent = TRUE);
				if(!identical(class(out), "try-error"))
					break;
				j <- j + 1;
			}
			
		}
		result$traits[[i]]$outcomes <- out;
	} #--- end of statement of for(i in 1:length(resp.var)) ---#
	options(old.options);
	detach("package:hdlm", unload = TRUE);
	class(result) <- "MultiMarkerAnalysis";
	return(result);
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

print.MultiMarkerAnalysis <- function(data, p.value = 0.05, level = 1)
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
		cat("The phenotype is ", data$traits[[i]]$name, sep = "");
		cat("\n")
		cat(rep("-", times = 50));
		cat("\n");
		outcomes <- data$traits[[i]]$outcomes;
		if(identical(class(outcomes), "try-error"))
		{
			cat("It could not fit the model with multiple markers on this trait.\n");
		}
		else
		{
			out.summary <- summary(data$traits[[i]]$outcomes, level = level);
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
	detach("package:hdlm", unload = TRUE);
}