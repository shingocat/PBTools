###############################################################################
# TODO: restricted the phenotypic data basing on the specified requirements.
#		And by far it will restrict missing rate, the coding label of pyramided 
#		line illegal or not (DEFAULT), use grand mean as the phenotypic value of
#		level that missed all observations
#		
#
# ARGUMENT: 
#	phenodata - The object of class "PhenotypicData"
#	missing.rate - missing rate of the whole phenotypic data;
#	check.label - checking wheather the genotypic labels are legal or not for 
#					pyramided line only, otherwise it will be discarded.
#	use.grand.mean - using the grand mean as the phenotypic value of the genotypic
#						level which missed all value;
# 
# Author: mqin
# Date: Nov 13, 2014
# FileName: restrict.pheno.data.R
###############################################################################


restrict.pheno.data <- function(
		phenodata,
		missing.rate = 0.2
		#use.grand.mean = TRUE
		#check.label = TRUE
	) 
	UseMethod("restrict.pheno.data");

restrict.pheno.data.PhenotypicData <- function(
		phenodata,
		missing.rate = 0.2
		#use.grand.mean = TRUE
		#check.label = TRUE
)
{
	if(!inherits(phenodata, "PhenotypicData"))
		stop("The argument of phenodata must be of class PhenotypicData");
	if(!is.numeric(missing.rate))
		stop("\tError: The argument of missing.rate should of value between 0 and 1!");
	if(missing.rate > 1 || missing.rate < 0)
		stop("\tError: The argument of missing.rate should of value between 0 and 1!");
	for(i in 1:length(phenodata$traits))
	{
		for(j in 1:length(phenodata$traits[[i]]$envs))
		{
			phenodata$traits[[i]]$envs[[j]]$restricted <- list();
			isRestricted <- FALSE;
			if(phenodata$traits[[i]]$envs[[j]]$missing.rate > missing.rate)
				isRestricted <- TRUE;
			phenodata$traits[[i]]$envs[[j]]$restricted$isTRUE <- isRestricted;
			phenodata$traits[[i]]$envs[[j]]$restricted$missing.rate.cond <- missing.rate;
			if(isRestricted)
			{
				phenodata$traits[[i]]$envs[[j]]$restricted$messages <- "Many observations are missing!";
			}else
			{
				phenodata$traits[[i]]$envs[[j]]$restricted$messages <- NA;
			}
			#--- checking whether missing rate by geno is meet---#
			missing.rate.by.geno <- phenodata$traits[[i]]$envs[[j]]$missing.rate.by.geno;
			old.phenodata <- phenodata$traits[[i]]$envs[[j]]$data;
			new.phenodata <- old.phenodata;
			trait.name <- phenodata$traits[[i]]$name;
			geno <- phenodata$traits[[i]]$envs[[j]]$design$geno;
			if(any(missing.rate.by.geno[,"RATE"] == 1))
			{
				for(k in which(missing.rate.by.geno[, "RATE"] == 1))
				{
					geno.level.name <- missing.rate.by.geno[k,1];
					warning(paste("\tWARNING: The genotypic level ", geno.level.name, 
									" would be deleted on trait ", trait.name, "!\n",
									"\tBecause there are no observed values on this levels!", sep =""));
					new.phenodata <- new.phenodata[new.phenodata[[geno]] != geno.level.name ,];
				}
				new.phenodata[ , geno] <- as.factor(as.character(new.phenodata[ , geno]));
			}
			phenodata$traits[[i]]$envs[[j]]$data <- new.phenodata;
			phenodata$traits[[i]]$envs[[j]]$old.data <- old.phenodata;
			#--- recomputed the observation read, observation used, missing.rate and missing.rate.by.geno after using grand mean insteated of missing value ---#
			obsread <- nrow(new.phenodata);
			obsused <- nrow(subset(new.phenodata, subset=(is.na(new.phenodata[,trait.name]) == FALSE)));
			missing.rate.used <- (obsread - obsused) / obsread;
			phenodata$traits[[i]]$envs[[j]]$obsread <- obsread;
			phenodata$traits[[i]]$envs[[j]]$obsused <- obsused;
			phenodata$traits[[i]]$envs[[j]]$missing.rate <- missing.rate.used;
			missing.rate.by.geno <- aggregate(new.phenodata[[trait.name]], by = list(new.phenodata[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
			names(missing.rate.by.geno) <- c(geno, "RATE");
			phenodata$traits[[i]]$envs[[j]]$missing.rate.by.geno <- missing.rate.by.geno;
			#--- checking whether it is full level or not when the population is pyramided line---#
			if(identical(phenodata$pop.type, "PL"))
			{	
				new.levels <- levels(as.factor(as.character(new.phenodata[[geno]])));
				gene.num <- phenodata$gene.num;
				isFullLevel <- TRUE;
				if(gene.num == 2)
				{
					if(!all(levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T)) %in% new.levels))
						isFullLevel <- FALSE;
				} else if(gene.num == 3)
				{
					if(!all(levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C=c("CC","Cc","cc")), sep="", lex.order = T)) %in% new.levels))
						isFullLevel <- FALSE;
				} else if(gene.num == 4)
				{
					if(!all(levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C=c("CC","Cc","cc"), D = c("DD","Dd","dd")), sep="", lex.order = T)) %in% new.levels))
						isFullLevel <- FALSE;		
				}
				phenodata$traits[[i]]$envs[[j]]$FullLevel <- isFullLevel;
			}
#			if(use.grand.mean)
#			{
#				missing.rate.by.geno <- phenodata$traits[[i]]$envs[[j]]$missing.rate.by.geno;
#				if(any(missing.rate.by.geno[,"RATE"] == 1))
#				{
#					old.phenodata <- phenodata$traits[[i]]$envs[[j]]$data;
#					new.phenodata <- old.phenodata;
#					trait.name <- phenodata$traits[[i]]$name;
#					grand.mean <- mean(old.phenodata[,trait.name], na.rm=TRUE);
#					geno <- phenodata$traits[[i]]$envs[[j]]$design$geno;
#					for(k in which(missing.rate.by.geno[,"RATE"] == 1))
#					{
#						new.phenodata[new.phenodata[,geno] == missing.rate.by.geno[k,1], trait.name] <- grand.mean;
#					}
#					phenodata$traits[[i]]$envs[[j]]$data <- new.phenodata;
#					phenodata$traits[[i]]$envs[[j]]$old.data <- old.phenodata;
#					#--- recompute the observation readed, observation used, missing.rate and missing.rate.by.geno after using grand mean insteated of missing value ---#
#					obsread <- nrow(new.phenodata);
#					obsused <- nrow(subset(new.phenodata, subset=(is.na(new.phenodata[,trait.name]) == FALSE)));
#					missing.rate <- (obsread - obsused) / obsread;
#					phenodata$traits[[i]]$envs[[j]]$obsread <- obsread;
#					phenodata$traits[[i]]$envs[[j]]$obsused <- obsused;
#					phenodata$traits[[i]]$envs[[j]]$missing.rate <- missing.rate;
#					missing.rate.by.geno <- aggregate(new.phenodata[[trait.name]], by = list(new.phenodata[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
#					names(missing.rate.by.geno) <- c("GENO", "RATE");
#					phenodata$traits[[i]]$envs[[j]]$missing.rate.by.geno <- missing.rate.by.geno;
#				}
#			}
		}
	}
	phenodata$isRestricted <- TRUE;
	return(phenodata);
}