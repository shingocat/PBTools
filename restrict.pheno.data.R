###############################################################################
# TODO: restricted the phenotypic data basing on the specified requirements.
#		And by far it will restrict missing rate, the coding label of pyramided 
#		line illegal or not (DEFAULT), use grand mean as the phenotypic value of
#		level that missed all observations
#		
#
# ARGUMENT: 
#	data - The object of class "PhenotypicData"
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
		data,
		missing.rate = 0.2
		#use.grand.mean = TRUE
		#check.label = TRUE
	) 
	UseMethod("restrict.pheno.data");

restrict.pheno.data.PhenotypicData <- function(
		data,
		missing.rate = 0.2
		#use.grand.mean = TRUE
		#check.label = TRUE
)
{
	if(!inherits(data, "PhenotypicData"))
		stop("The argument of data must be of class PhenotypicData");
	if(!is.numeric(missing.rate))
		stop("\tError: The argument of missing.rate should of value between 0 and 1!");
	if(missing.rate > 1 || missing.rate < 0)
		stop("\tError: The argument of missing.rate should of value between 0 and 1!");
	for(i in 1:length(data$traits))
	{
		for(j in 1:length(data$traits[[i]]$sites))
		{
			data$traits[[i]]$sites[[j]]$restricted <- list();
			isRestricted <- FALSE;
			if(data$traits[[i]]$sites[[j]]$missing.rate > missing.rate)
				isRestricted <- TRUE;
			data$traits[[i]]$sites[[j]]$restricted$isTRUE <- isRestricted;
			data$traits[[i]]$sites[[j]]$restricted$missing.rate <- missing.rate;
			if(isRestricted)
			{
				data$traits[[i]]$sites[[j]]$restricted$manyNAWarning <- "Many observations are missing!";
			}else
			{
				data$traits[[i]]$sites[[j]]$restricted$manyNAWarning <- NA;
			}
			missing.rate.by.geno <- data$traits[[i]]$sites[[j]]$missing.rate.by.geno;
			old.data <- data$traits[[i]]$sites[[j]]$data;
			new.data <- old.data;
			trait.name <- data$traits[[i]]$name;
			geno <- data$traits[[i]]$sites[[j]]$design$geno;
			if(any(missing.rate.by.geno[,"RATE"] == 1))
			{
				for(k in which(missing.rate.by.geno[, "RATE"] ==1))
				{
					geno.level.name <- missing.rate.by.geno[k,1];
					warning(paste("\tWARNING: The genotypic level ", geno.level.name, 
									" would be deleted on trait ", trait.name, "!\n",
									"\tBecause there are no observed values on this levels!", sep =""));
					new.data <- new.data[new.data[[geno]] != geno.level.name ,];
				}
				new.data[ , geno] <- as.factor(as.character(new.data[ , geno]));
			}
			data$traits[[i]]$sites[[j]]$data <- new.data;
			data$traits[[i]]$sites[[j]]$old.data <- old.data;
			#--- recomputed the observation read, observation used, missing.rate and missing.rate.by.geno after using grand mean insteated of missing value ---#
			obsread <- nrow(new.data);
			obsused <- nrow(subset(new.data, subset=(is.na(new.data[,trait.name]) == FALSE)));
			missing.rate <- (obsread - obsused) / obsread;
			data$traits[[i]]$sites[[j]]$obsread <- obsread;
			data$traits[[i]]$sites[[j]]$obsused <- obsused;
			data$traits[[i]]$sites[[j]]$missing.rate <- missing.rate;
			missing.rate.by.geno <- aggregate(new.data[[trait.name]], by = list(new.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
			names(missing.rate.by.geno) <- c(geno, "RATE");
			data$traits[[i]]$sites[[j]]$missing.rate.by.geno <- missing.rate.by.geno;
			#--- checking whether it is full level or not when the population is pyramided line---#
			if(identical(data$pop.type, "PL"))
			{	
				new.levels <- levels(as.factor(as.character(new.data[[geno]])));
				gene.num <- data$gene.num;
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
				data$traits[[i]]$sites[[j]]$FullLevel <- isFullLevel;
			}
#			if(use.grand.mean)
#			{
#				missing.rate.by.geno <- data$traits[[i]]$sites[[j]]$missing.rate.by.geno;
#				if(any(missing.rate.by.geno[,"RATE"] == 1))
#				{
#					old.data <- data$traits[[i]]$sites[[j]]$data;
#					new.data <- old.data;
#					trait.name <- data$traits[[i]]$name;
#					grand.mean <- mean(old.data[,trait.name], na.rm=TRUE);
#					geno <- data$traits[[i]]$sites[[j]]$design$geno;
#					for(k in which(missing.rate.by.geno[,"RATE"] == 1))
#					{
#						new.data[new.data[,geno] == missing.rate.by.geno[k,1], trait.name] <- grand.mean;
#					}
#					data$traits[[i]]$sites[[j]]$data <- new.data;
#					data$traits[[i]]$sites[[j]]$old.data <- old.data;
#					#--- recompute the observation readed, observation used, missing.rate and missing.rate.by.geno after using grand mean insteated of missing value ---#
#					obsread <- nrow(new.data);
#					obsused <- nrow(subset(new.data, subset=(is.na(new.data[,trait.name]) == FALSE)));
#					missing.rate <- (obsread - obsused) / obsread;
#					data$traits[[i]]$sites[[j]]$obsread <- obsread;
#					data$traits[[i]]$sites[[j]]$obsused <- obsused;
#					data$traits[[i]]$sites[[j]]$missing.rate <- missing.rate;
#					missing.rate.by.geno <- aggregate(new.data[[trait.name]], by = list(new.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
#					names(missing.rate.by.geno) <- c("GENO", "RATE");
#					data$traits[[i]]$sites[[j]]$missing.rate.by.geno <- missing.rate.by.geno;
#				}
#			}
		}
	}
	data$doRestricted <- TRUE;
	return(data);
}