###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Nov 13, 2014
# FileName: restrict.pheno.data.R
###############################################################################


restrict.pheno.data <- function(
		data,
		missing.rate = 0.8
	) 
	UseMethod("restrict.pheno.data");

restrict.pheno.data.PhenotypicData <- function(
		data,
		missing.rate = 0.8
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
			data$traits[[i]]$sites[[j]]$restricted$manyNAWarning <- "Many observations are missing!";
		}
	}
	
	return(data);
}