###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Sep 29, 2014
# FileName: combindGenoAndPhenoOnIL.R
###############################################################################



combineGenoAndPhenoOnIL <- function(
		geno.data, 
		pheno.data,
		geno.factor
	)
	UseMethod(combineGenoAndPhenoOnIL);
combineGenoAndPhenoOnIL <- function(
		geno.data, 
		pheno.data,
		geno.factor
)
{
	if(is.null(geno.data))
	{
		stop("\tError: The genotypic data should not be null!\n");
	}
	if(is.null(pheno.data))
	{
		stop("\tError: The phenotypic data should not be null!\n");
	}
	if(!is.data.frame(geno.data))
	{
		stop("\tError: The genotypic data should be data frame!\n");
	}
	if(!is.data.frame(pheno.data))
	{
		stop("\tError: The phenotypic data should be data frame!\n");
	}
	
	
}