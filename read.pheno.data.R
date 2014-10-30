###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Oct 15, 2014
# FileName: read.pheno.data.R
###############################################################################


read.pheno.data <- function(
		file,
		pop.type,
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),	
		resp.var,
		geno,
		block,
		row = NULL,
		column = NULL,
		rep = NULL,
		env = NULL,
		na.code = NA
)
{
	UseMethod("read.pheno.data");
}
read.pheno.data.default <- function(
		file,
		pop.type,
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),	
		resp.var,
		geno,
		block,
		row = NULL,
		column = NULL,
		rep = NULL,
		env = NULL,
		na.code = NA
)
{
	if(is.null(file))
	{
		stop("\tError: The argument of file could not be null!\n");
	}
	if(!is.character(file))
	{
		stop("\tError: The argument of file should be character!\n");
	}
	if(!file.exists(file))
	{
		stop("\tError: The file does not exist! Please checking again!\n");
	}
	if(is.null(pop.type))
	{
		stop("\tError: The argument of pop.type could not be null!\n");
	}
	if(!is.character(pop.type))
		stop("\tError: The argument of pop.type should be of character type!\n");
	if(length(pop.type) != 1)
	{
		stop("\tError: The argument of pop.type should be of length on 1!\n");
	}
	if(!(pop.type %in% c("SSSL", "PL", "IL")))
	{
		stop("\tError: The argument of pop.type should be of value \"SSSL\",\"PL\",\"IL\".\n");
	}
	if(is.null(exptl.design))
	{
		stop("\tError: The argument of exptl.design could not be null!\n");
	}
	if(!is.character(exptl.design))
		stop("\tError: The argument of exptl.design should be of character type!\n");
	if(length(exptl.design) != 1)
	{
		stop("\tError: The argument of exptl.design should be of length on 1!\n");
	}
	if(!(exptl.design %in% c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol")))
	{
		stop("\tError: The argument of exptl.design should be of value \"RCB\", \"AugRCB\", \"AugLS\", \"Alpha\", \"RowCol\", \"LatinAlpha\", \"LatinRowCol\".\n");
	}
	if(is.null(respvar))
	{
		stop("\tError: The argument of respvar could not be null!\n");
	}
	if(is.null(geno))
	{
		stop("\tError: The argument of geno could not be null!\n");
	}
	if(length(geno) != 1)
	{
		stop("\tError: The argument of geno should be of length of 1!\n");
	}
	if(exptl.design == "RCB" || exptl.design == "AugRCB")
	{
		if(is.null(block))
			stop("\tError: The argument of block could not be null when the exptl.design is RCB or AugRCB!\n");
		if(length(block) != 1)
			stop("\tError: The argument of block should be of length of 1!\n");
	} else if(exptl.design == "AugLS")
	{
		if(is.null(row))
			stop("\tError: The argument of row could not be null when the exptl.design is AugLS!\n");
		if(length(row) != 1)
			stop("\tError: The argument of row should be of length of 1!\n");
		if(is.null(column))
			stop("\tError: The argument of column could not be null when the exptl.design is AugLS!\n");
		if(length(column) != 1)
			stop("\tError: The argument of column should be of length of 1!\n");
	} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
	{
		if(is.null(block))
			stop("\tError: The argument of block could not be null when the exptl.design is Alpha or LattinAlpha!\n");
		if(length(block) != 1)
			stop("\tError: The argument of block should be of length of 1!\n");
		if(is.null(rep))
			stop("\tError: The argument of rep could not be null when the exptl.design is Alpha or LattinAlpha!\n");
		if(length(rep) != 1)
			stop("\tError: The argument of rep should be of length of 1!\n");
	} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
	{
		if(is.null(row))
			stop("\tError: The argument of row could not be null when the exptl.design is RowCol or LatinRowCol!\n");
		if(length(row) != 1)
			stop("\tError: The argument of row should be of length of 1!\n");
		if(is.null(column))
			stop("\tError: The argument of column could not be null when the exptl.design is RowCol or LatinRowCol!\n");
		if(length(column) != 1)
			stop("\tError: The argument of column should be of length of 1!\n");
		if(is.null(rep))
			stop("\tError: The argument of rep could not be null when the exptl.design is RowCol or LatinRowCol!\n");
		if(length(rep) != 1)
			stop("\tError: The argument of rep should be of length of 1!\n");
	}
		
	# reading data from file
	data <-  try(read.csv(file = file, header = T, na.strings=na.code), silent = TRUE);
	if(identical(class(data), "try-error"))
	{
		stop("\tError: There are some problems on in the file!\n");
	}
	
	# checking whether all design factor names exists 
	# if true, make all the design factor to be factor
	col.names <- colnames(data);
	if(!(geno %in% col.names))
	{
		stop("\tError: The argument of geno does not match any column in the file!\n");
	} else
	{
		geno <- as.character(geno);
		data[, geno] <- factor(trimStrings(data[, geno]));
	}
	if(all(resp.var %in% col.names))
		stop("\tError: The argument of resp.var does not match columns in the file!\n");
	if(exptl.design == "RCB" || exptl.design == "AugRCB")
	{
		if(!(block %in% col.names))
		{
			stop("\tError: The argument of block does not match any column in the file!\n");
		} else
		{
			block <- as.character(block);
			data[, block] <- factor(trimStrings(data[,block]));
		}
	} else if(exptl.design == "AugLS")
	{
		if(!(row %in% col.names))
		{
			stop("\tError: The argument of row does not match any column in the file!\n");
		} else
		{
			row <- as.character(row);
			
		}
		if(!(column %in% col.names))
			stop("\tError: The argument of column does not match any column in the file!\n");
	} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
	{
		if(!(block %in% col.names))
			stop("\tError: The argument of block does not match any column in the file!\n");
		if(!(rep %in% col.names))
			stop("\tError: The argument of rep does not match any column in the file!\n");
	} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
	{
		if(!(row %in% col.names))
			stop("\tError: The argument of row does not match any column in the file!\n");
		if(!(column %in% col.names))
			stop("\tError: The argument of column does not match any column in the file!\n");
		if(!(rep %in% col.names))
			stop("\tError: The argument of rep does not match any column in the file!\n");
	}
	
	if(is.null(env))
	{
		env = "EnvLevel";
		data = cbind(data, EnvLevel = 1);
	} else
	{
		if(length(env) != 1)
			stop("\tError: The argument of env should of length of 1!\n");
		if(!(env %in% col.names))
			stop("\tError: The argument of env does not match any column in the file!\n");
	}
	
	## checking all the design factor not include any missing value
	##
	## not implemented right now
	
	
	
}
