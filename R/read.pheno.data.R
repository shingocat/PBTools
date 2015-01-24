###############################################################################
# TODO: reading phenotypic data from phenodata and formating data to a list for analyzing later;
# 
# ARGUMENT:
# phenodata - the phenodata of phenotypic data, should be of csv, txt or data.frame of R format with header;
# type - the phenotypic data type, whether it is raw data or means
# pop.type - a character string descrided population type, only accepted "SSSL", "PY" and "IL" value;
# gene.num - the number of gene for pyramided line design, it must be specified when population type is pyramieded lines
#			and only accepted 2, 3, and 4 genes by far.
# resp.var - a vector of strings; The trait names are in the phenodata; 
# exptl.design - the experimental design, should be of value on RCB, AugRCB, AugLS, Alpha, RowCol, LatinAlpha, LatinRowCol;
# geno - a string; variable name of the treatment/genotype variable;
# block - a string; variable name of the block variable
# row - a string; variable name of the row variable
# column - a string; variable name of the column variable; NULL, if design is RCB, Alpha, LatinAlpha
# rep - a string; variable name of the replication variable; NULL, if design is RCB
# env - a string; variable name of the environment variable
# na.code - a string code of missing value, default is NA;
#
# Author: mqin
# Date: Oct 15, 2014
# phenodataName: read.pheno.data.R
###############################################################################


read.pheno.data <- function(
		phenodata,
		type = c("RAW", "MEAN"),
		pop.type = c("SSSL", "PL", "IL"),
		gene.num,
		resp.var,
		geno,
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),	
		block = NULL,
		row = NULL,
		column = NULL,
		rep = NULL,
		env = NULL,
		na.code = NA,
		sep,
		...
)
{
	#--- mandatory arguments---#
	
	if(missing(phenodata))
		stop("\tError: The argument of phenodata could not be null!\n");
	if(missing(type))
		stop("\tError: The argument of type could not be null!\n");
	type <- match.arg(type);	
#	if(!is.character(type))
#		stop("\tError: The argument of type should be of character type!\n");
#	if(length(type) != 1)
#		stop("\tError: The argument of type should be of length one.\n");
	if(missing(pop.type))
		stop("\tError: The argument of pop.type could not be null!\n");
	pop.type <- match.arg(pop.type);	
#	if(!is.character(pop.type))
#		stop("\tError: The argument of pop.type should be of character type!\n");
#	if(length(pop.type) != 1)
#		stop("\tError: The argument of pop.type should be of length one.\n");
	if(missing(resp.var))
		stop("\tError: The argument of resp.var could not be null!\n");
	if(missing(geno))
		stop("\tError: The argument of geno could not be null!\n");
	if(!is.character(geno))
		stop("\tError: The argument of geno should be of character!\n");
	if(length(geno) != 1)
		stop("\tError: The argument of geno should be of length of 1!\n");
	if(missing(na.code))
		stop("\tError: The argument of na.code could not be null!\n");
	#--- checking whether the phenotypic data is mean or raw data.---#
	if(identical(type, "RAW"))
	{
		#--- if it is raw phenotypic data, the experimental design could not be null.---#
		if(missing(exptl.design))
			stop("\tError: The argument of exptl.design could not be null when phenotypic data is 'RAW' type!\n");
		exptl.design <- match.arg(exptl.design);
		if(exptl.design == "RCB" || exptl.design == "AugRCB")
		{
			if(missing(block))
				stop("\tError: The argument of block could not be null when the exptl.design is RCB or AugRCB!\n");
			if(length(block) != 1)
				stop("\tError: The argument of block should be of length of 1!\n");
		} else if(exptl.design == "AugLS")
		{
			if(missing(row))
				stop("\tError: The argument of row could not be null when the exptl.design is AugLS!\n");
			if(length(row) != 1)
				stop("\tError: The argument of row should be of length of 1!\n");
			if(missing(column))
				stop("\tError: The argument of column could not be null when the exptl.design is AugLS!\n");
			if(length(column) != 1)
				stop("\tError: The argument of column should be of length of 1!\n");
		} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
		{
			if(missing(block))
				stop("\tError: The argument of block could not be null when the exptl.design is Alpha or LattinAlpha!\n");
			if(length(block) != 1)
				stop("\tError: The argument of block should be of length of 1!\n");
			if(missing(rep))
				stop("\tError: The argument of rep could not be null when the exptl.design is Alpha or LattinAlpha!\n");
			if(length(rep) != 1)
				stop("\tError: The argument of rep should be of length of 1!\n");
		} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
		{
			if(missing(row))
				stop("\tError: The argument of row could not be null when the exptl.design is RowCol or LatinRowCol!\n");
			if(length(row) != 1)
				stop("\tError: The argument of row should be of length of 1!\n");
			if(missing(column))
				stop("\tError: The argument of column could not be null when the exptl.design is RowCol or LatinRowCol!\n");
			if(length(column) != 1)
				stop("\tError: The argument of column should be of length of 1!\n");
			if(missing(rep))
				stop("\tError: The argument of rep could not be null when the exptl.design is RowCol or LatinRowCol!\n");
			if(length(rep) != 1)
				stop("\tError: The argument of rep should be of length of 1!\n");
		}
		
	} 
	
	#--- checking the population type ---#
	#--- if the population type is pyramided line, then the gene.num argument could not be null---#
	if(identical(pop.type, "PL"))
	{
		if(missing(gene.num))
			stop("\tError: The argument of gen.num could not be null when pop.type is \"PL\" and only accept one of integer 2, 3 and 4.\n");
		if(!is.numeric(gene.num))
			stop("\tError: The argument of gen.num only accepted one of integer 2, 3 and 4.\n");
		if(length(gene.num) != 1)
			stop("\tError: The argument of gen.num only accepted one of integer 2, 3 and 4.\n");
	}
	
	UseMethod("read.pheno.data");
}

#--- This S3 method is for character class, is used to read data from file---#
read.pheno.data.character <- function(
		phenodata,
		type = c("RAW", "MEAN"),
		pop.type = c("SSSL", "PL", "IL"),
		gene.num,
		resp.var,
		geno,
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),	
		block = NULL,
		row = NULL,
		column = NULL,
		rep = NULL,
		env = NULL,
		na.code = NA,
		sep,
		...
)
{	
	if(!file.exists(phenodata))
		stop("\tError: The phenodata file does not exist! Please checking again!\n");
	if(missing(sep))	
	{
		phenodata <- try(read.csv(phenodata, header = TRUE, na.strings = na.code, check.names = FALSE), silent = TRUE);
	}else
	{
		phenodata <- try(read.csv(phenodata, header = TRUE, na.strings = na.code, check.names = FALSE, sep = sep), silent = TRUE);
	}
		
	if(inherits(phenodata, "try-error"))
	{
		stop(paste("\tError: There are some problems in reading phenotypic data file!\n\t\t", cat(phenodata), "\n", sep = ""));		
	} else
	{
		#--- checking whether the arguments is specified---#
		if(type == "RAW")
		{
			if(pop.type == "PL")
			{
				if(exptl.design == "RCB" || exptl.design == "AugRCB")
				{
					if(missing(env))
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										block = block, na.code = na.code), silent = TRUE);
					} else
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										block = block, env = env, na.code = na.code), silent = TRUE);
					}
				} else if (exptl.design == "AugLS")
				{
					if(missing(env))
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										row = row, column = column, na.code = na.code), silent = TRUE);
					} else
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										row = row, column = column, env = env, na.code = na.code), silent = TRUE);
					}
				} else if (exptl.design == "Alpha" || exptl.design == "LatinAlpha")
				{
					if(missing(env))
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										block = block, rep = rep, na.code = na.code), silent = TRUE);
					} else
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										block = block, rep = rep, env = env, na.code = na.code), silent = TRUE);
					}
				} else if (exptl.design == "RowCol" || exptl.design == "LatinRowCol")
				{
					if(missing(env))
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										row = row, column = column, rep = rep, na.code = na.code), silent = TRUE);
					} else
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										row = row, column = column, rep = rep, env = env, na.code = na.code), silent = TRUE);
					}
				}
			} else
			{
				if(exptl.design == "RCB" || exptl.design == "AugRCB")
				{
					if(missing(env))
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										block = block, na.code = na.code), silent = TRUE);
					} else
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										block = block, env = env, na.code = na.code), silent = TRUE);
					}
				} else if (exptl.design == "AugLS")
				{
					if(missing(env))
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										row = row, column = column, na.code = na.code), silent = TRUE);
					} else
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										row = row, column = column, env = env, na.code = na.code), silent = TRUE);
					}
				} else if (exptl.design == "Alpha" || exptl.design == "LatinAlpha")
				{
					if(missing(env))
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										block = block, rep = rep, na.code = na.code), silent = TRUE);
					} else
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										block = block, rep = rep, env = env, na.code = na.code), silent = TRUE);
					}
				} else if (exptl.design == "RowCol" || exptl.design == "LatinRowCol")
				{
					if(missing(env))
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										row = row, column = column, rep = rep, na.code = na.code), silent = TRUE);
					} else
					{
						phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
										resp.var = resp.var, geno = geno, exptl.design = exptl.design,
										row = row, column = column, rep = rep, env = env, na.code = na.code), silent = TRUE);
					}
				}
			} #--- end stmt of if(pop.type == "PL") ---#
		} else
		{
			if(pop.type == "PL")
			{
				if(missing(env))
				{	
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num, 
									resp.var = resp.var, geno = geno, na.code = na.code), silent = TRUE);
				} else 
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num, 
									resp.var = resp.var, geno = geno, env = env, na.code = na.code), silent = TRUE);
				}
			} else
			{
				if(missing(env))
				{	
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, na.code = na.code), silent = TRUE);
				} else 
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type,
									resp.var = resp.var, geno = geno, env = env, na.code = na.code), silent = TRUE);
				}
			}
		} #--- end stmt of if(type == "RAW") ---#
#		data <- try(innerFunc.read.pheno.data(data, type = type, pop.type = pop.type, gene.num = gene.num,
#						exptl.design = exptl.design, resp.var = resp.var, geno = geno, block = block, 
#						row = row, column = column, rep = rep, env = env, na.code = na.code), silent = TRUE);
		if(inherits(phenodata, "try-error"))
			stop(paste("\tError: There are some problems in reading phenotypic data on inner function!\n\t\t", cat(phenodata), "\n", sep = ""));		
	} #--- end stmt of if(inherits(phenodata, "try-error")) ---#
	return(phenodata);
}

#--- This S3 method is for data.frame, is used to read data from data.frame in R---#
read.pheno.data.data.frame <- function(
		phenodata,
		type = c("RAW", "MEAN"),
		pop.type = c("SSSL", "PL", "IL"),
		gene.num,
		resp.var,
		geno,
		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),	
		block = NULL,
		row = NULL,
		column = NULL,
		rep = NULL,
		env = NULL,
		na.code = NA
)
{
	#--- checking whether the arguments is specified---#
	if(type == "RAW")
	{
		if(pop.type == "PL")
		{
			if(exptl.design == "RCB" || exptl.design == "AugRCB")
			{
				if(missing(env))
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									block = block, na.code = na.code), silent = TRUE);
				} else
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									block = block, env = env, na.code = na.code), silent = TRUE);
				}
			} else if (exptl.design == "AugLS")
			{
				if(missing(env))
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									row = row, column = column, na.code = na.code), silent = TRUE);
				} else
				{
					data <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									row = row, column = column, env = env, na.code = na.code), silent = TRUE);
				}
			} else if (exptl.design == "Alpha" || exptl.design == "LatinAlpha")
			{
				if(missing(env))
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									block = block, rep = rep, na.code = na.code), silent = TRUE);
				} else
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									block = block, rep = rep, env = env, na.code = na.code), silent = TRUE);
				}
			} else if (exptl.design == "RowCol" || exptl.design == "LatinRowCol")
			{
				if(missing(env))
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									row = row, column = column, rep = rep, na.code = na.code), silent = TRUE);
				} else
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num,
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									row = row, column = column, rep = rep, env = env, na.code = na.code), silent = TRUE);
				}
			}
		} else
		{
			if(exptl.design == "RCB" || exptl.design == "AugRCB")
			{
				if(missing(env))
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									block = block, na.code = na.code), silent = TRUE);
				} else
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									block = block, env = env, na.code = na.code), silent = TRUE);
				}
			} else if (exptl.design == "AugLS")
			{
				if(missing(env))
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									row = row, column = column, na.code = na.code), silent = TRUE);
				} else
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									row = row, column = column, env = env, na.code = na.code), silent = TRUE);
				}
			} else if (exptl.design == "Alpha" || exptl.design == "LatinAlpha")
			{
				if(missing(env))
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									block = block, rep = rep, na.code = na.code), silent = TRUE);
				} else
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									block = block, rep = rep, env = env, na.code = na.code), silent = TRUE);
				}
			} else if (exptl.design == "RowCol" || exptl.design == "LatinRowCol")
			{
				if(missing(env))
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									row = row, column = column, rep = rep, na.code = na.code), silent = TRUE);
				} else
				{
					phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
									resp.var = resp.var, geno = geno, exptl.design = exptl.design,
									row = row, column = column, rep = rep, env = env, na.code = na.code), silent = TRUE);
				}
			}
		} #--- end stmt of if(pop.type == "PL") ---#
	} else
	{
		if(pop.type == "PL")
		{
			if(missing(env))
			{	
				phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num, 
								resp.var = resp.var, geno = geno, na.code = na.code), silent = TRUE);
			} else 
			{
				phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, gene.num = gene.num, 
								resp.var = resp.var, geno = geno, env = env, na.code = na.code), silent = TRUE);
			}
		} else
		{
			if(missing(env))
			{	
				phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type, 
								resp.var = resp.var, geno = geno, na.code = na.code), silent = TRUE);
			} else 
			{
				phenodata <- try(innerFunc.read.pheno.data(phenodata, type = type, pop.type = pop.type,
								resp.var = resp.var, geno = geno, env = env, na.code = na.code), silent = TRUE);
			}
		}
	} #--- end stmt of if(type == "RAW") ---#
	if(inherits(phenodata, "try-error"))
		stop(paste("\tError: There are some problems in reading phenotypic data on inner function!\n\t\t", cat(phenodata), "\n", sep = ""));		
	return(phenodata);
}

#--- inner function to read pheno data ---#
innerFunc.read.pheno.data <- function(
		phenodata,
		type,
		pop.type,
		gene.num,
		exptl.design,	
		resp.var,
		geno,
		block,
		row,
		column,
		rep,
		env,
		na.code
)
{
	#--- trim space of all value in phenodata
	phenodata <- apply(phenodata, 2, trimStrings);
	phenodata <- as.data.frame(phenodata);
	# if true, make all the design factor to be factor
	col.names <- colnames(phenodata);
	design.factor <- c();
	if(!(geno %in% col.names))
	{
		stop("\tError: The argument of geno does not match any column in the phenodata!\n");
	} else
	{
		geno <- as.character(geno);
		phenodata[, geno] <- factor(phenodata[, geno]);
		design.factor <- c(design.factor, geno);
		#--- checking the geno level should be more than two levels ---#
		if(nlevels(phenodata[, geno]) <= 1)
			stop("\tThe genotypic variable should be more than two levels!;\n");
	}
	# --- checking geno factor should have correct coded labels and levels when population type is pyramided line--- #
	if(identical(pop.type, "PL"))
	{
		genoLevels <- levels(phenodata[, geno]);
		
		# -- levels of genes combinations should be depend on the parameters of input values.
		if(all(genoLevels %in% levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T)))
				|| all(genoLevels %in% levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C=c("CC","Cc","cc")), sep="", lex.order = T)))
				|| all(genoLevels %in% levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C=c("CC","Cc","cc"), D = c("DD","Dd","dd")), sep="", lex.order = T))))
		{
			# do noting right now
		} else
		{
			stop(
					paste("\tError: The geno argument should be coded as", 
							"aabb", "Aabb", "AAbb", "aaBb", "AaBb", "AABb", "aaBB", "AaBB", "AABB,", 
							" if it is bigenes design,\n",
							"\tor ",
							"aabbcc", "aabbCc", "aabbCC", "aaBbcc", "aaBbCc", "aaBbCC", "aaBBcc", "aaBBCc",
							"aaBBCC", "Aabbcc", "AabbCc", "AabbCC", "AaBbcc", "AaBbCc", "AaBbCC", "AaBBcc",
							"AaBBCc", "AaBBCC", "AAbbcc" ,"AAbbCc", "AAbbCC", "AABbcc", "AABbCc", "AABbCC",
							"AABBcc", "AABBCc", "AABBCC,",
							" if it is trigenes desing, \n",
							"\tor ",
							"aabbccdd", "aabbccDd", "aabbccDD", "aabbCcdd", "aabbCcDd", "aabbCcDD",
							"aabbCCdd", "aabbCCDd", "aabbCCDD", "aaBbccdd", "aaBbccDd", "aaBbccDD",
							"aaBbCcdd", "aaBbCcDd", "aaBbCcDD", "aaBbCCdd", "aaBbCCDd", "aaBbCCDD",
							"aaBBccdd", "aaBBccDd", "aaBBccDD", "aaBBCcdd", "aaBBCcDd", "aaBBCcDD",
							"aaBBCCdd", "aaBBCCDd", "aaBBCCDD", "Aabbccdd", "AabbccDd", "AabbccDD",
							"AabbCcdd", "AabbCcDd", "AabbCcDD", "AabbCCdd", "AabbCCDd", "AabbCCDD",
							"AaBbccdd", "AaBbccDd", "AaBbccDD", "AaBbCcdd", "AaBbCcDd", "AaBbCcDD",
							"AaBbCCdd", "AaBbCCDd", "AaBbCCDD", "AaBBccdd", "AaBBccDd", "AaBBccDD",
							"AaBBCcdd", "AaBBCcDd", "AaBBCcDD", "AaBBCCdd", "AaBBCCDd", "AaBBCCDD",
							"AAbbccdd", "AAbbccDd", "AAbbccDD", "AAbbCcdd", "AAbbCcDd", "AAbbCcDD",
							"AAbbCCdd", "AAbbCCDd", "AAbbCCDD", "AABbccdd", "AABbccDd", "AABbccDD",
							"AABbCcdd", "AABbCcDd", "AABbCcDD", "AABbCCdd", "AABbCCDd", "AABbCCDD",
							"AABBccdd", "AABBccDd", "AABBccDD", "AABBCcdd", "AABBCcDd", "AABBCcDD",
							"AABBCCdd", "AABBCCDd", "AABBCCDD,",
							" if it is quadragenes design\n",
							
							sep = " ") );
		}
	}
	
	if(!all(resp.var %in% col.names))
		stop("\tError: The specified trait name(s) in resp.var does not match column in the phenodata!\n");
	if(type == "RAW"){
		if(exptl.design == "RCB" || exptl.design == "AugRCB")
		{
			if(!(block %in% col.names))
			{
				stop("\tError: The argument of block does not match any column in the phenodata!\n");
			} else
			{
				block <- as.character(block);
				phenodata[, block] <- factor(phenodata[,block]);
				design.factor <- c(design.factor, block);
				#--- checking the block variable should be more than two levels! ---#
				if(nlevels(phenodata[, block]) <= 1)
					stop("\tError: The block variable should be more than two levels!\n");
			}
		} else if(exptl.design == "AugLS")
		{
			if(!(row %in% col.names))
			{
				stop("\tError: The argument of row does not match any column in the phenodata!\n");
			} else
			{
				row <- as.character(row);
				phenodata[, row] <- factor(phenodata[,row]);
				design.factor <- c(design.factor, row);
				#--- checking the row variable should be more than two levels!---#
				if(nlevels(phenodata[,row]) <= 1)
					stop("\tError: The row variable should not be less than two levels!\n");
			}
			if(!(column %in% col.names))
			{
				stop("\tError: The argument of column does not match any column in the phenodata!\n");
			} else
			{
				column <- as.character(column);
				phenodata[, column] <- factor(phenodata[,column]);
				design.factor <- c(design.factor, column);
				#---checking the column variable should be more than two levels!---#
				if(nlevels(phenodata[,column]) <= 1)
					stop("\tError: The column vairable should not be less than two levels!\n");
			}
		} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
		{
			if(!(block %in% col.names))
			{
				stop("\tError: The argument of block does not match any column in the phenodata!\n");
			} else
			{
				block <- as.character(block);
				phenodata[ , block] <- factor(phenodata[ , block]);
				design.factor <- c(design.factor, block);
				#--- checking block variable should not be less than two levels!---#
				if(nlevels(phenodata[ , block]) <= 1)
					stop("\tError: The block variable should not be less than two levels!\n");
			}
			if(!(rep %in% col.names))
			{
				stop("\tError: The argument of rep does not match any column in the phenodata!\n");
			} else
			{
				rep <- as.character(rep);
				phenodata[ , rep] <- factor(phenodata[ , rep]);
				design.factor <- c(design.factor, rep);
				#--- checking rep variable should not be less than two levels!---#
				if(nlevels(phenodata[ , rep]) <= 1)
					stop("\tError: The rep variable should not be less than two levels!\n");
			}
		} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
		{
			if(!(row %in% col.names))
			{
				stop("\tError: The argument of row does not match any column in the phenodata!\n");
			} else
			{
				row <- as.charachter(row);
				phenodata[,row] <- factor(phenodata[ , row]);
				design.factor <- c(design.factor, row);
				#--- checking row variable should not be less than two levels!---#
				if(nlevels(phenodata[,row]) <= 1)
					stop("\tError: The row variable should not be less than two levels!\n");
			}
			if(!(column %in% col.names))
			{
				stop("\tError: The argument of column does not match any column in the phenodata!\n");
			} else
			{
				column <- as.character(column);
				phenodata[ , column] <- factor(phenodata[ , column]);
				design.factor <- c(design.factor, column);
				#--- checking column variable should not be less than two levels!---#
				if(nlevels(phenodata[, column]) <= 1)
					stop("\tError: The column variable should not be less than two levels!\n");
			}
			if(!(rep %in% col.names))
			{
				stop("\tError: The argument of rep does not match any column in the phenodata!\n");
			} else
			{
				rep <- as.character(rep);
				phenodata[ , rep] <- factor(phenodata[ , rep]);
				design.factor <- c(design.factor, rep);
				#--- checking rep variable should not be less than two levels!---#
				if(nlevels(phenodata[,rep]) <= 1)
					stop("\tError: The rep variable should not be less than two levels!\n");
			}
		}
		
		# --- if design is Latinized Row-Column, check if the data follow case1 or case3 labeling --- #
		if(exptl.design == "LatinRowCol")
		{
			lengthPerCross <- tapply(phenodata[,respvar[i]], phenodata[ ,c(row, column)], length);
			if(all(lengthPerCross <= 1, na.rm = TRUE))
			{
				if(nlevels(phenodata[ , row]) > nlevels(phenodata[ , column]))
				{
					longerRow <- TRUE;
				} else
				{
					longerRow <- FALSE;
				}
			} else
			{
				stop("The levels of the row/column variable should be continuous across replicates.");
			}
		}
	}
	#--- modify on 2014-11-20 ---#
	#--- if user does not specify env argument all phenotypic data treating as on the same env---#
	isSingleEnv <- TRUE;
	
	if(missing(env))
	{
		env = "EnvLevel";
		env.name = "DefaultEnv";
		phenodata = cbind(phenodata, EnvLevel = env.name);
		phenodata[[env]] <- factor(phenodata[[env]]);
		isSingleEnv <- TRUE;
	} else
	{
		if(!is.character(env))
			stop("\tError: The argument of env should of character type!\n");
		if(length(env) != 1)
			stop("\tError: The argument of env should of length of 1!\n");
		if(!(env %in% col.names))
			stop("\tError: The argument of env does not match any column in the phenodata!\n");
		env <- as.character(env);
		phenodata[ , env] <- factor(phenodata[ , env]);
		env.number <- nlevels(phenodata[ , env]);
		if(env.number > 1)
			isSingleEnv <- FALSE;
	}
	design.factor <- c(design.factor, env);
	
	## checking all the design factor not include any missing value
	for( i in 1:length(design.factor))
	{
		if(na.code %in% phenodata[ , design.factor[i]] || NA %in% phenodata[ , design.factor[i]])
		{
			sotp(paste("\tError: The design factor ", design.factor[i], " could not contain any missing value!\n", sep = ""));
		}
	}
	
	#--- format the na.code of response variable to R default code "NA"---#
	if(!(identical(na.code, NA) || identical(na.code, "NA")))
	{
		for(i in 1:length(resp.var))
		{
			phenodata[ , resp.var[i]][which(phenodata[ , resp.var[i]] == na.code)] <- NA;
		}
	}
	
	#--- building and reformating to a list structer ---#
	outcomes <- list();
	outcomes$isRestricted <- FALSE;
	outcomes$isMean <- ifelse(type == "MEAN", TRUE, FALSE);
	outcomes$isSingleEnv <- isSingleEnv;
	outcomes$na.code <- na.code;
	outcomes$pop.type <- pop.type;
	if(identical(pop.type, "PL"))
	{
		outcomes$gene.num <- gene.num;
	} else
	{
		outcomes$gene.num <- nlevels(phenodata[,geno]);
	}
	outcomes$gene.code <- levels(phenodata[,geno]);
	outcomes$trait.names <- resp.var;
	outcomes$trait.number <- length(resp.var);
	outcomes$raw.data <- phenodata;
	outcomes$traits <- list();
	for(i in 1:length(resp.var))
	{	
		temp.data <- phenodata[ , c(design.factor, resp.var[i])];
		outcomes$traits[[i]] <- list();
		outcomes$traits[[i]]$name <- resp.var[i];
		outcomes$traits[[i]]$env.number <- nlevels(temp.data[,env]);
		outcomes$traits[[i]]$env.names <- levels(temp.data[,env]);
		outcomes$traits[[i]]$envs <- list();
		if(isSingleEnv)
		{	
			outcomes$traits[[i]]$envs[[1]] <- list();
			outcomes$traits[[i]]$envs[[1]]$name <- levels(temp.data[,env]);
			outcomes$traits[[i]]$envs[[1]]$obsread <- nrow(temp.data);
			outcomes$traits[[i]]$envs[[1]]$obsused <- nrow(subset(temp.data, subset= (is.na(temp.data[ , resp.var[i]]) == FALSE)));
			outcomes$traits[[i]]$envs[[1]]$missing.rate <- 1 - outcomes$traits[[i]]$envs[[1]]$obsused / outcomes$traits[[i]]$envs[[1]]$obsread;
			missing.rate.by.geno <- aggregate(temp.data[[resp.var[i]]], by = list(temp.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
			names(missing.rate.by.geno) <- c(geno, "RATE");
			outcomes$traits[[i]]$envs[[1]]$missing.rate.by.geno <- missing.rate.by.geno;
			outcomes$traits[[i]]$envs[[1]]$data <- temp.data;
			if(type == "RAW")
			{
				outcomes$traits[[i]]$envs[[1]]$design <- list();
				outcomes$traits[[i]]$envs[[1]]$design$exptl.design <- exptl.design;
				outcomes$traits[[i]]$envs[[1]]$design$geno <- geno;
				outcomes$traits[[i]]$envs[[1]]$design$geno.levels <- levels(temp.data[ , geno]);
				if(!missing(block))
				{
					outcomes$traits[[i]]$envs[[1]]$design$block <- block;
					outcomes$traits[[i]]$envs[[1]]$design$block.levels <- levels(temp.data[ , block]);
				}
				if(!missing(row))
				{
					outcomes$traits[[i]]$envs[[1]]$design$row <- row;
					outcomes$traits[[i]]$envs[[1]]$design$row.levels <- levels(temp.data[ , row]);
				}
				if(!missing(column))
				{
					outcomes$traits[[i]]$envs[[1]]$design$column <- column;
					outcomes$traits[[i]]$envs[[1]]$design$column.levels <- levels(temp.data[ , column]);
				}
				if(exptl.design == "LatinRowCol")
					outcomes$traits[[i]]$sites[[1]]$longerRow <- longerRow;
				if(!missing(rep))
				{
					outcomes$traits[[i]]$envs[[1]]$design$rep <- rep;
					outcomes$traits[[i]]$envs[[1]]$design$rep.levels <- levels(temp.data[ , rep]);
				}
				
				outcomes$traits[[i]]$envs[[1]]$design$env <- env;
				outcomes$traits[[i]]$envs[[1]]$design$env.levels <- levels(temp.data[ , env]);
			} else
			{
				outcomes$traits[[i]]$envs[[1]]$design <- list();
				outcomes$traits[[i]]$envs[[1]]$design$geno <- geno;
				outcomes$traits[[i]]$envs[[1]]$design$geno.levels <- levels(temp.data[ , geno]);
				outcomes$traits[[i]]$envs[[1]]$design$env <- env;
				outcomes$traits[[i]]$envs[[1]]$design$env.levels <- levels(temp.data[ , env]);
			}
		}else
		{
			env.levels <- levels(temp.data[[env]]);
			for(j in 1:length(env.levels))
			{
				# --- create temp.data with one environmental level only --- #
				env.data <- subset(temp.data, subset = (temp.data[[env]] == env.levels[j]));
				outcomes$traits[[i]]$envs[[j]] <- list();
				outcomes$traits[[i]]$envs[[j]]$name <- env.levels[j];
				outcomes$traits[[i]]$envs[[j]]$obsread <- nrow(env.data);
				outcomes$traits[[i]]$envs[[j]]$obsused <- nrow(subset(env.data, subset= (is.na(env.data[ , resp.var[i]]) == FALSE)));
				outcomes$traits[[i]]$envs[[j]]$missing.rate <- 1 - outcomes$traits[[i]]$envs[[j]]$obsused / outcomes$traits[[i]]$envs[[j]]$obsread;
				#outcomes$traits[[i]]$sites[[j]]$missing.rate.by.geno <- aggregate(temp.data[[resp.var[i]]], by = list(temp.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
				missing.rate.by.geno <- aggregate(env.data[[resp.var[i]]], by = list(env.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
				names(missing.rate.by.geno) <- c(geno, "RATE");
				outcomes$traits[[i]]$envs[[j]]$missing.rate.by.geno <- missing.rate.by.geno;
				outcomes$traits[[i]]$envs[[j]]$data <-env.data;
				if(type == "RAW")
				{
					outcomes$traits[[i]]$envs[[j]]$design <- list();
					outcomes$traits[[i]]$envs[[j]]$design$exptl.design <- exptl.design;
					outcomes$traits[[i]]$envs[[j]]$design$geno <- geno;
					outcomes$traits[[i]]$envs[[j]]$design$geno.levels <- levels(env.data[ , geno]);
					if(!missing(block))
					{
						outcomes$traits[[i]]$envs[[j]]$design$block <- block;
						outcomes$traits[[i]]$envs[[j]]$design$block.levels <- levels(env.data[ , block]);
					}
					if(!missing(row))
					{
						outcomes$traits[[i]]$envs[[j]]$design$row <- row;
						outcomes$traits[[i]]$envs[[j]]$design$row.levels <- levels(env.data[ , row]);
					}
					if(!missing(column))
					{
						outcomes$traits[[i]]$envs[[j]]$design$column <- column;
						outcomes$traits[[i]]$envs[[j]]$design$column.levels <- levels(env.data[ , column]);
					}
					if(exptl.design == "LatinRowCol")
						outcomes$traits[[i]]$envs[[j]]$longerRow <- longerRow;
					if(!missing(rep))
					{
						outcomes$traits[[i]]$envs[[j]]$design$rep <- rep;
						outcomes$traits[[i]]$envs[[j]]$design$rep.levels <- levels(env.data[ , rep]);
					}
					outcomes$traits[[i]]$envs[[j]]$design$env <- env;
					outcomes$traits[[i]]$envs[[j]]$design$env.levels <- levels(env.data[ , env]);
				} else
				{
					outcomes$traits[[i]]$envs[[j]]$design <- list();
					outcomes$traits[[i]]$envs[[j]]$design$geno <- geno;
					outcomes$traits[[i]]$envs[[j]]$design$geno.levels <- levels(temp.data[ , geno]);
					outcomes$traits[[i]]$envs[[j]]$design$env <- env;
					outcomes$traits[[i]]$envs[[j]]$design$env.levels <- levels(temp.data[ , env]);
				}
			}
		} ## end of statement if(isSingleEnv)
	} ## end of statement of for(i in 1:length(resp.var))
	class(outcomes) <- "PhenotypicData";
	return(outcomes);
}

#--- This S3 method is for default, is used to read data from other class in R ---#
#read.pheno.data.default <- function(
#		phenodata,
#		type = c("RAW", "MEAN"),
#		pop.type = c("SSSL","PL", "IL"),
#		gene.num,
#		exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"),	
#		resp.var,
#		geno,
#		block,
#		row = NULL,
#		column = NULL,
#		rep = NULL,
#		env = NULL,
#		na.code = NA
#)
#{
#	if(missing(phenodata))
#		stop("\tError: The argument of phenodata could not be null!\n");
#	if(inherits(phenodata, "character"))
#	{
#		if(!file.exists(phenodata))
#			stop("\tError: The phenodata file does not exist! Please checking again!\n");
#		isFile <- TRUE;
#	} else if(inherits(phenodata, "data.frame"))
#	{
#		isFile <- FALSE;
#	} else
#	{
#		stop("\tError: The argument of phenodata should be character of file path or data.frame format!\n");
#	}
#	if(missing(type))
#		stop("\tError: The argument of type could not be null.\n");
#	type <- match.arg(type);
#	if(missing(pop.type))
#		stop("\tError: The argument of pop.type could not be null!\n");
#	pop.type <- match.arg(pop.type);
##	if(!is.character(pop.type))
##		stop("\tError: The argument of pop.type should be of character type!\n");
##	if(length(pop.type) != 1)
##	{
##		stop("\tError: The argument of pop.type should be of length on 1!\n");
##	}
##	if(!(pop.type %in% c("SSSL", "PL", "IL")))
##	{
##		stop("\tError: The argument of pop.type should be of value \"SSSL\",\"PL\",\"IL\".\n");
##	}
#	#-- checking the argument of gen.num when the pop.type is PL --#
#	if(identical(pop.type, "PL"))
#	{
#		if(missing(gene.num))
#			stop("\tError: The argument of gen.num could not be null when pop.type is \"PL\" and only accept one of integer 2, 3 and 4.\n");
#		if(!is.numeric(gene.num))
#			stop("\tError: The argument of gen.num only accepted one of integer 2, 3 and 4.\n");
#		if(length(gene.num) != 1)
#			stop("\tError: The argument of gen.num only accepted one of integer 2, 3 and 4.\n");
#	}
#	
#	if(type == "RAW")
#	{
#		if(missing(exptl.design))
#		{
#			stop("\tError: The argument of exptl.design could not be null!\n");
#		}
#		exptl.design <- match.arg(exptl.design);
##		if(!is.character(exptl.design))
##			stop("\tError: The argument of exptl.design should be of character type!\n");
##		if(length(exptl.design) != 1)
##		{
##			stop("\tError: The argument of exptl.design should be of length on 1!\n");
##		}
##		if(!(exptl.design %in% c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol")))
##		{
##			stop("\tError: The argument of exptl.design should be of value \"RCB\", \"AugRCB\", \"AugLS\", \"Alpha\", \"RowCol\", \"LatinAlpha\", \"LatinRowCol\".\n");
##		}
#		if(exptl.design == "RCB" || exptl.design == "AugRCB")
#		{
#			if(missing(block))
#				stop("\tError: The argument of block could not be null when the exptl.design is RCB or AugRCB!\n");
#			if(!is.character(block))
#				stop("\tError: The argument of block should be of type character.\n");
#			if(length(block) != 1)
#				stop("\tError: The argument of block should be of length of 1!\n");
#		} else if(exptl.design == "AugLS")
#		{
#			if(missing(row))
#				stop("\tError: The argument of row could not be null when the exptl.design is AugLS!\n");
#			if(!is.character(row))
#				stop("\tError: The argument of row should be of character type!\n");
#			if(length(row) != 1)
#				stop("\tError: The argument of row should be of length of 1!\n");
#			if(missing(column))
#				stop("\tError: The argument of column could not be null when the exptl.design is AugLS!\n");
#			if(!is.character(column))
#				stop("\tError: The argument of column should be of character type!\n");
#			if(length(column) != 1)
#				stop("\tError: The argument of column should be of length of 1!\n");
#		} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
#		{
#			if(missing(block))
#				stop("\tError: The argument of block could not be null when the exptl.design is Alpha or LattinAlpha!\n");
#			if(!is.character(block))
#				stop("\tError: The argument of block should be of character type!\n");
#			if(length(block) != 1)
#				stop("\tError: The argument of block should be of length of 1!\n");
#			if(missing(rep))
#				stop("\tError: The argument of rep could not be null when the exptl.design is Alpha or LattinAlpha!\n");
#			if(!is.character(rep))
#				stop("\tError: The argument of rep should be of character type!\n");
#			if(length(rep) != 1)
#				stop("\tError: The argument of rep should be of length of 1!\n");
#		} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
#		{
#			if(missing(row))
#				stop("\tError: The argument of row could not be null when the exptl.design is RowCol or LatinRowCol!\n");
#			if(!is.character(row))
#				stop("\tError: The argument of row should be of character type!\n");
#			if(length(row) != 1)
#				stop("\tError: The argument of row should be of length of 1!\n");
#			if(missing(column))
#				stop("\tError: The argument of column could not be null when the exptl.design is RowCol or LatinRowCol!\n");
#			if(!is.character(column))
#				stop("\tError: The argument of column should be of character type!\n");
#			if(length(column) != 1)
#				stop("\tError: The argument of column should be of length of 1!\n");
#			if(missing(rep))
#				stop("\tError: The argument of rep could not be null when the exptl.design is RowCol or LatinRowCol!\n");
#			if(!is.character(rep))
#				stop("\tError: The argument of rep should be of character type!\n");
#			if(length(rep) != 1)
#				stop("\tError: The argument of rep should be of length of 1!\n");
#		}
#	}
#	if(missing(resp.var))
#		stop("\tError: The argument of resp.var could not be null!\n");
#	if(missing(geno))
#		stop("\tError: The argument of geno could not be null!\n");
#	if(!is.character(geno))
#		stop("\tError: The argument of geno should be of character type!\n");
#	if(length(geno) != 1)
#		stop("\tError: The argument of geno should be of length of 1!\n");
#	
#	
#	# reading data from file if the phenodata is a file
#	if(isFile)
#	{
#		data <-  try(read.csv(file = phenodata, header = T, na.strings=na.code), silent = TRUE);
#		if(identical(class(data), "try-error"))
#			stop("\tError: There are some problems on in the phenodata!\n");
#	} else
#	{
#		data <- as.data.frame(phenodata);
#	}
#	
#	#--- trim space of all value in the data ---#
#	data <- apply(data, 2, trimStrings);
#	
#	# checking whether all design factor names exists 
#	# if true, make all the design factor to be factor
#	col.names <- colnames(data);
#	design.factor <- c();
#	if(!(geno %in% col.names))
#	{
#		stop("\tError: The argument of geno does not match any column in the phenodata!\n");
#	} else
#	{
#		geno <- as.character(geno);
#		data[, geno] <- factor(data[, geno]);
#		design.factor <- c(design.factor, geno);
#		#--- checking the geno level should be more than two levels ---#
#		if(nlevels(data[, geno]) <= 1)
#			stop("\tThe genotypic variable should be more than two levels!;\n");
#	}
#	# --- checking geno factor should have correct coded labels and levels when population type is pyramided line--- #
#	if(identical(pop.type, "PL"))
#	{
#		genoLevels <- levels(data[, geno]);
#		
#		# -- levels of genes combinations should be depend on the parameters of input values.
#		if(all(genoLevels %in% levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T)))
#				|| all(genoLevels %in% levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C=c("CC","Cc","cc")), sep="", lex.order = T)))
#				|| all(genoLevels %in% levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C=c("CC","Cc","cc"), D = c("DD","Dd","dd")), sep="", lex.order = T))))
#		{
#			# do noting right now
#		} else
#		{
#			stop(
#					paste("\tError: The geno argument should be coded as", 
#							"aabb", "Aabb", "AAbb", "aaBb", "AaBb", "AABb", "aaBB", "AaBB", "AABB,", 
#							" if it is bigenes design,\n",
#							"\tor ",
#							"aabbcc", "aabbCc", "aabbCC", "aaBbcc", "aaBbCc", "aaBbCC", "aaBBcc", "aaBBCc",
#							"aaBBCC", "Aabbcc", "AabbCc", "AabbCC", "AaBbcc", "AaBbCc", "AaBbCC", "AaBBcc",
#							"AaBBCc", "AaBBCC", "AAbbcc" ,"AAbbCc", "AAbbCC", "AABbcc", "AABbCc", "AABbCC",
#							"AABBcc", "AABBCc", "AABBCC,",
#							" if it is trigenes desing, \n",
#							"\tor ",
#							"aabbccdd", "aabbccDd", "aabbccDD", "aabbCcdd", "aabbCcDd", "aabbCcDD",
#							"aabbCCdd", "aabbCCDd", "aabbCCDD", "aaBbccdd", "aaBbccDd", "aaBbccDD",
#							"aaBbCcdd", "aaBbCcDd", "aaBbCcDD", "aaBbCCdd", "aaBbCCDd", "aaBbCCDD",
#							"aaBBccdd", "aaBBccDd", "aaBBccDD", "aaBBCcdd", "aaBBCcDd", "aaBBCcDD",
#							"aaBBCCdd", "aaBBCCDd", "aaBBCCDD", "Aabbccdd", "AabbccDd", "AabbccDD",
#							"AabbCcdd", "AabbCcDd", "AabbCcDD", "AabbCCdd", "AabbCCDd", "AabbCCDD",
#							"AaBbccdd", "AaBbccDd", "AaBbccDD", "AaBbCcdd", "AaBbCcDd", "AaBbCcDD",
#							"AaBbCCdd", "AaBbCCDd", "AaBbCCDD", "AaBBccdd", "AaBBccDd", "AaBBccDD",
#							"AaBBCcdd", "AaBBCcDd", "AaBBCcDD", "AaBBCCdd", "AaBBCCDd", "AaBBCCDD",
#							"AAbbccdd", "AAbbccDd", "AAbbccDD", "AAbbCcdd", "AAbbCcDd", "AAbbCcDD",
#							"AAbbCCdd", "AAbbCCDd", "AAbbCCDD", "AABbccdd", "AABbccDd", "AABbccDD",
#							"AABbCcdd", "AABbCcDd", "AABbCcDD", "AABbCCdd", "AABbCCDd", "AABbCCDD",
#							"AABBccdd", "AABBccDd", "AABBccDD", "AABBCcdd", "AABBCcDd", "AABBCcDD",
#							"AABBCCdd", "AABBCCDd", "AABBCCDD,",
#							" if it is quadragenes design\n",
#							
#							sep = " ") );
#		}
#	}
#	
#	if(!all(resp.var %in% col.names))
#		stop("\tError: The specified trait name(s) in resp.var does not match column in the phenodata!\n");
#	if(type == "RAW")
#	{
#		if(exptl.design == "RCB" || exptl.design == "AugRCB")
#		{
#			if(!(block %in% col.names))
#			{
#				stop("\tError: The argument of block does not match any column in the phenodata!\n");
#			} else
#			{
#				block <- as.character(block);
#				data[, block] <- factor(data[,block]);
#				design.factor <- c(design.factor, block);
#				#--- checking the block variable should be more than two levels! ---#
#				if(nlevels(data[, block]) <= 1)
#					stop("\tError: The block variable should be more than two levels!\n");
#			}
#		} else if(exptl.design == "AugLS")
#		{
#			if(!(row %in% col.names))
#			{
#				stop("\tError: The argument of row does not match any column in the phenodata!\n");
#			} else
#			{
#				row <- as.character(row);
#				data[, row] <- factor(data[,row]);
#				design.factor <- c(design.factor, row);
#				#--- checking the row variable should be more than two levels!---#
#				if(nlevels(data[,row]) <= 1)
#					stop("\tError: The row variable should not be less than two levels!\n");
#			}
#			if(!(column %in% col.names))
#			{
#				stop("\tError: The argument of column does not match any column in the phenodata!\n");
#			} else
#			{
#				column <- as.character(column);
#				data[, column] <- factor(data[,column]);
#				design.factor <- c(design.factor, column);
#				#---checking the column variable should be more than two levels!---#
#				if(nlevels(data[,column]) <= 1)
#					stop("\tError: The column vairable should not be less than two levels!\n");
#			}
#		} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
#		{
#			if(!(block %in% col.names))
#			{
#				stop("\tError: The argument of block does not match any column in the phenodata!\n");
#			} else
#			{
#				block <- as.character(block);
#				data[ , block] <- factor(data[ , block]);
#				design.factor <- c(design.factor, block);
#				#--- checking block variable should not be less than two levels!---#
#				if(nlevels(data[ , block]) <= 1)
#					stop("\tError: The block variable should not be less than two levels!\n");
#			}
#			if(!(rep %in% col.names))
#			{
#				stop("\tError: The argument of rep does not match any column in the phenodata!\n");
#			} else
#			{
#				rep <- as.character(rep);
#				data[ , rep] <- factor(data[ , rep]);
#				design.factor <- c(design.factor, rep);
#				#--- checking rep variable should not be less than two levels!---#
#				if(nlevels(data[ , rep]) <= 1)
#					stop("\tError: The rep variable should not be less than two levels!\n");
#			}
#		} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
#		{
#			if(!(row %in% col.names))
#			{
#				stop("\tError: The argument of row does not match any column in the phenodata!\n");
#			} else
#			{
#				row <- as.charachter(row);
#				data[,row] <- factor(data[ , row]);
#				design.factor <- c(design.factor, row);
#				#--- checking row variable should not be less than two levels!---#
#				if(nlevels(data[,row]) <= 1)
#					stop("\tError: The row variable should not be less than two levels!\n");
#			}
#			if(!(column %in% col.names))
#			{
#				stop("\tError: The argument of column does not match any column in the phenodata!\n");
#			} else
#			{
#				column <- as.character(column);
#				data[ , column] <- factor(data[ , column]);
#				design.factor <- c(design.factor, column);
#				#--- checking column variable should not be less than two levels!---#
#				if(nlevels(data[, column]) <= 1)
#					stop("\tError: The column variable should not be less than two levels!\n");
#			}
#			if(!(rep %in% col.names))
#			{
#				stop("\tError: The argument of rep does not match any column in the phenodata!\n");
#			} else
#			{
#				rep <- as.character(rep);
#				data[ , rep] <- factor(data[ , rep]);
#				design.factor <- c(design.factor, rep);
#				#--- checking rep variable should not be less than two levels!---#
#				if(nlevels(data[,rep]) <= 1)
#					stop("\tError: The rep variable should not be less than two levels!\n");
#			}
#		}
#		
#		# --- if design is Latinized Row-Column, check if the data follow case1 or case3 labeling --- #
#		if(exptl.design == "LatinRowCol")
#		{
#			lengthPerCross <- tapply(temp.data[,respvar[i]], temp.data[ ,c(row, column)], length);
#			if(all(lengthPerCross <= 1, na.rm = TRUE))
#			{
#				if(nlevels(temp.data[ , row]) > nlevels(temp.data[ , column]))
#				{
#					longerRow <- TRUE;
#				} else
#				{
#					longerRow <- FALSE;
#				}
#			} else
#			{
#				stop("The levels of the row/column variable should be continuous across replicates.");
#			}
#		}
#	}
#	#--- modify on 2014-11-20 ---#
#	#--- if user does not specify env argument all phenotypic data treating as on the same env---#
#	isSingleEnv <- TRUE;
#	
#	if(is.null(env))
#	{
#		env = "EnvLevel";
#		env.name = "DefaultEnv";
#		data = cbind(data, EnvLevel = env.name);
#		data[[env]] <- factor(data[[env]]);
#	} else
#	{
#		if(!is.character(env))
#			stop("\tError: The argument of env should be of character type!\n");
#		if(length(env) != 1)
#			stop("\tError: The argument of env should be of length of one!\n");
#		if(!(env %in% col.names))
#			stop("\tError: The argument of env does not match any column in the phenodata!\n");
#		env <- as.character(env);
#		data[ , env] <- factor(data[ , env]);
#		env.number <- nlevels(data[, env]);
#		if(env.number > 1)
#			isSingleEnv <- FALSE;
#	}
#	design.factor <- c(design.factor, env);
#	
#	
#	## checking all the design factor not include any missing value
#	for( i in 1:length(design.factor))
#	{
#		if(na.code %in% data[ , design.factor[i]] || NA %in% data[ , design.factor[i]])
#		{
#			sotp(paste("\tError: The design factor ", design.factor[i], " could not contain any missing value!\n", sep = ""));
#		}
#	}
#	
#	#--- format the na.code of response variable to R default code "NA"---#
#	if(!(identical(na.code, NA) || identical(na.code, "NA")))
#	{
#		for(i in 1:length(resp.var))
#		{
#			data[ , resp.var[i]][which(data[ , resp.var[i]] == na.code)] <- NA;
#		}
#	}
#	
#	#--- building and reformating to a list structer ---#
#	outcomes <- list();
#	outcomes$isRestricted <- FALSE;
#	outcomes$isMean <- ifelse(type == "MEAN", TRUE, FALSE);
#	outcomes$isSingleEnv <- isSingleEnv;
#	outcomes$na.code <- na.code;
#	outcomes$pop.type <- pop.type;
#	if(identical(pop.type, "PL"))
#	{
#		outcomes$gene.num <- gene.num;
#	}
#	else
#	{
#		outcomes$gene.num <- nlevels(data[,geno]);
#	}
#	outcomes$trait.names <- resp.var;
#	outcomes$raw.data <- data;
##	outcomes$design <- list();
##	outcomes$design$exptl.design <- exptl.design;
##	outcomes$design$geno <- geno;
##	outcomes$design$geno.levels <- levels(data[ , geno]);
##	outcomes$design$block <- block;
##	outcomes$design$block.levels <- levels(data[ , block]);
##	outcomes$design$row <- row;
##	outcomes$design$row.levels <- levels(data[ , row]);
##	outcomes$design$column <- column;
##	outcomes$design$column.levels <- levels(data[ , column]);
##	outcomes$design$rep <- rep;
##	outcomes$design$rep.levels <- levels(data[ , rep]);
##	outcomes$design$env <- env;
##	outcomes$design$env.levels <- levels(data[ , env]);
#	## each trait build one list
#	outcomes$traits <- list();
#	for(i in 1:length(resp.var))
#	{	
#		temp.data <- data[ , c(design.factor, resp.var[i])];
#		outcomes$traits[[i]] <- list();
#		outcomes$traits[[i]]$name <- resp.var[i];
#		outcomes$traits[[i]]$env.number <- nlevels(temp.data[,env]);
#		outcomes$traits[[i]]$env.names <- levels(temp.data[,env]);
#		outcomes$traits[[i]]$envs <- list();
#		if(isSingleEnv)
#		{	
#			outcomes$traits[[i]]$envs[[1]] <- list();
#			outcomes$traits[[i]]$envs[[1]]$name <- levels(temp.data[,env]);
#			outcomes$traits[[i]]$envs[[1]]$obsread <- nrow(temp.data);
#			outcomes$traits[[i]]$envs[[1]]$obsused <- nrow(subset(temp.data, subset= (is.na(temp.data[ , resp.var[i]]) == FALSE)));
#			outcomes$traits[[i]]$envs[[1]]$missing.rate <- 1 - outcomes$traits[[i]]$sites[[1]]$obsused / outcomes$traits[[i]]$sites[[1]]$obsread;
#			missing.rate.by.geno <- aggregate(temp.data[[resp.var[i]]], by = list(temp.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
#			names(missing.rate.by.geno) <- c(geno, "RATE");
#			outcomes$traits[[i]]$envs[[1]]$missing.rate.by.geno <- missing.rate.by.geno;
#			outcomes$traits[[i]]$envs[[1]]$data <- temp.data;
#			if(type == "RAW")
#			{
#				outcomes$traits[[i]]$envs[[1]]$design <- list();
#				outcomes$traits[[i]]$envs[[1]]$design$exptl.design <- exptl.design;
#				outcomes$traits[[i]]$envs[[1]]$design$geno <- geno;
#				outcomes$traits[[i]]$envs[[1]]$design$geno.levels <- levels(temp.data[ , geno]);
#				outcomes$traits[[i]]$envs[[1]]$design$block <- block;
#				outcomes$traits[[i]]$envs[[1]]$design$block.levels <- levels(temp.data[ , block]);
#				outcomes$traits[[i]]$envs[[1]]$design$row <- row;
#				outcomes$traits[[i]]$envs[[1]]$design$row.levels <- levels(temp.data[ , row]);
#				outcomes$traits[[i]]$envs[[1]]$design$column <- column;
#				outcomes$traits[[i]]$envs[[1]]$design$column.levels <- levels(temp.data[ , column]);
#				if(exptl.design == "LatinRowCol")
#					outcomes$traits[[i]]$sites[[1]]$longerRow <- longerRow;
#				outcomes$traits[[i]]$envs[[1]]$design$rep <- rep;
#				outcomes$traits[[i]]$envs[[1]]$design$rep.levels <- levels(temp.data[ , rep]);
#				outcomes$traits[[i]]$envs[[1]]$design$env <- env;
#				outcomes$traits[[i]]$envs[[1]]$design$env.levels <- levels(temp.data[ , env]);
#			}
#		}else
#		{
#			env.levels <- levels(temp.data[[env]]);
#			for(j in 1:length(env.levels))
#			{
#				# --- create temp.data with one environmental level only --- #
#				env.data <- subset(temp.data, subset = (temp.data[[env]] == env.levels[j]));
#				outcomes$traits[[i]]$envs[[j]] <- list();
#				outcomes$traits[[i]]$envs[[j]]$name <- env.levels[j];
#				outcomes$traits[[i]]$envs[[j]]$obsread <- nrow(env.data);
#				outcomes$traits[[i]]$envs[[j]]$obsused <- nrow(subset(env.data, subset= (is.na(env.data[ , resp.var[i]]) == FALSE)));
#				outcomes$traits[[i]]$envs[[j]]$missing.rate <- 1 - outcomes$traits[[i]]$sites[[j]]$obsused / outcomes$traits[[i]]$sites[[j]]$obsread;
#				#outcomes$traits[[i]]$sites[[j]]$missing.rate.by.geno <- aggregate(temp.data[[resp.var[i]]], by = list(temp.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
#				missing.rate.by.geno <- aggregate(env.data[[resp.var[i]]], by = list(env.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
#				names(missing.rate.by.geno) <- c(geno, "RATE");
#				outcomes$traits[[i]]$envs[[j]]$missing.rate.by.geno <- missing.rate.by.geno;
#				outcomes$traits[[i]]$envs[[j]]$data <-env.data;
#				if(type == "RAW")
#				{
#					outcomes$traits[[i]]$envs[[j]]$design <- list();
#					outcomes$traits[[i]]$envs[[j]]$design$exptl.design <- exptl.design;
#					outcomes$traits[[i]]$envs[[j]]$design$geno <- geno;
#					outcomes$traits[[i]]$envs[[j]]$design$geno.levels <- levels(env.data[ , geno]);
#					outcomes$traits[[i]]$envs[[j]]$design$block <- block;
#					outcomes$traits[[i]]$envs[[j]]$design$block.levels <- levels(env.data[ , block]);
#					outcomes$traits[[i]]$envs[[j]]$design$row <- row;
#					outcomes$traits[[i]]$envs[[j]]$design$row.levels <- levels(env.data[ , row]);
#					outcomes$traits[[i]]$envs[[j]]$design$column <- column;
#					outcomes$traits[[i]]$envs[[j]]$design$column.levels <- levels(env.data[ , column]);
#					if(exptl.design == "LatinRowCol")
#						outcomes$traits[[i]]$envs[[j]]$longerRow <- longerRow;
#					outcomes$traits[[i]]$envs[[j]]$design$rep <- rep;
#					outcomes$traits[[i]]$envs[[j]]$design$rep.levels <- levels(env.data[ , rep]);
#					outcomes$traits[[i]]$envs[[j]]$design$env <- env;
#					outcomes$traits[[i]]$envs[[j]]$design$env.levels <- levels(env.data[ , env]);
#				}
#			}
#		} ## end of statement if(isSingleEnv)
#	} ## end of statement of for(i in 1:length(resp.var))
#	class(outcomes) <- "PhenotypicData";
#	return(outcomes);
#}

##################################################################################
# Todo: print the description of phenotypic data
#
# ARGUMENT:
#	data - the instance of class PhenotypicData
#	levels - the level of information to print, accepted only 1 or 2, i.e. 1 for concise, 2 for precise;
#

print.PhenotypicData <- function(data, level = 1)
{
	if(missing(level))
		level <- 1;
	if(!is.numeric(level))
		stop("\tError: The argument of level should be of value 1 or 2 where 1 is for concise details, and 2 is for precise details.\n");
	if(length(level) != 1)
		stop("\tError: The argument of level should be of length 1.\n");
	if(!any(level %in% c(1,2)))
		stop("\tError: The argument of level only accepted the integer 1 or 2. \n");
	cat("Phenotypic data Summary");
	if(data$isRestricted)
		cat("(After Restricted)");
	cat(":\n");	
	cat("\tThe population type is : ", if(data$pop.type == "IL") "Introgression Lines" else if (data$pop.type == "PL") "Pyramided Lines" else "Single Segment Substitution Lines", ".\n");
	cat("\tThe number of response variable is:", data$trait.number, ".\n");
	cat("\tSingle environment experiment :", if(data$isSingleEnv) "YES" else "NO", ".\n");
	cat("\tPhenotypic data is mean data :", if(data$isMean) "YES" else "NO", ".\n");
	if(data$isRestricted)
	{
		cat("\tRestricted conditions:\n");
		cat("\t\tMissing Rate = ", data$traits[[1]]$envs[[1]]$restricted$missing.rate.cond, ".\n", sep = "");
	} 
	cat(rep("-", 40), sep="");
	cat("\n");
	
	#else
	#{
		for(i in 1:data$trait.number)
		{
			cat("\t\t", "The trait name is:", data$traits[[i]]$name, ".\n");
			for(j in 1:data$traits[[i]]$env.number)
			{	
				cat("\t\t\t", "The environment is: ", data$traits[[i]]$envs[[j]]$name, ".\n", sep = " ");
				if(data$isRestricted)
					cat("\t\t\t", "Is restricted: ", ifelse(data$traits[[i]]$envs[[j]]$restricted$isTRUE, "YES", "NO"), ".\n", sep = " ")
				if(!data$isMean){					
					exptl.design <- data$traits[[i]]$envs[[j]]$design$exptl.design;
					cat("\t\t\t", "The experiment design is ", exptl.design, ".\n");
					if(level == 2)
					{					
						if(exptl.design == "RCB" || exptl.design == "AugRCB")
						{
							cat("\t\t\t", "The block variable is ", data$traits[[i]]$envs[[j]]$design$block, ".\n");
							cat("\t\t\t", "The levels of block variable are ", data$traits[[i]]$envs[[j]]$design$block.levels, ".\n");
						} else if(exptl.design == "AugLS")
						{
							cat("\t\t\t", "The row variable is ", data$traits[[i]]$envs[[j]]$design$row, ".\n");
							cat("\t\t\t", "The levels of row variable are ", data$traits[[i]]$envs[[j]]$design$row.levels, ".\n");
							cat("\t\t\t", "The column variable is ", data$traits[[i]]$envs[[j]]$design$column, ".\n");
							cat("\t\t\t", "The levels of column variable are ", data$traits[[i]]$envs[[j]]$design$column.levels, ".\n");
						} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
						{
							cat("\t\t\t", "The block variable is ", data$traits[[i]]$envs[[j]]$design$block, ".\n");
							cat("\t\t\t", "The levels of block variable are ", data$traits[[i]]$envs[[j]]$design$block.levels, ".\n");
							cat("\t\t\t", "The rep variable is ", data$traits[[i]]$envs[[j]]$design$rep, ".\n");
							cat("\t\t\t", "The levels of rep variable are ", data$traits[[i]]$envs[[j]]$design$rep.levels, ".\n");
						} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
						{
							cat("\t\t\t", "The row variable is ", data$traits[[i]]$envs[[j]]$design$row, ".\n");
							cat("\t\t\t", "The levels of row variable are ", data$traits[[i]]$envs[[j]]$design$row.levels, ".\n");
							cat("\t\t\t", "The column variable is ", data$traits[[i]]$envs[[j]]$design$column, ".\n");
							cat("\t\t\t", "The levels of column variable are ", data$traits[[i]]$envs[[j]]$design$column.levels, ".\n");
							cat("\t\t\t", "The rep variable is ", data$traits[[i]]$envs[[j]]$design$rep, ".\n");
							cat("\t\t\t", "The levels of rep variable are ", data$traits[[i]]$envs[[j]]$design$rep.levels, ".\n");
						}
					}
				}
				cat("\t\t\t", "The Genotypic variable is ", data$traits[[i]]$envs[[j]]$design$geno, ".\n");
				if(level == 2)
				{
					cat("\t\t\t", "The level of Genotypic variable is ", data$traits[[i]]$envs[[j]]$design$geno.levels, ".\n");
				}
				cat("\t\t\t", "The read records are ", data$traits[[i]]$envs[[j]]$obsread, ".\n");
				cat("\t\t\t", "The used records are ", data$traits[[i]]$envs[[j]]$obsused, ".\n");
				cat("\t\t\t", "The missing rate are ", data$traits[[i]]$envs[[j]]$missing.rate, ".\n");
				cat("\n");
			}
			cat(rep("-", 40), "\n", sep="");
			
		} #--- end stmt of for(j in 1:data$traits[[i]]$envs.number)---#
	#} #--- end stmt of if(data$isRestricted)---#
	
}
