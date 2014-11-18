###############################################################################
# TODO: reading phenotypic data from file and formating data to a list for analyzing later;
# 
# ARGUMENT:
# file - the file of phenotypic data, should be of csv or txt format with header;
# pop.type - a character string descrided population type, only accepted "SSSL", "PY" and "IL" value;
# gen.num - the number of gene for pyramided line design, it must be specified when population type is pyramieded lines
#			and only accepted 2, 3, and 4 genes by far.
# resp.var - a vector of strings; variable names of the response variables; 
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
# FileName: read.pheno.data.R
###############################################################################


read.pheno.data <- function(
		file,
		pop.type,
		gen.num,
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
		gen.num,
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
	if(missing(file))
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
	if(missing(pop.type))
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
	#-- checking the argument of gen.num when the pop.type is PL --#
	if(identical(pop.type, "PL"))
	{
		if(missing(gen.num))
			stop("\tError: The argument of gen.num could not be null when pop.type is \"PL\" and only accept one of integer 2, 3 and 4.\n");
		if(!is.numeric(gen.num))
			stop("\tError: The argument of gen.num only accepted one of integer 2, 3 and 4.\n");
		if(length(gen.num) != 1)
			stop("\tError: The argument of gen.num only accepted one of integer 2, 3 and 4.\n");
	}
	
	if(missing(exptl.design))
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
	if(missing(resp.var))
	{
		stop("\tError: The argument of resp.var could not be null!\n");
	}
	if(missing(geno))
	{
		stop("\tError: The argument of geno could not be null!\n");
	}
	if(length(geno) != 1)
	{
		stop("\tError: The argument of geno should be of length of 1!\n");
	}
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
	
	# reading data from file
	data <-  try(read.csv(file = file, header = T, na.strings=na.code), silent = TRUE);
	if(identical(class(data), "try-error"))
	{
		stop("\tError: There are some problems on in the file!\n");
	}
	
	# checking whether all design factor names exists 
	# if true, make all the design factor to be factor
	col.names <- colnames(data);
	design.factor <- c();
	if(!(geno %in% col.names))
	{
		stop("\tError: The argument of geno does not match any column in the file!\n");
	} else
	{
		geno <- as.character(geno);
		data[, geno] <- factor(trimStrings(data[, geno]));
		design.factor <- c(design.factor, geno);
	}
	# --- checking geno factor should have correct coded labels and levels --- #
	if(identical(pop.type, "PL"))
	{
		genoLevels <- levels(data[,match(geno, names(data))]);
		if(length(genoLevels) <= 1)
		{
			stop("\tError: The levels of geno argument cannot be less than two levels.\n");
		}
		
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
			design.factor <- c(design.factor, block);
		}
	} else if(exptl.design == "AugLS")
	{
		if(!(row %in% col.names))
		{
			stop("\tError: The argument of row does not match any column in the file!\n");
		} else
		{
			row <- as.character(row);
			data[, row] <- factor(trimStrings(data[,row]));
			design.factor <- c(design.factor, row);
		}
		if(!(column %in% col.names))
		{
			stop("\tError: The argument of column does not match any column in the file!\n");
		} else
		{
			column <- as.character(column);
			data[, column] <- factor(trimStrings(data[,column]));
			design.factor <- c(design.factor, column);
		}
	} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
	{
		if(!(block %in% col.names))
		{
			stop("\tError: The argument of block does not match any column in the file!\n");
		} else
		{
			block <- as.character(block);
			data[ , block] <- factor(trimStrings(data[ , block]));
			design.factor <- c(design.factor, block);
		}
		if(!(rep %in% col.names))
		{
			stop("\tError: The argument of rep does not match any column in the file!\n");
		} else
		{
			rep <- as.character(rep);
			data[ , rep] <- factor(trimStrings(data[ , rep]));
			design.factor <- c(design.factor, rep);
		}
	} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
	{
		if(!(row %in% col.names))
		{
			stop("\tError: The argument of row does not match any column in the file!\n");
		} else
		{
			row <- as.charachter(row);
			data[,row] <- factor(trimStrings(data[ , row]));
			design.factor <- c(design.factor, row);
		}
		if(!(column %in% col.names))
		{
			stop("\tError: The argument of column does not match any column in the file!\n");
		} else
		{
			column <- as.character(column);
			data[ , column] <- factor(trimStrings(data[ , column]));
			design.factor <- c(design.factor, column);
		}
		if(!(rep %in% col.names))
		{
			stop("\tError: The argument of rep does not match any column in the file!\n");
		} else
		{
			rep <- as.character(rep);
			data[ , rep] <- factor(trimStrings(data[ , rep]));
			design.factor <- c(design.factor, rep);
		}
	}
	
	isEnvNull = FALSE;
	
	if(is.null(env))
	{
		env = "EnvLevel";
		data = cbind(data, EnvLevel = 1);
		data[[env]] <- factor(data[[env]]);
		isEnvNull = TRUE;
		design.factor <- c(design.factor, env);
	} else
	{
		if(length(env) != 1)
			stop("\tError: The argument of env should of length of 1!\n");
		if(!(env %in% col.names))
			stop("\tError: The argument of env does not match any column in the file!\n");
		isEnvNULL = FALSE;
		env <- as.character(env);
		data[ , env] <- factor(trimStrings(data[ , env]));
		design.factor <- c(design.factor, env);
	}
	
	## checking all the design factor not include any missing value
	for( i in 1:length(design.factor))
	{
		if(na.code %in% data[ , design.factor[i]] || NA %in% data[ , design.factor[i]])
		{
			sotp(paste("\tError: The design factor ", design.factor[i], " could not contain any missing value!\n", sep = ""));
		}
	}
	
	## format the na.code of response variable to R default code "NA"
	if(!(identical(na.code, NA) || identical(na.code, "NA")))
	{
		for(i in 1:length(resp.var))
		{
			data[ , resp.var[i]][which(data[ , resp.var[i]] == na.code)] <- NA;
		}
	}
	
	## building and reformating to a list structer 
	outcomes <- list();
	outcomes$na.code <- na.code;
	outcomes$pop.type <- pop.type;
	if(identical(pop.type, "PL"))
		outcomes$gen.num <- gen.num;
	outcomes$resp.var <- resp.var;
	outcomes$raw.data <- data;
#	outcomes$design <- list();
#	outcomes$design$exptl.design <- exptl.design;
#	outcomes$design$geno <- geno;
#	outcomes$design$geno.levels <- levels(data[ , geno]);
#	outcomes$design$block <- block;
#	outcomes$design$block.levels <- levels(data[ , block]);
#	outcomes$design$row <- row;
#	outcomes$design$row.levels <- levels(data[ , row]);
#	outcomes$design$column <- column;
#	outcomes$design$column.levels <- levels(data[ , column]);
#	outcomes$design$rep <- rep;
#	outcomes$design$rep.levels <- levels(data[ , rep]);
#	outcomes$design$env <- env;
#	outcomes$design$env.levels <- levels(data[ , env]);
	## each trait build one list
	outcomes$traits <- list();
	for(i in 1:length(resp.var))
	{
		outcomes$traits[[i]] <- list();
		outcomes$traits[[i]]$name <- resp.var[i];
		outcomes$traits[[i]]$sites <- list();
		if(isEnvNull)
		{	
			outcomes$single.site <- TRUE;
			outcomes$traits[[i]]$sites[[1]] <- list();
			outcomes$traits[[i]]$sites[[1]]$name <- levels(data[[env]]);
			outcomes$traits[[i]]$sites[[1]]$obsread <- nrow(data);
			outcomes$traits[[i]]$sites[[1]]$obsused <- nrow(subset(data, subset= (is.na(data[ , resp.var[i]]) == FALSE)));
			outcomes$traits[[i]]$sites[[1]]$missing.rate <- 1 - outcomes$traits[[i]]$sites[[1]]$obsused / outcomes$traits[[i]]$sites[[1]]$obsread;
			outcomes$traits[[i]]$sites[[1]]$missing.rate.by.geno <- aggregate(data[[resp.var[i]]], by = list(data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
			outcomes$traits[[i]]$sites[[1]]$data <- data;
			outcomes$traits[[i]]$sites[[1]]$design <- list();
			outcomes$traits[[i]]$sites[[1]]$design$exptl.design <- exptl.design;
			outcomes$traits[[i]]$sites[[1]]$design$geno <- geno;
			outcomes$traits[[i]]$sites[[1]]$design$geno.levels <- levels(data[ , geno]);
			outcomes$traits[[i]]$sites[[1]]$design$block <- block;
			outcomes$traits[[i]]$sites[[1]]$design$block.levels <- levels(data[ , block]);
			outcomes$traits[[i]]$sites[[1]]$design$row <- row;
			outcomes$traits[[i]]$sites[[1]]$design$row.levels <- levels(data[ , row]);
			outcomes$traits[[i]]$sites[[1]]$design$column <- column;
			outcomes$traits[[i]]$sites[[1]]$design$column.levels <- levels(data[ , column]);
			outcomes$traits[[i]]$sites[[1]]$design$rep <- rep;
			outcomes$traits[[i]]$sites[[1]]$design$env <- env;
			outcomes$traits[[i]]$sites[[1]]$design$env.levels <- levels(data[ , env]);
		}else
		{
			outcomes$single.site <- FALSE;
			env.levels <- levels(data[[env]]);
			for(j in 1:length(env.levels))
			{
				# --- create temp.data with one environment level only --- #
				temp.data <- subset(data, subset = (data[[env]] == env.levels[j]));
				outcomes$traits[[i]]$sites[[j]] <- list();
				outcomes$traits[[i]]$sites[[j]]$name <- env.levels[j];
				outcomes$traits[[i]]$sites[[j]]$obsread <- nrow(temp.data);
				outcomes$traits[[i]]$sites[[j]]$obsused <- nrow(subset(temp.data, subset= (is.na(temp.data[ , resp.var[i]]) == FALSE)));
				outcomes$traits[[i]]$sites[[j]]$missing.rate <- 1 - outcomes$traits[[i]]$sites[[j]]$obsused / outcomes$traits[[i]]$sites[[j]]$obsread;
				outcomes$traits[[i]]$sites[[j]]$missing.rate.by.geno <- aggregate(temp.data[[resp.var[i]]], by = list(temp.data[[geno]]), FUN = function(x){sum(is.na(x)) / length(x)});
				outcomes$traits[[i]]$sites[[j]]$data <- temp.data;
				outcomes$traits[[i]]$sites[[j]]$design <- list();
				outcomes$traits[[i]]$sites[[j]]$design$exptl.design <- exptl.design;
				outcomes$traits[[i]]$sites[[j]]$design$geno <- geno;
				outcomes$traits[[i]]$sites[[j]]$design$geno.levels <- levels(temp.data[ , geno]);
				outcomes$traits[[i]]$sites[[j]]$design$block <- block;
				outcomes$traits[[i]]$sites[[j]]$design$block.levels <- levels(temp.data[ , block]);
				outcomes$traits[[i]]$sites[[j]]$design$row <- row;
				outcomes$traits[[i]]$sites[[j]]$design$row.levels <- levels(temp.data[ , row]);
				outcomes$traits[[i]]$sites[[j]]$design$column <- column;
				outcomes$traits[[i]]$sites[[j]]$design$column.levels <- levels(temp.data[ , column]);
				outcomes$traits[[i]]$sites[[j]]$design$rep <- rep;
				outcomes$traits[[i]]$sites[[j]]$design$env <- env;
				outcomes$traits[[i]]$sites[[j]]$design$env.levels <- levels(temp.data[ , env]);
			}
		} ## end of statement if(isEnvNull)
	} ## end of statement of for(i in 1:length(resp.var))
	class(outcomes) <- "PhenotypicData";
	return(outcomes);
}

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
	cat("Phenotypic data details:\n");
	cat("\tThe population type is : ", if(data$pop.type == "IL") "Introgression Lines" else if (data$pop.type == "PY") "Pyramided Lines" else "Single Segment Substitution Lines", ".\n");
	cat("\tThe number of response variable is:", length(data$resp.var), ".\n");
	cat("\tIs single site experiment :", if(data$single.site) "Yes" else "NO", ".\n");
	cat(rep("-.-", 30));
	cat("\n");
	for(i in 1:length(data$resp.var))
	{
		cat("\t\t", "The response variable is: ", data$resp.var[i], ".\n");
		for(j in 1:length(data$traits[[i]]$sites))
		{
			cat("\t\t\t", "The experiment design is ", data$traits[[i]]$sites[[j]]$design$exptl.design, ".\n");
			if(level == 2)
			{
				exptl.design <- data$traits[[i]]$sites[[j]]$design$exptl.design;
				if(exptl.design == "RCB" || exptl.design == "AugRCB")
				{
					cat("\t\t\t", "The block variable is ", data$traits[[i]]$sites[[j]]$design$block, ".\n");
					cat("\t\t\t", "The levels of block variable are ", data$traits[[i]]$sites[[j]]$design$block.levels, ".\n");
				} else if(exptl.design == "AugLS")
				{
					cat("\t\t\t", "The row variable is ", data$traits[[i]]$sites[[j]]$design$row, ".\n");
					cat("\t\t\t", "The levels of row variable are ", data$traits[[i]]$sites[[j]]$design$row.levels, ".\n");
					cat("\t\t\t", "The column variable is ", data$traits[[i]]$sites[[j]]$design$column, ".\n");
					cat("\t\t\t", "The levels of column variable are ", data$traits[[i]]$sites[[j]]$design$column.levels, ".\n");
				} else if(exptl.design == "Alpha" || exptl.design == "LatinAlpha")
				{
					cat("\t\t\t", "The block variable is ", data$traits[[i]]$sites[[j]]$design$block, ".\n");
					cat("\t\t\t", "The levels of block variable are ", data$traits[[i]]$sites[[j]]$design$block.levels, ".\n");
					cat("\t\t\t", "The rep variable is ", data$traits[[i]]$sites[[j]]$design$rep, ".\n");
					cat("\t\t\t", "The levels of rep variable are ", data$traits[[i]]$sites[[j]]$design$rep.levels, ".\n");
				} else if(exptl.design == "RowCol" || exptl.design == "LatinRowCol")
				{
					cat("\t\t\t", "The row variable is ", data$traits[[i]]$sites[[j]]$design$row, ".\n");
					cat("\t\t\t", "The levels of row variable are ", data$traits[[i]]$sites[[j]]$design$row.levels, ".\n");
					cat("\t\t\t", "The column variable is ", data$traits[[i]]$sites[[j]]$design$column, ".\n");
					cat("\t\t\t", "The levels of column variable are ", data$traits[[i]]$sites[[j]]$design$column.levels, ".\n");
					cat("\t\t\t", "The rep variable is ", data$traits[[i]]$sites[[j]]$design$rep, ".\n");
					cat("\t\t\t", "The levels of rep variable are ", data$traits[[i]]$sites[[j]]$design$rep.levels, ".\n");
				}
			}
			cat("\t\t\t", "The environment is ", data$traits[[i]]$sites[[j]]$name, ".\n");
			cat("\t\t\t", "The Genotypic variable is ", data$traits[[i]]$sites[[j]]$design$geno, ".\n");
			if(level == 2)
			{
				cat("\t\t\t", "The level of Genotypic variable is ", data$traits[[i]]$sites[[j]]$design$geno.levels, ".\n");
			}
			cat("\t\t\t", "The read records are ", data$traits[[i]]$sites[[j]]$obsread, ".\n");
			cat("\t\t\t", "The used records are ", data$traits[[i]]$sites[[j]]$obsused, ".\n");
			cat("\t\t\t", "The missing rate are ", data$traits[[i]]$sites[[j]]$missing.rate, ".\n");
			cat("\n");
		}
		cat(rep("-.-", 30), "\n");
		
	}
	
}
