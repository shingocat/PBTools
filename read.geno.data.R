###############################################################################
# TODO: Add comment
# 
#	[Arguments]
#		source: The CSV source containts genotypic data where the first column must be the individuals names and the first rows must be the header;
#		dp.code: The marker code of donor parent;
#		rp.code: The marker code of recurrent parent;
#		ht.code: The marker code of heterozygous;
#		na.code: The marker code of missing value;
#		BCn: The backcross generation;
#		Fn: The intercross generation;
#		replace.illegal.code: replace.illegal code to NA if true, else it will stop execution 
# Author: mqin
# Date: Oct 7, 2014
# FileName: read.geno.data.R
###############################################################################


read.geno.data <- function(
		source,
		dp.code,
		rp.code,
		ht.code, 
		na.code,
		BCn = 3,
		Fn = 2,
		replace.illegal.code = TRUE,
#		doRestrict = TRUE
		...
)
{
	if(missing(source))
		stop("\tError: The source argument could not be null!\n");
	if(missing(dp.code))
	{
		stop("\tError: The dp.code argument could not be null!\n");
	} else
	{
		dp.code <- as.character(dp.code);
	}
	if(missing(rp.code))
	{
		stop("\tError: The rp.code argument could not be null!\n");
	} else
	{
		rp.code <- as.character(rp.code);
	}
	if(missing(ht.code))
	{
		stop("\tError: The ht.code argument could not be null!\n");
	} else
	{
		ht.code <- as.character(ht.code);
	}
	if(missing(na.code))
	{
		stop("\tError: The na.code argument could not be null!\n");
	} else
	{
		na.code <- as.character(na.code);
	}
	if(missing(BCn))
	{
		warning("\tWarning: The BCn argument is not specified and it will be set to default value 3!\n");
		BCn <- 3;
	}
	if(!is.numeric(BCn))
		stop("\tError: The BCn argument should be of integer type! Default value is 3.\n");
	if( !((BCn %% 1) == 0))
		stop("\tError: The BCn argument should be of integer type! Default value is 3.\n");
	if(missing(Fn))
	{
		warning("\tWarning: The Fn argument is not specified and it will be set to default value 2!\n");
		Fn <- 2;
	}
	if(!is.numeric(Fn))
		stop("\tError: The Fn argument should be of integer type! Default value is 2.\n");
	if(!((Fn %% 1) == 0))
		stop("\tError: The Fn argument should be of integer type! Default value is 2.\n");	
	if(missing(replace.illegal.code))
		replace.illegal.code <- TRUE;
	if(!is.logical(replace.illegal.code))
		stop("\tErrow: The argument of replace.illegal.code should be of logical type!\n");
	UseMethod("read.geno.data");
}
#--- s3 method for data type is data.frame---#
read.geno.data.data.frame <- function(
		source,
		dp.code,
		rp.code,
		ht.code, 
		na.code,
		BCn = 3,
		Fn = 2,
		replace.illegal.code = TRUE,
		...
)
{
	result <- innerFunc.read.geno(source, dp.code, rp.code, ht.code, na.code, BCn, Fn);
	return(result);
}


#--- s3 method for data type is character---#
read.geno.data.character <- function(
		source,
		dp.code,
		rp.code,
		ht.code, 
		na.code,
		BCn = 3,
		Fn = 2,
		replace.illegal.code = TRUE,
#		doRestrict = TRUE
		...
)
{
	if(!file.exists(source))
		stop("\tError: The file does not exist! Please checking again!\n");

#	if(!is.logical(doRestrict))
#	{
#		stop("\tError: The doRestrict argument should be of logical type! Default value is TRUE.\n");
#	}
#	
	# read the file 
	data <- try(read.csv(file = source, header = T, na.strings=na.code), silent = FALSE);
	if(identical(class(data), "try-error"))
	{
		stop("\tError: There are some problems on in the file!\n");
	}
	
	result <- innerFunc.read.geno(data, dp.code, rp.code, ht.code, na.code, BCn, Fn);
	return(result);
}


innerFunc.read.geno <- function(
		source,
		dp.code,
		rp.code,
		ht.code,
		na.code,
		BCn,
		Fn,
		replace.illegal.code = TRUE
)
{
	data <- source;
	# checking whether there are illegal marker code in the file, only accept four type of marker code
	data.gen <- data[,1];
	geno.name <- colnames(data)[1];
	data.without.gen <- data[,-1];
	marker.names <- names(data.without.gen);
	marker.code <- c(dp.code, rp.code, ht.code, na.code);
	checking.result <- apply(data.without.gen,2, function(x) all(x %in% marker.code));
	if(!all(checking.result))
	{	
		if(replace.illegal.code)
		{
			data.without.gen[!checking.result] <- apply(data.without.gen[!checking.result], 2, function(x){
						x[!(x %in% marker.code)] <- na.code; 
						x;
					} );
		} else
		{
			stop(paste("\tError: There are illegal code in the marker : ", names(checking.result[!checking.result]), " !\n", sep = " " ));
	
		}
	}
	# built a list for storing all value of original
	result <- list();
	result$BCn <- BCn;
	result$Fn <- Fn;
	result$geno.name <- geno.name;
	result$original <- list();
	result$original$data <- data;
	result$original$dp.code <- dp.code;
	result$original$rp.code <- rp.code;
	result$original$ht.code <- ht.code;
	result$original$na.code <- na.code;
	
	# storing all value of processed data
	result$processed <- list();
	result$processed$dp.code <- 2;
	result$processed$rp.code <- 0;
	result$processed$ht.code <- 1;
	result$processed$na.code <- NA;
	# recoding original data
	library(car);
	recodes.str <- paste("'",dp.code , "'=2; '", rp.code,"'=0; '", ht.code, "'=1;'", na.code, "'=NA;",  sep = "");
	processed.geno.data <- apply(data.without.gen, 2, Recode, recodes.str);
	processed.geno.data <- as.data.frame(processed.geno.data);
	processed.data <- cbind(data.gen, processed.geno.data);
	colnames(processed.data) <- c(geno.name, marker.names);
	result$processed$data <- processed.data;
	
	# compute each marker levels
	marker.levels <- apply(processed.geno.data, 2, table);
	if(is.matrix(marker.levels))
		marker.levels <- as.data.frame(marker.levels)
	marker.levels <- unlist(lapply(marker.levels, length));
	result$processed$marker.levels <- marker.levels;
	
#	summary for each marker
#	marker.missing.rates <- apply(processed.geno.data, 2, function(x) length(which(is.na(x)))/length(x));
#	result$marker.missing.rate <- marker.missing.rates;	
	marker.geno.summary <- apply(processed.geno.data, 2, function(x)
			{
				dp.count <- length(which(x==2));
				rp.count <- length(which(x==0));
				ht.count <- length(which(x==1));
				na.count <- length(which(is.na(x)));
				total <- length(x);
				missing.rate <- na.count / total;
				result <- c(dp.count, rp.count, ht.count, na.count, total, missing.rate);
				names(result) <- c("dp.count", "rp.count", "ht.count", "na.count", "total", "missing.rate");
				return(result)
			});
	result$processed$marker.summary <- marker.geno.summary;
	
	# expected genotypic frequence 
	exp.genoFreq <- getExpectGenotypicFreq(BCn = BCn , Fn = Fn, BC.parent = 2);
	result$exp.genoFreq <- exp.genoFreq;
	
	# setting the result to be class GenotypicData
	class(result) <- "GenotypicData";
	detach(package:car, unload = TRUE);
	return(result);
}

print.GenotypicData <- function(
		data
)
{
	cat("Since modeified the class structure, needed to reimplemented\n");
}