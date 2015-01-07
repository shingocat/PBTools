###############################################################################
# TODO: Add comment
# 
# [Arguments]
# 	marker.geno: [GenotypicData] marker coding data.
#	missing.rate: [double: 0.80] the rate of missing value;
#	cor.rate: [double: 0.90] the rate of correlation value;
#	mono.reduced: [logical: TRUE] indicating whether deleted monomorphysm marker;
#	donor.minicount: [integer: 5] The minimum count number of donor number;
# Author: mqin
# Date: Oct 13, 2014
# FileName: filter.geno.data.R
###############################################################################


restrict.geno.data <- function(
		geno.data,
		missing.rate = 0.80,
		cor.rate = 0.90,
		mono.reduced = TRUE,
		donor.minicount = 5
)
{
	if(is.null(geno.data$processed))
	{
		stop("\tError: The geno.data argument is not the instance of class GenotypicData!\n");
	}
	if(!is.double(missing.rate))
	{
		stop("\tError: The missing.rate argument must be double type!\n");
	}
	if(missing.rate < 0 | missing.rate > 1)
	{
		stop("\tError: The missing.rate argument must be of value between 0 and 1!\n");
	}
	if(!is.double(cor.rate))
	{
		stop("\tError: The cor.rate argument must be double type!\n");
	}
	if(cor.rate < 0 | cor.rate > 1)
	{
		stop("\tError: The cor.rate argument must be of value between 0 and 1!\n");
	}
	if(!is.logical(mono.reduced))
	{
		stop("\tError: The mono.reduced argument must be logical type!\n");
	}
	if(!((donor.minicount %% 1) == 0))
	{
		stop("\tError: The donor.minicount argument must be integer type!\n");
	}
	if(donor.minicount < 0)
	{
		stop("\tError: The donor.minicount argument must be larger than zero!\n");
	}
	UseMethod("restrict.geno.data");
}
restrict.geno.data.GenotypicData <- function(
		geno.data,
		missing.rate = 0.80,
		cor.rate = 0.90,
		mono.reduced = TRUE,
		donor.minicount = 5
)
{
	
	geno.data$restricted <- list();
	
	geno.data$restricted$missing.rate <- missing.rate;
	geno.data$restricted$cor.rate <- cor.rate;
	geno.data$restricted$mono.reduced <- mono.reduced;
	geno.data$restricted$donor.minicount <- donor.minicount;
	
	processed.data <- as.data.frame(geno.data$processed$data);
	# restrict missing rate large than 0.80
	marker.match.missing.rate <- which(geno.data$processed$marker.summary[6,] > missing.rate);
	marker.match.missing.rate.index <- which(names(processed.data) %in% names(marker.match.missing.rate));
	if(length(marker.match.missing.rate.index) != 0)
	{
		geno.data$restricted$missing.rate.marker <- names(marker.match.missing.rate);
		processed.data <- processed.data[ , -marker.match.missing.rate.index];
	}
	# restrict monomorphysm marker
	marker.match.mono <- which(geno.data$processed$marker.levels <= 1);
	marker.match.mono.index <- which(names(processed.data) %in% names(marker.match.mono));
	if(length(marker.match.mono.index) != 0)
	{
		geno.data$restricted$mono.reduced.marker <- names(marker.match.mono);
		processed.data <- processed.data[ , -marker.match.mono.index];
	}
	# restric minimum donor count
	marker.match.mini <- which(geno.data$processed$marker.summary[1,] <= donor.minicount);
	marker.match.mini.index <- which(names(processed.data) %in% names(marker.match.mini));
	if(length(marker.match.mini.index) != 0 )
	{
		geno.data$restricted$donor.minicount.marker <- names(marker.match.mini);
		processed.data <- processed.data[ , - marker.match.mini.index];
	}
	# restric correlation rate, method by Dr.Wang, modify by me
	marker.names <- colnames(processed.data)[-1];
	processed.data[ , -1] <- apply(processed.data[ , -1], 2, as.character);
	processed.data[ , -1] <- apply(processed.data[ ,-1], 2, as.numeric);
	test <- cor(processed.data[,-1],, use="pairwise.complete.obs");
	diag(test) = 0;
	iden = apply(test, 1, function(x) which(abs(x) > cor.rate)); # considering negative correlation
	
	if(is.list(iden))
	{
		lens = unlist(lapply(iden, length));
	}
	else if(is.vector(iden))
	{
		iden = as.list(iden);
		iden = lapply(iden, function(x){names(x) <- marker.names[x]; return(x)});
		lens = unlist(lapply(iden, length));
	}
	else
	{
		iden = as.data.frame(iden);
		iden = as.list(iden);
		iden = lapply(iden, function(x){names(x) <- marker.names[x]; return(x)});
		lens = unlist(lapply(iden, length));
	}
		
	# remove high_correlated markers step by step
	idx_max = which.max(lens);
	dis = NULL;
	while(length(iden))
	{
		if(lens[idx_max] > 1)
		{
			dis = c(dis, idx_max);
			for(i in iden[[idx_max]])
			{
				iden[[i]] <- setdiff(iden[[i]], idx_max);
			}
			iden[[idx_max]] <- setdiff(iden[[idx_max]],iden[[idx_max]]);
			lens <- unlist(lapply(iden,length));
			idx_max <- which.max(lens);
		}
		else{
			rest=which(lens==1);
			pairs=data.frame(marker_1 = rest, marker_2 = 0,stringsAsFactors=FALSE)
			for(i in 1:length(rest)){
				marker_2 <-iden[[rest[i]]];
				pairs[i,2] <- ifelse(marker_2 %in% pairs[1:(i-1),1],NA, marker_2)
			}
			pairs=na.omit(pairs)
			for(i in 1:nrow(pairs)){
				#deletion more NA marker
				inf1 <- geno.data$processed$marker.summary[c(1,4),marker.names[pairs[i,1]]];
				inf2 <- geno.data$processed$marker.summary[c(1,4),marker.names[pairs[i,2]]];
#				inf1=round(report[report$SNPname==snps[pairs[i,1]],4:5],digits=3)
#				inf2=round(report[report$SNPname==snps[pairs[i,2]],4:5],digits=3)
				dis1=ifelse(inf1[2]>inf2[2],pairs[i,1],
						ifelse(inf1[1] > inf2[1],pairs[i,2], pairs[i,1]));
				names(dis1) <- NULL;
				dis=c(dis,dis1)
			}
			iden=NULL
		}
	}
	marker.match.cor.index <- which(names(processed.data) %in% marker.names[dis]);
	if(length(marker.match.cor.index) != 0 )
	{
		processed.data <- processed.data[ , - marker.match.cor.index];
		geno.data$restricted$cor.marker <- marker.names[dis];
	}
	#--- adding imputation after restricted data ---#
	processed.data[,-1] <- imputation.GenotypicData(processed.data[,-1]);
	geno.data$restricted$data <- processed.data;
	geno.data$isRestricted <- TRUE;
	return(geno.data);
}
restrict.geno.data.default <- function(
		geno.data,
		missing.rate = 0.80,
		cor.rate = 0.90,
		mono.reduced = TRUE,
		donor.minicount = 5
)
{
	stop("\tError: Only accepted the read.geno.data return object!\n");
}

imputation.GenotypicData <- function(data, method = c("MajorGene"))
{
	method <- match.arg(method);
	if(missing(data))
		stop("\tError: The argument of data could not be null!\n");
	outcomes <- c();
	if(method == "MajorGene")
	{
		outcomes <- apply(data, 2, replace.by.major.gene);
	}
	return(outcomes);
}

replace.by.major.gene <- function(x)
{
	n <- length(x);
	x <- factor(x);
	y <- table(x, useNA = "always");
	prob <- y / n;
	major.gene <- names(which(max(prob) == prob));
	#--- if there are more than one genotype fitting for our standar---#
	#--- use the first one at now, but it is not good way ---#
	#--- modified this later if I have time---#
	if(length(major.gene) > 1)
		major.gene <- major.gene[1]
	x[is.na(x)] <- major.gene;
	return(x);
}
