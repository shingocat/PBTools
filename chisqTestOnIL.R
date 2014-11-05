###############################################################################
# TODO: Chisq Test on each marker on introgression line population which are 
# 		selected by phenotype!
# 
#
# [Arguments]
#	marker.matrix: Header is marker names and each row is genotypic value of each marker!.
#	BCn: The generation of backcross;
#	Fn:	The generation of selfing or intercross;
#	dp.code: The marker code of donor parent;
#	rp.code: The marker code of recurrent parent;
#	ht.code: The marker code of heterozygous;
#	na.code: The marker code of missing data;
# 	odd.code: The marker code of abnormal genotype which is not donor or recurrent genotype;
#	inc.ht: Logic value indicated whether heterozygous observation should be treated as donor observations.
#	pop.min: The minimum size of population, if the population size is less than this size, the routine will be stop.
#	do.ref: Logic value indicated whether marker expected frequency distribution depended on referenced marker matrix;
#	ref.matrix: referenced marker matrix should be used to compute marker frequency distribution when it is not null,
#				and frequecy distribution of marker if it is not in this referenced matrix, will be used theory frequency.
#				referenced marker matrix is same as marker.matrix format.
#	byCol: The marker is list on column else on row of marker.matrix;
#	simulate.p.value: a logical indicating whether to compute p-values by Monte Carlo simulation.
#	B: an integer specifying the number of replicates used in the Monte Carlo test.
#
# Author: mqin
# Date: Sep 25, 2014
# FileName: chisq.il.R
###############################################################################

chisqTestOnIL <- function(
		marker.matrix, 
		BCn, 
		Fn, 
		dp.code = "B", 
		rp.code = "A", 
		ht.code = "H", 
		na.code = NA, 
		odd.code = c("C", "D"),
		inc.ht = TRUE,
		pop.min = 20, 
#		do.ref = FALSE, 
		ref.matrix = NULL,  
		byCol = TRUE,
		simulate.p.value = TRUE,
		B = 1000
)
{
	if(is.null(marker.matrix))
	{
		stop("\tError: The argument of marker.matrix should not be null!\n");
	}
	if(!is.matrix(marker.matrix))
	{
		stop("\tError: The argument of marker.matrix should be matrix!\n");
	}
	if(is.null(BCn))
	{
		stop("\tError: The argument of BCn should not be null!\n");
	}
	if(is.null(Fn))
	{
		stop("\tError: The argument of Fn should not be null!\n");
	}
	if(!is.numeric(BCn))
	{
		stop("\tError: The argument of BCn should be numeric type!\n");
	}
	if(!is.numeric(Fn))
	{
		stop("\tError: The argument of Fn should be numeric type!\n");
	}
	if(length(BCn) != 1)
	{
		stop("\tError: The argument of BCn should be of length one!\n");
	}
	if(length(Fn) != 1)
	{
		stop("\tError: The argument of Fn should be of length one!\n");
	}
	if(BCn < 0)
	{
		stop("\tError: The argument of BCn should not be negative!\n");
	}
	if(Fn < 0)
	{
		stop("\tError: The argument of Fn should not be negative!\n");
	}
	if(is.null(dp.code))
	{
		stop("\tError: The argument of dp.code should not be null!\n");
	}
	if(!is.character(dp.code))
	{
		stop("\tError: The argument of dp.code should be characteric type!\n");
	}
	if(length(dp.code) != 1)
	{
		stop("\tError: The argument of dp.code should be of length one!\n");
	}
	if(is.null(rp.code))
	{
		stop("\tError: The argument of rp.code should not be null!\n");
	}
	if(!is.character(rp.code))
	{
		stop("\tError: The argument of rp.code should be characteric type!\n");
	}
	if(length(rp.code) != 1)
	{
		stop("\tError: The argument of rp.code should be of length one!\n");
	}
	
	if(!byCol)
	{
		marker.matrix <- t(marker.matrix);
	}
	markers.num <- ncol(marker.matrix);
	markers.names <- colnames(marker.matrix);
	if(is.null(markers.names))
		markers.names <- paste("M", 1:markers.num, sep="")
	lines.num <- nrow(marker.matrix);
	marker.matrix <- apply(marker.matrix,2,factor);	
	
	if(lines.num < pop.min)
	{
		stop(paste("\tError: There are minimum requriement of populaiton size is ", pop.min, " using chisq test!\n", sep=""));
	}
	
	# Cheacking whether there are unaccepted code in marker matrix
	marker.code <- c(dp.code, rp.code, ht.code, na.code, odd.code);
	checking.result <- apply(marker.matrix, 2, function(x) all(x %in% marker.code));
	if(!all(checking.result))
	{
		stop(paste("\tError: There are illegal code in the marker : ", names(checking.result[!checking.result]), "!\n", sep = ""));
	}
	
	do.ref <- FALSE;
	refFreq <- function(x, dp.code, ht.code, rp.code)
	{
		result <- c();
		x <- factor(x);
		x.table <- table(x, useNA = "no");
		rp.num <- x.table[as.character(rp.code)];
		dp.num <- x.table[as.character(dp.code)];
		ht.num <- x.table[as.character(ht.code)];
		if(is.na(rp.num))
			rp.num <- 0;
		if(is.na(dp.num))
			dp.num <- 0;
		if(is.na(ht.num))
			ht.num <- 0;
		sample.size <- sum(c(dp.num, ht.num, rp.num));
		freq <- c(dp.num, ht.num, rp.num)/sample.size;
		result <- freq;
		names(result) <- c("Donor", "Heterozygous", "Recurrent");
		return(result);
	}
	# If the referenced marker matrix is not null, computed the referenced frequence of each marker.
	if(!is.null(ref.matrix))
	{
		do.ref <- TRUE;
		ref.markers.names <- colnames(ref.matrix);
		# checking whether there are illegal code in referenced marker matrix
		checking.result <- apply(ref.matrix, 2, function(x) all(x %in% marker.code));
		if(!all(checking.result))
		{
			stop(paste("\tError: There are illegal code in the referenced matrix : ", names(checking.result[!checking.result]), "!\n", sep=""));
		}
		ref.markers.freq <- t(apply(ref.matrix, 2, refFreq, dp.code, ht.code, rp.code));
	}
	
# using apply method, it is faster!
#
#	innerFunc <- function(x, donor.marker.code, recurrent.marker.code)
#	{
#		result <- c();
#		x <- factor(x);
#		x.table <- table(x, useNA = "no");
#		sample.size <- sum(x.table);
#		
#		recurrent.genotype.num <- x.table[as.character(recurrent.marker.code)];
#		if(is.na(recurrent.genotype.num))
#		{
#			recurrent.genotype.num <- 0;
#		}
#		names(recurrent.genotype.num) <- c();
#
#		donor.genotype.num <- x.table[as.character(donor.marker.code)];
#		if(is.na(donor.genotype.num))
#		{
#			donor.genotype.num <- 0;
#		}
#		names(donor.genotype.num) <- c();
#		
#		# if there are other genotype levels, it will be add the other genotype levels individuals
#		# to the donor genotype individual when donor genotype is less than 5 observated individuals;
#		other.genotype.num <- x.table[-c(which(names(x.table) == as.character(donor.marker.code)), 
#						which(names(x.table) == as.character(recurrent.marker.code)))];
#		if(length(other.genotype.num) == 0)
#		{
#			other.genotype.num <- 0;
#		}
#		other.genotype.num <- sum(other.genotype.num);
#		names(other.genotype.num) <- c();
#		
#		if(donor.genotype.num < 3)
#		{
#			donor.genotype.num <- donor.genotype.num + other.genotype.num;
#		}
#		
#		if((length(x.table) == 1 && recurrent.genotype.num == sample.size) ||
#				((recurrent.genotype.num > donor.genotype.num) && (donor.genotype.num < 3)))
#		{
#			# TypeI denoted as the genotype of marker are always of recurrent genotype, and would not compute this kind of marker;
#			# treated as no effect marker.
#			result <- c("TypeI", donor.genotype.num + recurrent.genotype.num, donor.genotype.num, recurrent.genotype.num, NA, NA);
#			return(result);
#		}
#		
#		if((length(x.table) == 1 && donor.genotype.num == sample.size) ||
#			((donor.genotype.num > recurrent.genotype.num) && (recurrent.genotype.num < 3)))
#		{
#			# TypeII denoted as the genotype of marker are always of donor genotype or the num of donor marker genotype is large than
#			# the num of recurrent marker genotype which are less than 3 individuals;
#			# , and would not compute this kind of marker;
#			# treated as potential linkage marker.
#			result <- c("TypeII", sample.size, donor.genotype.num, recurrent.genotype.num, NA, NA);
#			return(result);
#		}
#		
#		expect.freq <- getExpectGenotypicFreq(BCn, Fn, 2);
#		expect.freq <- expect.freq[c(1,3)]/sum(expect.freq[c(1,3)]);
#		
#		if(floor(sample.size * expect.freq[1]) < 1)
#		{
#			# TypeIII denoted as the expected least genotype, i.e. donor genotype, are less than 1 on this sample size
#			# could not process the chisq testing!
#			warning("\tWarning: the minimum of sample size for the num of donor genotype are large than 1 are not reached, could not process chisq testing!\n");
#			result <- c("TypeIII",sample.size, donor.genotype.num, recurrent.genotype.num, NA, NA);
#			return(result);
#		}
#		
#		
#		data <- c(donor.genotype.num, recurrent.genotype.num);
#		outcomes <- chisq.test(data, p = expect.freq, simulate.p.value = T, B = 1000);
#		statistic <- round(outcomes$statistic, digits = 2);
#		p.value <- round(outcomes$p.value, digits = 5);
#		names(p.value) <- NA;
#		# TypeIV denoted as the normal case for chisq testing!
#		result <- c("TypeIV", sample.size, donor.genotype.num, recurrent.genotype.num,
#				statistic, p.value);
#		return(result);
#		
#	}
#	
#	temp <- apply(marker.matrix, 2, innerFunc, donor.marker.code = donor.marker.code, recurrent.marker.code = recurrent.marker.code);
#	temp <- as.data.frame(t(temp));
#	temp <- cbind(names.markers, temp);
#	rownames(temp) <- c();
#	colnames(temp) <- c("Marker","Type", "SampleSize", "DonorType", "RecurrentType", "Chisq", "P.Value");
#	#print(temp);
#	result <- list();
#	result$TypeI <- list(); # for all the genotype of marker are recurrent type or the num of donor genotype are less than 3;
#	result$TypeI$Markers <- temp[temp$Type == "TypeI", "Marker"];
#	result$TypeI$outcomes <- temp[temp$Type == "TypeI", c("Marker", "SampleSize", "DonorType", "RecurrentType")];
#	result$TypeI$counts <- length(result$TypeI$Markers);
#	result$TypeI$Messages <- c("For all the genotype of marker are recurrent type or the num of donor genotype are less than 3;");
#	result$TypeII <- list(); # for all the genotype of marker are donor type or the num of recurrent genotype are less than 3, potential linkage marker!;
#	result$TypeII$Markers <- temp[temp$Type == "TypeII", "Marker"];
#	result$TypeII$outcomes <- temp[temp$Type == "TypeII", c("Marker", "SampleSize", "DonorType", "RecurrentType")];
#	result$TypeII$counts <- length(result$TypeII$Markers);
#	result$TypeII$Messages <- c("For all the genotype of marker are donor type or the num of recurrent genotype are less than 3, potential linkage marekr!");
#	result$TypeIII <- list(); # for the minimum of sample size for the num of donor genotype are large than 1 are not reached;
#	result$TypeIII$Markers <- temp[temp$Type == "TypeIII", "Marker"];
#	result$TypeIII$outcomes <- temp[temp$Type == "TypeIII", c("Marker", "SampleSize", "DonorType", "RecurrentType")];
#	result$TypeIII$counts <- length(result$TypeIII$Markers);
#	result$TypeIII$Messages <- c("For the minimum of sample size for the num of donor genotype are large than 1 are not reached;");
#	result$TypeIV <- list(); # for the normal chisq case;
#	result$TypeIV$Markers <- temp[temp$Type == "TypeIV", "Marker"];
#	result$TypeIV$outcomes <- temp[temp$Type == "TypeIV", c("Marker", "SampleSize", "DonorType", "RecurrentType", "Chisq", "P.Value")];
#	result$TypeIV$counts <- length(result$TypeIV$Markers);
#	result$TypeIV$Messages <- c("For the normal chisq case;");
#	return(result);

# using for loop, it is more slow than using apply method!
#
	result <- c();
	temp <- c(); 
	for(i in 1:markers.num)
	{
		geno <- marker.matrix[,i];
		geno <- factor(geno);
		geno.table <- table(geno, useNA = "no");
		
		dp.num <- geno.table[as.character(dp.code)];	
		if(is.na(dp.num))
			dp.num <- 0;
		rp.num <- geno.table[as.character(rp.code)];
		if(is.na(rp.num))
			rp.num <- 0;
		ht.num <- geno.table[as.character(ht.code)];
		if(is.na(ht.num))
			ht.num <- 0;
		
		sample.size <- sum(c(dp.num, rp.num, ht.num));
		expect.freq <- getExpectGenotypicFreq(BCn, Fn, 2);
		
		
		if(dp.num == 0 && ht.num == 0)
		{
			# TypeI denoted as the genotypic value of this marker are always recurrent, and would not compute this kind of marker;
			# treated as no effect marker. In the other words, donor and heterozygous number is zero!	
			if(do.ref)
			{
				if(markers.names[i] %in% ref.markers.names)
				{
					freq <- ref.markers.freq[which(markers.names[i] == ref.markers.names),];
					temp <- c(markers.names[i], "TypeI", sample.size, dp.num, ht.num, rp.num, "N", paste(freq, collapse=":"), NA, NA);
				} else
				{
					temp <- c(markers.names[i], "TypeI", sample.size, dp.num, ht.num, rp.num, "Y", paste(expect.freq, collapse=":"), NA, NA);
				}
			}else
			{
				temp <- c(markers.names[i], "TypeI", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
			}
			names(temp) <- c();
			result <- rbind(result, temp);
			next;
		} else if(rp.num == 0 || (dp.num > rp.num && rp.num < 3))
		{
			# TypeII denoted as the genotype of marker are always of donor and heterozygous genotype or the num of donor marker genotype is large than
			# the num of recurrent marker genotype which are less than 3 individuals;
			# , and would not compute this kind of marker;
			# treated as potential linkage marker.
			if(do.ref)
			{
				if(markers.names[i] %in% ref.markers.names)
				{
					freq <- ref.markers.freq[which(markers.names[i] == ref.markers.names),];
					temp <- c(markers.names[i], "TypeII", sample.size, dp.num, ht.num, rp.num, "N", paste(freq, collapse=":"), NA, NA);
				} else
				{
					temp <- c(markers.names[i], "TypeII", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
				}
			}else
			{
				temp <- c(markers.names[i], "TypeII", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
			}
			names(temp) <- c();
			result <- rbind(result, temp);
			next;
		} else if(floor(sample.size * expect.freq[1]) < 1)
		{
			# TypeIII denoted as the expected least genotype, i.e. donor genotype, are less than 1 on this sample size
			# could not process the chisq testing!
			if(do.ref)
			{
				if(markers.names[i] %in% ref.markers.names)
				{
					freq <- ref.markers.freq[which(markers.names[i] == ref.markers.names),];
					temp <- c(markers.names[i], "TypeIII", sample.size, dp.num, ht.num, rp.num, "N", paste(freq, collapse=":"), NA, NA);
				} else
				{
					temp <- c(markers.names[i], "TypeIII", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
				}
			}else
			{
				temp <- c(markers.names[i], "TypeIII", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
			}
			names(temp) <- c();
			result <- rbind(result, temp);
			next;
		} else
		{
			# ThypeIV  denoted as the normal case for chisq testing!
			#
			# if the argument of inc.ht is TRUE, than the heterozygous will add to donor type;
			if(inc.ht)
			{
				data <- c(dp.num + ht.num, rp.num);
			} else
			{
				data <- c(dp.num, rp.num);
			}

			if(do.ref)
			{
				if(markers.names[i] %in% ref.markers.names)
				{
					freq <- ref.markers.freq[which(markers.names[i] == ref.markers.names), ];
					if(inc.ht)
					{
						prob <- c(freq[1] + freq[2], freq[3]) / sum(freq);
					} else
					{
						prob <- c(freq[1], freq[3]) / sum(freq[c(1,3)]);
					}
					
					# if one of the referenced marker freq is equal to zero, and then the prob is used expected.freq.
					if(any(prob == 0))
					{
						if(inc.ht)
						{
							prob <- c(expect.freq[1] + expect.freq[2], expect.freq[3]) / sum(expect.freq);
						} else
						{
							prob <- c(expect.freq[1], expect.freq[3]) / sum(expect.freq[c(1,3)]);
						}
						outcomes <- chisq.test(data, p = prob, simulate.p.value = simulate.p.value, B = B);
						statistic <- round(outcomes$statistic, digits = 2);
						p.value <- round(outcomes$p.value, digits = 4);
						names(p.value) <- c();
						temp <- c(markers.names[i], "TypeV", sample.size, dp.num, ht.num, rp.num, "Y", paste(expect.freq, collapse=":"), statistic, p.value);
						names(temp) <- c();
						result <- rbind(result, temp);
						next;
					}
					outcomes <- chisq.test(data, p = prob, simulate.p.value = simulate.p.value, B = B);
					statistic <- round(outcomes$statistic, digits = 2);
					p.value <- round(outcomes$p.value, digits = 4);
					names(p.value) <- c();
					temp <- c(markers.names[i], "TypeIV", sample.size, dp.num, ht.num, rp.num, "N", paste(freq, collapse=":"), statistic, p.value);
				} else
				{
					if(inc.ht)
					{
						prob <- c(expect.freq[1] + expect.freq[2], expect.freq[3]) / sum(expect.freq);
					} else
					{
						prob <- c(expect.freq[1], expect.freq[3]) / sum(expect.freq[c(1,3)]);
					}
					outcomes <- chisq.test(data, p = prob, simulate.p.value = simulate.p.value, B = B);
					statistic <- round(outcomes$statistic, digits = 2);
					p.value <- round(outcomes$p.value, digits = 4);
					names(p.value) <- c();
					temp <- c(markers.names[i], "TypeIV", sample.size, dp.num, ht.num, rp.num, "Y", paste(expect.freq, collapse=":"), statistic, p.value);
				}
			}else
			{
				if(inc.ht)
				{
					prob <- c(expect.freq[1] + expect.freq[2], expect.freq[3]) / sum(expect.freq);
				} else
				{
					prob <- c(expect.freq[1], expect.freq[3]) / sum(expect.freq[c(1,3)]);
				}
				outcomes <- chisq.test(data, p = expect.freq[c(1,3)]/sum(expect.freq[c(1,3)]), simulate.p.value = simulate.p.value, B = B);
				statistic <- round(outcomes$statistic, digits = 2);
				p.value <- round(outcomes$p.value, digits = 4);
				names(p.value) <- c();
				temp <- c(markers.names[i], "TypeIV", sample.size, dp.num, ht.num, rp.num, "Y", paste(expect.freq, collapse=":"), statistic, p.value);
			}
			names(temp) <- c();
			result <- rbind(result, temp);
			next;
		}
	}
	result <- as.data.frame(result);
	colnames(result) <- c("Marker", "Type", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp","ExptFreq","Freq(DP:HT:RP)","Chisq", "P.Value");
	rownames(result) <- c();
	
	
	outcomes <- list();
	outcomes$result <- result;
	outcomes$TypeI <- list(); # for all the genotype of marker are recurrent type or the num of donor genotype are less than 3;
	outcomes$TypeI$Markers <- result[result$Type == "TypeI", "Marker"];
	names(outcomes$TypeI$Markers) <- c();
	outcomes$TypeI$outcomes <- result[result$Type == "TypeI", c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp")];
	outcomes$TypeI$counts <- length(outcomes$TypeI$Markers);
	outcomes$TypeI$Messages <- c("For all the genotype of marker are recurrent type or the num of donor genotype are less than 3;");
	outcomes$TypeII <- list(); # for all the genotype of marker are donor type or the num of recurrent genotype are less than 3, potential linkage marker!;
	outcomes$TypeII$Markers <- result[result$Type == "TypeII", "Marker"];
	names(outcomes$TypeII$Markers) <- c();
	outcomes$TypeII$outcomes <- result[result$Type == "TypeII", c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp")];
	outcomes$TypeII$counts <- length(outcomes$TypeII$Markers);
	outcomes$TypeII$Messages <- c("For all the genotype of marker are donor type or the num of recurrent genotype are less than 3, potential linkage marekr!");
	outcomes$TypeIII <- list(); # for the minimum of sample size for the num of donor genotype are large than 1 are not reached;
	outcomes$TypeIII$Markers <- result[result$Type == "TypeIII", "Marker"];
	names(outcomes$TypeIII$Markers) <- c();
	outcomes$TypeIII$outcomes <- result[result$Type == "TypeIII", c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp")];
	outcomes$TypeIII$counts <- length(outcomes$TypeIII$Markers);
	outcomes$TypeIII$Messages <- c("For the minimum of sample size for the num of donor genotype are large than 1 are not reached;");
	outcomes$TypeIV <- list(); # for the normal chisq case;
	outcomes$TypeIV$Markers <- result[result$Type == "TypeIV", "Marker"];
	names(outcomes$TypeIV$Markers) <- c();
	outcomes$TypeIV$outcomes <- result[result$Type == "TypeIV",  c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp", "Chisq", "P.Value")];
	outcomes$TypeIV$counts <- length(outcomes$TypeIV$Markers);
	outcomes$TypeIV$Messages <- c("For the normal chisq case;");
	outcomes$TypeV <- list();
	outcomes$TypeV$Markers <- result[result$Type == "TypeV", "Marker"];
	names(outcomes$TypeV$Markers) <- c();
	outcomes$TypeV$outcomes <- result[result$Type == "TypeV", c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp", "Chisq", "P.Value")];
	outcomes$TypeV$counts <- length(outcomes$TypeV$Markers);
	outcomes$TypeV$Messages <- c("For genotypic frequence of donor and/or heterozygous is zero on the referenced marker, it will be compute chisq test using expected frequency!");
	class(outcomes) <- "ChisqOnIL";
	return(outcomes);
}

