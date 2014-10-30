###############################################################################
# TODO: Add comment
# 
#
# [Arguments]
#	marker.matrix:
#	BCn: The generation of backcross;
#	Fn:	The generation of selfing or intercross;
#	donor.marker.code: The marker code of donor parent;
#	recurrent.marker.code: The marker code of recurrent.marker.code;
#	byCol: The marker is list on column else on row of marker.matrix;
#
# Author: mqin
# Date: Sep 25, 2014
# FileName: chisq.il.R
###############################################################################

chisqTestOnIL <- function(
		marker.matrix, BCn, Fn, donor.marker.code = 3, recurrent.marker.code = 1, byCol = TRUE)
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
		stop("\tError: The argument of BCn should be numeric!\n");
	}
	if(!is.numeric(Fn))
	{
		stop("\tError: The argument of Fn should be numeric!\n");
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
	
	if(!byCol)
	{
		marker.matrix <- t(marker.matrix);
	}
	num.markers <- ncol(marker.matrix);
	names.markers <- colnames(marker.matrix);
	num.lines <- nrow(marker.matrix);
	marker.matrix <- apply(marker.matrix,2,factor);
	
	if(num.lines < 20)
	{
		stop("\tError: There are at least 20 records requriements of using chisq test!\n");
	}
	
	innerFunc <- function(x, donor.marker.code, recurrent.marker.code)
	{
		result <- c();
		x <- factor(x);
		x.table <- table(x, useNA = "no");
		sample.size <- sum(x.table);
		
		recurrent.genotype.num <- x.table[as.character(recurrent.marker.code)];
		if(is.na(recurrent.genotype.num))
		{
			recurrent.genotype.num <- 0;
		}
		names(recurrent.genotype.num) <- c();

		donor.genotype.num <- x.table[as.character(donor.marker.code)];
		if(is.na(donor.genotype.num))
		{
			donor.genotype.num <- 0;
		}
		names(donor.genotype.num) <- c();
		
		# if there are other genotype levels, it will be add the other genotype levels individuals
		# to the donor genotype individual when donor genotype is less than 5 observated individuals;
		other.genotype.num <- x.table[-c(which(names(x.table) == as.character(donor.marker.code)), 
						which(names(x.table) == as.character(recurrent.marker.code)))];
		if(length(other.genotype.num) == 0)
		{
			other.genotype.num <- 0;
		}
		other.genotype.num <- sum(other.genotype.num);
		names(other.genotype.num) <- c();
		
		if(donor.genotype.num < 3)
		{
			donor.genotype.num <- donor.genotype.num + other.genotype.num;
		}
		
		if((length(x.table) == 1 && recurrent.genotype.num == sample.size) ||
				((recurrent.genotype.num > donor.genotype.num) && (donor.genotype.num < 3)))
		{
			# TypeI denoted as the genotype of marker are always of recurrent genotype, and would not compute this kind of marker;
			# treated as no effect marker.
			result <- c("TypeI", donor.genotype.num + recurrent.genotype.num, donor.genotype.num, recurrent.genotype.num, NA, NA);
			return(result);
		}
		
		if((length(x.table) == 1 && donor.genotype.num == sample.size) ||
			((donor.genotype.num > recurrent.genotype.num) && (recurrent.genotype.num < 3)))
		{
			# TypeII denoted as the genotype of marker are always of donor genotype or the num of donor marker genotype is large than
			# the num of recurrent marker genotype which are less than 3 individuals;
			# , and would not compute this kind of marker;
			# treated as potential linkage marker.
			result <- c("TypeII", sample.size, donor.genotype.num, recurrent.genotype.num, NA, NA);
			return(result);
		}
		
		expect.freq <- getExpectGenotypicFreq(BCn, Fn, 2);
		expect.freq <- expect.freq[c(1,3)]/sum(expect.freq[c(1,3)]);
		
		if(floor(sample.size * expect.freq[1]) < 1)
		{
			# TypeIII denoted as the expected least genotype, i.e. donor genotype, are less than 1 on this sample size
			# could not process the chisq testing!
			warning("\tWarning: the minimum of sample size for the num of donor genotype are large than 1 are not reached, could not process chisq testing!\n");
			result <- c("TypeIII",sample.size, donor.genotype.num, recurrent.genotype.num, NA, NA);
			return(result);
		}
		
		
		data <- c(donor.genotype.num, recurrent.genotype.num);
		outcomes <- chisq.test(data, p = expect.freq, simulate.p.value = T, B = 1000);
		statistic <- round(outcomes$statistic, digits = 2);
		p.value <- round(outcomes$p.value, digits = 5);
		names(p.value) <- NA;
		# TypeIV denoted as the normal case for chisq testing!
		result <- c("TypeIV", sample.size, donor.genotype.num, recurrent.genotype.num,
				statistic, p.value);
		return(result);
		
	}
	
	temp <- apply(marker.matrix, 2, innerFunc, donor.marker.code = donor.marker.code, recurrent.marker.code = recurrent.marker.code);
	temp <- as.data.frame(t(temp));
	temp <- cbind(names.markers, temp);
	rownames(temp) <- c();
	colnames(temp) <- c("Marker","Type", "SampleSize", "DonorType", "RecurrentType", "Chisq", "P.Value");
	#print(temp);
	result <- list();
	result$TypeI <- list(); # for all the genotype of marker are recurrent type or the num of donor genotype are less than 3;
	result$TypeI$Markers <- temp[temp$Type == "TypeI", "Marker"];
	result$TypeI$outcomes <- temp[temp$Type == "TypeI", c("Marker", "SampleSize", "DonorType", "RecurrentType")];
	result$TypeI$counts <- length(result$TypeI$Markers);
	result$TypeI$Messages <- c("For all the genotype of marker are recurrent type or the num of donor genotype are less than 3;");
	result$TypeII <- list(); # for all the genotype of marker are donor type or the num of recurrent genotype are less than 3, potential linkage marker!;
	result$TypeII$Markers <- temp[temp$Type == "TypeII", "Marker"];
	result$TypeII$outcomes <- temp[temp$Type == "TypeII", c("Marker", "SampleSize", "DonorType", "RecurrentType")];
	result$TypeII$counts <- length(result$TypeII$Markers);
	result$TypeII$Messages <- c("For all the genotype of marker are donor type or the num of recurrent genotype are less than 3, potential linkage marekr!");
	result$TypeIII <- list(); # for the minimum of sample size for the num of donor genotype are large than 1 are not reached;
	result$TypeIII$Markers <- temp[temp$Type == "TypeIII", "Marker"];
	result$TypeIII$outcomes <- temp[temp$Type == "TypeIII", c("Marker", "SampleSize", "DonorType", "RecurrentType")];
	result$TypeIII$counts <- length(result$TypeIII$Markers);
	result$TypeIII$Messages <- c("For the minimum of sample size for the num of donor genotype are large than 1 are not reached;");
	result$TypeIV <- list(); # for the normal chisq case;
	result$TypeIV$Markers <- temp[temp$Type == "TypeIV", "Marker"];
	result$TypeIV$outcomes <- temp[temp$Type == "TypeIV", c("Marker", "SampleSize", "DonorType", "RecurrentType", "Chisq", "P.Value")];
	result$TypeIV$counts <- length(result$TypeIV$Markers);
	result$TypeIV$Messages <- c("For the normal chisq case;");
	return(result);

#	result <- c();
#	for(i in 1:num.markers)
#	{
#		marker.genotype <- marker.matrix[,i];
#		marker.genotype <- factor(marker.genotype);
#		marker.genotype.table <- table(marker.genotype);
#		donor.genotype.num <- marker.genotype.table[as.character(donor.marker.code)];	
#		if(is.na(donor.genotype.num))
#		{
#			donor.genotype.num <- 0;
#			next;
#		}
#		recurrent.genotype.num <- marker.genotype.table[as.character(recurrent.marker.code)];
#		if(is.na(recurrent.genotype.num))
#		{
#			recurrent.genotype.num <- 0;
#		}
#		expect.freq <- getExpectGenotypicFreq(BCn, Fn, 2);
#		expect.freq <- expect.freq[c(1,3)]/sum(expect.freq[c(1,3)]);
#		data <- c(donor.genotype.num, recurrent.genotype.num);
#		outcomes <- chisq.test(data, p = expect.freq, simulate.p.value = T, B = 1000);
#		expected.donor.num <- outcomes$expected[as.character(donor.marker.code)];
#		if(is.na(expected.donor.num))
#		{
#			expected.donor.num <- 0;
#		}
#		expected.donor.num <- round(expected.donor.num);
#		expected.recurrent.num <- outcomes$expected[as.character(recurrent.marker.code)];
#		if(is.na(expected.recurrent.num))
#		{
#			expected.recurrent.num <- 0;
#		}
#		expected.recurrent.num <- round(expected.recurrent.num);
#		statistic <- round(outcomes$statistic, digits = 2);
#		p.value <- round(outcomes$p.value, digits = 5);
#		temp <- c(names.markers[i], donor.genotype.num, recurrent.genotype.num, expected.donor.num,
#				expected.recurrent.num, statistic, p.value);
#		names(temp) <- c();
#		result <- rbind(result, temp);
#		rownames(result) <- c();
#	}
#	result <- as.data.frame(result);
#	colnames(result) <- c("Marker", "Obs.Don.Num", "Obs.Rec.Num", "Exp.Don.Num", "Exp.Rec.Num",
#			"Chisq", "P.Value");
#	rownames(result) <- c();
#	
#	return(result);
	
}
