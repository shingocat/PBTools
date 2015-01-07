###############################################################################
# TODO: Getting the default expect genotypic frequency of each marker on 
#		the specified generation.
# [Argument]
# BCn: denoted the backcross generation, default is 5;
# Fn: denoted the selfing or intercrossing generation, default is 1:
# BC.parent: denoted the backcross parent, if it is 1, then will return "AA", "Aa", if it is 2, then will return "aa", "Aa"
# codominant : denoted whether the marker is codominant or not! 
#
# [return]
# 	return a vector of expect genotypic frequency, if there are only 
#		backcross, it will be on two values, if including selfing or 
#		intercrossing, it will be on three values.
# 
# Author: mqin
# Date: Sep 25, 2014
# FileName: getExpectGenotypicFreq.R
###############################################################################


getExpectGenotypicFreq <- function(BCn=5, Fn=1, BC.parent = 1)
#getExpectGenotypicFreq <- function(BCn=5, Fn=1, BC.parent = 1, codominant = TRUE)	# not implemented right now	
{
	if(!is.numeric(BCn))
	{
		stop("\tError: The argument of BCn should be numeric!\n");
	}
	if(!is.numeric(Fn))
	{
		stop("\tError: The argument of Fn should be numeric!\n");
	}
	if(!is.numeric(BC.parent))
	{
		stop("\tError: The argument of BC.parnet should be numeric!\n");
	}
	if(length(BCn) != 1)
	{
		stop("\tError: The argument of BCn should be of length one!\n");
	}
	if(length(Fn) != 1)
	{
		stop("\tError: The argument of Fn should be of length one!\n");
	}
	if(length(BC.parent) != 1)
	{
		stop("\tError: The argument of BC.parent should be of length one!\n");
	}
	if(is.na(match(BC.parent, c(1,2))))
	{
		stop("\tError: The argument of BC.parent should be 1 or 2!\n");
	}
	if(BCn < 0)
	{
		stop("\tError: The argument of BCn should not be negative!\n");
	}
	if(Fn < 0)
	{
		stop("\tError: The argument of Fn should not be negative!\n");
	}
	
	geno.labels <- c("AA", "Aa", "aa");
	geno.freq <- c();
	if(BCn != 0 && Fn == 0)
	{
		heter.freq = (1/2)^BCn;
		homo.freq = 1 - heter.freq;
		geno.freq = c(homo.freq, heter.freq);
		if(BC.parent == 1)
		{
			names(geno.freq) <- c(geno.labels[1:2]);
		} else
		{
			names(geno.freq) <- c(geno.labels[3:2]);
		}
		
	} else if(BCn == 0 && Fn != 0)
	{
		heter.freq = 2^(Fn -1) / 4^(Fn-1);
		homo.freq = 1 - heter.freq;
		geno.freq <- c(homo.freq/2, heter.freq, homo.freq/2);
		names(geno.freq) <- geno.labels;
	} else if(BCn != 0 && Fn != 0)
	{
		heter.freq = (1/2)^BCn;
		if(BC.parent == 1)
		{
			homo.freq.A = 1 - heter.freq;
			homo.freq.a = 0;
			for( i in 1:Fn)
			{
				heter.freq = heter.freq * (2/4);
				homo.freq.A = homo.freq.A + (heter.freq)/2;
				homo.freq.a = homo.freq.a + (heter.freq)/2;
			}
			geno.freq <- c(homo.freq.A, heter.freq, homo.freq.a);
			names(geno.freq) <- geno.labels;
		} else 
		{
			homo.freq.A = 0;
			homo.freq.a = 1 - heter.freq;
			for( i in 1:Fn)
			{
				heter.freq = heter.freq * (2/4);
				homo.freq.A = homo.freq.A + (heter.freq)/2;
				homo.freq.a = homo.freq.a + (heter.freq)/2;
			}
			geno.freq <- c(homo.freq.A, heter.freq, homo.freq.a);
			names(geno.freq) <- geno.labels;
		}
	} else
	{
		warning("\tWarning: There are no specified generation!\n");
		geno.freq <- c(0,0,0);
		names(geno.freq) <- geno.labels;
	}
	
	return(geno.freq);
}
