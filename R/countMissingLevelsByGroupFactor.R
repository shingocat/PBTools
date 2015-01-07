###############################################################################
# TODO: Count whether missing whole value of response variable by group factor 
# 		or missing whole levels of value of response variable, if Yes return 
#		TRUE else FALSE
#
##		respvar: response variable be count
##		groupfactor: dividing response variable by these group factor
## 		data: containning the original data
#
#
# Author: mqin
# Date: Sep 20, 2014
# FileName: countByGroupFactor.R
###############################################################################


countMissingLevelsByGroupFactor <- function(respvar, groupfactor, data)
{
	if(missing(respvar))
	{
		stop("\tError: The respvar must not be null!\n");
	}
	if(missing(groupfactor))
	{
		stop("\tError: The groupfactor must not be null!\n");
	}
	if(missing(data))
	{
		stop("\tError: The data must not be null!\n");
	}
	if(!is.data.frame(data))
	{
		stop("\tError: The data must be a data frame object in R!\n");
	}
	
	# checking whether the groupfactor name(s) exist in the data;
	if(!all(!is.na(match(groupfactor, names(data)))))
	{
		stop("\tError: Some of groupfactor names are not in the data!\n");
	}
	
	# factoring the groupfactor variable of data
	for(i in 1 : length(groupfactor))
	{
		data[[groupfactor[i]]] <- factor(data[[groupfactor[i]]]);
	}
	is.Exist = FALSE;
	
	# Checking the groupfactor interaction levels is larger than particular factor levels of group factor
	temp <- list();
	for(i in 1:length(groupfactor))
	{
		temp[[groupfactor[i]]] <- levels(data[[groupfactor[i]]]);
	}
	#--- modify this line because it would show warning message In ans * length(l) + if1 : ---#
    #--- longer object length is not a multiple of shorter object length ---#
	#--- replace by a inner function ---#
	# inter.levels <- levels(interaction(temp, sep =""));
	int <- function(temp)
	{
		int.temp <- list()
		if(length(temp) >= 1)
			int.temp[[1]] <- factor(temp[[1]]);
		if(length(temp) == 1)
		{
			return(int.temp);
		} else
		{
			for(i in 2:length(temp))
			{
				x <- NULL;
				for(j in 1:length(temp[[i]]))
				{
					x <- c(x, paste(int.temp[[i - 1]], temp[[i]][j], sep=""));
				}
				int.temp[[i]] <- x;
			}
			int.temp[[length(int.temp)]];
			return(int.temp[[length(int.temp)]]);
		}
	}
	inter.levels <- int(temp);
	
	data.levels <- levels(as.factor(do.call(paste0, data[match(groupfactor,names(data))])));
	if(!all(inter.levels %in% data.levels))
	{
		is.Exist = TRUE
		return (is.Exist);
	}
	
	# Checking whether the responvar vairable are missing all the values on level of groupfactor combination.
	# it means only with record but without values;
	temp <- aggregate(data[respvar], by =data[groupfactor],FUN = function(x){all(is.na(x))});
	if(any(temp[[respvar]]))
	{
		is.Exist = TRUE;
		return (is.Exist);
	}
	
	return(is.Exist);
}
