###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Sep 25, 2014
# FileName: Haldane.R
###############################################################################


getMarkerGeneDist <- function(data, method = c("Haldane"), byCol = TRUE)
{
	if(!is.matrix(data))
	{
		stop("\tError: only acceptting matrix type of data argument!\n");
	}
	if(byCol)
	{
		if(ncol(data) <=1)
		{
			stop("\tError: There must be at least two markers on column!\n");
		}
	} else
	{
		if(nrow(data) <= 1)
		{
			stop("\tError: There must be at least two markers on row!\n");
		}
	}
}
