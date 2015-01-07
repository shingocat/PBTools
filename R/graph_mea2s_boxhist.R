# --------------------------------------------------------
# ARGUMENTS:
# data - a dataframe
# respvar - a string; variable name of the response variable
# result - output of single environment analysis
# box - logical; indicates if boxplot(s) is(are) to be created
# hist - logical; indicates if histogram(s) is(are) to be created
#
# Author: Alaine Gulles
# --------------------------------------------------------

graph.mea2s.boxhist <- function(data, respvar, result, box = FALSE, hist = FALSE) UseMethod("graph.mea2s.boxhist")

graph.mea2s.boxhist.default <- function(data, respvar, result, box = FALSE, hist = FALSE) {

  #dir.create("plots")
  
  if (box == TRUE) {
    #create boxplot of raw data (1 file for each respvar)
    for (i in (1:length(respvar))) {
  		boxfile = paste(getwd(),"/boxplotMea2S_",respvar[i],".png",sep = "")
  		if (!all(is.na(data[,match(respvar[i], names(data))]))) {
  		  png(filename = boxfile) #par(mfrow = n2mfrow(length(respvar)));
  		  xlabel = respvar[i]
  		  boxplot((data[,match(respvar[i], names(data))]), data = data,xlab = xlabel, main = paste("Boxplot of ", respvar[i], sep=""));
  		  dev.off()
  		}
  	}
  }
	
  if (hist == TRUE) {
    #create histogram of the raw data (1 file for each envt in each respvar)
  	for (i in (1:length(respvar))) {
  	  if (!is.null(result[[i]]$data)) {
  	    histfile = paste(getwd(),"/histMea2S_",respvar[i],".png",sep="")
  	    png(filename = histfile) # par(mfrow = n2mfrow(length(respvar)*nlevels(data[,match(env, names(data))])));
  	    xlabel = respvar[i]
  	    hist(result[[i]]$data[,match(respvar[i], names(result[[i]]$data))], xlab = xlabel, main = paste("Histogram of ", respvar[i], sep=""));
  	    dev.off()
  	  }
  	}	
  }	
  
}
