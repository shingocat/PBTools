# --------------------------------------------------------
# ARGUMENTS:
# data - a dataframe
# respvar - a string; variable name of the response variable
# env - a string; variable name of the environment variable
# result - output of single environment analysis
# box - logical; indicates if boxplot(s) is(are) to be created
# hist - logical; indicates if histogram(s) is(are) to be created
#
# Author: Alaine Gulles
# --------------------------------------------------------

graph.sea.boxhist <- function(data, respvar, env, result, box = FALSE, hist = FALSE) UseMethod("graph.sea.boxhist")

graph.sea.boxhist.default <- function(data, respvar, env, result, box = FALSE, hist = FALSE) {

  #dir.create("plots")
  if (is.null(env)) {
    env = "EnvLevel"
    data <- cbind(data, EnvLevel=1)
  }
  
  if (box == TRUE) {
    #create boxplot of raw data (1 file for each respvar)
    for (i in (1:length(respvar))) {
  		boxfile = paste(getwd(),"/boxplotSEA_",respvar[i],".png",sep = "")
  		if (!all(is.na(data[,match(respvar[i], names(data))]))) {
  		  png(boxfile) #par(mfrow = n2mfrow(length(respvar)));
  		  xlabel = respvar[i]
  		  boxplot(data[,match(respvar[i], names(data))]~data[,match(env, names(data))],data = data,xlab = xlabel, main = paste("Boxplot of ", respvar[i], sep=""));
  		  dev.off()
  		}
  	}
  }
	
  if (hist == TRUE) {
    #create histogram of the raw data (1 file for each envt in each respvar)
  	for (i in (1:length(respvar))) {
  		for (j in (1:nlevels(factor(data[,match(env, names(data))])))) {
        if (!is.null(result$output[[i]]$site[[j]]$data)) {
          siteLabel = result$output[[i]]$site[[j]]$env[[1]]
          histfile = paste(getwd(),"/histSEA_",respvar[i],"_",siteLabel,".png",sep="")
          png(histfile) # par(mfrow = n2mfrow(length(respvar)*nlevels(data[,match(env, names(data))])));
          xlabel = respvar[i]
          hist(result$output[[i]]$site[[j]]$data[,match(respvar[i], names(result$output[[i]]$site[[j]]$data))], xlab = xlabel,main = paste("Histogram of ", respvar[i], " for ", env, " = ", siteLabel, sep=""));
          dev.off()
        }
  		}	
  	}	
  }	
  
}
