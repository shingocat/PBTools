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

graph.mea2s.predint <- function(data, respvar, geno, result) UseMethod("graph.mea2s.predint")

graph.mea2s.predint.default <- function(data, respvar, geno, result) {

  #dir.create("plots")
  
  # Prediction intervals of Genotype means 
  for (i in (1:length(respvar))) {
    filename = paste(getwd(), "/predIntMea2S_",respvar[i],".png",sep="");
	png(filename = filename); #par(mfrow = n2mfrow(length(respvar)*nlevels(result[[i]]$site[[j]]$data[,match(env, names(result[[i]]$site[[j]]$data))]))); 
	if (length(levels(result[[i]]$data[,match(geno, names(result[[i]]$data))])) <= 50) {
		print(dotplot(ranef(result[[i]]$model, postVar=TRUE))[[1]]) } else { (qqmath(ranef(result[[i]]$model, postVar=TRUE))[1]) }; #,xlab=xlabel
	dev.off();
  }
  
}
