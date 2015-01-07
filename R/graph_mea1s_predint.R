# --------------------------------------------------------
# ARGUMENTS:
# data - a dataframe
# respVar - a string; variable name of the response variable
# env - a string; variable name of the environment variable
# result - output of single environment analysis
# box - logical; indicates if boxplot(s) is(are) to be created
# hist - logical; indicates if histogram(s) is(are) to be created
#
# Author: Alaine Gulles
# --------------------------------------------------------

graph.mea1s.predint <- function(data, respVar, geno, result) UseMethod("graph.mea1s.predint")

graph.mea1s.predint.default <- function(data, respVar, geno, result) {

  #dir.create("plots")
  
  # Prediction intervals of Genotype means 
  for (i in (1:length(respVar))) {
    filename = paste(getwd(), "/predIntMea1S_",respVar[i],".png",sep="");
	png(filename = filename); #par(mfrow = n2mfrow(length(respVar)*nlevels(result[[i]]$site[[j]]$data[,match(env, names(result[[i]]$site[[j]]$data))]))); 
	if (length(levels(result$output[[i]]$data[,match(geno, names(result$output[[i]]$data))])) <= 50) {
		print(dotplot(ranef(result$output[[i]]$model, postVar=TRUE))[[1]])} else { (qqmath(ranef(result$output[[i]]$model, postVar=TRUE))[1]) } #,xlab=xlabel
  	dev.off()
  }
  
}
