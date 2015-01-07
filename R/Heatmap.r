#-------------------------------------------------
# This function outputs heatmap for single-environment analysis
# Created by:  Paul Eilers, Gerrit Gort, Sabine Schnabel for R Analytical Pipeline (RAP) 
# Modified by: Alaine A. Gulles
#              Rose Imee Zhella A. Morantte
#-------------------------------------------------

Heatmap <- function(R, genAs = "fixed", row, column, respvar, model, env) UseMethod("Heatmap")

Heatmap.default <- function(R, genAs = "fixed", row, column, respvar, model, env) {
  
  if (is.null(env)) {
    origEnv<-NULL
    env = "EnvLevel"
    R <- cbind(R, EnvLevel=1)
  } else {
    origEnv<-"not NULL"
  }
  
  R[,match(env, names(R))] <- factor(trimStrings(R[,match(env, names(R))]))
  
  if (any(is.na(as.numeric(trimStrings(R[,match(row, names(R))]))))) {
    stop("The specified row variable is not numeric.")
  } else {
    R[,match(row, names(R))] <- as.numeric(trimStrings(R[,match(row, names(R))]))
  }
  if (any(is.na(as.numeric(trimStrings(R[,match(column, names(R))]))))) {
    stop("The specified column variable is not numeric.")
  } else {
    R[,match(column, names(R))] <- as.numeric(trimStrings(R[,match(column, names(R))]))
  }
  
  envLevels <- levels(R[,match(env, names(R))])
  
  warningList<-list()
  for (i in 1:length(respvar)) {
    warningList[[i]]<-list()
  	for (j in 1:length(envLevels)) {
  	  warningList[[i]]$site[[j]]<-list()
      resid = paste(respvar[i],"_","resid_",genAs,sep = "")
    	if (resid %in% names(R)) {
    	  envlevel = levels(R[,match(env, names(R))])[j]
    	  Rdata1 = subset(R, R[match(env,names(R))] == envlevel)
        
    	  #check first if each row-column combination contains one residual
    	  lengthPerCross<-tapply(Rdata1[,respvar[i]], Rdata1[,c(row,column)], length)
        if (all(lengthPerCross==1, na.rm=TRUE)) {
          warningList[[i]]$site[[j]]<-"unique"
          
          if (is.null(origEnv)) {
            lp = levelplot(Rdata1[[match(resid,names(Rdata1))]] ~ Rdata1[[match(row,names(Rdata1))]] * Rdata1[[match(column,names(Rdata1))]],
                           xlab = row, ylab = column,
                           main = paste("Heatmap of Residuals (", resid, ")", sep = ""),
                           col.regions = colorRampPalette(c("yellow","red"))(50))
          } else {
            lp = levelplot(Rdata1[[match(resid,names(Rdata1))]] ~ Rdata1[[match(row,names(Rdata1))]] * Rdata1[[match(column,names(Rdata1))]],
                           xlab = row, ylab = column,
                           main = paste("Heatmap of Residuals (", resid, ")\n", env, "=", envlevel, sep = ""),
                           col.regions = colorRampPalette(c("yellow","red"))(50))
          }
          png(paste(getwd(), "/heatmap_", resid, "_", envlevel,".png", sep=""))
          print(lp)
          dev.off()
        } else {
          if (is.null(origEnv)) {
            warningList[[i]]$site[[j]]<-paste("Heatmap of residuals (genotype = ", genAs, ") for response variable = ", respvar[i], " is not done.\n       Each residual should have a corresponding unique row-column label.", sep="")
          } else {
            warningList[[i]]$site[[j]]<-paste("Heatmap of residuals (genotype = ", genAs, ") for response response variable = ", respvar[i], " and ", env, " = ", envlevel, " is not done.\n       Each residual should have a corresponding unique row-column label.", sep="")
          }
          
        }
    	} else {
    	  warningList[[i]]$site[[j]]<- "empty"
    	}
  	}
	}
  return(warningList)
}