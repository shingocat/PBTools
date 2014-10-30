# -------------------------------------------------------
# GENOTYPIC AND PHENOTYPIC CORRELATION 
# File Created by: Alaine A. Gulles 
# Script Modified by: Rose Imee Zhella A. Morantte
#                     Nellwyn L. Sales
# --------------------------------------------------------

# --------------------------------------------------------
# ARGUMENTS:
# exptl.design = RCB, AugRCB, AugLS, Alpha, RowCol, Latinized Alpha, Latinized RowCol
# data - a dataframe
# respvar - a vector of strings containing variable names of the response variable
# geno - a string; variable name of the treatment/genotype variable
# row - a string; variable name of the blocking variable or row variable
# column - a string; variable name of the column variable
#        - NULL, if design is RCB, Aug RCB, Alpha Lattice, Latinized Alpha
# rep - a string; variable name of the replication variable
#       - NULL, if design is RCB, Aug RCB, Aug LS, 
# env - a string; variable name of the environment variable
# excludeLevels - logical; TRUE if some genotype levels will be excluded
# excludeList - vector of genotype levels to exclude
# --------------------------------------------------------

genoNpheno.corr <- function(exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"), data, respvar, geno, row, column = NULL, rep = NULL, env, excludeLevels=FALSE, excludeList = NULL) UseMethod("genoNpheno.corr") 

genoNpheno.corr.default <- function(exptl.design = c("RCB", "AugRCB", "AugLS", "Alpha", "RowCol", "LatinAlpha", "LatinRowCol"), data, respvar, geno, row, column = NULL, rep = NULL, env, excludeLevels=FALSE, excludeList = NULL){

  library(lme4) 
  
	if (is.null(env)) {
	  env = "EnvLevel"
	  data <- cbind(data, EnvLevel=1)
	}
  
  if (!is.null(excludeList)) {
    excludeList<-trimStrings(excludeList)
  }
  
	if (length(respvar) <= 1) stop("Correlation cannot be performed. At least two response variables are required.")
	if (exptl.design == "RCB" || exptl.design == "AugRCB") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }} 
	if (exptl.design == "AugLS") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(column, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}
	if (exptl.design == "Alpha" || exptl.design == "LatinAlpha") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(rep, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}
	if (exptl.design == "RowCol" || exptl.design == "LatinRowCol") { if ( is.na(match(respvar, names(data))) ||  is.na(match(geno, names(data))) || is.na(match(row, names(data))) || is.na(match(column, names(data))) || is.na(match(rep, names(data))) || is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }}

	data[,match(geno, names(data))] <- factor(trimStrings(data[,match(geno, names(data))]))
	data[,match(env, names(data))] <- factor(trimStrings(data[,match(env, names(data))]))
	data[,match(row, names(data))] <- factor(trimStrings(data[,match(row, names(data))]))
	if (!is.null(column)) data[,match(column, names(data))] <- factor(trimStrings(data[,match(column, names(data))]))
	if (!is.null(rep)) data[,match(rep, names(data))] <- factor(trimStrings(data[,match(rep, names(data))]))
	
  EnvLevels <- list()
	GenoCorr <- list()
	PhenoCorr <- list()
	NumObs <- list()

	for (i in (1:nlevels(data[,match(env, names(data))]))) {
    EnvLevels[[i]] <- levels(data[,match(env, names(data))])[i]
		GenoCorr[[i]] <- matrix(NA, nrow = length(respvar), ncol = length(respvar))
		PhenoCorr[[i]] <- matrix(NA, nrow = length(respvar), ncol = length(respvar))
		NumObs[[i]] <- matrix(NA, nrow = length(respvar), ncol = length(respvar))
		temp.data <- subset(data, data[,match(env, names(data))] == levels(data[,match(env, names(data))])[i])
    
		#if design is Latinized Row-Column, check if the data follow case1 or case3 labeling --- #
		if (exptl.design == "LatinRowCol") {
		  lengthPerCross<-tapply(temp.data[,respvar[1]], temp.data[,c(row,column)], length)
		  if (all(lengthPerCross<=1, na.rm=TRUE)) {
		    if (nlevels(data[, row])>nlevels(data[, column])) {
		      longerRow<-TRUE
		    } else {
		      longerRow<-FALSE
		    }
		  } else {
		    stop("The levels of the row/column variable should be continuous across replicates.")
		  }
		}
		
		if (!excludeLevels) {
		  if (exptl.design == "RCB" | exptl.design == "AugRCB") { myformula3 <- paste("Ysum ~ 1 + (1|", row, ") + (1|", geno,")", sep = "") }
		  if (exptl.design == "AugLS") { myformula3 <- paste("Ysum ~ 1 + (1|", row, ") + (1|", column, ") + (1|", geno,")", sep = "") }
		  if (exptl.design == "Alpha") { myformula3 <- paste("Ysum ~ 1 + (1|", rep,"/", row,") + (1|", geno,")", sep = "") }
		  if (exptl.design == "RowCol") { myformula3 <- paste("Ysum ~ 1 + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,") + (1|", geno,")", sep = "") }
		  if (exptl.design == "LatinAlpha") { myformula3 <- paste("Ysum ~ 1 + (1|", rep,") + (1|", row,") + (1|", rep, ":", row, ") + (1|", geno,")", sep = "") } 
		  if (exptl.design == "LatinRowCol") { 
		    if (longerRow) {
		      myformula3 <- paste("Ysum ~ 1 + (1|", rep,") + (1|", column, ") + (1|", rep,":", column,") + (1|", row, ") + (1|", geno,")", sep = "") 
		    } else {
		      myformula3 <- paste("Ysum ~ 1 + (1|", rep,") + (1|", row, ") + (1|", rep,":", row,") + (1|", column, ") + (1|", geno,")", sep = "")
		    }
		  }
		}else {
		  if (exptl.design == "RCB" | exptl.design == "AugRCB") { myformula3 <- paste("Ysum ~ 1 + Check + (1|", row, ") + (1|Test:Check)", sep = "") }
		  if (exptl.design == "AugLS") { myformula3 <- paste("Ysum ~ 1 + Check + (1|", row, ") + (1|", column, ") + (1|Test:Check)", sep = "") }
		  if (exptl.design == "Alpha") { myformula3 <- paste("Ysum ~ 1 + Check + (1|", rep,"/", row,") + (1|Test:Check)", sep = "") }
		  if (exptl.design == "RowCol") { myformula3 <- paste("Ysum ~ 1 + Check + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,") + (1|Test:Check)", sep = "") }
		  if (exptl.design == "LatinAlpha") { myformula3 <- paste("Ysum ~ 1 + Check + (1|", rep,") + (1|", row,") + (1|", rep, ":", row, ") + (1|Test:Check)", sep = "") } 
		  if (exptl.design == "LatinRowCol") { 
		    if (longerRow) {
		      myformula3 <- paste("Ysum ~ 1 + Check + (1|", rep,") + (1|", column, ") + (1|", rep,":", column,") + (1|", row, ") + (1|Test:Check)", sep = "") 
		    } else {
		      myformula3 <- paste("Ysum ~ 1 + Check + (1|", rep,") + (1|", row, ") + (1|", rep,":", row,") + (1|", column, ") + (1|Test:Check)", sep = "") 
		    }
		  }
		}
    
		for (j in (1:(length(respvar) - 1))) {
			for (k in (1:length(respvar))) {
				if (k > j) {
          
				  #check if respvar[[j]] and respvar[[k]] are all NAs
				  checkData1 <- subset(temp.data, subset = (is.na(temp.data[,match(respvar[[j]], names(temp.data))]) == F))
				  checkData2 <- subset(temp.data, subset = (is.na(temp.data[,match(respvar[[k]], names(temp.data))]) == F))
          
				  if (nrow(checkData1)==0 | nrow(checkData2)==0) {
				    GenoCorr[[i]][j,k] <- "error"
				    PhenoCorr[[i]][j,k] <- "error"
				    GenoCorr[[i]][k,j] <- "error"
				    PhenoCorr[[i]][k,j] <- "error"
				    NumObs[[i]][k,j] <- 0 
            NumObs[[i]][j,k] <- 0
				    next
            
				  } else {
				    temp.data$Ysum <- temp.data[,respvar[[j]]] + temp.data[,respvar[[k]]]
				    mydata <- subset(temp.data, subset = (is.na(Ysum) == F))
				    
				    # --- if design if AugRCB or AugLS, automatically define the vector newCheckList
				    if (exptl.design == "AugRCB" | exptl.design == "AugLS") {
				      library(doBy)
				      nobs <- summaryBy(formula(paste(respvar[j], " ~ ", geno, sep = "")), data=mydata, FUN=length)
				      excludeList <- as.vector(nobs[nobs[,match(paste(respvar[j],".length", sep=""), names(nobs))]>1, geno])
				      detach("package:doBy")
				    }
				    
				    if (excludeLevels) {
              #check if all items in excludeList are in genoLevels
				      excludeList<-trimStrings(excludeList)
				      levelsGenoTemp<-levels(mydata[,match(geno, names(temp.data))])
				      isExcludePresent<-excludeList %in% levelsGenoTemp
				      deletedCheck<-excludeList[!isExcludePresent]
              
              if (length(deletedCheck)==(length(excludeList))) {
                trmt.label <- geno 
                if (exptl.design == "RCB" | exptl.design == "AugRCB") { myformula3 <- paste("Ysum ~ 1 + (1|", row, ") + (1|", geno,")", sep = "") }
                if (exptl.design == "AugLS") { myformula3 <- paste("Ysum ~ 1 + (1|", row, ") + (1|", column, ") + (1|", geno,")", sep = "") }
                if (exptl.design == "Alpha") { myformula3 <- paste("Ysum ~ 1 + (1|", rep,"/", row,") + (1|", geno,")", sep = "") }
                if (exptl.design == "RowCol") { myformula3 <- paste("Ysum ~ 1 + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,") + (1|", geno,")", sep = "") }
                if (exptl.design == "LatinAlpha") { myformula3 <- paste("Ysum ~ 1 + (1|", rep,") + (1|", row,") + (1|", rep, ":", row, ") + (1|", geno,")", sep = "") } 
                if (exptl.design == "LatinRowCol") { 
                  if (longerRow) {
                    myformula3 <- paste("Ysum ~ 1 + (1|", rep,") + (1|", column, ") + (1|", rep,":", column,") + (1|", row, ") + (1|", geno,")", sep = "") 
                  } else {
                    myformula3 <- paste("Ysum ~ 1 + (1|", rep,") + (1|", row, ") + (1|", rep,":", row,") + (1|", column, ") + (1|", geno,")", sep = "")
                  }
                }
              } else {
                trmt.label <- "Test:Check"
                for (l in (1:nrow(mydata))) {
                  if (is.na(match(mydata[l,match(geno, names(mydata))], excludeList))) {
                    mydata$Test[l] <- levels(mydata[,match(geno, names(mydata))])[mydata[l,match(geno, names(mydata))]]
                    mydata$Check[l] <- "0"
                  }
                  else {
                    mydata$Test[l] <- "0"
                    mydata$Check[l] <- excludeList[match(mydata[l,match(geno, names(mydata))], excludeList)]
                  }
                }
                mydata$Test <- factor(mydata$Test)
                mydata$Check <- factor(mydata$Check)
              }
				    } else { trmt.label <- geno }
				    myformula1 <- gsub("Ysum ~", paste(respvar[j]," ~", sep=""), myformula3, fixed = TRUE)
				    myformula2 <- gsub("Ysum ~", paste(respvar[k]," ~", sep=""), myformula3, fixed = TRUE)
				    
				    # --- DEFINE MODEL --- #
				    model1 <- lmer(formula(myformula1), data = mydata)
				    model2 <- lmer(formula(myformula2), data = mydata)
				    model3 <- lmer(formula(myformula3), data = mydata)
				    
				    # --- COMPUTE HARMONIC MEANS --- #
				    no.reps <- data.frame(n = tapply(eval(parse(text = paste("mydata$Ysum"))), eval(parse(text = paste("mydata$", geno, sep = ""))), FUN = length))
				    no.reps <- as.numeric(1/mean(1/no.reps$n, na.rm=TRUE)) 
				    
				    # --- COMPUTE GENOTYPIC AND PHENOTYPIC VARIANCES --- #
				    varcomp1 <- NULL
				    for (l in (1:length(VarCorr(model1)))) { varcomp1 <- rbind(varcomp1, data.frame(Groups = names(VarCorr(model1))[l], Variance = VarCorr(model1)[[l]][1], Std.Dev. = attr(VarCorr(model1)[[l]], "stddev")[[1]])) }
				    varcomp1 <- rbind(varcomp1, data.frame(Groups = "Residual", Variance = attr(VarCorr(model1), "sc")**2, Std.Dev. = attr(VarCorr(model1), "sc")))
				    genetic.var1 <- varcomp1[varcomp1[,1] == trmt.label, "Variance"]
				    pheno.var1 <- varcomp1[varcomp1[,1] == trmt.label, "Variance"] + (varcomp1[varcomp1[,1] == "Residual", "Variance"]/no.reps)
				    
				    varcomp2 <- NULL
				    for (l in (1:length(VarCorr(model2)))) { varcomp2 <- rbind(varcomp2, data.frame(Groups = names(VarCorr(model2))[l], Variance = VarCorr(model2)[[l]][1], Std.Dev. = attr(VarCorr(model2)[[l]], "stddev")[[1]])) }
				    varcomp2 <- rbind(varcomp2, data.frame(Groups = "Residual", Variance = attr(VarCorr(model2), "sc")**2, Std.Dev. = attr(VarCorr(model2), "sc")))
				    genetic.var2 <- varcomp2[varcomp2[,1] == trmt.label, "Variance"]
				    pheno.var2 <- varcomp2[varcomp2[,1] == trmt.label, "Variance"] + (varcomp2[varcomp2[,1] == "Residual", "Variance"]/no.reps)
				    
				    varcomp3 <- NULL
				    for (l in (1:length(VarCorr(model3)))) { varcomp3 <- rbind(varcomp3, data.frame(Groups = names(VarCorr(model3))[l], Variance = VarCorr(model3)[[l]][1], Std.Dev. = attr(VarCorr(model3)[[l]], "stddev")[[1]])) }
				    varcomp3 <- rbind(varcomp3, data.frame(Groups = "Residual", Variance = attr(VarCorr(model3), "sc")**2, Std.Dev. = attr(VarCorr(model3), "sc")))
				    genetic.var <- (varcomp3[varcomp3[,1] == trmt.label, "Variance"] - genetic.var1 - genetic.var2)/2
				    pheno.var <- genetic.var + (((varcomp3[varcomp3[,1] == "Residual", "Variance"] - varcomp1[varcomp1[,1] == "Residual", "Variance"] - varcomp2[varcomp2[,1] == "Residual", "Variance"])/2)/no.reps)
				    
            genoCorrValue<-genetic.var/sqrt(genetic.var1*genetic.var2)
            phenoCorrValue<-pheno.var/sqrt(pheno.var1*pheno.var2)
            
				    lowerLimit<--1
				    
				    # --- check if geno and pheno correlation values are within -1 to 1, if not, set it to -1 or 1
            if (!is.nan(genoCorrValue)) {
              if (genoCorrValue<lowerLimit) genoCorrValue<--1
              if (genoCorrValue>1) genoCorrValue<-1
            }
				    
				    if (!is.nan(phenoCorrValue)) {
				      if (phenoCorrValue<lowerLimit) phenoCorrValue<--1
				      if (phenoCorrValue>1) phenoCorrValue<-1
				    }
				                
				    GenoCorr[[i]][j,k] <- formatC(genoCorrValue, format="f")
				    PhenoCorr[[i]][j,k] <- formatC(phenoCorrValue, format="f")
				    NumObs[[i]][j,k] <- length(mydata$Ysum)
				   
				    # --- assign values to the lower triangular part of the matrix
				    GenoCorr[[i]][k,j] <- GenoCorr[[i]][[j,k]]
				    PhenoCorr[[i]][k,j] <- PhenoCorr[[i]][j,k]
				    NumObs[[i]][k,j] <- NumObs[[i]][j,k]
				  }
				}  ## end stmt -- if (k > j)
			} ## end stmt -- for (k in (1:length(respvar)))
		} ## end stmt -- for (j in (1:(length(respvar) - 1)))
		colnames(GenoCorr[[i]]) <- rownames(GenoCorr[[i]]) <- respvar
		GenoCorr[[i]] <- format(GenoCorr[[i]], justify="right")
		GenoCorr[[i]] <- noquote(gsub("NA", "", GenoCorr[[i]]))
		
		colnames(PhenoCorr[[i]]) <- rownames(PhenoCorr[[i]]) <- respvar
		PhenoCorr[[i]] <- format(PhenoCorr[[i]], justify="right")
		PhenoCorr[[i]] <- noquote(gsub("NA", "", PhenoCorr[[i]]))
		
		NumObs[[i]] <- as.table(NumObs[[i]])
		colnames(NumObs[[i]]) <- rownames(NumObs[[i]]) <- respvar
	} ## end stmt -- for (i in (1:nlevels(data[,match(env, names(data))])))
	detach("package:lme4")
	return(list(EnvLevels=EnvLevels, GenoCorr = GenoCorr, PhenoCorr = PhenoCorr, NumObs = NumObs))
}
