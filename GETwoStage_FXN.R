# -------------------------------------------------------
# This function persforms multi-environment two stage test
#
# ARGUMENTS:
# data - data frame
# respvar - vector of strings; variable names of the reponse variables
# stderr - string; column name corresponding to the standard error of the means
# sigma2 - string; column name corresponding to the residual variances
# numrep - string; column name corresponding to the number of replicates
# geno - string; name of genotype factor
# env - string; name of environment factor
# weight - weight options are: "none", "stderr" and "stdmean"
# is.genoRandom - logical; TRUE if genotype factor is random
#
# File Created by: Alaine A. Gulles 
# File Modified by: Alaine A. Gulles 
# Script Created by: Violeta Bartolome
# Script Modified by: Violeta Bartolome
#                     Alaine A. Gulles
#                     Rose Imee Zhella Morantte
#                     Nellwyn Sales
# --------------------------------------------------------

GETwoStage.test <- function(data, respvar, stderr = NULL, sigma2, numrep, geno, env, weight = c("none", "stderr", "stdmean"), is.genoRandom = FALSE) UseMethod("GETwoStage.test")

GETwoStage.test.default <- function(data, respvar, stderr = NULL, sigma2, numrep, geno, env, weight = c("none", "stderr", "stdmean"), is.genoRandom = FALSE) {
	                                     
     library(lme4) 
     options(show.signif.stars=FALSE)
 
     # --- TRIM THE STRINGS --- #
     respvar<-trimStrings(respvar)
     if (!is.null(stderr)) { stderr <-trimStrings(stderr) }
     sigma2 <-trimStrings(sigma2)
     numrep <-trimStrings(numrep)
     geno <-trimStrings(geno)
     env <-trimStrings(env)
 
	# --- CHECK INPUT --- #
	if (is.na(match(respvar, names(data))) ||  
      #is.na(match(stderr , names(data))) || 
	    is.na(match(sigma2 , names(data))) || 
	    is.na(match(numrep, names(data))) || 
	    is.na(match(geno, names(data))) || 
	    is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }
	
	weight <- match.arg(weight)
  
     # --- SET GENO AND ENV TO FACTORS --- #
	data[,match(geno, names(data))] <- factor(data[,match(geno, names(data))])
	data[,match(env, names(data))] <- factor(data[,match(env, names(data))])
	result <- list()

	for (i in (1:length(respvar))) {
    
	     # --- CREATE TEMP.DATA WHICH CONTAINS ALL NON-MISSING OBSERVATIONS --- #
		temp.data <- subset(data, subset = (is.na(data[,match(respvar[i], names(data))]) == F))
    
		# --- CREATE COMPUTEDWEIGHT COLUMN --- #
		if (weight == "stderr")  {  temp.data$computedWeight <- 1/(temp.data[, match(stderr[i], names(temp.data))]^2)
          } else {  temp.data$computedWeight <- NULL }
    
		# --- CREATE Y COLUMN --- #
		#if (weight == "stdmean") {
          #  temp.data$Y <- temp.data[,match(respvar[i], names(temp.data))]/temp.data[,match(stderr[i], names(temp.data))]
          #} else {
          #  temp.data$Y <- temp.data[,match(respvar[i], names(temp.data))] 
          #}
    
		result[[i]] <- list()
		result[[i]]$respvar <- respvar[i]
    
		# --- COMPUTE RESPONSE RATE --- #
		obsread <- nrow(data)
		obsused <- nrow(temp.data)
          responseRate <- obsused/obsread
		result[[i]]$obsread <- obsread
		result[[i]]$obsused <- obsused
		result[[i]]$responseRate <- responseRate
    
          if (responseRate < 0.80) {
               result[[i]]$manyNAWarning <- "Too many missing observations. Cannot proceed with the analysis."
               next
          } else {
               temp.data[,match(geno, names(temp.data))] <- factor(trimStrings(temp.data[,match(geno, names(temp.data))]))
               temp.data[,match(env, names(temp.data))] <- factor(trimStrings(temp.data[,match(env, names(temp.data))]))
      
               # --- get levels of genotype and environment --- #
               levelsGeno<-levels(temp.data[,match(geno, names(temp.data))])
               levelsEnv<-levels(temp.data[,match(env, names(temp.data))])
      
               result[[i]]$nlevelsGeno <- length(levelsGeno)
               result[[i]]$nlevelsEnv <- length(levelsEnv)
      
               # --- if max length of the characters of levelsGeno or levelsEnv greater than 4, recode the levels
               if (max(nchar(levelsGeno))>4 || max(nchar(levelsEnv))>4) {
        
                    # --- recode genotype and environment levels --- #
                    newCodingGeno<-data.frame(Genotype=levelsGeno, Code=paste("G",seq(1:length(levelsGeno)), sep=""))
                    newCodingEnv<-data.frame(Environment=levelsEnv, Code=paste("E",seq(1:length(levelsEnv)), sep=""))
        
                    result[[i]]$newCodingGeno <- newCodingGeno
                    result[[i]]$newCodingEnv <- newCodingEnv
                    recodedLevels <- TRUE
                    result[[i]]$recodedLevels <- recodedLevels
        
                    # --- attach the new labels to temp.data --- #
                    temp.data$CodedGeno <- newCodingGeno$Code[match(temp.data[,geno], newCodingGeno$Genotype)]
                    temp.data$CodedEnv <- newCodingEnv$Code[match(temp.data[,env], newCodingEnv$Environment)]
        
               } else {
        
                    temp.data$CodedGeno <- temp.data[,match(geno, names(temp.data))]
                    temp.data$CodedEnv <- temp.data[,match(env, names(temp.data))]
        
                    recodedLevels <- FALSE
                    result[[i]]$recodedLevels <- recodedLevels
               }
                  
               # --- CONSTRUCT THE MODEL --- #
               if (is.genoRandom) { trt.stmt <- paste("(1|", geno,")", sep = "")
               } else {  trt.stmt <- paste(geno, sep = "")	}
      
               #myformula1 <- paste("Y ~ 1 + ", trt.stmt, " + (1|", env,")", sep = "")
               myformula1 <- paste(respvar[i], " ~ 1 + ", trt.stmt, " + (1|", env,")", sep = "")
               model <- lmer(formula(myformula1), weights = temp.data$computedWeight, data = temp.data)
               result[[i]]$formula1 <- myformula1
               result[[i]]$model <- model
      
               # --- SUMMARY STATISTICS PER SITE --- #
               sumStat.Env <- data.frame(tapply(temp.data[, match(respvar[i], names(temp.data))], temp.data[, match(env, names(temp.data))], mean))
               sumStat.Env <- data.frame(rownames(sumStat.Env), sumStat.Env)
               colnames(sumStat.Env) <- c(env, "means")
               rownames(sumStat.Env) <- NULL
      
               numRep <- 1/mean(1/tapply(temp.data[, match(numrep[i], names(temp.data))] , temp.data[, match(env, names(temp.data))], mean))
               numEnv <- nlevels(temp.data[, match(env, names(temp.data))])
               error.variance <- sum(tapply(temp.data[, match(sigma2[i], names(temp.data))] , temp.data[, match(env, names(temp.data))], mean) * tapply(temp.data[, match(numrep[i], names(temp.data))] , temp.data[, match(env, names(temp.data))], mean))/sum(tapply(temp.data[, match(numrep[i], names(temp.data))] , temp.data[, match(env, names(temp.data))], mean))
      
               # the next line of code was revised for R version 3.0.2 by AAGulles 07.28.2014     
               # varcomp_her <- summary(model)@REmat
               
               varcomp_her <- summary(model)$sigma  # added by AAGulles 07.28.2014
               
               # the next line of code was revised for R version 3.0.2 by AAGulles 07.28.2014     
               # error.gxe.var <- as.numeric(varcomp_her[varcomp_her[,1] == "Residual", "Variance"])    #added by NSales 
               
               error.gxe.var <-  varcomp_her**2   # added by AAGulles 07.28.2014    
               
               #if (weight == "none")    sigma2.ge <- attr(VarCorr(model), "sc")**2 - error.variance
               #if (weight == "stderr")  sigma2.ge <- 1
               #if (weight == "stdmean") sigma2.ge <- attr(VarCorr(model), "sc")**2 * error.variance
      
               if (weight == "none")    sigma2.ge <- error.gxe.var - error.variance
               if (weight == "stderr")  sigma2.ge <- 1
               #if (weight == "stdmean") sigma2.ge <- error.gxe.var * error.variance
      
               # --- VARIANCE COMPONENTS --- #
               varcomp <- NULL
               for (j in (1:length(VarCorr(model)))) {
                    varcomp <- rbind(varcomp, data.frame(Groups = names(VarCorr(model))[j], Variance = VarCorr(model)[[j]][1], Std.Dev. = attr(VarCorr(model)[[j]], "stddev")[[1]])) 
               }
      
               # if the value of sigma2.ge is negative, set it to zero
               if (sigma2.ge < 0) {  sigma2.ge <- 0 }
      
               #varcomp <- rbind(varcomp, data.frame(Groups = paste(geno, ":", env, sep = ""), Variance = sigma2.ge, Std.Dev. = sqrt(sigma2.ge)))  #disabled by NSALES 
               varcomp <- rbind(varcomp, data.frame(Groups = "Residual", Variance = attr(VarCorr(model), "sc")**2, Std.Dev. = attr(VarCorr(model), "sc")))
               result[[i]]$varcomp.table <- varcomp	
      
               # --- TEST OF SIGNIFICANCE OF GENO EFFECT --- #
               myformula2 <- gsub(paste(" + ", trt.stmt, sep = ""), "", myformula1, fixed = TRUE)
      
               if (!is.genoRandom) {
                    #model1 <- lmer(formula(myformula1), weights = temp.data$computedWeight, data = temp.data, REML = F)
                    #model2 <- lmer(formula(myformula2), weights = temp.data$computedWeight, data = temp.data, REML = F)
                    #model.comp <- anova(model2, model1)
                    #attr(model.comp, "heading")[3] <- paste("model2: ", myformula2, sep = "")
                    #attr(model.comp, "heading")[4] <- paste("model1: ", myformula1, "\n", sep = "")
                    #attr(model.comp, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF GENOTYPIC EFFECT:\n", sep = "")
                    #result[[i]]$testsig.Geno <- model.comp
        
               } else {
                    # --- compare the two models --- #
                    model2 <- lmer(formula(myformula2), weights = temp.data$computedWeight, data = temp.data, REML = T)
                    models.table2<-modelComparisonTable(model, model2)
                    result[[i]]$formula2 <- myformula2
                    result[[i]]$testsig.Geno <- models.table2
               }
      
               if (is.genoRandom) {
                    # --- TEST OF SIGNIFICANCE OF ENVIRONMENT EFFECT USING LRT --- #
                    myformula3 <- gsub(paste(" + (1|", env, ")", sep = ""), "", myformula1, fixed = TRUE)
                    model3 <- lmer(formula(myformula3), weights = temp.data$computedWeight, data = temp.data, REML = T)
                    models.table3<-modelComparisonTable(model, model3)
                    result[[i]]$formula3 <- myformula3
                    result[[i]]$testsig.Env <- models.table3
        
                    # --- ESTIMATE HERITABILITY --- #
                    genetic.var <- varcomp[varcomp[,1] == geno, "Variance"]
                    heritability <- genetic.var/(genetic.var + (sigma2.ge/numEnv) + (error.variance/(numRep*numEnv)))
                    #heritability <- as.matrix(round(heritability,digits = 2))         #RIZAM 091811
                    #rownames(heritability) <- ""                    #RIZAM 091811
                    heritability <- format(round(heritability,digits = 2), digits=2, nsmall=2, scientific=FALSE)
                    result[[i]]$heritability <- heritability
        
                    # --- PREDICTED MEANS/MEANS OF GENOTYPE --- #
                    sumStat.Geno <- eval(parse(text = paste("coef(model)$", geno, sep = ""))); 
                    sumStat.Geno <- cbind(rownames(sumStat.Geno), sumStat.Geno) 
                    colnames(sumStat.Geno) <- c(geno, "Mean")
               } else {
                    myformula4 <- gsub("~ 1", "~ 0", myformula1, fixed = TRUE)
                    model.noint <- lmer(formula(myformula4), weights = temp.data$computedWeight, data = temp.data)
                    sumStat.Geno <- data.frame(summary(model.noint)$coefficients)[,1:2]
                    rownames(sumStat.Geno) <- gsub(geno,"", rownames(sumStat.Geno))
                    sumStat.Geno <- cbind(rownames(sumStat.Geno), sumStat.Geno)
                    colnames(sumStat.Geno) <- c(geno, "LSMean", "StdErrMean")
        
                    # --- COMPUTE STANDARD ERROR OF THE DIFFERENCE --- #
                    noEntries<-nlevels(temp.data[,match(geno, names(temp.data))])
                    covs <- as.matrix(vcov(model.noint)[1:noEntries, 1:noEntries])
                    vars <- diag(covs)
                    vdiff <- outer(vars, vars, "+") - 2 * covs
                    sed <- sqrt(vdiff[upper.tri(vdiff)])
        
                    # --- DISPLAY SED TABLE --- #
                    minSed<-formatC(as.numeric(format(min(sed), scientific=FALSE)), format="f")
                    meanSed<-formatC(as.numeric(format(mean(sed), scientific=FALSE)), format="f")
                    maxSed<-formatC(as.numeric(format(max(sed), scientific=FALSE)), format="f")
                    sedCol<-rbind(minSed, meanSed, maxSed)
                    rowNames<-rbind("Minimum  ", "Average  ", "Maximum  ")
                    sedTable<-as.table(cbind(rowNames, sedCol))
                    rownames(sedTable)<-c("","","")
                    colnames(sedTable)<-c("","Estimate")
                    result[[i]]$sedTable <- sedTable
               }
               rownames(sumStat.Geno) <- NULL
               result[[i]]$means.Geno <- sumStat.Geno
               result[[i]]$means.Env  <- sumStat.Env
               result[[i]]$residuals  <- resid(model)
               result[[i]]$fitted.values  <- fitted(model)
               result[[i]]$data  <- subset(temp.data, select = c(env, geno, respvar[i], stderr[i], sigma2[i], numrep[i], "CodedGeno", "CodedEnv"))
          }
	}
	detach("package:lme4")
	return(result)
}