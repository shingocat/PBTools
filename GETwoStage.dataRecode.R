# -------------------------------------------------------
# This function performs data recoding for genotype and environment columns for Multi-Env two stage analysis
#
# ARGUMENTs:
# data - data frame}
# respvar - vector of strings; variable names of the response variables}
# geno - string; name of genotype factor}
# env - string; name of environment factor
#
# File Created by: Nellwyn L. Sales 
# --------------------------------------------------------

GETwoStage.dataRecode <- function(data, respvar, geno, env) UseMethod("GETwoStage.dataRecode")

GETwoStage.dataRecode.default <- function(data, respvar, geno, env) {
  
  # --- TRIM THE STRINGS --- #
  respvar<-trimStrings(respvar)
  geno <-trimStrings(geno)
  env <-trimStrings(env)
  
  # --- CHECK INPUT --- #
  if (is.na(match(respvar, names(data))) ||  
        is.na(match(geno, names(data))) || 
        is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }
  
  result <- list()
  
  for (i in (1:length(respvar))) {
    
    # --- CREATE TEMP.DATA WHICH CONTAINS ALL NON-MISSING OBSERVATIONS --- #
    temp.data <- subset(data, subset = (is.na(data[,match(respvar[i], names(data))]) == F))
    
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
      
      result[[i]]$data<-temp.data
    }
  }
  return(result)
}