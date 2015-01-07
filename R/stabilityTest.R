#--------------------------------------------------------------------------------
# This function performs stability analysis and writes the results to a text file
# Created by: Nellwyn L. Sales
#
# PARAMETERS:
# data - data frame
# outFileName - path and filename of the text file to be created (example: "E:/testingFolder/output.txt")
# respvar - vector of the names of response variable
# environment - string; name of environment factor
# genotype - string; name of genotype factor
# useFinlay - logical; set TRUE is Finlay-Wilkinson model will be used
# useShukla - logical; set TRUE is Shukla's model will be used
#
# FUNCTIONS NEEDED: stability.analysis, trimStrings, reshape {stats}
#--------------------------------------------------------------------------------

stabilityTest <- function (data, outFileName, respvar, environment, genotype, useFinlay=TRUE, useShukla=FALSE) UseMethod("stabilityTest")

stabilityTest.default <- function (data, outFileName, respvar, environment, genotype, useFinlay=TRUE, useShukla=FALSE) {

  #trim the strings
  outFileName<-trimStrings(outFileName)
  respvar<-trimStrings(respvar)
  environment<-trimStrings(environment)
  genotype<-trimStrings(genotype)
  
  #add title
  capture.output(cat("\nSTABILITY ANALYSIS\n\n\n"),file=outFileName,append = FALSE)
  
  result<-list()
  
  for (i in 1:length(respvar)) {
    temp.data <- data
    
    result$respvar[[i]]<-list()
    
    #set environment and genotype to factors
    temp.data[,environment]<-factor(temp.data[,environment])
    temp.data[,genotype]<-factor(temp.data[,genotype])
    
    capture.output(cat("------------------------------\n"),file=outFileName,append = TRUE)
    capture.output(cat("RESPONSE VARIABLE: ", respvar[i], "\n"),file=outFileName,append = TRUE)
    capture.output(cat("------------------------------\n\n\n"),file=outFileName,append = TRUE)
    result$respvar[[i]]<-respvar[i]
    
    if (nlevels(temp.data[,environment])>4) {
      
      #check if data contains one value for each env-geno combination
      rep<-tapply(temp.data[, respvar[i]] , temp.data[,c(environment, genotype)], function(x) length(which(!is.na(x))))
      
      #compute response rate
      responseRate<-1-(sum(rep==0, na.rm=TRUE)/(nlevels(temp.data[,environment])*nlevels(temp.data[,genotype])))
      
      if (responseRate < 0.80) {
        
        capture.output(cat("*** \nERROR: Too many missing observations. Cannot proceed with the analysis.\n***\n\n\n"),file=outFileName,append = TRUE)
        suppressWarnings(result$respvar[[i]]$Error<-"ERROR: Too many missing observations. Cannot proceed with the analysis.")
        
      } else {
        
        #if data contains more than 1 observation per env-geno combination, compute for the means
        if (any(rep>1, na.rm=TRUE)) {
          meansWide<-as.data.frame(tapply(temp.data[, respvar[i]] , data[,c(environment, genotype)], function(x) mean(x,na.rm=TRUE)))
          meansWide<-cbind(envTemp=rownames(meansWide), meansWide)
          colnames(meansWide)[1]<-environment
          temp.data<-data.frame(reshape(meansWide, direction="long", varying=colnames(meansWide)[2:ncol(meansWide)], v.names=respvar[i], timevar=genotype, idvar=environment, times=colnames(meansWide)[2:ncol(meansWide)]), row.names = NULL)
        }
        
        if (useFinlay) {
          capture.output(cat("MODEL: FINLAY-WILKINSON\n\n"),file=outFileName,append = TRUE)
          
          result.temp<-try(stability.analysis(temp.data, respvar[i], genotype, environment, method = "regression"), silent=TRUE)
          
          if (!is.null(class(result.temp)) & class(result.temp)=="try-error") {  
            msg <- trimStrings(strsplit(result.temp, ":")[[1]])
            msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
            msg <- gsub("\"", "", msg)
            capture.output(cat("*** \nERROR in stability.analysis function:\n  ",msg, "\n***\n\n\n", sep = ""), file=outFileName,append = TRUE)
            suppressWarnings(result$respvar[[i]]$Error<-paste("ERROR: ", msg, sep=""))
          } else {
            capture.output(result.temp[[1]][[1]]$stability,file=outFileName,append = TRUE)
            capture.output(cat("\n\n"),file=outFileName,append = TRUE)
            suppressWarnings(result$respvar[[i]]$Finlay<-result.temp[[1]][[1]]$stability)
          }
        } 
        
        if (useShukla) {
          capture.output(cat("MODEL: SHUKLA\n\n"),file=outFileName,append = TRUE)
          
          result.temp2<-try(stability.analysis(temp.data, respvar[i], genotype, environment, method = "shukla"), silent=TRUE)
          
          if (!is.null(class(result.temp2)) & class(result.temp2)=="try-error") {  
            msg <- trimStrings(strsplit(result.temp2, ":")[[1]])
            msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "))
            msg <- gsub("\"", "", msg)
            capture.output(cat("*** \nERROR in stability.analysis function:\n  ",msg, "\n***\n\n\n", sep = ""), file=outFileName,append = TRUE)
            suppressWarnings(result$respvar[[i]]$Error<-paste("ERROR: ", msg, sep=""))
          } else {
            capture.output(result.temp2[[1]][[1]]$stability,file=outFileName,append = TRUE)
            capture.output(cat("\n\n"),file=outFileName,append = TRUE)
            suppressWarnings(result$respvar[[i]]$Shukla<-result.temp2[[1]][[1]]$stability)
          }
        } 
      }
      
    } else {
      capture.output(cat("*** \nERROR: The environment factor should have at least five levels.\n***\n\n\n"),file=outFileName,append = TRUE)
      suppressWarnings(result$respvar[[i]]$Error<-"ERROR: The environment factor should have at least five levels.")
    } #end of if (nlevel(temp.data[,environment])>4)
    
    capture.output(cat("\n"),file=outFileName,append = TRUE)
    
  } #end of for (i in 1:length(respvar))
  
  return(result)
} 
