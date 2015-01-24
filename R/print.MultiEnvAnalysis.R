print.MultiEnvAnalysis <- function(data, level = 1)
{
  if(!inherits(data, "MultiEnvAnalysis"))
    stop("\tError: The argument of data must be of class MultiSiteAnalysis!\n");
  if(missing(level))
    level <- 1;
  if(!is.numeric(level))
    stop("\tError: The argument of level should be of value 1 or 2 where 1 is for concise details, and 2 is for precise details.\n");
  if(length(level) != 1)
    stop("\tError: The argument of level should be of length 1.\n");
  if(!any(level %in% c(1,2)))
    stop("\tError: The argument of level only accepted the integer 1 or 2. \n");
  if(level == 1)
  {
    cat("Multiple Environment Analyis:");
    cat("\n");
    for(i in 1:data$trait.number)
    {
      trait.name <- data$traits[[i]]$name;
      cat("On the trait: ", trait.name, ".\n", sep ="");
      cat("Environment is fixed: ", ifelse(data$traits[[i]]$analysis$mea$envFixed, "YES", "NO"),".\n", sep = "");		
      cat("Variance Table:\n");
      print(data$traits[[i]]$analysis$mea$varcomp.table, row.names = FALSE);
      cat("\n");
      cat("Adjusted Means Across Environment:\n");
      print(data$traits[[i]]$analysis$mea$wide.GenoEnv, row.names = FALSE);
      cat(rep("-", times=50), sep = "");
      cat("\n")
    }
  } else
  {
    for(i in 1:length(data$traits))
    {
      if(is.null(data$traits[[i]]$analysis$mea))
        next;
      cat(rep("=",40), sep="");
      cat("\n");
      cat("RESPONSE VARIABLE: ");
      cat(data$traits[[i]]$name);
      cat("\n");
      cat(rep("=",40), sep="");
      cat("\n");
      cat(rep("-",40), sep="");
      cat("\n");
      cat("ENVIRONMENT AS: ");
      cat(ifelse(data$traits[[i]]$analysis$mea$envFixed, "FIXED", "RANDOM"));
      cat("\n");
      cat(rep("-",40), sep="");
      cat("\n\n");
      cat("DATA SUMMARY:\n");
      cat("Number of observations read:", data$traits[[i]]$analysis$mea$obsread);
      cat("\n");
      cat("Number of observations used:", data$traits[[i]]$analysis$mea$obsused);
      cat("\n");
      print(data$traits[[i]]$analysis$mea$factor.summary, row.names=FALSE);
      cat("\n\n");
      cat("VARIANCE COMPONENT TABLE:")
      cat("\n");
      print(data$traits[[i]]$analysis$mea$varcomp.table, row.names=FALSE);
      cat("\n\n");
      print(data$traits[[i]]$analysis$mea$testsig.Geno, row.names=FALSE);
      cat("\n\n");
      print(data$traits[[i]]$analysis$mea$testsig.Env,row.names=FALSE);
      cat("\n\n");
      print(data$traits[[i]]$analysis$mea$testsig.GenoEnv, row.names=FALSE);
      cat("\n\n");
      cat("GENOTYPE BY ENVIRONMENT MEANS:\n")
      print(data$traits[[i]]$analysis$mea$wide.GenoEnv, row.names=FALSE);
      cat("\n");
      cat("GENOTYPE LSMEANS AND STANDARD ERRORS:\n");
      print(data$traits[[i]]$analysis$mea$means.Geno, row.names=FALSE);
      cat("\n");      
      cat("STANDARD ERROR OF THE DIFFERENCE (SED):\n")
      print(data$traits[[1]]$analysis$mea$sedTable, row.names=FALSE);
      cat("\n\n");
    }
  }
}