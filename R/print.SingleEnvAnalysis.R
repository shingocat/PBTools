print.SingleEnvAnalysis <- function
(
  data,
  level = 1
)
{
  if(!inherits(data, "SingleEnvAnalysis"))
    stop("\tError: The argument of data must be of class SingleSiteAnalysis!\n");
  if(missing(level))
    level <- 1;
  if(!is.numeric(level))
    stop("\tError: The argument of level should be of value 1 or 2 where 1 is for concise details, and 2 is for precise details.\n");
  if(length(level) != 1)
    stop("\tError: The argument of level should be of length 1.\n");
  if(!any(level %in% c(1,2)))
    stop("\tError: The argument of level only accepted the integer 1 or 2. \n");
  #cat(rep("-", times = 50), sep = "");
  if(level == 1){
    cat("Single Environment Analysis:");
    cat("\n");
    for(i in 1:data$trait.number)
    {
      trait.name <- data$traits[[i]]$name;
      cat("The phenotypic trait is ", trait.name, ".\n", sep = "");
      for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
      {
        env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
        cat("*******************Start of ", env.name, " environment********************\n", sep = "");
        #cat("\t\tOn the environment:", env.name, ".\n", sep = "");
        cat("The variance table:\n");
        print(data$traits[[i]]$analysis$sea$envs[[j]]$varcomp.table, row.names = FALSE);
        cat("\n");
        cat("Adjust Means of each genotype:\n");
        print(data$traits[[i]]$analysis$sea$envs[[j]]$sum.out, row.names = FALSE);
        cat("*******************End of ", env.name, " environment********************\n", sep = "");
        cat("\n");
      }#--- end stmt of for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
      cat("\n\n");
    } #--- end stmt of for(i in 1:data$trait.number) ---#
  }else
  {
    for(i in 1:data$trait.number)
    {
      if(is.null(data$traits[[i]]$analysis$sea))
      {
        next;
      } else
      {
        trait.name <- data$traits[[i]]$name;
        cat(rep("=",times=40), sep="");
        cat("\n");
        cat("RESPONSE VARIABLE: ", trait.name, "\n", sep="");
        cat(rep("=",times=40), sep="");
        cat("\n");
        cat("\n");
        for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
        {
          env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
          cat(rep("-",times=40), sep="");
          cat("\n");
          cat("ANALYSIS FOR: Env = ", env.name, "\n", sep="");
          cat(rep("-",times=40), sep="");
          cat("\n");
          cat("\n");
          cat("Data Summary:\n");
          cat("\n");
          cat("Number of observations read: ", data$traits[[i]]$envs[[j]]$obsread, "\n", sep = "");
          cat("Number of observations used: ", data$traits[[i]]$envs[[j]]$obsused, "\n", sep = "");
          print(data$traits[[i]]$analysis$sea$envs[[j]]$factor.summary, row.names=FALSE);
          cat("\n");
          cat("Variance Components Table: \n");
          cat("\n");
          print(data$traits[[i]]$analysis$sea$envs[[j]]$varcomp.table, row.names = FALSE);
          cat("\n\n");
          cat("Testing for the Significance of Genotypic Effect: \n");
          cat("\n");
          print(data$traits[[i]]$analysis$sea$envs[[j]]$m2Vm1, row.names = FALSE);
          cat("\n")
          print(data$traits[[i]]$analysis$sea$envs[[j]]$geno.test,row.names = FALSE);
          cat("\n");
          cat("Genotoype LSMean and Standard Errors:\n");
          cat("\n");
          print(data$traits[[i]]$analysis$sea$envs[[j]]$summary.statistic, row.names = FALSE);
          cat("\n");
          cat("Standard Error of The Difference (SED):\n");
          print(data$traits[[i]]$analysis$sea$envs[[j]]$sedTable, row.names = FALSE);
          cat("\n");
#           #--- contrast outcomes---#
#           if(!is.null(data$traits[[i]]$analysis$sea$envs[[j]]$contrast))
#           {
#             cat(rep("-", times = 40), sep = "");
#             cat("\n");
#             cat("Contrast Analysis\n")
#             cat(rep("-", times = 40), sep = "");
#             cat("\n");
#             if(nrow(data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome) == 0)
#               cat("There are no significant contrasts!\n")
#             else
#               print(data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome, row.names = FALSE);
#             cat("\n");
#           }
        } # end stmt of for(j in 1:1:length(data$traits[[i]]$analysis$sea$envs))
      }# end stmt of if(is.null(data$traits[[i]]$analysis$sea))
    }# end stmt of for(i in 1:data$trait.number)
  }
}