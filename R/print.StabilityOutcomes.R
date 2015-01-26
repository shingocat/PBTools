# This function is used for webapp to print output.
# data-should be of class MultiEnvAnalysis
print.StabilityOutcomes <- function
(
  data,
  ...
)
{
  if(!inherits(data,"MultiEnvAnalysis"))
    stop("\tError: The data argument should be of class MultiEnvAnalysis!\n");
  cat(rep("=",40),sep="");
  cat("\n");
  cat("Stability Analysis:\n")
  cat(rep("=",40),sep="");
  cat("\n");
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    env.number <- data$traits[[i]]$analysis$mea$nlevelsEnv;
    cat(rep("-",40),sep="");
    cat("\n");
    cat("Trait Name: ", trait.name, "\n", sep="");
    cat("Number of Environment: ", env.number, "\n", sep="");
    cat(rep("-",40),sep="");
    cat("\n");
    if(is.null(data$traits[[i]]$analysis$mea$stability))
    {
      cat("***There are no any outcomes of stability analysis!\n");
      next;
    } else
    {
      if(!is.null(data$traits[[i]]$analysis$mea$stability$finlay.wilkinson))
      {
        cat(data$traits[[i]]$analysis$mea$stability$finlay.wilkinson);
        cat("\n");
        print(data$traits[[i]]$analysis$mea$stability$slope, row.names=TRUE);
        cat("\n\n")
      }
      if(!is.null(data$traits[[i]]$analysis$mea$stability$shukla))
      {
        cat(data$traits[[i]]$analysis$mea$stability$shukla);
        cat("\n");
        print(data$traits[[i]]$analysis$mea$stability$par, row.names=TRUE);
        cat("\n\n");
      }
    }
  }
}