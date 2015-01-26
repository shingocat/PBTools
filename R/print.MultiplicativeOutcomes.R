#This is funciton is design for webapp used
print.MultiplicativeOutcomes <- function(
  data,
  ...
)
{
  cat(rep("=",40),sep="");
  cat("\n");
  cat("Multiplicative Analysis:\n")
  cat(rep("=",40),sep="");
  cat("\n");
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(is.null(data$traits[[i]]$analysis$mea))
    {
      warning(paste("\tWarnings: There are no Multi-Environment Analysis on "
                    ," this trait ", trait.name, ".\n", sep=""));
      next;
    }
    env.number <- data$traits[[i]]$analysis$mea$nlevelsEnv;
    geno.number <-  data$traits[[i]]$analysis$mea$nlevelsGeno;
    if(is.null(data$traits[[i]]$analysis$mea$ammi))
    {
      warning(paste("\tWarnings: There are no AMMI Analysis on "
                    ," this trait ", trait.name, ".\n", sep=""));
      next;
    } else
    {
      cat(rep("-",40),sep="");
      cat("\n");
      cat("AMMI Analysis:\n")
      cat("Trait Name: ", trait.name, "\n", sep="");
      cat("Number of Genotype: ", geno.number, "\n", sep="");
      cat("Number of Environment: ", env.number, "\n", sep="");
      cat(rep("-",40),sep="");
      cat("\n");
      if(is.null(data$traits[[i]]$analysis$mea$ammi$analysis))
      { cat(paste("\tWarnings: There are no AMMI Analysis outcomes on "
                  ," this trait ", trait.name, ".\n", sep=""));
      }else
      {
        print(data$traits[[i]]$analysis$mea$ammi$analysis);
      }
      cat("\n");
    }
    if(is.null(data$traits[[i]]$analysis$mea$gge))
    {
      warning(paste("\tWarnings: There are no GGE Analysis on "
                    ," this trait ", trait.name, ".\n", sep=""));
      next;
    } else
    {
      cat(rep("-",40),sep="");
      cat("\n");
      cat("GGE Analysis:\n")
      cat("Trait Name: ", trait.name, "\n", sep="");
      cat("Number of Genotype: ", geno.number, "\n", sep="");
      cat("Number of Environment: ", env.number, "\n", sep="");
      cat(rep("-",40),sep="");
      cat("\n");
      if(is.null(data$traits[[i]]$analysis$mea$gge$analysis))
      {  
        cat(paste("\tWarnings: There are no GGE Analysis Outcomes on "
                  ," this trait ", trait.name, ".\n", sep=""));
      }else
      {  
        print(data$traits[[i]]$analysis$mea$gge$analysis);
      }
      cat("\n");
    }
    
  }#end stmt  for(i in 1:length(data))
}