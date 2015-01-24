# this print method is to design to using pbtools on webapp

print.ContrastOutcomes <- function(
  data,
  ...
)
{
  UseMethod("print.ContrastOutcomes");
}
print.ContrastOutcomes.SingleEnvAnalysis <- function
(
  data,
  ...
)
{
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(!is.null(data$traits[[i]]$analysis$sea))
    {
      cat(rep("=", times = 40), sep = "");
      cat("\n");
      cat("Contrast Analysis\n")
      cat(rep("=", times = 40), sep = "");
      cat("\n");
      for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
      {
        env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
        #--- contrast outcomes---#
        if(!is.null(data$traits[[i]]$analysis$sea$envs[[j]]$contrast))
        {
          cat("Environment: ", env.name, "\n",sep ="");
          contrast.type = data$traits[[i]]$analysis$sea$envs[[j]]$contrast$type;
          if(contrast.type == "RecurrentParent")
          {
            cat("CONTRAST TYPE : Comparing With Recurrent Parent.\n" );
            cat("Recurrent: ",data$traits[[i]]$analysis$sea$envs[[j]]$contrast$recurrent, sep="");
            cat("\n");
            cat("Alpha: ",data$traits[[i]]$analysis$sea$envs[[j]]$contrast$alpha, sep="");
            cat("\n");
          } else if(contrast.type == "Custom")
          {
            cat("CONTRAST TYPE : User Specified Contrast.\n" );
            cat("\n");
          } else if(contrast.type == "Default")
          {
            cat("CONTRAST TYPE : Default Pyramided Lines Contrast.\n" );
            cat("Gene Number: ",data$traits[[i]]$analysis$sea$envs[[j]]$contrast$gene.number,sep="");
            cat("\n");
          }
          if(nrow(data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome) == 0)
            cat("There are no significant contrasts!\n")
          else
            print(data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome, row.names = FALSE);
          cat("\n");
        } else
        {
          cat("NO CONTRAST OUTCOMES ON THIS TRAIT ", trait.name ,"on ", env.name," Environment.", sep="");
          cat("\n");
        }
        cat(rep("-", times = 40), sep = "");
        cat("\n");
      } # end stmt of for(j in 1:1:length(data$traits[[i]]$analysis$sea$envs))
    }# end stmt of if(is.null(data$traits[[i]]$analysis$sea))
  }# end stmt of for(i in 1:data$trait.number)
}

print.ContrastOutcomes.MultiEnvAnalysis <- function
(
  data,
  ...
)
{
}







