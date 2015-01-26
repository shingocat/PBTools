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
  cat(rep("=", times = 40), sep = "");
  cat("\n");
  cat("Contrast Analysis\n")
  cat(rep("=", times = 40), sep = "");
  cat("\n");
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(!is.null(data$traits[[i]]$analysis$sea))
    {
      for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
      {
        env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
        #--- contrast outcomes---#
        if(!is.null(data$traits[[i]]$analysis$sea$envs[[j]]$contrast))
        {
          cat(rep("-", times = 40), sep = "");
          cat("\n");
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
          cat(rep("-", times = 40), sep = "");
          cat("\n");
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
  cat(rep("=", times = 40), sep = "");
  cat("\n");
  cat("Contrast Analysis\n")
  cat(rep("=", times = 40), sep = "");
  cat("\n");
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(!is.null(data$traits[[i]]$analysis$mea))
    {
        #--- contrast outcomes---#
        if(!is.null(data$traits[[i]]$analysis$mea$contrast))
        {
          cat(rep("-", times = 40), sep = "");
          cat("\n");
          contrast.type = data$traits[[i]]$analysis$mea$contrast$type;
          if(contrast.type == "RecurrentParent")
          {
            cat("CONTRAST TYPE : Comparing With Recurrent Parent.\n" );
            cat("Recurrent: ",data$tratis[[i]]$analysis$mea$contrast$recurrent, sep="");
            cat("\n");
          } else if(contrast.type == "Custom")
          {
            cat("CONTRAST TYPE : User Specified Contrast.\n" );
            cat("\n");
          } else if(contrast.type == "Default")
          {
            cat("CONTRAST TYPE : Default Pyramided Lines Contrast.\n" );
            cat("Gene Number: ",data$traits[[i]]$analysis$mea$contrast$gene.number,sep="");
            cat("\n");
          }
          cat(rep("-", times = 40), sep = "");
          cat("\n");
          cat("Contrast On Genotype:\n");
          if(nrow(data$traits[[i]]$analysis$mea$contrast$contrastOnGeno) == 0)
            cat("There are no significant contrasts!\n")
          else
            print(data$traits[[i]]$analysis$mea$contrast$contrastOnGeno, row.names = FALSE);
          cat("\n");
          cat("Contrast On Environment:\n");
          if(nrow(data$traits[[i]]$analysis$mea$contrast$contrastAcrossEnv) == 0)
            cat("There are no significant contrasts!\n")
          else
            print(data$traits[[i]]$analysis$mea$contrast$contrastAcrossEnv, row.names = FALSE);
          cat("\n");
        } else
        {
          cat("NO CONTRAST OUTCOMES ON THIS TRAIT ", trait.name ,"on ", env.name," Environment.", sep="");
          cat("\n");
        }
        cat(rep("-", times = 40), sep = "");
        cat("\n");
    }# end stmt of if(is.null(data$traits[[i]]$analysis$mea))
  }# end stmt of for(i in 1:data$trait.number)
}







