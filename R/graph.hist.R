# graph histogram plot
#
# [Arguments]
#   data - data of SingleEnvAnalysis or MultiEnvAnalysis
#   path - path to create plot;
#   single.env - logical, whether plot hist on single environment.


graph.hist <- function
(
  data,
  path,
  single.env = TRUE,
  ...
)
{
  UseMethod("graph.hist");
}

graph.hist.SingleEnvAnalysis <- function
(
  data,
  path,
  single.env = TRUE,
  ...
)
{
  if(missing(path))
    path <- getwd();
  #creat hist plot of trait after SingleEnvAnalysis
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(is.null(data$traits[[i]]$analysis$sea))
    {
      warning(cat("\tSkip the ", trait.name, " histogram plot.\n",sep = ""));
      next;
    } else
    {
      if(single.env)
      {
        for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
        {
          env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
          histfile <- paste(path,"/histplot_",trait.name,"_",env.name,".png",sep="");
          if(!all(is.na(data$traits[[i]]$envs[[j]]$data[,trait.name])))
          {
            png(histfile);
            xlabel = trait.name;
            hist(as.numeric(data$traits[[i]]$envs[[j]]$data[,trait.name]), xlab = xlabel, 
                    main = paste("Histogram of ", trait.name, sep=""));
            dev.off();
          }
        }
      } else
      {
        env.label <- data$traits[[i]]$envs[[1]]$design$env;
        histfile <- paste(path,"/histplot_", trait.name,"_ALL_Env.png", sep= ""); 
        if(!all(is.na(data$raw.data[,trait.name])))
        {
          png(histfile);
          xlabel = trait.name;
          hist(as.numeric(data$raw.data[,trait.name]), 
                 xlab = xlabel, main = paste("Histogram of ", trait.name, sep=""));
          dev.off();
        }
      }
    }
  }
}