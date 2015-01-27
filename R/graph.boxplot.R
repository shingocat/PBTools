# for graphing the SingleEnvAnalysis and MultiEnvAnalysis boxplot
#
# [Arguments]
#   data - SingleEnvAnalysis or MultiEnvAnalysis Outcome;
#   path - the path to create boxplot file
#   single.env - logical, whether include all environment under this trait.
#
graph.boxplot <- function
(
  data,
  path,
  single.env = FALSE,
  ...
)
{
  UseMethod("graph.boxplot");
}
graph.boxplot.SingleEnvAnalysis <- function
(
  data,
  path,
  single.env = FALSE,
  ...
)
{
  if(missing(path))
    path <- getwd();
  #create boxplot of traits after SingleEnvAnalysis on each environment.
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(is.null(data$traits[[i]]$analysis$sea))
    {
      warning(cat("\tSkip the ", trait.name, " boxplot\n",sep = ""));
      next;
    } else
    {
      if(single.env)
      {
        for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
        {
          env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
          boxfile <- paste(path,"/boxplot_",trait.name,"_",env.name,".png",sep="");
          if(!all(is.na(data$traits[[i]]$envs[[j]]$data[,trait.name])))
          {
            png(boxfile);
            xlabel = trait.name;
            boxplot(as.numeric(as.character(data$traits[[i]]$envs[[j]]$data[,trait.name])),
                    xlab = xlabel, main = paste("Boxplot of ", trait.name, sep=""));
            dev.off();
          }
        }
      } else
      {
        env.label <- data$traits[[i]]$envs[[1]]$design$env;
        boxfile <- paste(path,"/boxplot_", trait.name,"_ALL_Env.png", sep= ""); 
        if(!all(is.na(data$raw.data[,trait.name])))
        {
          png(boxfile);
          xlabel = trait.name;
          boxplot(as.numeric(as.character(data$raw.data[,trait.name])) ~ as.factor(data$raw.data[,env.label]), 
                  data = data$raw.data, xlab = xlabel, main = paste("Boxplot of ", trait.name, sep=""));
          dev.off();
        }
      }
    }
  }
}

graph.boxplot.MultiEnvAnalysis <- function
(
  data,
  path,
  single.env = FALSE,
  ...
)
{
  if(missing(path))
    path <- getwd();
  #create boxplot of traits after SingleEnvAnalysis on each environment.
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(is.null(data$traits[[i]]$analysis$mea))
    {
      warning(cat("\tSkip the ", trait.name, " boxplot\n",sep = ""));
      next;
    } else
    {
      boxfile = paste(getwd(),"/boxplotMea1S_",trait.name,".png",sep = "");
      if (!all(is.na(data$traits[[i]]$analysis$mea$data[,trait.name]))) {
        png(filename = boxfile); #par(mfrow = n2mfrow(length(respvar)));
        xlabel = trait.name;
        boxplot((data$traits[[i]]$analysis$mea$data[,trait.name]), data = data,
                xlab = xlabel, main = paste("Boxplot of ", trait.name, sep=""));
        dev.off()
      }
    }
  }#end stmt  for(i in 1:length(data$traits))
}
