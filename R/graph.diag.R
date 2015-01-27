# graph diagnostic plot
#
#[Arguments]
# data - SingleEnvAnalysis or MultiEnvAnalysis outcomes;
# path - the path for create plot
# single.env - logical, whether plot diagnostic of single environment

graph.diag <- function
(
  data,
  path,
  ...
)
{
  UseMethod("graph.diag");
}

graph.diag.SingleEnvAnalysis <- function
(
  data,
  path,
  ...
)
{
  if(missing(path))
    path <- getwd();
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(is.null(data$traits[[i]]$analysis$sea))
    {
      warning(cat("\tSkip the ", trait.name, " diagnostic plot.\n",sep = ""));
      next;
    } else
    {
      for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
      {
        env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
        xlabel <- trait.name;
        cc = rgb(0.3, 0.3, 0.7, 0.2);
        filename = paste(path, "/diagPlot_",trait.name,"_", env.name, ".png", sep="");
        residValues <- data$traits[[i]]$analysis$sea$envs[[j]]$residuals;
        fittedValues <- data$traits[[i]]$analysis$sea$envs[[j]]$fitted.values;
        if(!is.null(residValues) & !is.null(fittedValues))
        {
          png(filename);
          old.pars <- par();
          par(mfrow = c(2,2));
          #scatterplot of residuals against fitted values;
          plot(residValues ~ fittedValues, xlab = "Predicted Values", pch = 15,
               cex = 0.7, col = cc, ylab = "Residuals", 
               main = "Scatterplot of Residuals \nagainst Fitted values");
          #qqplot of residuals
          qqnorm(residValues);
          qqline(residValues, col = 2, main = title, sub = xlabel);
          #freq dist of residuals
          hist(residValues, main = "Histogram of Residuals", col = cc, 
               xlab = "Residual", ylab = "Frequency");
          #creat blank plot
          plot(seq(1:10) ~ seq(1:10), type = "n", boxed = FALSE, 
               axes = FALSE, xlab = "", ylab = "");
          notestring <- paste("Note: Residuals plotted are taken from the model, response variable ",
                              trait.name, " and environment = ", env.name, ".", sep = "");
          text(5, 7, paste(strwrap(notestring, width = 30), sep = "", collapse = "\n"));
          par(old.pars);
          dev.off();
        }
      }
    }
  }
}

graph.diag.MultiEnvAnalysis <- function
(
  data,
  path,
  ...
)
{
  if(missing(path))
    path <- getwd();
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(is.null(data$traits[[i]]$analysis$mea))
    {
      warning(cat("\tSkip the ", trait.name, " diagnostic plot.\n",sep = ""));
      next;
    } else
    {
      is.random <- FALSE; #genotype must be fixed
      xlabel = trait.name;
      cc = rgb(0.3, 0.3, 0.7, 0.2);
      
      if (is.random) {
        filename1 = paste(getwd(),"/diagPlotsMea1s_",trait.name, "_random.png",sep="");
      } else
      {
        filename1 = paste(getwd(),"/diagPlotsMea1s_",trait.name, "_fixed.png",sep="");
      }
      residValues<-data$traits[[i]]$analysis$mea$residuals;
      fittedValues<-data$traits[[1]]$analysis$mea$fitted.values;
      
      if (!is.null(residValues) & (!is.null(fittedValues))) {
        png(filename = filename1);
        par(mfrow = c(2,2)) #par(mfrow = n2mfrow(length(respvar)*nlevels(data[,match(env, names(data))])));
        
        #scatterplot of residuals against fitted values
        plot(residValues~fittedValues, xlab = "Predicted Values", pch = 15, 
             cex = 0.7, col = cc, ylab = "Residuals", 
             main = "Scatterplot of Residuals \nagainst Fitted Values");
        
        #qqplot of residuals
        qqnorm(residValues); 
        qqline(residValues, col=2, main=title,sub=xlabel);
        
        #freq dist of residuals
        hist(residValues, main = "Histogram of Residuals",  
             col = cc, xlab = "Residual", ylab = 'Frequency' );
        
        #create blank plot
        #plot(seq(1:10)~seq(1:10), type="n", boxed = FALSE, axes=FALSE, xlab="", ylab="")
        plot(seq(1:10)~seq(1:10), type="n", bty = "n", axes=FALSE, xlab="", ylab="");
        if (is.random) {
          noteString <- paste("NOTE: Residuals plotted are taken from the model where genotype is random and response variable = ", 
                              trait.name, ".", sep="")
        } else {
          noteString <- paste("NOTE: Residuals plotted are taken from the model where genotype is fixed and response variable = ", 
                              trait.name, ".", sep="")
        }
        text(5,7,paste(strwrap(noteString,width=30), sep="", collapse="\n"))
        
        dev.off();
      }
    }
  }#end stmt for(i in 1:length(data$traits))
}
