print.MultiMarkerAnalysis <- function(
  data, 
  p.value = 0.05, 
  level = 1)
{
  if(!inherits(data, "MultiMarkerAnalysis"))
    stop("\tError: The argument of data should be of class MultiMarkerAnalysis!\n");
  if(!is.numeric(p.value))
    stop("\tError: The argument of p.value should be of value between 0 and 1.\n");
  if(p.value > 1 || p.value < 0)
    stop("\tError: The argument of p.value should be of value between 0 and 1.\n");
  if(!level %in% c(1,2,3))
    stop("\tError: The argument of level should be one of value 1, 2 and 3.\t");
  cat("Multiple Marker Analysis (p.value < ", p.value, "):\n", sep = "");
  
  #library(hdlm);
  for(i in 1:length(data$traits))
  {
    if(is.null(data$traits[[i]]$analysis$mma))
      next;
    for(j in 1:length(data$traits[[i]]$analysis$mma$envs))
    {
      cat(rep("=", times=40), sep="");
      cat("\n");
      cat("Trait ", data$traits[[i]]$name, sep = "");
      cat("\n")
      cat("Environment ", data$traits[[i]]$analysis$mma$envs[[j]]$name, sep = "");
      cat("\n");
      cat(rep("=", times = 40), sep="");
      cat("\n");
      outcomes <- data$traits[[i]]$analysis$mma$envs[[j]]$outcomes;
      if(identical(class(outcomes), "try-error"))
      {
        cat("It could not fit the model with multiple markers on this trait.\n");
      }else
      {
        out.summary <- summary(outcomes, level = level);
        coef <- out.summary$coefficients;
        lower.bound <- out.summary$lower.bound;
        upper.bound <- out.summary$upper.bound;
        PValue <- out.summary$p.value;
        temp <- data.frame(Marker = names(coef), Estimate = coef, CI.Lower = lower.bound, CI.Upper = upper.bound, P.value = PValue);
        #--- remove intercept in the temp data frame ---#
        temp <- temp[-1,];
        temp <- temp[which(temp$P.value <= p.value),];
        if(nrow(temp) == 0)
        {
          cat("There are no significant markers!");
        } else
        {
          print(, row.names = FALSE);
        }
      }
      cat("\n");
      cat(rep("-", times=40), sep="");
      cat("\n");
    }
  }
  #detach("package:hdlm", unload = TRUE);
}