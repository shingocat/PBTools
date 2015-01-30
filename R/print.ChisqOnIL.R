print.ChisqOnIL <- function
(
  data, 
  level = 1,
  ...
)
{
  if(level == 1)
  {
    cat(rep("=", 40), sep="");
    cat("\n");
    cat("Chisq Test:\n");
    cat(rep("=", 40), sep="");
    cat("\n");
    cat(rep("-",40), sep="");
    cat("\n");
    if(!is.null(data$TypeI))
    {
      cat(data$TypeI$Messages);
      cat("\n");
      cat("Type I Markers: ", data$TypeI$counts, "\n", sep="");
      if(nrow(data$TypeI$outcomes) == 0)
      {
        cat("There are no records!\n");
      }
      else
      {
        print(data$TypeI$outcomes, row.names=FALSE);
      }
    }
    cat(rep("-",40), sep="");
    cat("\n");
    if(!is.null(data$TypeII))
    {
      cat(data$TypeII$Messages);
      cat("\n");
      cat("Type II Markers: ", data$TypeII$counts, "\n", sep="");
      if(nrow(data$TypeII$outcomes) == 0)
      {cat("There are no records!\n");}
      else
      {print(data$TypeII$outcomes, row.names=FALSE);}
    }
    cat(rep("-",40), sep="");
    cat("\n");
    if(!is.null(data$TypeIII))
    {
      cat(data$TypeIII$Messages);
      cat("\n");
      cat("Type III Markers: ", data$TypeIII$counts, "\n", sep="");
      if(nrow(data$TypeIII$outcomes) == 0)
      {cat("There are no records!\n");}
      else
      {print(data$TypeIII$outcomes, row.names=FALSE);}
    }
    cat(rep("-",40), sep="");
    cat("\n");
    if(!is.null(data$TypeIV))
    {
      cat(data$TypeIV$Messages);
      cat("\n");
      cat("Type IV Markers: ", data$TypeIV$counts, "\n", sep="");
      if(nrow(data$TypeIV$outcomes) == 0)
      {cat("There are no records!\n");}
      else
      {print(data$TypeIV$outcomes, row.names=FALSE);}
    }
    cat(rep("-",40), sep="");
    cat("\n");
    if(!is.null(data$TypeV))
    {
      cat(data$TypeV$Messages);
      cat("\n");
      cat("Type V Markers: ", data$TypeV$counts, "\n", sep="");
      if(nrow(data$TypeV$outcomes) == 0)
      {cat("There are no records!\n");}
      else
      {print(data$TypeV$outcomes, row.names=FALSE);}
    }
    cat(rep("-",40), sep="");
    cat("\n");
  } else
  {
    
  }
}