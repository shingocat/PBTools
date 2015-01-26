# This is designed for Webapp to read contrast files and build a list 
#
#[Arguments]
# files - a string vector containing all the file absoluted path
#
#

getContrastListFromFiles <- function
(
  files,
  ...
)
{
  UseMethod("getContrastListFromFiles");
}

getContrastListFromFiles.default <- function
(
files,
...
)
{
  if(missing(files))
    stop("The argument files could not be null!");
  outcomes <- list();
  for(i in 1:length(files))
  {
    contrast <- read.csv(file=files[i], header=TRUE, row.names =1);
    if(all(is.na(contrast[ncol(contrast)])))
      contrast <- contrast[, -ncol(contrast)];
    labels <- rownames(contrast);
    if(nrow(contrast) == 1)
    {
      contrast <- apply(contrast, 2, as.character);
      contrast <- t(as.data.frame(contrast));
      contrast <- apply(contrast, 2, as.numeric);
      contrast <- t(as.data.frame(contrast))
      contrast <- as.matrix(contrast);
    } else
    {
      contrast <- apply(contrast, 2, as.character);
      contrast <- apply(contrast, 2, as.numeric);
      contrast <- as.matrix(contrast);
    }
    rownames(contrast) <- labels;
    outcomes[[i]] <- contrast;
  }
  return(outcomes);
}  
