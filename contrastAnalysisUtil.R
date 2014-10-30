# ------------------------------------------------------------------------------
# File and Script Created by: Alaine A. Gulles 08.07.2014
#                             for International Rice Research Institute
# ------------------------------------------------------------------------------
# Description: Utilities for contrast analysis
# Arguments: 
# model - a lmerMod class object
# contrast - a contrast matrix
# control - control level
# alpha - level of significance
# Returned Value: a data frame
# ------------------------------------------------------------------------------

contrastCompute <- function(model, contrast, alpha) UseMethod("contrastCompute")

contrastCompute.default <- function(model, contrast, alpha) {
     mc1 <- eval(parse(text = paste("glht(model, linfct = mcp(", names(model@frame)[2]," = contrast))", sep = "")))
     interval <- confint(mc1, level = 1 - alpha)
     intervalConfint <- as.data.frame(interval$confint)
     result <- data.frame(rownames(intervalConfint), intervalConfint)
     colnames(result) <- c("Contrast", "Difference", "Lower", "Upper")
     rownames(result) <- NULL
     return(result)    
}

userDefContrast <- function(model, contrast, alpha = 0.05) UseMethod("userDefContrast")

userDefContrast.default <- function(model, contrast, alpha = 0.05) {
     trmtLevels <- levels(model@frame[,2])
     contrMatName <- colnames(contrast)
     numLevels <- length(trmtLevels)
     numContrast <- nrow(contrast)
     
     # determine if the number of columns in the contrast matrix is equal to the levels of the treatment
     if (numLevels != ncol(contrast)) { stop("Error: The number of columns in the contrast matrix is not equal to the levels of the treatment.") }
     # determine if the labels in the contrast matrix is the same as the data set
     if (any(contrMatName %in% trmtLevels == FALSE)) { stop("Error: At least one of the column names in the contrast matrix does not match the levels in the data set.") }
     # determine if the number of contrast in the contrast matrix is at most treatment df
     if (numContrast >= numLevels) { stop("Error: Too many rows in the contrast matrix.") }
     # determine if the matrix is a contrast matrix
     if (!isContrast(contrast)) { stop("Error: At least one of the rows in the matrix is not a contrast.") }
     
     # re-arrange the column of the contrast to follow the levels and order of the dataset
     newContrast <- contrast[,match(trmtLevels,contrMatName)]
     
     # determine if the matrix is an orthogonal contrast
     #if (numContrast <= 2) {
     #     if (!isOrthogonal(contrast)) { stop("Error: At least one pair of contrast is not orthogonal.") }
     #}
     
     # 
     result <- contrastCompute(model, newContrast, alpha)
     return(result)     
}

compareControlContrast <- function(model, control, alpha = 0.05) UseMethod("compareControlContrast")

compareControlContrast.default <- function(model, control, alpha = 0.05) {
     temp <- NULL
     trmtLevels <- levels(model@frame[,2])
     tempLevels <- 1:length(trmtLevels)
     names(tempLevels) <- trmtLevels
     numLevels <- length(trmtLevels)
     
     contrast1 <- contrMat(tempLevels, type = "Dunnett", base = as.numeric(match(control, names(tempLevels))))     
     result <- contrastCompute(model, contrast1, alpha)
     return(result)
}

