# ------------------------------------------------------------------------------
# File and Script Created by: Alaine A. Gulles 08.07.2014
#                             for International Rice Research Institute
# ------------------------------------------------------------------------------
# Description: Perform contrast analysis
# Arguments: 
# model - a lmerMod class object
# contrastOpt - contrast option; c("control", "user")
# contrast - NULL or a contrast matrix; required if contrastOpt == "user"
# controlLevel - NULL or control level
# alpha - level of significance
# Returned Value: a data frame
# ------------------------------------------------------------------------------

contrastAnalysis <- function(model, contrastOpt = c("control", "user"), controlLevel = NULL, contrast = NULL, alpha = 0.05) UseMethod("contrastAnalysis")

contrastAnalysis.default <- function(model, contrastOpt = c("control", "user"), controlLevel = NULL, contrast = NULL, alpha = 0.05) {
     contrastOpt <- match.arg(contrastOpt)
     temp <- NULL
     trmtLevels <- levels(model@frame[,2])
     
     library(multcomp)
     if (contrastOpt == "control") {
          if (is.null(controlLevel)) { controlLevel <- trmtLevels[1]
          } else {
               if (any(controlLevel %in% trmtLevels)) { controlLevel <- trmtLevels[na.exclude(match(controlLevel, trmtLevels))]
               } else { controlLevel <- trmtLevels[1] }
          }
          for (i in (1:length(controlLevel))) {
               temp <- rbind(temp, compareControlContrast(model, control = controlLevel[i], alpha))
          }
     } else {
          # perform user specified contrast
          if (is.null(contrast)) { stop("Error: Missing contrast matrix.") }
          if (is.data.frame(contrast)) {
               if (ncol(contrast) == length(trmtLevels) + 1) {
                    contrastLabel <- contrast[,1]
                    contrast <- as.matrix(contrast[,2:ncol(contrast)])
                    rownames(contrast) <- contrastLabel
               } else {
                    if (ncol(contrast) < length(trmtLevels)) { stop("Error: Contrast matrix contain too few columns.") 
                    } else {
                         if (ncol(contrast) > length(trmtLevels)) { stop("Error: Contrast matrix contain too fee columns.") 
                         } else { contrast <- as.matrix(contrast) }
                    }
               } 
          }
          if (!is.matrix(contrast)) { stop("Error: contrast should be a matrix.") }
          temp <- userDefContrast(model, contrast, alpha)     
     }
     
     signif <- temp[temp[,"Lower"] <= 0 & temp[,"Upper"] >= 0,]
     
     if (nrow(signif) != 0) {
          rownames(signif) <- 1:nrow(signif)
          cat("SIGNIFICANT LINEAR CONTRAST AT ALPHA = ", alpha, ":\n", sep = "")
          print(signif)
     } else { cat("At alpha = ", alpha, ", no significant contrast.\n", sep = "") }
     
     return(signif)
     detach("package:multcomp")
}
