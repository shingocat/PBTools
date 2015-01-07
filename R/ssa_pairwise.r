# -------------------------------------------------------
# PAIRWISE MEAN COMPARISONS AFTER SINGLE SITE ANALYSIS 
# File Created by: Alaine A. Gulles 04.04.2013
# Script Modified by: Rose Imee Zhella A. Morantte
#                     Nellwyn L. Sales
# --------------------------------------------------------

# --------------------------------------------------------
# ARGUMENTS:
# model - model from ssa.test
# type - "Tukey" or "Dunnett"
# alpha - level of significance; default is 0.05
# control.level - vector of genotype levels considered as controls;
#               - if not NULL, comparison with control is performed
# --------------------------------------------------------

ssa.pairwise <- function(model, type = "Tukey", alpha = 0.05, control.level = NULL) UseMethod("ssa.pairwise")

ssa.pairwise.default <- function(model, type = "Tukey", alpha = 0.05, control.level = NULL) {
  
	library(lme4)
	library(multcomp) 
	temp <- NULL
	trmt.levels <- levels(model@frame[,2])
	n <- c(1:length(trmt.levels))
	names(n) <- trmt.levels
  
	# --- check if items in control.levels are in trmt.levels, if not adjust control.level and create warnings --- #
	controlTestWarning="NONE"
	if (!is.null(control.level)) {
	  control.level<-trimStrings(control.level)
	  canProceedToAnalysis<-all(control.level %in% trmt.levels)
	  if (!canProceedToAnalysis) {
	    isControlPresent<-control.level %in% trmt.levels
	    deletedControl<-control.level[!isControlPresent]
	    
	    if (length(deletedControl)==length(control.level)) {
	      stop("The specified control level(s) is(are) not present in this trial. Cannot perform comparison with control.") 
	    } else {
	      control.level<-control.level[isControlPresent]
	      deletedList=NULL
	      for (k in 1:length(deletedControl)) {
	        if (k == 1) {
	          deletedList<-paste(deletedControl[k])
	        } else {
	          deletedList<-paste(deletedList, ", ", deletedControl[k], sep="")
	        }
	      }
	      controlTestWarning <-paste("The following control level(s) is(are) not present in this trial and deleted from the list of control levels: ", deletedList, sep="")
	    }
	  }
	}
  
   if (type == "Dunnett") {
        for (z in (1:length(control.level))) {
             contrast1 <- contrMat(n, type = "Dunnett", base = as.numeric(match(control.level[z], names(n))))
             mc1 <- eval(parse(text = paste("glht(model, linfct = mcp(", names(model@frame)[2]," = contrast1))", sep = "")))
             #mc1 <- glht(model, linfct = contrast1)
             interval <- confint(mc1, level = 1 - alpha)
             interval.confint <- as.data.frame(interval$confint)
             signif2 <- subset(interval.confint, as.logical(lwr <= 0 & 0 <= upr) == F)
             newtemp <- NULL
             if (nrow(signif2) != 0) {
                  trmt.comp <- strsplit(rownames(signif2), split = " - ")
                  temp <- rbind(temp,data.frame(t(data.frame(trmt.comp)), signif2, row.names = NULL))
             } 
        }
        
        if (!is.null(temp)) {
            rownames(temp) <- 1:nrow(temp) 
            colnames(temp)[1:5] <- c("Trmt[i]", "Trmt[j]", "Difference", "Lower", "Upper")
            return(list(result=temp, controlTestWarning=controlTestWarning))
        } else {
             temp <- as.matrix("(No significant pairwise comparisons.)")
             row.names(temp) <- ""
             return(list(result=noquote(temp[1,1]), controlTestWarning=controlTestWarning))
        }
   } else {
        # type == "Tukey"
        contrast1 <- contrMat(n, type = "Tukey")
        mc1 <- eval(parse(text = paste("glht(model, linfct = mcp(", names(model@frame)[2]," = contrast1))", sep = "")))
        #mc1 <- glht(model, linfct = contrast1)
        interval <- confint(mc1, level = 1 - alpha)
        interval.confint <- as.data.frame(interval$confint)
        signif2 <- subset(interval.confint, as.logical(lwr <= 0 & 0 <= upr) == F)
        if (nrow(signif2) != 0) {
             trmt.comp <- strsplit(rownames(signif2), split = " - ")
             for (i in (1:nrow(signif2))) temp <- rbind(temp, trmt.comp[[i]])
             temp <- data.frame(temp, signif2, row.names = NULL)
             colnames(temp)[1:3] <- c("Trmt[i]", "Trmt[j]", "Difference", "Lower", "Upper")
             return(list(result=temp, controlTestWarning=controlTestWarning))
             #	} else {return(temp <- signif2)}
        } else {
             temp <- as.matrix("(No significant pairwise comparisons.)")
             row.names(temp) <- ""
             return(list(result=noquote(temp[1,1]), controlTestWarning=controlTestWarning))
        }
  }
  detach("package:multcomp")
	detach("package:lme4")	
}
