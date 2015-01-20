###############################################################################
# TODO: Do single envrironment analysis for raw phenotypic data.
# [Arguments]
#	data - an object of PhenotypicData and it should be the raw data type!
# 
# Author: mqin
# Date: Nov 14, 2014
# FileName: doSSA.R
###############################################################################


doSEA <- function(
  data
)
{
  UseMethod("doSEA");
}
doSEA.PhenotypicData <- function(
  data
)
{
  if(!inherits(data, "PhenotypicData"))
    stop("\tError: The argument of data must be of class PhenotypicData!\n");
  if(data$isMean)
    stop("\tError: The phenotypic data is mean type! It could not process single site analysis!\n");
  if(!data$isRestricted)
  {
    warning(paste("\tWarning: It will use the default parameter to restrict phenotypic data;\n",
                  "\tMissing rate should be not larger than 0.2 of each trait on each environment!\n",
                  #"\tAnd used the grand mean of each trait on each environment to instead of missing all obseravtion!\n",
                  sep=""));
    data <- restrict.pheno.data(data);
  }
  
  options(show.signif.stars=FALSE);
  #	library(lme4);
  
  for(i in (1:length(data$traits)))
  {
    #--- saving all analysis output to this list.---#
    if(is.null(data$traits[[i]]$analysis))
      data$traits[[i]]$analysis <- list();
    #--- saving single environment analysis output on this list---#
    #--- if there are former outcomes, deleted these!---#
    if(is.null(data$traits[[i]]$analysis$sea))
    {
      data$traits[[i]]$analysis$sea <- list();
    } else
    {
      data$traits[[i]]$analysis$sea <- NULL;
      data$traits[[i]]$analysis$sea <- list();
    }
    data$traits[[i]]$analysis$sea$envs <- list();
    trait.name <- data$traits[[i]]$name;
    #--- reformating this after all conditions are met, it will be kept the outcomes---#
    for(j in (1:length(data$traits[[i]]$envs)))
    {
      #--- checking whether this site is restricted
      if(data$traits[[i]]$envs[[j]]$restricted$isTRUE)
      {
        warning(paste("\tWarning: Too many missing observations. Cannot proceed with analysis on ", env.name, " of ",
                      trait.name, ".\n", sep = ""));
        next;
      }
      
      #--- retrive all design factor name and experimental design---#
      exptl.design <- data$traits[[i]]$envs[[j]]$design$exptl.design;
      geno <- data$traits[[i]]$envs[[j]]$design$geno;
      env <- data$traits[[i]]$envs[[j]]$design$env;
      env.name <- data$traits[[i]]$envs[[j]]$name;
      temp.data <- data$traits[[i]]$envs[[j]]$data;		
      
      #--- construct model, geno as fixed factor ---#
      if(exptl.design == "RCB" || exptl.design == "AugRCB")
      {
        block <- data$traits[[i]]$envs[[j]]$design$block;
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + (1|",block, ")", sep="");
      } else if(exptl.design == "AugLS")
      {
        row <- data$traits[[i]]$envs[[j]]$design$row;
        column <- data$traits[[i]]$envs[[j]]$design$column;
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + (1|", row , "), + (1|", column, ")", sep ="");
      } else if(exptl.design == "Alpha")
      {
        block <- data$traits[[i]]$envs[[j]]$design$block;
        rep <- data$traits[[i]]$envs[[j]]$design$rep;
        myformula1 <- paste(trait.name, " ~ 1 + ", geno," + (1|", rep,"/", block,")", sep = "");
      } else if(exptl.design == "RowCol")
      {
        row <- data$traits[[i]]$envs[[j]]$design$row;
        column <- data$traits[[i]]$envs[[j]]$design$column;
        rep <- data$traits[[i]]$envs[[j]]$design$rep;
        myformula1 <- paste(trait.name, " ~ 1 + ", geno," + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,")", sep = "");
      } else if(exptl.design == "LatinAlpha")
      {
        block <- data$traits[[i]]$envs[[j]]$design$block;
        rep <- data$traits[[i]]$envs[[j]]$design$rep;
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + (1|", rep,") + (1|", block,") + (1|", rep, ":", block, ")", sep = "");
      } else if(exptl.design == "LatinRowCol")
      {
        row <- data$traits[[i]]$envs[[j]]$design$row;
        column <- data$traits[[i]]$envs[[j]]$design$column;
        rep <- data$traits[[i]]$envs[[j]]$design$rep;
        longerRow <- data$traits[[i]]$envs[[j]]$design$longerRow;
        if(longerRow)
        {
          myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + (1|", rep,") + (1|", column, ") + (1|", rep,":", column,") + (1|", row, ")", sep = ""); 
        } else
        {
          myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + (1|", rep,") + (1|", row, ") + (1|", rep,":", row,") + (1|", column, ")", sep = "") ;
        }	
      }
      
      #--- call lmer function using myformula1 to get variance components table ---#
      #--- model <- lmer(formula(myformula1), data = temp.data) ---#
      model <- try(lmer(formula(myformula1), data = temp.data), silent=TRUE);
      
      #--- if there is some problem ocurred, do next site---#
      if (!is.null(model) && class(model)=="try-error") {  
        msg <- trimStrings(strsplit(model, ":")[[1]]);
        msg <- trimStrings(paste(strsplit(msg, "\n")[[length(msg)]], collapse = " "));
        msg <- gsub("\"", "", msg);
        
        #data$traits[[i]]$analysis$sea$envs[[j]]$lmerRun <- "ERROR";
        #data$traits[[i]]$analysis$sea$envs[[j]]$lmerError <- msg;
        warning(paste("\tWarning: There are error occured on lmer model computed. Cannot proceed with analysis on ", env.name, " of ",
                      trait.name, ".\n", sep = ""));
        next;
      } 
      
      #--- kept all legal single site analysis outcomes---#
      data$traits[[i]]$analysis$sea$envs[[j]] <- list();
      data$traits[[i]]$analysis$sea$envs[[j]]$name <- env.name;				
      
      #--- compute harmonic mean---#
      no.reps <- data.frame(n = tapply(eval(parse(text = paste("temp.data$", trait.name, sep = ""))), eval(parse(text = paste("temp.data$", geno, sep = ""))), FUN = length));
      no.reps <- as.numeric(1/sapply(1/no.reps, mean));
      data$traits[[i]]$analysis$sea$envs[[j]]$numreps <- no.reps;
      
      #data$traits[[i]]$analysis$sea$envs[[j]]$lmerRun <- "NO ERROR";
      data$traits[[i]]$analysis$sea$envs[[j]]$formula1 <- myformula1;
      data$traits[[i]]$analysis$sea$envs[[j]]$model <- model;
      
      #--- get variance components ---#
      varcomp <- NULL;
      for(k in (1:length(VarCorr(model))))
      {
        varcomp <- rbind(varcomp, data.frame(
          Groups = names(VarCorr(model))[k], 
          Variance = VarCorr(model)[[k]][1],
          Std.Dev. = attr(VarCorr(model)[[k]],"stddev")[[1]]
        ));
      }
      varcomp <- rbind(varcomp, data.frame(
        Groups = "Residual",
        Variance = attr(VarCorr(model),"sc")^2,
        Std.Dev. = attr(VarCorr(model), "sc")
      ));
      attr(varcomp, "heading") <- "Variance Components for Random Effects\n";
      data$traits[[i]]$analysis$sea$envs[[j]]$varcomp.table <- varcomp;
      
      #--- for saving variance and sum of reps---#
      data$traits[[i]]$analysis$sea$envs[[j]]$varcompnRep <- as.data.frame(attr(VarCorr(model), "sc")^2);
      data$traits[[i]]$analysis$sea$envs[[j]]$varcompnRep$No.Rep <- no.reps;
      data$traits[[i]]$analysis$sea$envs[[j]]$varcompnRep$Env <- env.name;
      colnames(data$traits[[i]]$analysis$sea$envs[[j]]$varcompnRep) <- c(paste(trait.name, "sigma2", sep="_"),
                                                                         paste(trait.name, "No.Rep", sep="_"),
                                                                         env);
      #--- test significance of genotype effect using maximum likelihood ratio test ---#
      myformula2 <- gsub(paste(" + ", geno, sep=""), "", myformula1, fixed = TRUE);
      model1 <- lmer(formula(myformula1), data = temp.data, REML = T);
      model2 <- lmer(formula(myformula2), data = temp.data, REML = T);
      
      data$traits[[i]]$analysis$sea$envs[[j]]$formula2 <- myformula2;
      m2Vm1 <- anova(model2, model1);
      attr(m2Vm1, "heading") <- attr(m2Vm1, "heading")[-1];
      attr(m2Vm1, "heading") <- attr(m2Vm1, "heading")
      data$traits[[i]]$analysis$sea$envs[[j]]$m2Vm1 <- m2Vm1;
      
      libary(lmerTest);
      model1b <- lmer(formula(myformula1), data = temp.data, REML = T);
      a.table <- anova(model1b);
      data$traits[[i]]$analysis$sea$envs[[j]]$geno.test <- a.table;
      detach("package:lmerTest");
      
      #--- conmpute genotype means---#
      myformula3 <- gsub("~ 1", "~ 0", myformula1);
      model3 <- lmer(formula(myformula3), data = temp.data);
      sumStat.table <- data.frame(summary(model3)$coefficients)[ , 1:2];
      rownames(sumStat.table) <- gsub(geno, "", rownames(sumStat.table));
      sumStat.table <- cbind(rownames(sumStat.table), sumStat.table);
      rownames(sumStat.table) <- NULL;
      colnames(sumStat.table) <- c(geno, "LSMean", "StdErrMean");
      data$traits[[i]]$analysis$sea$envs[[j]]$summary.statistic <- sumStat.table;
      
      #--- display standard error of the differents ---#
      noEntries <- nlevels(temp.data[,geno]);
      covs <- as.matrix(vcov(model3)[1:noEntries, 1:noEntries]);
      vars <- diag(covs);
      vdiff <- outer(vars, vars, "+") - 2 * covs;
      sed <- sqrt(vdiff[upper.tri(vdiff)]);
      
      #--- display sed table ---#
      minSed <- formatC(as.numeric(format(min(sed), scientific = FALSE)), format = "f");
      meanSed <- formatC(as.numeric(format(mean(sed), scientific = FALSE)), format = "f");
      maxSed <- formatC(as.numeric(format(max(sed), scientific = FALSE)), format = "f");
      sedCol <- rbind(minSed, meanSed, maxSed);
      rowNames <- rbind("Minimu", "Average", "Maximum");
      sedTable <- as.table(cbind(rowNames, sedCol));
      rownames(sedTable) <- c("", "", "");
      colnames(sedTable) <- c("", "Estimate");
      data$traits[[i]]$analysis$sea$envs[[j]]$sedTable <- sedTable;
      
      #--- for saving to file ---#
      data$traits[[i]]$analysis$sea$envs[[j]]$sum.out <- sumStat.table;
      data$traits[[i]]$analysis$sea$envs[[j]]$sum.out$Env <- env.name;
      colnames(data$traits[[i]]$analysis$sea$envs[[j]]$sum.out) <- c(geno, paste(trait.name, "Mean", sep = "_"), paste(trait.name, "StdErrMean", sep = "_"), env);
      
      data$traits[[i]]$analysis$sea$envs[[j]]$residuals <- resid(model1);
      data$traits[[i]]$analysis$sea$envs[[j]]$fitted.values <- fitted(model1);
    }#--- end statement of for(j in (1:length(data$triatis[[i]]$envs))) ---#
    
  }#--- end statement of for(i in (1:length(data$traits)))---#
  
  class(data) <- c("SingleEnvAnalysis", class(data));
  #	detach("package:lme4", unload = TRUE);
  return(data);
}


print.SingleEnvAnalysis <- function
(
  data,
  level = 1
)
{
  if(!inherits(data, "SingleEnvAnalysis"))
    stop("\tError: The argument of data must be of class SingleSiteAnalysis!\n");
  if(missing(level))
    level <- 1;
  if(!is.numeric(level))
    stop("\tError: The argument of level should be of value 1 or 2 where 1 is for concise details, and 2 is for precise details.\n");
  if(length(level) != 1)
    stop("\tError: The argument of level should be of length 1.\n");
  if(!any(level %in% c(1,2)))
    stop("\tError: The argument of level only accepted the integer 1 or 2. \n");
  #cat(rep("-", times = 50), sep = "");
  if(level == 1){
    cat("Single Environment Analysis:");
    cat("\n");
    for(i in 1:data$trait.number)
    {
      trait.name <- data$traits[[i]]$name;
      cat("The phenotypic trait is ", trait.name, ".\n", sep = "");
      for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
      {
        env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
        cat("*******************Start of ", env.name, " environment********************\n", sep = "");
        #cat("\t\tOn the environment:", env.name, ".\n", sep = "");
        cat("The variance table:\n");
        print(data$traits[[i]]$analysis$sea$envs[[j]]$varcomp.table, row.names = FALSE);
        cat("\n");
        cat("Adjust Means of each genotype:\n");
        print(data$traits[[i]]$analysis$sea$envs[[j]]$sum.out, row.names = FALSE);
        cat("*******************End of ", env.name, " environment********************\n", sep = "");
        cat("\n");
      }#--- end stmt of for(j in 1:length(data$traits[[i]]$analysis$sea$envs))
      cat("\n\n");
    } #--- end stmt of for(i in 1:data$trait.number) ---#
  }else
  {
    for(i in 1:data$trait.number)
    {
      if(is.null(data$traits[[i]]$analysis$sea))
      {
        next;
      } else
      {
        trait.name <- data$traits[[i]]$name;
        cat(rep("-",times=30), sep="");
        cat("\n");
        cat("RESPONSE VARIABLE: ", trait.name, "\n", sep="");
        cat(rep("-",times=30), sep="");
        cat("\n");
        cat("\n");
        for(j in 1:1:length(data$traits[[i]]$analysis$sea$envs))
        {
          env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
          cat(rep("-",times=30), sep="");
          cat("\n");
          cat("ANALYSIS FOR: Env = ", env.name, "\n", sep="");
          cat(rep("-",times=30), sep="");
          cat("\n");
          cat("\n");
          cat("Data Summary:\n");
          cat("\n");
          cat("Number of observations read: ", data$traits[[i]]$envs[[j]]$obsread, "\n", sep = "");
          cat("Number of observations used: ", data$traits[[i]]$envs[[j]]$obsused, "\n", sep = "");
          cat(" Factors\t","Number of Levels\t", "Levels\n", sep="");
          if(data$traits[[i]]$envs[[j]]$design$exptl.design == "RCB" ||
              data$traits[[i]]$envs[[j]]$design$exptl.design == "AugRCB")
          {
            cat(" ", data$traits[[i]]$envs[[j]]$design$geno, "\t", length(data$traits[[i]]$envs[[j]]$design$geno),"\t", data$traits[[i]]$envs[[j]]$design$geno.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$block, "\t", length(data$traits[[i]]$envs[[j]]$design$block),"\t", data$traits[[i]]$envs[[j]]$design$block.levels,"\n",sep="");
          } else if(data$traits[[i]]$envs[[j]]$design$exptl.design == "AugLS")
          {
            cat(" ", data$traits[[i]]$envs[[j]]$design$geno, "\t", length(data$traits[[i]]$envs[[j]]$design$geno),"\t", data$traits[[i]]$envs[[j]]$design$geno.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$row, "\t", length(data$traits[[i]]$envs[[j]]$design$row),"\t", data$traits[[i]]$envs[[j]]$design$row.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$column, "\t", length(data$traits[[i]]$envs[[j]]$design$column),"\t", data$traits[[i]]$envs[[j]]$design$column.levels,"\n",sep="");
          } else if(data$traits[[i]]$envs[[j]]$design$exptl.design == "Alpha")
          {
            cat(" ", data$traits[[i]]$envs[[j]]$design$geno, "\t", length(data$traits[[i]]$envs[[j]]$design$geno),"\t", data$traits[[i]]$envs[[j]]$design$geno.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$block, "\t", length(data$traits[[i]]$envs[[j]]$design$block),"\t", data$traits[[i]]$envs[[j]]$design$block.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$rep, "\t", length(data$traits[[i]]$envs[[j]]$design$rep),"\t", data$traits[[i]]$envs[[j]]$design$rep.levels,"\n",sep="");
          } else if(data$traits[[i]]$envs[[j]]$design$exptl.design == "RowCol")
          {
            cat(" ", data$traits[[i]]$envs[[j]]$design$geno, "\t", length(data$traits[[i]]$envs[[j]]$design$geno),"\t", data$traits[[i]]$envs[[j]]$design$geno.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$row, "\t", length(data$traits[[i]]$envs[[j]]$design$row),"\t", data$traits[[i]]$envs[[j]]$design$row.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$column, "\t", length(data$traits[[i]]$envs[[j]]$design$column),"\t", data$traits[[i]]$envs[[j]]$design$column.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$rep, "\t", length(data$traits[[i]]$envs[[j]]$design$rep),"\t", data$traits[[i]]$envs[[j]]$design$rep.levels,"\n",sep="");
          } else if(data$traits[[i]]$envs[[j]]$design$exptl.design == "LatinAlpha")
          {
            cat(" ", data$traits[[i]]$envs[[j]]$design$geno, "\t", length(data$traits[[i]]$envs[[j]]$design$geno),"\t", data$traits[[i]]$envs[[j]]$design$geno.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$block, "\t", length(data$traits[[i]]$envs[[j]]$design$block),"\t", data$traits[[i]]$envs[[j]]$design$block.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$rep, "\t", length(data$traits[[i]]$envs[[j]]$design$rep),"\t", data$traits[[i]]$envs[[j]]$design$rep.levels,"\n",sep="");
          } else if(data$traits[[i]]$envs[[j]]$design$exptl.design == "LatinRowCol")
          {
            cat(" ", data$traits[[i]]$envs[[j]]$design$geno, "\t", length(data$traits[[i]]$envs[[j]]$design$geno),"\t", data$traits[[i]]$envs[[j]]$design$geno.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$row, "\t", length(data$traits[[i]]$envs[[j]]$design$row),"\t", data$traits[[i]]$envs[[j]]$design$row.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$column, "\t", length(data$traits[[i]]$envs[[j]]$design$column),"\t", data$traits[[i]]$envs[[j]]$design$column.levels,"\n",sep="");
            cat(" ", data$traits[[i]]$envs[[j]]$design$rep, "\t", length(data$traits[[i]]$envs[[j]]$design$rep),"\t", data$traits[[i]]$envs[[j]]$design$rep.levels,"\n",sep="");
          }
          cat("\n");
          cat("Variance Components Table: \n");
          cat("\n");
          print(data$traits[[i]]$analysis$sea$envs[[j]]$varcomp.table, row.names = FALSE);
          cat("\n\n");
          cat("Testing for the Significance of Genotypic Effect: \n");
          cat("\n");
          
          
        } # end stmt of for(j in 1:1:length(data$traits[[i]]$analysis$sea$envs))
      }# end stmt of if(is.null(data$traits[[i]]$analysis$sea))
    }# end stmt of for(i in 1:data$trait.number)
  }
}