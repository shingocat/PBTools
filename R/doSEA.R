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
      
      factor.summary <- c();
      env.levels <- data$traits[[i]]$envs[[j]]$design$env.levels;
      geno.levels <- data$traits[[i]]$envs[[j]]$design$geno.levels;
      factor.summary <- rbind(c(env, length(env.levels), 
                                ifelse(length(env.levels) <= 4, env.levels, 
                                       paste(env.levels[1:3],"...",env.levels[length(env.levels)])
                                )));
      factor.summary <- rbind(factor.summary,c(geno, length(geno.levels), 
                                               ifelse(length(geno.levels) <= 4, geno.levels, 
                                                      paste(geno.levels[1:3],"...",geno.levels[length(geno.levels)])
                                               )));
      
      #--- construct model, geno as fixed factor ---#
      if(exptl.design == "RCB" || exptl.design == "AugRCB")
      {
        block <- data$traits[[i]]$envs[[j]]$design$block;
        block.levels <- data$traits[[i]]$envs[[j]]$design$block.levels;
        factor.summary <- rbind(factor.summary,c(block, length(block.levels), 
                                                 ifelse(length(block.levels) <= 4, block.levels, 
                                                        paste(block.levels[1:3],"...",block.levels[length(block.levels)])
                                                 )));
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + (1|",block, ")", sep="");
      } else if(exptl.design == "AugLS")
      {
        row <- data$traits[[i]]$envs[[j]]$design$row;
        column <- data$traits[[i]]$envs[[j]]$design$column;
        row.levels <- data$traits[[i]]$envs[[j]]$design$row.levels;
        factor.summary <- rbind(factor.summary,c(row, length(row.levels), 
                                                 ifelse(length(row.levels) <= 4, row.levels, 
                                                        paste(row.levels[1:3],"...",row.levels[length(row.levels)])
                                                 )));
        column.levels <- data$traits[[i]]$envs[[j]]$design$column.levels;
        factor.summary <- rbind(factor.summary,c(column, length(column.levels), 
                                                 ifelse(length(column.levels) <= 4, column.levels, 
                                                        paste(column.levels[1:3],"...",column.levels[length(column.levels)])
                                                 )));
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + (1|", row , "), + (1|", column, ")", sep ="");
      } else if(exptl.design == "Alpha")
      {
        block <- data$traits[[i]]$envs[[j]]$design$block;
        rep <- data$traits[[i]]$envs[[j]]$design$rep;
        block.levels <- data$traits[[i]]$envs[[j]]$design$block.levels;
        factor.summary <- rbind(factor.summary,c(block, length(block.levels), 
                                                 ifelse(length(block.levels) <= 4, block.levels, 
                                                        paste(block.levels[1:3],"...",block.levels[length(block.levels)])
                                                 )));
        rep.levels <- data$traits[[i]]$envs[[j]]$design$rep.levels;
        factor.summary <- rbind(factor.summary,c(rep, length(rep.levels), 
                                                 ifelse(length(rep.levels) <= 4, rep.levels, 
                                                        paste(rep.levels[1:3],"...",rep.levels[length(rep.levels)])
                                                 )));
        myformula1 <- paste(trait.name, " ~ 1 + ", geno," + (1|", rep,"/", block,")", sep = "");
      } else if(exptl.design == "RowCol")
      {
        row <- data$traits[[i]]$envs[[j]]$design$row;
        column <- data$traits[[i]]$envs[[j]]$design$column;
        rep <- data$traits[[i]]$envs[[j]]$design$rep;
        row.levels <- data$traits[[i]]$envs[[j]]$design$row.levels;
        factor.summary <- rbind(factor.summary,c(row, length(row.levels), 
                                                 ifelse(length(row.levels) <= 4, row.levels, 
                                                        paste(row.levels[1:3],"...",row.levels[length(row.levels)])
                                                 )));
        column.levels <- data$traits[[i]]$envs[[j]]$design$column.levels;
        factor.summary <- rbind(factor.summary,c(column, length(column.levels), 
                                                 ifelse(length(column.levels) <= 4, column.levels, 
                                                        paste(column.levels[1:3],"...",column.levels[length(column.levels)])
                                                 )));
        rep.levels <- data$traits[[i]]$envs[[j]]$design$rep.levels;
        factor.summary <- rbind(factor.summary,c(rep, length(rep.levels), 
                                                 ifelse(length(rep.levels) <= 4, rep.levels, 
                                                        paste(rep.levels[1:3],"...",rep.levels[length(rep.levels)])
                                                 )));
        myformula1 <- paste(trait.name, " ~ 1 + ", geno," + (1|", rep,") + (1|", rep,":", row,") + (1|", rep, ":", column,")", sep = "");
      } else if(exptl.design == "LatinAlpha")
      {
        block <- data$traits[[i]]$envs[[j]]$design$block;
        rep <- data$traits[[i]]$envs[[j]]$design$rep;
        block.levels <- data$traits[[i]]$envs[[j]]$design$block.levels;
        factor.summary <- rbind(factor.summary,c(block, length(block.levels), 
                                                 ifelse(length(block.levels) <= 4, block.levels, 
                                                        paste(block.levels[1:3],"...",block.levels[length(block.levels)])
                                                 )));
        rep.levels <- data$traits[[i]]$envs[[j]]$design$rep.levels;
        factor.summary <- rbind(factor.summary,c(rep, length(rep.levels), 
                                                 ifelse(length(rep.levels) <= 4, rep.levels, 
                                                        paste(rep.levels[1:3],"...",rep.levels[length(rep.levels)])
                                                 )));
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + (1|", rep,") + (1|", block,") + (1|", rep, ":", block, ")", sep = "");
      } else if(exptl.design == "LatinRowCol")
      {
        row <- data$traits[[i]]$envs[[j]]$design$row;
        column <- data$traits[[i]]$envs[[j]]$design$column;
        rep <- data$traits[[i]]$envs[[j]]$design$rep;
        row.levels <- data$traits[[i]]$envs[[1]]$design$row.levels;
        factor.summary <- rbind(factor.summary,c(row, length(row.levels), 
                                                 ifelse(length(row.levels) <= 4, row.levels, 
                                                        paste(row.levels[1:3],"...",row.levels[length(row.levels)])
                                                 )));
        column.levels <- data$traits[[i]]$envs[[1]]$design$column.levels;
        factor.summary <- rbind(factor.summary,c(column, length(column.levels), 
                                                 ifelse(length(column.levels) <= 4, column.levels, 
                                                        paste(column.levels[1:3],"...",column.levels[length(column.levels)])
                                                 )));
        rep.levels <- data$traits[[i]]$envs[[1]]$design$rep.levels;
        factor.summary <- rbind(factor.summary,c(rep, length(rep.levels), 
                                                 ifelse(length(rep.levels) <= 4, rep.levels, 
                                                        paste(rep.levels[1:3],"...",rep.levels[length(rep.levels)])
                                                 )));
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
	  #--- modify by QINMAO on 2015-6-8, bug fix that when the trait value sometimes is treated as
	  #--- factor nor numeric ---#
	  temp.data[[trait.name]] <- as.numeric(as.character(temp.data[[trait.name]]));
	  
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
      factor.summary <- as.data.frame(factor.summary);
      colnames(factor.summary) <- c("Factor", "No. of Levels", "Levels");
      data$traits[[i]]$analysis$sea$envs[[j]]$factor.summary <- factor.summary;
      
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
      
      library(lmerTest);
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


