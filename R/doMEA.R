###############################################################################
# TODO: Do Multiple Environment Analysis
#
# [Arguments]
#	data - an object of PhenotypicData and It could not be of mean type!
#	is.EnvFixed - specified whether environment is fixed or random
# 
# Author: mqin
# Date: Nov 14, 2014
# FileName: doMSA.R
###############################################################################


doMEA <- function(
  data,
  is.EnvFixed = TRUE
)
{
  UseMethod("doMEA");
}
doMEA.PhenotypicData <- function(
  data,
  is.EnvFixed = TRUE
)
{
  #library(lme4);
  #library(lsmeans);
  options(show.signif.stars = FALSE);
  pre.opt <- options()$warn;
  #--- comment out this option ---#
  #option(warn = -1) 
  
  if(!inherits(data, "PhenotypicData"))
    stop("\tError: The argument of data must be of class PhenotypicData!\n");
  if(data$isMean)
    stop("\tError: The mean's type phenotypic data could not process multiple environment analysis!\n");
  if(!data$isRestricted)
  {
    warning(paste("\tWarning: It will use the default parameter to restrict phenotypic data;\n",
                  "\tMissing rate should be not larger than 0.2 of each trait on each environment!\n", 
                  #"\tAnd used the grand mean of each trait on each environment to instead of missing all obseravtion!\n", 
                  sep=""));
    data <- restrict.pheno.data(data);
  }
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    #--- combinded all envs data under this traits---#
    temp.data <- NULL;
    if(length(data$traits[[i]]$envs) <= 1)
    {
      warning(paste("\tWarning: It could not implement mult-site analysis on ", trait.name, 
                    " trait,\n\tBecause of less than two environment! \n\tTry next trait!", sep =""));
      next;
    }
    for(j in 1:length(data$traits[[i]]$envs))
    {
      #--- whether used original data ??---#
      #--- or according to restricted information---#
      #--- after combind all the data, checking missing rate again---#
      temp.data <- rbind(temp.data, data$traits[[i]]$envs[[j]]$data);
    }
    temp.data[ , trait.name] <- as.numeric(as.character(temp.data[,trait.name]));
    
    #--- checking missing rate again not meet by the specified value---#
    missing.rate <- data$traits[[i]]$envs[[1]]$restricted$missing.rate.cond;
    if(sum(is.na(temp.data[,trait.name])) / nrow(temp.data) > missing.rate)
    {
      warnings(paste("\tWarning: It could not process multiple environment analysis because of", 
                     "missing rate is large than specified ", missing.rate, ". Try next trait!\n", sep = ""));
      next;
    }
    #---retrived all desigen factor variable---#
    geno <- data$traits[[i]]$envs[[1]]$design$geno;
    env <- data$traits[[i]]$envs[[1]]$design$env;
    exptl.design <- data$traits[[i]]$envs[[1]]$design$exptl.design;
    if(exptl.design %in% c("RCB", "AugRCB"))
    {
      block <- data$traits[[i]]$envs[[1]]$design$block;
    } else if(exptl.design == "AugLS")
    {
      row <- data$traits[[i]]$envs[[1]]$design$row;
      column <- data$traits[[i]]$envs[[1]]$design$column;
    } else if (exptl.design %in% c("Alpha", "LatinAlpha"))
    {
      block <- data$traits[[i]]$envs[[1]]$design$block;
      rep <- data$traits[[i]]$envs[[1]]$design$rep;
    } else if(exptl.design %in% c("RowCol", "LatinRowCol"))
    {
      row <- data$traits[[i]]$envs[[1]]$design$row;
      column <- data$traits[[i]]$envs[[1]]$design$column;
      rep <- data$traits[[i]]$envs[[1]]$design$rep;
    }
    
    #--- get levels of genotype and environment---#
    levelsGeno <- levels(temp.data[ ,geno]);
    levelsEnv <- levels(temp.data[,env]);
    
    #---checking whether missing whole levels of this trait---#
    if(countMissingLevelsByGroupFactor(trait.name, groupfactor = c(geno, env), data = temp.data) && is.EnvFixed)
    {
      warning(paste("\tWarning: There are missing one levels of combination between Genotype and Environment on \n" 
                    , "\t", trait.name , " when the env factor is setted to be fixed.",
                    " Please checking your data or setting env factor to be random!\n\tTry next trait!", sep = ""));
      next;
    }
    
    # --- if design is Latinized Row-Column, check if the data follow case1 or case3 labeling --- #
    if(exptl.design == "LatinRowCol")
    {
      lengthPerCross <- tapply(temp.data[ , trait.name], temp.data[ ,c(row, column)], length);
      if(all(lengthPerCross <= nlevels(temp.data[,env]), na.rm=TRUE))
      {
        if(nlevels(temp.data[,row]) > nlevels(temp.data[,column]))
        {
          longerRow <- TRUE;
        } else
        {
          longerRow <- FALSE;
        }
      } else
      {
        stop("The levels of the row/column variable should be continuous across replicates.");
      }	
    }
    
    # --- create temp.data without missing observations --- #
    used.data <- subset(temp.data, subset = (is.na(temp.data[,trait.name]) == FALSE));
    
    #---compute number of observation read, used and response rate---#
    obsread <- nrow(temp.data);
    obsused <- nrow(used.data);
    responseRate <- (obsused/obsread);
    if(responseRate < 0.8)
    {
      warning("Too many missing observation. Could not process multiple environment analysis. Try next trait!");
      next;
    } else
    {
      if(is.null(data$traits[[i]]$analysis))
        data$traits[[i]]$analysis <- list();
      if(is.null(data$traits[[i]]$analysis$mea))
      {
        data$traits[[i]]$analysis$mea <- list();
      } else
      {
        data$traits[[i]]$analysis$mea <- NULL;
        data$traits[[i]]$analysis$mea <- list();
      }
      data$traits[[i]]$analysis$mea$obsread <- obsread;
      data$traits[[i]]$analysis$mea$obsused <- obsused;
      data$traits[[i]]$analysis$mea$nlevelsGeno <- length(levelsGeno);
      data$traits[[i]]$analysis$mea$nlevelsEnv <- length(levelsEnv);
      
      #--- compute summary statistics per environment---#
      sumStat.Env <- DescriptiveStatistics(used.data, trait.name, env, c("min", "mean", "max", "var", "sd"));
      sumStat.Env <- sumStat.Env[ , c(2:ncol(sumStat.Env))];
      
      #--- construct the model---#
      if(is.EnvFixed) 
      {
        env.stmt <- paste(env," + ", geno, ":", env, sep ="");
      }else
      {
        env.stmt <- paste("(1|", env, ") + (1|",geno, ":", env, ")", sep="");
      }
      factor.summary <- c();
      env.levels <- data$traits[[i]]$envs[[1]]$design$env.levels;
      geno.levels <- data$traits[[i]]$envs[[1]]$design$geno.levels;
      factor.summary <- rbind(c(env, length(env.levels), 
                              ifelse(length(env.levels) <= 4, env.levels, 
                                     paste(env.levels[1:3],"...",env.levels[length(env.levels)])
                              )));
      factor.summary <- rbind(factor.summary,c(geno, length(geno.levels), 
                              ifelse(length(geno.levels) <= 4, geno.levels, 
                                     paste(geno.levels[1:3],"...",geno.levels[length(geno.levels)])
                              )));
      if(exptl.design == "RCB" || exptl.design =="AugRCB")
      {
        block.levels <- data$traits[[i]]$envs[[1]]$design$block.levels;
        factor.summary <- rbind(factor.summary,c(block, length(block.levels), 
                                ifelse(length(block.levels) <= 4, block.levels, 
                                       paste(block.levels[1:3],"...",block.levels[length(block.levels)])
                                )));
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , block , ":", env,")", sep="");
      } else if(exptl.design == "AugLS")
      {
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
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , row , ":", env, ") + (1|", column, ":", env, ")", sep="");
      } else if(exptl.design == "Alpha")
      {
        block.levels <- data$traits[[i]]$envs[[1]]$design$block.levels;
        factor.summary <- rbind(factor.summary,c(block, length(block.levels), 
                                ifelse(length(block.levels) <= 4, block.levels, 
                                       paste(block.levels[1:3],"...",block.levels[length(block.levels)])
                                )));
        rep.levels <- data$traits[[i]]$envs[[1]]$design$rep.levels;
        factor.summary <- rbind(factor.summary,c(rep, length(rep.levels), 
                                ifelse(length(rep.levels) <= 4, rep.levels, 
                                       paste(rep.levels[1:3],"...",rep.levels[length(rep.levels)])
                                )));
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , rep , ":", env,") + (1|", rep, ":", block,":",env,")", sep = "" );
      } else if(exptl.design == "RowCol")
      {
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
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + " , env.stmt , " + (1|", rep , ":", env,") + (1|", rep,":", row,":", env,") + (1|", rep,":", column,":", env,")", sep = "");
      } else if(exptl.design == "LatinAlpha")
      {
        block.levels <- data$traits[[i]]$envs[[1]]$design$block.levels;
        factor.summary <- rbind(factor.summary,c(block, length(block.levels), 
                                ifelse(length(block.levels) <= 4, block.levels, 
                                       paste(block.levels[1:3],"...",block.levels[length(block.levels)])
                                )));
        rep.levels <- data$traits[[i]]$envs[[1]]$design$rep.levels;
        factor.summary <- rbind(factor.summary,c(rep, length(rep.levels), 
                                ifelse(length(rep.levels) <= 4, rep.levels, 
                                       paste(rep.levels[1:3],"...",rep.levels[length(rep.levels)])
                                )));
        myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + " , env.stmt , " + (1|" , rep , ":", env, ") + (1|", block, ":", env,") + (1|", rep, ":", block, ":", env,")", sep = "");
      } else if(exptl.design == "LatinRowCol")
      {
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
        if(longerRow)
        {
          myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + ", env.stmt , " + (1|" , rep, ":", env, ") + (1|", column, ":", env, ") + (1|", rep, ":", column, ":", env, ") + (1|", row, ":", env,")", sep = "");
        } else
        {
          myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + ", env.stmt , " + (1|" , rep, ":", env, ") + (1|", row, ":", env, ") + (1|", rep, ":" + row, ":", env,") + (1|", column, ":" , env, ")", sep = "");
        }
      }
      
      model <- lmer(formula(myformula1), data = used.data);
      data$traits[[i]]$analysis$mea$formula1 <- myformula1;
      data$traits[[i]]$analysis$mea$model <- model;
      data$traits[[i]]$analysis$mea$envFixed <- is.EnvFixed;
      factor.summary <- as.data.frame(factor.summary);
      colnames(factor.summary) <- c("Factor", "No. of Levels", "Levels");
      data$traits[[i]]$analysis$mea$factor.summary <- factor.summary;
      
      #--- variance components ---#
      varcomp <- NULL;
      for(k in (1:length(VarCorr(model))))
      {
        varcomp <- rbind(varcomp, 
                         data.frame(
                           Groups = names(VarCorr(model))[k], 
                           Variance = VarCorr(model)[[k]][1], 
                           Std.Dev. = attr(VarCorr(model)[[k]], "stddev")[[1]]
                         )
        );
      }
      varcomp <- rbind(varcomp, 
                       data.frame(
                         Groups = "Residual", 
                         Variance = attr(VarCorr(model), "sc")^2, 
                         Std.Dev. = attr(VarCorr(model), "sc")
                       )
      );
      
      data$traits[[i]]$analysis$mea$varcomp.table <- varcomp;
      
      #--- test for significance of genotypic effect using LRT && ANOVA table if geno is fixed ---#
      #--- myformula2 is full model minus the genotype term---#
      myformula2 <- sub(paste(" + ", geno, sep=""), "", myformula1, fixed = TRUE);
      model1 <- lmer(formula(myformula1), data = used.data, REML = T);
      model2 <- lmer(formula(myformula2), data = used.data, REML = T);
      
      data$traits[[i]]$analysis$mea$formula2 <- myformula2;
      
      #--- This is for genotype as fixed---#
      anova.table1 <- anova(model2, model1);
      rownames(anova.table1) <- c("Model2", "Model1");
      attr(anova.table1, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF GENOTYPIC EFFECT USING LIKELIHOOD RATIO TEST:\n", sep = "");
      attr(anova.table1, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
      attr(anova.table1, "heading")[3] <- paste("Formula for Model2: ", myformula2, sep = "");
      attr(anova.table1, "heading")[4] <- paste("", sep = "");
      data$traits[[i]]$analysis$mea$testsig.Geno <- anova.table1;
      
      #--- test for significance of environment effect using LRT if env is fixed or random---#
      #--- myformula3 is full model minus the environment term ---#
      if(is.EnvFixed)
      {
        myformula3 <- gsub(paste(" + ", env, sep = ""), "" , myformula1, fixed = TRUE);
        model3 <- lmer(formula(myformula3), data = used.data, REML = TRUE);
        data$traits[[i]]$analysis$mea$myformula3;
        
        anova.table2 <- anova(model3, model1);
        rownames(anova.table2) <- c("Model3", "Model1");
        attr(anova.table2, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF ENVIRONMENTAL EFFECT USING LIKELIHOOD RATIO TEST:\n", sep ="");
        attr(anova.table2, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
        attr(anova.table2, "heading")[3] <- paste("Formula for Model3: ", myformula3, sep = "");
        attr(anova.table2, "heading")[4] <- paste("", sep="");
        data$traits[[i]]$analysis$mea$testsig.Env <- anova.table2;
      } else
      {
        myformula3 <- gsub(paste(" + (1|", env, ")", sep = ""), "", myformula1, fixed = TRUE);
        model3 <- lmer(formula(myformula3), data = used.data, REML = TRUE);
        data$traits[[i]]$analysis$mea$formula3 <- myformula3;
        
        models.table2 <- modelComparisonTable(model1, model3);
        data$traits[[i]]$analysis$mea$testsig.Env <- models.table2;
      }
      
      #--- test of significance of genotype by environment effect using LRT---#
      #--- myformula4 is full model minus the geno by environment interaction term ---#
      if(is.EnvFixed)
      {
        myformula4 <- gsub(paste(" + ", geno, ":", env, sep = ""), "", myformula1, fixed = TRUE);
        model4 <- lmer(formula(myformula4), data = used.data, REML = TRUE);
        data$traits[[i]]$analysis$mea$formula4 <- myformula4;
        
        anova.table3 <- anova(model4, model1);
        rownames(anova.table3) <- c("Model4", "Model1");
        attr(anova.table3, "heading")[1] <- paste("TESTING FOR THE SIGNIFICANCE OF GENOTYPE BY ENVIRONMENT EFFECT USING LIKELIHOOD RATIO TEST:\n", sep="");
        attr(anova.table3, "heading")[2] <- paste("Formula for Model1: ", myformula1, sep = "");
        attr(anova.table3, "heading")[3] <- paste("Formula for Model4: ", myformula4, sep = "");
        attr(anova.table3, "heading")[4] <- paste("", sep = "");
        data$traits[[i]]$analysis$mea$testsig.GenoEnv <- anova.table3;
        
      } else 
      {
        myformula4 <- gsub(paste(" + (1|", geno, ":", env, ")", sep = ""), "", myformula1, fixed = TRUE);
        model4 <- lmer(formula(myformula4), data = used.data, REML = TRUE);
        
        models.table3<-modelComparisonTable(model1, model4);
        
        data$traits[[i]]$analysis$mea$formula4 <- myformula4;
        data$traits[[i]]$analysis$mea$testsig.GenoEnv <- models.table3;
      }
      
      #--- predicted means/lsmeans of genotype and SED (Standard Error of the Difference) statistics ---#
      #--- myformula5 is full model but without intercept ---#
      myformula5 <- gsub("~ 1", "~ 0", myformula1, fixed = TRUE);
      model.noint <- lmer(formula(myformula5), data = used.data);
      
      # --- Modify by QIN MAO FOR USING LSMENAS PACKAGE TO GET LSMEANS --- #
      sumStat.Geno <- lsmeans(model, geno, weights = "cells");
      sumStat.Geno <- summary(sumStat.Geno)[,1:3];
      colnames(sumStat.Geno) <- c(geno, "LSMean", "StdErrMean");
      rownames(sumStat.Geno) <- NULL;
      
      # --- display standard error of the differences --- #
      # ??? how to caculated standard error of the differences using with intercept model ??? #
      noEntries<-nlevels(used.data[,geno]);
      covs <- as.matrix(vcov(model.noint)[1:noEntries, 1:noEntries]);
      vars <- diag(covs);
      vdiff <- outer(vars, vars, "+") - 2 * covs;
      sed <- sqrt(vdiff[upper.tri(vdiff)]);
      
      # --- display SED Table --- #
      minSed<-formatC(as.numeric(format(min(sed), scientific=FALSE)), format="f");
      meanSed<-formatC(as.numeric(format(mean(sed), scientific=FALSE)), format="f");
      maxSed<-formatC(as.numeric(format(max(sed), scientific=FALSE)), format="f");
      sedCol<-rbind(minSed, meanSed, maxSed);
      rowNames<-rbind("Minimum  ", "Average  ", "Maximum  ");
      sedTable<-as.table(cbind(rowNames, sedCol));
      rownames(sedTable)<-c("","","");
      colnames(sedTable)<-c("","Estimate");
      data$traits[[i]]$analysis$mea$sedTable <- sedTable;
      
      # --- ESTIMATES OF EFFECTS --- #
      if(is.EnvFixed)
      {
        # --- ENVIRONMENT LSMEANS --- #
        sumStat.Env <- lsmeans(model1,env);
        sumStat.Env <- summary(sumStat.Env)[,1:3];
        colnames(sumStat.Env) <- c(env, "LSMean", "StdErrMean");
        rownames(sumStat.Env) <- NULL;
        
        # --- G X E MEANS --- #
        # IF ENVIRONMENT IS FIXED, GXE LSMEANS COULD BE RETRIEVED DIRECTLY BY LSMEANS;
        sumStat.GenoEnv <- lsmeans(model1, c(geno, env));
        sumStat.GenoEnv <- summary(sumStat.GenoEnv)[,1:4];
        colnames(sumStat.GenoEnv) <- c(geno, env, paste(trait.name,"LSMean",sep ="_"), "StdErrMean");
        rownames(sumStat.GenoEnv) <- NULL;
      } else
      {
        # -- GENOTYPE EFFECT --- #
        intercept <- fixef(model)[[1]];
        geno.effect <- as.data.frame(fixef(model)[-1]);
        geno.effect <- data.frame(gsub(geno, "", rownames(geno.effect)), geno.effect);
        colnames(geno.effect) <- c(geno, "geno_effect")
        rownames(geno.effect) <- NULL
        
        # -- ENVIRONMENT EFFECT--- #
        env.effect <- eval(parse(text = paste("ranef(model)$", env, sep = "")));
        env.effect <- data.frame(gsub(env, "", rownames(env.effect)), env.effect);
        colnames(env.effect) <- c(env, "env_effect");
        rownames(env.effect) <- NULL;
        
        # -- G X E EFFECT --- #
        GXE.effect <- as.data.frame(eval(parse(text = paste("ranef(model)$'", geno, ":", env, "'",sep = ""))));
        names <- t(as.data.frame(strsplit(rownames(GXE.effect), ":")));
        GXE.effect <- data.frame(names[,1], names[,2], GXE.effect+intercept);
        colnames(GXE.effect) <- c(geno, env, "ge_effect");
        rownames(GXE.effect) <- NULL;
        
        # -- G X E MEANS --- #
        sumStat.GenoEnv <- merge(GXE.effect, env.effect, by = env, all = TRUE);
        sumStat.GenoEnv <- merge(sumStat.GenoEnv, geno.effect, by = geno, all = TRUE);
        sumStat.GenoEnv <- data.frame(sumStat.GenoEnv[,match(geno, names(sumStat.GenoEnv))], sumStat.GenoEnv[,match(env,names(sumStat.GenoEnv))], rowSums(subset(sumStat.GenoEnv, select = c(ge_effect, env_effect, geno_effect)), na.rm = TRUE));
        colnames(sumStat.GenoEnv) <- c(geno, env, paste(trait.name, "LSMean", sep = "_"));
      }
      
      # -- create G x E MEANS with coded levels -- #
      sumStat.GenoEnvCode <- sumStat.GenoEnv;
      sumStat.GenoEnvCode$CodedGeno <- sumStat.GenoEnvCode[,match(geno, names(sumStat.GenoEnvCode))];
      sumStat.GenoEnvCode$CodedEnv <- sumStat.GenoEnvCode[,match(env, names(sumStat.GenoEnvCode))];
      
      # --- display G x E means in 2-way table --- #
      #wide.GenoEnv<-ToWide(sumStat.GenoEnv, paste(trait.name, "means", sep = "_"), env, geno)
      wide.GenoEnv <- reshape(
        sumStat.GenoEnv[,1:3], 
        v.names=paste(trait.name, "LSMean", sep = "_"), 
        timevar=env, 
        idv=geno, 
        direction="wide"
      );
      colnames(wide.GenoEnv) <- gsub(
        paste(trait.name, "LSMean.", sep = "_"), 
        "", 
        colnames(wide.GenoEnv));
      rownames(wide.GenoEnv) <- 1:nrow(wide.GenoEnv);
      
      data$traits[[i]]$analysis$mea$means.Geno <- sumStat.Geno;
      data$traits[[i]]$analysis$mea$means.Env <- sumStat.Env;
      data$traits[[i]]$analysis$mea$means.GenoEnv <- sumStat.GenoEnv;
      data$traits[[i]]$analysis$mea$wide.GenoEnv <- wide.GenoEnv;
      data$traits[[i]]$analysis$mea$means.GenoEnvCode <- sumStat.GenoEnvCode;
      data$traits[[i]]$analysis$mea$residuals <- resid(model1);
      data$traits[[i]]$analysis$mea$fitted.values <- fitted(model1);
      data$traits[[i]]$analysis$mea$data <- used.data;
      
      # --- if genotype is fixed, output MSE and harmonic mean for AMMI --- #
      # --- compute harmonic mean per environment level --- #
      envgenorep <- as.data.frame.table(tapply(temp.data[, trait.name], used.data[,c(env, geno)], length));
      envgenorep <- envgenorep[(is.na(envgenorep[,"Freq"]) == FALSE),];
      envgenorep$reciprocal <- 1/envgenorep$Freq;
      envgenorep2 <- as.data.frame.table(tapply(envgenorep[, "reciprocal"], envgenorep[,env], mean));
      envgenorep2$reciprocal2 <- 1/envgenorep2$Freq;
      envgenorep3 <- merge(envgenorep, envgenorep2, by.x=env, by.y="Var1");
      
      numrepEnv <- tapply(envgenorep3[,"reciprocal2"], envgenorep3[,env], mean);
      no.reps <- 1/mean(1/numrepEnv);
      
      data$traits[[i]]$analysis$mea$harmonicMean <- no.reps;
      data$traits[[i]]$analysis$mea$MSE <- varcomp[varcomp[,1] == "Residual", "Variance"];
      
    }
    
    
  } #--- end of statment for(i in 1:length(data$traits))
  
  class(data) <- c("MultiEnvAnalysis", class(data));
  #detach("package:lme4", unload = TRUE);
  #detach("package:lsmeans", unload = TRUE);
  return(data);
}
# 
# print.MultiEnvAnalysis <- function(data, level = 1)
# {
# 	if(!inherits(data, "MultiEnvAnalysis"))
# 		stop("\tError: The argument of data must be of class MultiSiteAnalysis!\n");
# 	if(missing(level))
# 		level <- 1;
# 	if(!is.numeric(level))
# 		stop("\tError: The argument of level should be of value 1 or 2 where 1 is for concise details, and 2 is for precise details.\n");
# 	if(length(level) != 1)
# 		stop("\tError: The argument of level should be of length 1.\n");
# 	if(!any(level %in% c(1,2)))
# 		stop("\tError: The argument of level only accepted the integer 1 or 2. \n");
# 	if(level == 1)
# 	{
#     cat("Multiple Environment Analyis:");
#   	cat("\n");
#   	for(i in 1:data$trait.number)
#   	{
#   		trait.name <- data$traits[[i]]$name;
#   		cat("On the trait: ", trait.name, ".\n", sep ="");
#   		cat("Environment is fixed: ", ifelse(data$traits[[i]]$analysis$mea$envFixed, "YES", "NO"),".\n", sep = "");		
#   		cat("Variance Table:\n");
#   		print(data$traits[[i]]$analysis$mea$varcomp.table, row.names = FALSE);
#   		cat("\n");
#   		cat("Adjusted Means Across Environment:\n");
#   		print(data$traits[[i]]$analysis$mea$wide.GenoEnv, row.names = FALSE);
#   		cat(rep("-", times=50), sep = "");
#   		cat("\n")
#   	}
# 	} else
# 	{
#     for(i in 1:length(data$traits))
#     {
#       if(is.null(data$traits[[i]]$analysis$mea))
#         next;
#       cat(rep("=",40), sep="");
#       cat("\n");
#       cat("RESPONSE VARIABLE: ");
#       cat(data$traits[[i]]$name);
#       cat("\n");
#       cat(rep("=",40), sep="");
#       cat("\n");
#       cat(rep("-",40), sep="");
#       cat("\n");
#       cat("ENVIRONMENT AS: ");
#       cat(data$traits[[i]]$analysis$mea$envFixed ? "FIXED" : "RANDOM");
#       cat("\n");
#       cat(rep("-",40), sep="");
#       cat("\n");
#     }
# 	}
# }
