stability.analysis.MultiEnvAnalysis <- function
(
  data, 
  method = c("regression", "shukla"),
  ...
)
{
  if(!inherits(data, "MultiEnvAnalysis"))
    stop("\tError: The data argument should be an object of MultiEnvAnalyssi!\n");
  method <- match.arg(method);
  
  for(i in 1:length(data$traits))
  {
    trait.name <- data$traits[[i]]$name;
    if(is.null(data$traits[[i]]$analysis$mea))
    {
      warning(cat("\tWarning: This is not Multi-Environment Analysis outcomes 
                    on this trait ",trait.name, ".\n", sep=""));
      next;
    } else
    {
      if(data$traits[[1]]$analysis$mea$nlevelsEnv <= 4)
      {
        warning(cat("\tWarning: The enivronment factor should have at least five levels  
                    on this trait ",trait.name, ".\n", sep=""));
        next;
      }
      
      if(is.null(data$traits[[i]]$analysis$mea$stability))
      {
        data$traits[[i]]$analysis$mea$stability <- list();
      } else
      {
        data$traits[[i]]$analysis$mea$stability <- NULL;
        data$traits[[i]]$analysis$mea$stability <- list();
      }
      
      trait.name <- paste(trait.name,"_LSMean",sep="");
      means.GenoEnvCode <- data$traits[[1]]$analysis$mea$means.GenoEnvCode;
      geno <- colnames(means.GenoEnvCode)[1];
      env <- colnames(means.GenoEnvCode)[2];
      
      if(method == "regression")
      {
        site.means <- data.frame(levels(means.GenoEnvCode[, match(env, names( means.GenoEnvCode))]),
                                 as.data.frame(tapply( means.GenoEnvCode[, match(trait.name, names( means.GenoEnvCode))], 
                                                       means.GenoEnvCode[, match(env, names( means.GenoEnvCode))], mean, na.rm = TRUE)));
        colnames(site.means) <- c(env, "site.index");
        rownames(site.means) <- NULL;
        temp.data <- subset( means.GenoEnvCode, select = c(geno, env, trait.name));
        temp.data <- subset(temp.data, subset = 
                              (is.na(temp.data[,match(trait.name, names(temp.data))]) == F));
        trt.nlevels <- nlevels(temp.data[,match(geno, names(temp.data))]);
        ge.means.wide <- reshape(temp.data, v.names = trait.name, idvar = env, 
                                 timevar = geno, direction = "wide");
        colnames(ge.means.wide) <- gsub(paste(trait.name, ".", sep = ""), "", colnames(ge.means.wide));
        data.all <- merge(site.means, ge.means.wide, by = env);
        slope <- as.matrix(c(1:(6*trt.nlevels)), nrow = trt.nlevels, ncol = 6);
        dim(slope) <- c(trt.nlevels, 6);
        
        if (var(data.all$site.index) != 0) {
          for (j in (1:trt.nlevels)) {
            if (length(na.omit(data.all[[j+2]])) < 3) { slope[j,] <- rep(NA, 6) } else {
              model <- lm(data.all[[j+2]] ~ site.index, data.all);
              temp <- summary(model);
              slope[j,] <- c(temp$coef[2,], anova(model)[1, "Mean Sq"], anova(model)[2, "Mean Sq"]);
            }
          }
          rownames(slope) <- levels(temp.data[, match(geno,names(temp.data))]);
          colnames(slope) <- c("Slope", "SE", "t-value", "Prob", "MSReg", "MSDev");
          slope <- data.frame(slope);
          slope <- subset(slope, subset = (is.na(SE)==F));
        } else { 
          slope <- NULL; 
        }
        data$traits[[i]]$analysis$mea$stability$slope <- slope;
        data$traits[[i]]$analysis$mea$stability$finlay.wilkinson <- "Stability Analysis using Finlay-Wilkinson Model";
      } else if(method == "shukla")
      {
        library(nlme);
        temp.data <- subset(means.GenoEnvCode, select = c(geno, env, trait.name));
        temp.data <- subset(temp.data, subset = (is.na(temp.data[,match(trait.name, names(temp.data))]) == F));
        #myformula1 <- paste(trait.name, " ~ 1 + ", geno, " + ", env, sep = "")
        #myformula2 <- paste("varIdent(form = ~ 1 |", geno,")", sep = "")
        command <- paste("gls(formula(", trait.name," ~ 1 + ",geno," + ", env,"), 
                           weights = varIdent(form = ~ 1 |", geno,"), data = temp.data)", sep = "");
        model1 <- eval(parse(text = command));
        data$traits[[i]]$analysis$mea$stability$shukla <- "Stability Analysis using Shukla's model"
        if (class(try(parse(text = "intervals(model1)"), TRUE)) == "try-error")
        { 
          data$traits[[i]]$analysis$mea$stability$par <- geterrmessage();  
        }else 
        {
          int <- intervals(model1);
          par <- as.data.frame(int$varStruct);
          par <- rbind(c(1,1,1), par);
          rownames(par) <- levels(temp.data[,match(geno, names(temp.data))]);
          sigma <- as.data.frame(int$sigma);
          par$lower <- par$lower*sigma[rownames(sigma) == "lower",];
          par$est.<- par$est.*sigma[rownames(sigma) == "est.",];
          par$upper <- par$upper*sigma[rownames(sigma) == "upper",];
          data$traits[[i]]$analysis$mea$stability$par <- par;
        }
        detach("package:nlme")
      } #end stmt else if(method == "shukla")
    } #end stmt if(is.null(data$traits[[i]]$analysis$mea))
  }#end stmt  for(i in 1:length(data$traits))
  return(data);
}