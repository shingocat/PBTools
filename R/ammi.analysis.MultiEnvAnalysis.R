ammi.analysis.MultiEnvAnalysis <- function
(
  data,
  number = TRUE,
  biplotPC12 = FALSE,
  biplotPC13 = FALSE,
  biplotPC23 = FALSE,
  ammi1 = FALSE,
  adaptMap = FALSE,
  ...
)
{
  if(!inherits(data, "MultiEnvAnalysis"))
    stop("\tError: The data argument should be of class MultiEnvAnalysis!\n");
  for(ii in 1:length(data$traits))
  {
    trait.name <- data$traits[[ii]]$name;
    if(is.null(data$traits[[ii]]$analysis$mea))
    {
      warning(paste("\tWarning: There are no Multi-Environment Analysis Outcomes on this trait", trait.name,".\n", sep=""));
      next;
    }
    if(is.null(data$traits[[ii]]$analysis$mea$ammi))
    {
      data$traits[[ii]]$analysis$mea$ammi <- list();
    } else
    {
      data$traits[[ii]]$analysis$mea$ammi <- NULL;
      data$traits[[ii]]$analysis$mea$ammi <- list();
    }
    genoEnvMeans <- data$traits[[ii]]$analysis$mea$means.GenoEnvCode
    ENV <- genoEnvMeans[,"CodedEnv"];
    GEN <- genoEnvMeans[,"CodedGeno"];
    REP <- data$traits[[ii]]$analysis$mea$harmonicMean;
    MSE <- data$traits[[ii]]$analysis$mea$MSE;
    trait.name <- paste(trait.name, "_LSMean", sep="");
    Y <- genoEnvMeans[,trait.name];
    yVar <- trait.name;
    ENV <- as.factor(ENV);
    GEN <- as.factor(GEN);
    nenv <- length(unique(ENV));
    ngen <- length(unique(GEN));
    minimo <- min(ngen, nenv);
    if (length(REP) > 1) {
      REP <- as.factor(REP);
      nrep <- length(unique(REP));
      modelo <- aov(Y ~ ENV + REP %in% ENV + GEN + ENV:GEN);
      mm <- anova(modelo);
      nn <- mm[2, ];
      mm[2, ] <- mm[3, ];
      mm[3, ] <- nn;
      row.names(mm)[2] <- "REP(ENV)";
      row.names(mm)[3] <- "GEN     ";
      mm[1, 4] <- mm[1, 3]/mm[2, 3];
      mm[1, 5] <- 1 - pf(mm[1, 4], mm[1, 1], mm[2, 1]);
      DFE <- df.residual(modelo);
      MSE <- deviance(modelo)/DFE;
      medy <- mean(Y,na.rm=TRUE);
    } else {
      DFE <- nenv * (ngen - 1) * (REP - 1);
      DFEa <- nenv * (REP - 1);
      nrep <- REP;
      modelo <- aov(Y ~ ENV + GEN + ENV:GEN);
      xx <- as.matrix(anova(modelo));
      xx <- rbind(xx[1, ], xx[1, ], xx[2:4, ]);
      xx[2, 1] <- DFEa;
      xx[2, 2:5] <- NA;
      xx[, 2] <- xx[, 2] * nrep;
      xx[, 3] <- xx[, 3] * nrep;
      xx[5, 1] <- DFE;
      xx[5, 3] <- MSE;
      xx[5, 2] <- MSE * DFE;
      xx[1, 4] <- NA;
      xx[3, 4] <- xx[3, 3]/MSE;
      xx[4, 4] <- xx[4, 3]/MSE;
      xx[1, 5] <- NA;
      xx[3, 5] <- 1 - pf(xx[3, 4], xx[3, 1], DFE);
      xx[4, 5] <- 1 - pf(xx[4, 4], xx[4, 1], DFE);
      row.names(xx)[1] <- "ENV     ";
      row.names(xx)[2] <- "REP(ENV)";
      medy <- mean(Y,na.rm=TRUE);
    }
    raw <- data.frame(ENV, GEN, Y);
    MEDIAS <- tapply(raw[, 3], raw[, c(1, 2)], mean);
    xx <- rownames(MEDIAS);
    yy <- colnames(MEDIAS);
    fila <- length(xx);
    col <- length(yy);
    total <- fila * col;
    x <- character(length = total);
    y <- character(length = total);
    z <- numeric(length = total);
    k <- 0;
    for (i in 1:fila) {
      for (j in 1:col) {
        k <- k + 1;
        x[k] <- xx[i];
        y[k] <- yy[j];
        z[k] <- MEDIAS[i, j];
      }
    }
    MEDIAS <- data.frame(ENV=x, GEN=y, Y=z);
    x <- MEDIAS[, 1];
    y <- MEDIAS[, 2];
    z <- MEDIAS[, 3];
    modelo2 <- lm(z ~ x + y);
    for (i in 1:length(z)) {
      if (is.na(z[i]))
        z[i] <- predict(modelo2, data.frame(x = MEDIAS[i,
                                                       1], y = MEDIAS[i, 2]));
    }
    MEDIAS <- data.frame(ENV = x, GEN = y, Y = z);
    modelo1 <- lm(Y ~ ENV + GEN, data = MEDIAS);
    residual <- modelo1$residuals;
    MEDIAS <- data.frame(MEDIAS, RESIDUAL = residual);
    mlabel <- names(MEDIAS);
    names(MEDIAS) <- c(mlabel[1:2], yVar, mlabel[4]);
    OUTRES <- MEDIAS[order(MEDIAS[, 1], MEDIAS[, 2]), ];
    OUTRES2 <- by(OUTRES[, 4], OUTRES[, c(2, 1)], function(x) sum(x));
    OUTMED <- by(OUTRES[, 3], OUTRES[, c(2, 1)], function(x) sum(x));
    s <- svd(OUTRES2);
    U <- s$u;
    L <- s$d;
    V <- s$v;
    L <- L[1:minimo];
    SS <- (L^2) * nrep;
    SUMA <- sum(SS);
    percent <- round(((1/SUMA) * SS) * 100, 1);
    minimo <- min(ngen, nenv);
    DFAMMI <- rep(0, minimo);
    acum <- DFAMMI;
    MSAMMI <- DFAMMI;
    F.AMMI <- DFAMMI;
    PROBF <- DFAMMI;
    acumula <- 0;
    for (i in 1:minimo) {
      DF <- (ngen - 1) + (nenv - 1) - (2 * i - 1);
      if (DF <= 0) break;
      DFAMMI[i] <- DF;
      acumula <- acumula + percent[i];
      acum[i] <- acum[i] + acumula;
      MSAMMI[i] <- SS[i]/DFAMMI[i];
      F.AMMI[i] <- round(MSAMMI[i]/MSE, 2);
      PROBF[i] <- round(1 - pf(F.AMMI[i], DFAMMI[i], DFE), 4);
    }
    SS <- round(SS, 6);
    MSAMMI <- round(MSAMMI, 6);
    SSAMMI <- data.frame(percent, acum, Df = DFAMMI, "Sum Sq" = SS,
                         "Mean Sq" = MSAMMI, "F value" = F.AMMI, Pr.F = PROBF);
    nssammi<-nrow(SSAMMI);
    SSAMMI<-SSAMMI[SSAMMI$Df>0,];
    nss<-nrow(SSAMMI);
    row.names(SSAMMI) <- paste("PC", 1:nss, sep = "");
    #cat("\nAnalysis\n")
    #print(SSAMMI)
    LL <- sqrt(diag(L));
    SCOREG <- U %*% LL;
    SCOREE <- V %*% LL;
    SCORES <- rbind(SCOREG, SCOREE);
    colnames(SCORES) <- paste("PC", 1:nssammi, sep = "");
    MSCORES <- SCORES[1:ngen, ];
    NSCORES <- SCORES[(ngen + 1):(ngen + nenv), ];
    MGEN <- data.frame(type = "GEN", Y = apply(OUTMED, 1, mean),
                       MSCORES);
    MENV <- data.frame(type = "ENV", Y = apply(OUTMED, 2, mean),
                       NSCORES);
    
    # added by NSALES
    commonRowNames<-intersect(row.names(MGEN),row.names(MENV));
    if (length(commonRowNames)>0) {
      row.names(MGEN)<-paste("G",row.names(MGEN),sep=""); #added by NSALES
      row.names(MENV)<-paste("E",row.names(MENV),sep=""); #added by NSALES
    }
    
    bplot <- rbind(MGEN, MENV);
    bplot<- bplot[,1:(nss+2)];
    mlabel <- names(bplot);
    names(bplot) <- c(mlabel[1], yVar, mlabel[c(-1, -2)]);
    maxy <- max(bplot[, 4]);
    miny <- min(bplot[, 4]);
    maxx <- max(bplot[, 3]);
    minx <- min(bplot[, 3]);
    row.names(bplot) <- c(row.names(MGEN), row.names(MENV));
    cp.name <- rownames(SSAMMI)[1:3];
    cp.per <- SSAMMI[1:3, 1];
    
    #create string to be used for filenames -added by NSALES
    forFilename<-"AMMI";
    
    # save PC scores to csv file - added by NSALES 
    pcScoreFilename<-paste(getwd(),"/", forFilename, "_PC_scores_",yVar,".csv",sep = "");
    bplot_print<-data.frame(Levels=rownames(bplot), bplot);
    write.csv(bplot_print, file=pcScoreFilename, row.names = FALSE);
    
    # Create new variable color_code
    Type.code <- as.numeric(bplot$type);
    
    # generate graphs
    if (biplotPC12) {
      if (ncol(bplot)>3) {
        png(filename = paste(getwd(),"/", forFilename, "_biplot_PC1_PC2_",yVar,".png",sep = ""));
        par(cex=0.8);
        plot(bplot[,3],bplot[,4],cex=0.8, xlab = "PC 1", ylab = "PC 2", frame = TRUE, 
             pch=" ",col=Type.code, main=paste(forFilename, " Biplot for ", gsub("_means", "", yVar), sep=""));
        if (number == TRUE) {
          text(MGEN[, 3], MGEN[, 4], cex = 0, text(MGEN[, 3], MGEN[, 4], 
                                                   labels = as.character(1:nrow(MGEN)), col = "blue"));
        }
        if (number == FALSE) {
          text(MGEN[, 3], MGEN[, 4], cex = 0, text(MGEN[, 3], MGEN[, 4], 
                                                   labels = row.names(MGEN), col = "blue"));
        }
        points(MENV[, 3], MENV[, 4], cex = 0, text(MENV[, 3], MENV[, 4], 
                                                   labels = row.names(MENV), col = "red"));
        abline(h = 0, v = 0, lty = 2.5, col = "green", lwd = 1);
        s <- seq(length(MENV[, 3]));
        arrows(0, 0, 0.9 * MENV[, 3][s], 0.9 * MENV[, 4][s], col = "brown", lwd = 1.8, 
               length = 0.1, code = 2);
        
        # overlay ellipse
        AMMI.contour(bplot, distance=.15, shape=20,col="black", lwd=2, lty=3);
        
        pc_percent<-paste("PC1=", cp.per[1], "%; PC2=", cp.per[2], "%", sep="");
        mtext(text = pc_percent, side=3, cex=0.8, adj=0);
        dev.off();
      }
    }
    
    if (biplotPC13) {
      if (ncol(bplot)>4) {
        png(filename = paste(getwd(),"/", forFilename, "_biplot_PC1_PC3_",yVar,".png",sep = ""));
        par(cex=0.8);
        plot(bplot[,3],bplot[,5],cex=0.8, xlab = "PC 1", ylab = "PC 3", 
             frame = TRUE, pch=" ",col=Type.code, main=paste(forFilename, " Biplot for ", 
                                                             gsub("_means", "", yVar), sep=""));
        if (number == TRUE) {
          text(MGEN[, 3], MGEN[, 5], cex = 0, text(MGEN[, 3], MGEN[, 5], 
                                                   labels = as.character(1:nrow(MGEN)), col = "blue"));
        }
        if (number == FALSE) {
          text(MGEN[, 3], MGEN[, 5], cex = 0, text(MGEN[, 3], MGEN[, 5], 
                                                   labels = row.names(MGEN), col = "blue"));
        }
        points(MENV[, 3], MENV[, 5], cex = 0, text(MENV[, 3], MENV[, 5], 
                                                   labels = row.names(MENV), col = "red"));
        abline(h = 0, v = 0, lty = 2.5, col = "green", lwd = 1);
        s <- seq(length(MENV[, 3]));
        arrows(0, 0, 0.9 * MENV[, 3][s], 0.9 * MENV[, 5][s], col = "brown", lwd = 1.8, 
               length = 0.1, code = 2);
        
        pc_percent<-paste("PC1=", cp.per[1], "%; PC3=", cp.per[3], "%", sep="");
        mtext(text = pc_percent, side=3, cex=0.8, adj=0);
        dev.off();
      }
    }
    
    if (biplotPC23) {
      if (ncol(bplot)>4) {
        png(filename = paste(getwd(),"/", forFilename, "_biplot_PC2_PC3_",yVar,".png",sep = ""));
        par(cex=0.8);
        plot(bplot[,4],bplot[,5],cex=0.8, xlab = "PC 2", ylab = "PC 3", frame = TRUE, pch=" ",
             col=Type.code, main=paste(forFilename, " Biplot for ", gsub("_means", "", yVar), sep=""));
        if (number == TRUE) {
          text(MGEN[, 4], MGEN[, 5], cex = 0, text(MGEN[, 4], MGEN[, 5], 
                                                   labels = as.character(1:nrow(MGEN)), col = "blue"));
        }
        if (number == FALSE) {
          text(MGEN[, 4], MGEN[, 5], cex = 0, text(MGEN[, 4], MGEN[, 5], 
                                                   labels = row.names(MGEN), col = "blue"));
        }
        points(MENV[, 4], MENV[, 5], cex = 0, text(MENV[, 4], MENV[, 5], 
                                                   labels = row.names(MENV), col = "red"));
        abline(h = 0, v = 0, lty = 2.5, col = "green", lwd = 1);
        s <- seq(length(MENV[, 3]));
        arrows(0, 0, 0.9 * MENV[, 4][s], 0.9 * MENV[, 5][s], col = "brown", lwd = 1.8, 
               length = 0.1, code = 2);
        
        pc_percent<-paste("PC2=", cp.per[2], "%; PC3=", cp.per[3], "%", sep="");
        mtext(text = pc_percent, side=3, cex=0.8, adj=0);
        dev.off();
      }
    }  
    
    if (ammi1) {
      # ammi1 biplot
      png(filename = paste(getwd(),"/AMMI1_biplot_",yVar,".png",sep = ""));
      par(cex=0.8);
      plot(bplot[,yVar],bplot$PC1,cex=0.8, main="AMMI1 Biplot", frame=TRUE,xlab=yVar, 
           ylab="PC1", pch=" ",col=Type.code);
      text(MGEN[,2], MGEN[,3], labels=row.names(MGEN),col="blue");
      text(MENV[,2], MENV[,3],labels=row.names(MENV),col=2);
      
      MEANS<-mean(bplot[,yVar]);
      abline(h=0,v= MEANS,lty=2,col="green");
      
      pc_percent2<-paste("PC1=", cp.per[1], "%", sep="");
      mtext(text = pc_percent2, side=3, cex=0.8, adj=0);
      dev.off();
    }
    
    if (adaptMap) {
      # adaptation map
      png(filename = paste(getwd(),"/adaptationMap_",yVar,".png",sep = ""));
      
      geno <- bplot[bplot$type=="GEN",c(2:3)];
      geno <- data.frame(name=rownames(geno),geno);
      env <- bplot[bplot$type=="ENV", c(2:3)];
      env <- data.frame(name=rownames(env),env);
      
      env.PC1 <- sort(bplot[bplot$type=="ENV",3]);
      
      pred.yld <- as.matrix(rep(0,nenv*ngen),nrow=nenv,ncol=ngen);
      dim(pred.yld) <- c(nenv,ngen);
      
      for (i in 1:ngen)   {
        for (j in 1:nenv)
          pred.yld[j,i] <- geno[i,2] + geno$PC1[i]*env.PC1[j];
      }
      pred.yld <- as.data.frame(pred.yld);
      names(pred.yld) <- geno$name;
      pred.yld <- cbind(pred.yld, env.PC1);
      
      MIN.Y <- min(pred.yld[1:ngen]);
      MAX.Y <- max(pred.yld[1:ngen]);
      
      MIN.X <- min(pred.yld[ngen+1]);
      MAX.X <- max(pred.yld[ngen+1]);
      
      par(mar=c(5,4,7,1)+0.1);
      Y3 <- data.frame(V1=env$PC1, env$name);
      Y3.sort <- Y3[order(Y3$V1),];
      
      par(cex=0.8);
      #graph the lines
      for (i in 1:ngen)  {
        if (i==1)  {
          yLabel<-paste("Predicted ", gsub("_means", "", yVar), sep="");
          XLabel<-paste("Environment IPCA1 ", "(PC1=", cp.per[1], "%)", sep="");
          plot(pred.yld[[ngen+1]],pred.yld[[i]], pch=" ", ylim=c(MIN.Y,MAX.Y), 
               xlim=c(MIN.X,MAX.X), main=, xlab=XLabel, ylab=yLabel);
          abline(lm(pred.yld[[i]]~pred.yld[[ngen+1]]),col=i);   
        }
        abline(lm(pred.yld[[i]]~pred.yld[[ngen+1]]),col=i);
      }
      
      title("Adaptation Map", line=4);
      text(rep(pred.yld[[nenv,ngen+1]]-abs(0.03*max(pred.yld[[nenv,ngen+1]])), ngen),
           pred.yld[nenv,c(1:ngen)], labels=geno$name,col=1);
      axis(side=3, at=Y3.sort[,1], labels=Y3.sort[,2],las=3, col=1);
      dev.off();
    }
    data$traits[[ii]]$analysis$mea$ammi$genXenv <- OUTRES2;
    data$traits[[ii]]$analysis$mea$ammi$analysis <- SSAMMI;
    data$traits[[ii]]$analysis$mea$ammi$means <- MEDIAS;
    data$traits[[ii]]$analysis$mea$ammi$biplot <- bplot;
    #   return(list(genXenv=OUTRES2, analysis=SSAMMI, means=MEDIAS, biplot=bplot))
  }#end stmt for(ii in 1:length(data$traits))
  return(data);
}