# 
#
#

doChisqTest <- function
(
  data,
  inc.ht = TRUE,
  pop.min = 20, 
  ref.matrix = NULL,  
  simulate.p.value = TRUE,
  B = 1000,
  ...
)
{
  UseMethod("doChisqTest");
}

doChisqTest.GenotypicData <- function
(
  data,
  inc.ht = TRUE,
  pop.min = 20, 
  ref.matrix = NULL,  
  simulate.p.value = TRUE,
  B = 1000,
  ...
)
{
  if(!inherits(data, "GenotypicData"))
    stop("\tError: The data argument should be of class GenotypicData!\n");
  if(missing(inc.ht))
    inc.ht <- TRUE;
  if(!is.logical(inc.ht))
    stop("\tError: The inc.ht should be logical value!\n");
  if(missing(pop.min))
    pop.min <- 20;
  if(!missing(ref.matrix) && !inherits(ref.matrix, "GenotypicData"))
    stop("\tError: The ref.matrix should be of class GenotypicData!\n");
  if(missing(simulate.p.value))
    simulate.p.value <- TRUE;
  if(!is.logical(simulate.p.value))
    stop("\tError: The simulate.p.value should be logical value!\n");
  if(missing(B))
    B <- 1000;
  dp.code <- data$processed$dp.code;
  rp.code <- data$processed$rp.code;
  ht.code <- data$processed$ht.code;
  na.code <- data$processed$na.code;
  BCn <- data$BCn;
  Fn <- data$Fn;
  geno.name <- data$geno.name
  marker.matrix <- data$processed$data;
  lineNames <- marker.matrix[,geno.name];
  marker.matrix <- marker.matrix[,-which(geno.name %in% colnames(marker.matrix))];
  rownames(marker.matrix) <- lineNames;
  marker.matrix <- as.matrix(marker.matrix);
  markers.num <- ncol(marker.matrix);
  markers.names <- colnames(marker.matrix);
  if(is.null(markers.names))
  {
    stop("\tError: The marker names could not be null!");
    #markers.names <- paste("M", 1:markers.num, sep="")
  } 
  lines.num <- nrow(marker.matrix);
  marker.matrix <- apply(marker.matrix,2,factor);	
  
  if(lines.num < pop.min)
  {
    stop(paste("\tError: There are minimum requriement of populaiton size is ", pop.min, " using chisq test!\n", sep=""));
  }
  
  # Cheacking whether there are unaccepted code in marker matrix
  marker.code <- c(dp.code, rp.code, ht.code, na.code);
  checking.result <- apply(marker.matrix, 2, function(x) all(x %in% marker.code));
  if(!all(checking.result))
  {
    stop(paste("\tError: There are illegal code in the marker : ", names(checking.result[!checking.result]), "!\n", sep = ""));
  }
  
  do.ref <- FALSE;
  refFreq <- function(x, dp.code, ht.code, rp.code)
  {
    result <- c();
    x <- factor(x);
    x.table <- table(x, useNA = "no");
    rp.num <- x.table[as.character(rp.code)];
    dp.num <- x.table[as.character(dp.code)];
    ht.num <- x.table[as.character(ht.code)];
    if(is.na(rp.num))
      rp.num <- 0;
    if(is.na(dp.num))
      dp.num <- 0;
    if(is.na(ht.num))
      ht.num <- 0;
    sample.size <- sum(c(dp.num, ht.num, rp.num));
    freq <- c(dp.num, ht.num, rp.num)/sample.size;
    result <- freq;
    names(result) <- c("Donor", "Heterozygous", "Recurrent");
    return(result);
  }
  # If the referenced marker matrix is not null, computed the referenced frequence of each marker.
  if(!missing(ref.matrix))
  {
    do.ref <- TRUE;
    ref.geno.name <- ref.matrix$geno.name;
    ref.matrix <- ref.matrix$processed$data;
    refLineNames <- ref.matrix[,ref.geno.name];
    ref.matrix <- ref.matrix[,-which(ref.geno.name %in% colnames(ref.matrix))];
    rownames(ref.matrix) <- refLineNames;
    ref.matrix <- as.matrix(ref.matrix);
    ref.markers.names <- colnames(ref.matrix);
    if(is.null(ref.markers.names))
      stop("\tError: The reference matrix marker names could not be empty!");
    # checking whether there are illegal code in referenced marker matrix
    checking.result <- apply(ref.matrix, 2, function(x) all(x %in% marker.code));
    if(!all(checking.result))
    {
      stop(paste("\tError: There are illegal code in the referenced matrix : ", names(checking.result[!checking.result]), "!\n", sep=""));
    }
    ref.markers.freq <- t(apply(ref.matrix, 2, refFreq, dp.code, ht.code, rp.code));
  }
  
  # using for loop, it is more slow than using apply method!
  #
  result <- c();
  temp <- c(); 
  for(i in 1:markers.num)
  {
    geno <- marker.matrix[,i];
    geno <- factor(geno);
    geno.table <- table(geno, useNA = "no");
    
    dp.num <- geno.table[as.character(dp.code)];	
    if(is.na(dp.num))
      dp.num <- 0;
    rp.num <- geno.table[as.character(rp.code)];
    if(is.na(rp.num))
      rp.num <- 0;
    ht.num <- geno.table[as.character(ht.code)];
    if(is.na(ht.num))
      ht.num <- 0;
    
    sample.size <- sum(c(dp.num, rp.num, ht.num));
    expect.freq <- getExpectGenotypicFreq(BCn, Fn, 2);
    
    
    if(dp.num == 0 && ht.num == 0)
    {
      # TypeI denoted as the genotypic value of this marker are always recurrent, and would not compute this kind of marker;
      # treated as no effect marker. In the other words, donor and heterozygous number is zero!	
      if(do.ref)
      {
        if(markers.names[i] %in% ref.markers.names)
        {
          freq <- ref.markers.freq[which(markers.names[i] == ref.markers.names),];
          temp <- c(markers.names[i], "TypeI", sample.size, dp.num, ht.num, rp.num, "N", paste(freq, collapse=":"), NA, NA);
        } else
        {
          temp <- c(markers.names[i], "TypeI", sample.size, dp.num, ht.num, rp.num, "Y", paste(expect.freq, collapse=":"), NA, NA);
        }
      }else
      {
        temp <- c(markers.names[i], "TypeI", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
      }
      names(temp) <- c();
      result <- rbind(result, temp);
      next;
    } else if(rp.num == 0 || (dp.num > rp.num && rp.num < 3))
    {
      # TypeII denoted as the genotype of marker are always of donor and heterozygous genotype or the num of donor marker genotype is large than
      # the num of recurrent marker genotype which are less than 3 individuals;
      # , and would not compute this kind of marker;
      # treated as potential linkage marker.
      if(do.ref)
      {
        if(markers.names[i] %in% ref.markers.names)
        {
          freq <- ref.markers.freq[which(markers.names[i] == ref.markers.names),];
          temp <- c(markers.names[i], "TypeII", sample.size, dp.num, ht.num, rp.num, "N", paste(freq, collapse=":"), NA, NA);
        } else
        {
          temp <- c(markers.names[i], "TypeII", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
        }
      }else
      {
        temp <- c(markers.names[i], "TypeII", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
      }
      names(temp) <- c();
      result <- rbind(result, temp);
      next;
    } else if(floor(sample.size * expect.freq[1]) < 1)
    {
      # TypeIII denoted as the expected least genotype, i.e. donor genotype, are less than 1 on this sample size
      # could not process the chisq testing!
      if(do.ref)
      {
        if(markers.names[i] %in% ref.markers.names)
        {
          freq <- ref.markers.freq[which(markers.names[i] == ref.markers.names),];
          temp <- c(markers.names[i], "TypeIII", sample.size, dp.num, ht.num, rp.num, "N", paste(freq, collapse=":"), NA, NA);
        } else
        {
          temp <- c(markers.names[i], "TypeIII", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
        }
      }else
      {
        temp <- c(markers.names[i], "TypeIII", sample.size, dp.num, ht.num, rp.num,"Y", paste(expect.freq, collapse=":"), NA, NA);
      }
      names(temp) <- c();
      result <- rbind(result, temp);
      next;
    } else
    {
      # ThypeIV  denoted as the normal case for chisq testing!
      #
      # if the argument of inc.ht is TRUE, than the heterozygous will add to donor type;
      if(inc.ht)
      {
        data <- c(dp.num + ht.num, rp.num);
      } else
      {
        data <- c(dp.num, rp.num);
      }
      
      if(do.ref)
      {
        if(markers.names[i] %in% ref.markers.names)
        {
          freq <- ref.markers.freq[which(markers.names[i] == ref.markers.names), ];
          if(inc.ht)
          {
            prob <- c(freq[1] + freq[2], freq[3]) / sum(freq);
          } else
          {
            prob <- c(freq[1], freq[3]) / sum(freq[c(1,3)]);
          }
          
          # if one of the referenced marker freq is equal to zero, and then the prob is used expected.freq.
          if(any(prob == 0))
          {
            if(inc.ht)
            {
              prob <- c(expect.freq[1] + expect.freq[2], expect.freq[3]) / sum(expect.freq);
            } else
            {
              prob <- c(expect.freq[1], expect.freq[3]) / sum(expect.freq[c(1,3)]);
            }
            outcomes <- chisq.test(data, p = prob, simulate.p.value = simulate.p.value, B = B);
            statistic <- round(outcomes$statistic, digits = 2);
            p.value <- round(outcomes$p.value, digits = 4);
            names(p.value) <- c();
            temp <- c(markers.names[i], "TypeV", sample.size, dp.num, ht.num, rp.num, "Y", paste(expect.freq, collapse=":"), statistic, p.value);
            names(temp) <- c();
            result <- rbind(result, temp);
            next;
          }
          outcomes <- chisq.test(data, p = prob, simulate.p.value = simulate.p.value, B = B);
          statistic <- round(outcomes$statistic, digits = 2);
          p.value <- round(outcomes$p.value, digits = 4);
          names(p.value) <- c();
          temp <- c(markers.names[i], "TypeIV", sample.size, dp.num, ht.num, rp.num, "N", paste(freq, collapse=":"), statistic, p.value);
        } else
        {
          if(inc.ht)
          {
            prob <- c(expect.freq[1] + expect.freq[2], expect.freq[3]) / sum(expect.freq);
          } else
          {
            prob <- c(expect.freq[1], expect.freq[3]) / sum(expect.freq[c(1,3)]);
          }
          outcomes <- chisq.test(data, p = prob, simulate.p.value = simulate.p.value, B = B);
          statistic <- round(outcomes$statistic, digits = 2);
          p.value <- round(outcomes$p.value, digits = 4);
          names(p.value) <- c();
          temp <- c(markers.names[i], "TypeIV", sample.size, dp.num, ht.num, rp.num, "Y", paste(expect.freq, collapse=":"), statistic, p.value);
        }
      }else
      {
        if(inc.ht)
        {
          prob <- c(expect.freq[1] + expect.freq[2], expect.freq[3]) / sum(expect.freq);
        } else
        {
          prob <- c(expect.freq[1], expect.freq[3]) / sum(expect.freq[c(1,3)]);
        }
        outcomes <- chisq.test(data, p = expect.freq[c(1,3)]/sum(expect.freq[c(1,3)]), simulate.p.value = simulate.p.value, B = B);
        statistic <- round(outcomes$statistic, digits = 2);
        p.value <- round(outcomes$p.value, digits = 4);
        names(p.value) <- c();
        temp <- c(markers.names[i], "TypeIV", sample.size, dp.num, ht.num, rp.num, "Y", paste(expect.freq, collapse=":"), statistic, p.value);
      }
      names(temp) <- c();
      result <- rbind(result, temp);
      next;
    }
  }
  result <- as.data.frame(result);
  colnames(result) <- c("Marker", "Type", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp","ExptFreq","Freq(DP:HT:RP)","Chisq", "P.Value");
  rownames(result) <- c();
  
  
  outcomes <- list();
  outcomes$result <- result;
  outcomes$TypeI <- list(); # for all the genotype of marker are recurrent type or the num of donor genotype are less than 3;
  outcomes$TypeI$Markers <- result[result$Type == "TypeI", "Marker"];
  names(outcomes$TypeI$Markers) <- c();
  outcomes$TypeI$outcomes <- result[result$Type == "TypeI", c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp")];
  outcomes$TypeI$counts <- length(outcomes$TypeI$Markers);
  outcomes$TypeI$Messages <- c("For all the genotype of marker are recurrent type or the num of donor genotype are less than 3;");
  outcomes$TypeII <- list(); # for all the genotype of marker are donor type or the num of recurrent genotype are less than 3, potential linkage marker!;
  outcomes$TypeII$Markers <- result[result$Type == "TypeII", "Marker"];
  names(outcomes$TypeII$Markers) <- c();
  outcomes$TypeII$outcomes <- result[result$Type == "TypeII", c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp")];
  outcomes$TypeII$counts <- length(outcomes$TypeII$Markers);
  outcomes$TypeII$Messages <- c("For all the genotype of marker are donor type or the num of recurrent genotype are less than 3, potential linkage marekr!");
  outcomes$TypeIII <- list(); # for the minimum of sample size for the num of donor genotype are large than 1 are not reached;
  outcomes$TypeIII$Markers <- result[result$Type == "TypeIII", "Marker"];
  names(outcomes$TypeIII$Markers) <- c();
  outcomes$TypeIII$outcomes <- result[result$Type == "TypeIII", c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp")];
  outcomes$TypeIII$counts <- length(outcomes$TypeIII$Markers);
  outcomes$TypeIII$Messages <- c("For the minimum of sample size for the num of donor genotype are large than 1 are not reached;");
  outcomes$TypeIV <- list(); # for the normal chisq case;
  outcomes$TypeIV$Markers <- result[result$Type == "TypeIV", "Marker"];
  names(outcomes$TypeIV$Markers) <- c();
  outcomes$TypeIV$outcomes <- result[result$Type == "TypeIV",  c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp", "Chisq", "P.Value")];
  outcomes$TypeIV$counts <- length(outcomes$TypeIV$Markers);
  outcomes$TypeIV$Messages <- c("For the normal chisq case;");
  outcomes$TypeV <- list();
  outcomes$TypeV$Markers <- result[result$Type == "TypeV", "Marker"];
  names(outcomes$TypeV$Markers) <- c();
  outcomes$TypeV$outcomes <- result[result$Type == "TypeV", c("Marker", "Samples", "Obs.Dp", "Obs.Ht", "Obs.Rp", "Chisq", "P.Value")];
  outcomes$TypeV$counts <- length(outcomes$TypeV$Markers);
  outcomes$TypeV$Messages <- c("For genotypic frequence of donor and/or heterozygous is zero on the referenced marker, it will be compute chisq test using expected frequency!");
  class(outcomes) <- "ChisqOnIL";
  return(outcomes);
}
