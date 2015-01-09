# Get the defualt Orthogonal Contrast Matrix 
# accroding to genes number
#
#
getDefGeneOrthogonalContrastMatrix <- function(
  gene.no = 2,
  ...
)
{
  UseMethod("getDefGeneOrthogonalContrastMatrix");
}

getDefGeneOrthogonalContrastMatrix.default <- function(
  gene.no = 2,
  ...
)
{
  #--- default is two genes ---#
  if(missing(gene.no))
    gene.no <- 2;
  if(!is.numeric(gene.no))
  {
    warning("\tWarnings: Only accept large than 2 positive integer value for gene.no! return two genes Orthogonal Contrast Matrix.\n");
    gene.no <- 2;
  }
  if(gene.no <= 1)
  {
    warning("\tWarnings: Only accept large than 2 positive integer value for gene.no! return two genes Orthogonal Contrast Matrix.\n");
    gene.no <- 2;
  }
  gene.no <- as.integer(gene.no);
  #--- uppercase gene label---#
  gene.ucl <- LETTERS[1:gene.no];
  #--- lowercase gen label---#
  gene.lcl <- letters[1:gene.no];
  
  #---build the coefficient dataframe---#
  coef.df <- matrix(NA, ncol = gene.no + 3^gene.no + 1, nrow = 3^gene.no);
  coef.df <- as.data.frame(coef.df);
  #---build the header names---#
  header.names <- c();
  #---additive codes---#
  a.codes <- paste("a", 1:gene.no, sep = "");
  #---dominance codes---#
  d.codes <- paste("d", 1:gene.no, sep = "");
  #---interaction codes---#
  int.eff.codes <- list();
  for(i in 1:gene.no)
  {
    int.eff.codes[[i]] <- c(paste("a",i,sep=""),paste("d",i,sep=""));
  }
  int.codes <- c();
  if(gene.no == 2)
  {
    int.codes <- levels(interaction(int.eff.codes, sep=" ", lex.order = T));
  } else
  {
    for(i in 2:gene.no)
    {
      combs <- combn(gene.no, i);
      for(j in 1:ncol(combs))
      {
        int.codes <- c(int.codes, levels(interaction(int.eff.codes[combs[,j]], sep=" ", lex.order =T)));
      }
    }
  } 
  gene.label <- paste("Gene", LETTERS[1:gene.no], sep = "");
  header.names <- c("Genotype", gene.label, "Mean", a.codes, d.codes, int.codes);
  colnames(coef.df) <- header.names; 
  #---get all genotype combination---#
  gene.comb <- list();
  for(i in 1:gene.no)
  {
    gene.comb[[i]] <- c(paste(gene.ucl[i], gene.ucl[i], sep = ""),paste(gene.ucl[i], gene.lcl[i], sep = ""),
                        paste(gene.lcl[i], gene.lcl[i], sep = ""));
  }
  geno.levels <- levels(interaction(gene.comb, sep=" ", lex.order = T));
  coef.df[[1]] <- geno.levels;
  #--- adding each gene label---#
  for(i in 1:gene.no)
  {
    coef.df[[i + 1]] <- unlist(lapply(strsplit(coef.df[[1]], " "), FUN = function(x){x[i];}))
  }
  #--- reformat the genotype column---#
  geno.levels <- levels(interaction(gene.comb, sep="", lex.order=T));
  coef.df[, "Genotype"] <- geno.levels;
  #---adding coefficients to data frame---#
  #---the mean coefficients---#
  coef.df[, "Mean"] <- rep(1, times = nrow(coef.df));
  #---the additive coefficients---#
  for(i in 1:length(a.codes))
  {
    #--- row by row to set the value---#
    for(j in 1:nrow(coef.df))
    {
      if(coef.df[j,gene.label[i]] == paste(gene.lcl[i], gene.lcl[i], sep = ""))
      {
        coef.df[j, a.codes[i]] <- -1;
      }
      else if(coef.df[j,gene.label[i]] == paste(gene.ucl[i], gene.ucl[i], sep = ""))
      {
        coef.df[j, a.codes[i]] <- 1;
      }
      else
      {
        coef.df[j, a.codes[i]] <- 0;
      }
    }
  }
  #---the dominance coefficients---#
  for(i in 1:length(d.codes))
  {
    #--- row by row set the value---#
    for(j in 1:nrow(coef.df))
    {
      if(coef.df[j,gene.label[i]] == paste(gene.lcl[i], gene.lcl[i], sep = "") ||
           coef.df[j,gene.label[i]] == paste(gene.ucl[i], gene.ucl[i], sep = "")
         )
      {
        coef.df[j, d.codes[i]] <- 0;
      }else
      {
        coef.df[j, d.codes[i]] <- 1;
      }
    }
  }
  #---the interaction coefficients---#
  for(i in 1:length(int.codes))
  {
    #---row by row set the value---#
    coef.df[,int.codes[i]] <- apply(coef.df[,unlist(strsplit(int.codes[i], " "))],1,prod);
  }
  #---rename the int.codes and assign the column name to coef.df---#
  #---interaction codes---#
  int.eff.codes <- list();
  for(i in 1:gene.no)
  {
    int.eff.codes[[i]] <- c(paste("a",i,sep=""),paste("d",i,sep=""));
  }
  int.codes <- c();
  if(gene.no == 2)
  {
    int.codes <- levels(interaction(int.eff.codes, sep="", lex.order = T));
  } else
  {
    for(i in 2:gene.no)
    {
      combs <- combn(gene.no, i);
      for(j in 1:ncol(combs))
      {
        int.codes <- c(int.codes, levels(interaction(int.eff.codes[combs[,j]], sep="", lex.order =T)));
      }
    }
  } 
  header.names <- c("Genotype", gene.label, "Mean", a.codes, d.codes, int.codes);
  colnames(coef.df) <- header.names; 
  
  #---get the reverse matrix of coefficients---#
  coef.matrix <- coef.df[ ,(length(gene.label) + 2):ncol(coef.df)];
  rownames(coef.matrix) <- coef.df[,"Genotype"];
  orthogonal.contrast.matrix <- solve(coef.matrix);
  #---remove Mean contrast coefficient---#
  orthogonal.contrast.matrix <- orthogonal.contrast.matrix[-which(rownames(orthogonal.contrast.matrix) %in% "Mean"),];

  outcomes <- list();
  outcomes$gene.number <- gene.no;
  outcomes$coeff.table <- coef.df;
  outcomes$orthogonal.contrast.matrix <- orthogonal.contrast.matrix;
  class(outcomes) <- "DefaultOrthogonalContrastMatrix";
  return(outcomes);
}


print.DefaultOrthogonalContrastMatrix <- function(
  x, 
  level = c(1,2),
  ...
  )
{
  if(missing(level))
    level <- 1;
  if(!is.numeric(level))
    level <- 1;
  if(!level %in% c(1,2))
    level <- 1;
  level <- as.integer(level);
  cat("The Orthogonal Contrast Matrix:\n");
  cat("Gene Number: ", x$gene.number,"\n");
  cat(rep("*",times=50), sep = "");
  cat("\n");
  if(level == 1)
  {
    print(x$orthogonal.contrast.matrix, row.names = FALSE);
  } else
  {
    cat("The coefficients table of all genotype:\n");
    print(x$coeff.table,row.names = FALSE);
    cat(rep("-",times=50),sep = "");
    cat("\n");
    cat("The orthogonal contrast matrix:\n")
    print(x$orthogonal.contrast.matrix, row.names = FALSE);
  }
}