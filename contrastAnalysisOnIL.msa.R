###############################################################################
# TODO: Perform contrast analysis
# 
# model - a lmerMod class object
# contrastOpt - contrast option; c("control", "user")
# genContrast - NULL or a contrast matrix; required if contrastOpt == "Custom"
# envContrast - NULL or a contrast 
# controlLevel - NULL or control level
# Returned Value: a data frame
#
# Author: mqin
# Date: Sep 10, 2014
# FileName: contrastAnalysis.pyramidedline.R
###############################################################################

contrastAnalysisOnIL.msa <- function(
		msaOutput, 
		contrastOpt = c("RecurrentParent", "Custom", "Default"), 
		recurrentParent = NULL, 
		genoContrast = NULL, 
		envContrast = NULL
	) UseMethod("contrastAnalysisOnIL.msa");
contrastAnalysisOnIL.msa.default <- function(
		msaOutput, 
		contrastOpt = c("RecurrentParent", "Custom", "Default"), 
		recurrentParent = NULL, 
		genoContrast = NULL, 
		envContrast = NULL
) 
{
	contrastOpt <- match.arg(contrastOpt);
	if(is.null(msaOutput))
	{
		stop("\tError: The argument of msaOutput could not be null!\n");
	}
	if(msaOutput$analysisType != "MultiSiteAnalysis")
	{
		stop("\tError: The argument of msaOutput should be the outcomes of MultiSiteAnalysis!\n");
	}
	
	result <- list();
	library(phia);
	library(lme4);
	
	if(contrastOpt == "RecurrentParent")
	{
		if(is.null(recurrentParent))
		{
			stop("\tError: The argument of recurrentParent should be null when contrastOpt is selected RecurrentParent!\n");
		}
		if(!is.vector(recurrentParent,mode = "character"))
		{
			stop("\tError: The argument of recurrentParent should be vector when contrastOpt is selected RecurrentParent!\n");
		}
		if(length(recurrentParent) != 1)
		{
			stop("\tError: The argument of recurrentParent should be of length one when contrastOpt is selected RecurrentParent!\n");
		}
		
		respvar.number  <- length(msaOutput$output);
		# -- checking conditions --- #
		for( i in 1:respvar.number)
		{
			is.envFixed <- msaOutput$output[[i]]$envfixed;
			model <- msaOutput$output[[i]]$model;
			trmtLevels <- levels(msaOutput$output[[i]]$model@frame[,2]);
			envLevels <- levels(msaOutput$output[[i]]$model@frame[,3]);
			# -- checking wheather the specified recurrentParent exist!
			if(!(recurrentParent %in% trmtLevels))
			{
				stop("\tError: The specified recurrentParent does not include in the model!\n");
			}
			if(is.envFixed)
			{
				if(!is.null(envContrast))
				{	
					if(!is.matrix(envContrast))
					{
						stop("\tError: The specified envContrast should be matrix fromat!\n");
					}
					envColContrastNames <- colnames(envContrast);
					if(length(envLevels) != length(envColContrastNames))
					{
						stop("\tError: The levels of specified envContrast is not equal to the levels of env factor on model!\n");
					}
					if(!(is.null(envColContrastNames)))
					{
						if(!all(envLevels %in% envColContrastNames))
						{
							stop("\tError: The column names of specified envContrast is not same as the names of env factor on model!\n");
						}
					} else
					{
						stop("\tError: The column names of specified envContrast could not be null and should be same as environment level names!\n");
					}
					# checking the env coefficients to zero!
					if(!all(rowSums(envContrast) == 0))
					{
						stop("\tError: Sum of each row should be equal to zero on envContrast!\n");
					}
				}
			} else
			{
				if(!is.null(envContrast))
				{
					warning("\tWarning: The specified envContrast could not be used when the env factor in the model is random!\n");
				}
			}	
		}
		
		for(i in 1:respvar.number)
		{
			result[[i]] <- list();
			result[[i]]$respvar <- msaOutput$output[[i]]$respvar;
			
			is.envFixed <- msaOutput$output[[i]]$envfixed;
			result[[i]]$envfixed <- is.envFixed;
			model <- msaOutput$output[[i]]$model;
			result[[i]]$model <- model;
			
			genoFactorName <- names(msaOutput$output[[i]]$model@frame)[2];
			envFactorName <- names(msaOutput$output[[i]]$model@frame)[3];
			trmtLevels <- levels(msaOutput$output[[i]]$model@frame[,2]);
			trmtLevels.without.rp <- trmtLevels[(trmtLevels != recurrentParent)];
			trmtContrast <- diag(1, nrow = length(trmtLevels.without.rp), ncol = length(trmtLevels.without.rp));
			trmtContrast <- cbind(trmtContrast, -1);
			colnames(trmtContrast) <- c(trmtLevels.without.rp, recurrentParent);
			rownames(trmtContrast) <- interaction(trmtLevels.without.rp, " - HJX", sep ="");
			trmtContrast <- as.data.frame(t(trmtContrast));
			trmtContrast <- list(genoFactorName = trmtContrast);
			names(trmtContrast) <- genoFactorName;
			if(is.envFixed)
			{
				if(is.null(envContrast))
				{	
					result[[i]]$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
					result[[i]]$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast, across = envFactorName);
				} else
				{
					result[[i]]$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
					tempEnvContrast = envContrast;
					tempEnvContrast = as.data.frame(t(tempEnvContrast));
					trmtContrast$tempEnvContrast <- tempEnvContrast;
					names(trmtContrast) <- c(genoFactorName, envFactorName);
					result[[i]]$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast);
				}
			} else
			{
				result[[i]]$contrastOnGeno <- testInteractions(model, custom=trmtContrast);
				result[[i]]$contrastAcrossEnv <- NA;
			}
		}
	} else if(contrastOpt == "Custom")
	{
		if(is.null(genoContrast))
		{
			stop("\tError: The argument of genoContrast could not be null when contrastOpt is Custom!\n");
		}
		if(!is.matrix(genoContrast))
		{
			stop("\tError: The genoContrast should be matrix format!\n");
		}
		respvar.number  <- length(msaOutput$output);
		# -- checking conditions --- #
		for( i in 1:respvar.number)
		{
			is.envFixed <- msaOutput$output[[i]]$envfixed;
			model <- msaOutput$output[[i]]$model;
			trmtLevels <- levels(msaOutput$output[[i]]$model@frame[,2]);
			envLevels <- levels(msaOutput$output[[i]]$model@frame[,3]);
			genoColContrastNames <- colnames(genoContrast);
			
			# -- checking wheather the specified genoContrast levels exist!
			if(!(is.null(genoColContrastNames)))
			{
				if(length(trmtLevels) != length(genoColContrastNames))
				{
					stop("\tError: The levels of specified genoContrast is not equal to the levels of geno factor on model!\n");
				}
				
				if(!all(genoColContrastNames %in% trmtLevels))
				{
					stop("\tError: Some of the specified genoContrast levels are not in geno factor on model!\n");
				}
			}else
			{
				stop("\tError: The column names of specified genoContrast could not be null and should be same as genotypic level names!\n");
			}
			
			# checking the env coefficients to zero!
			if(!all(rowSums(genoContrast) == 0))
			{
				stop("\tError: Sum of each row should be equal to zero on genoContrast!\n");
			}
			
			if(is.envFixed)
			{
				if(!is.null(envContrast))
				{	
					if(!is.matrix(envContrast))
					{
						stop("\tError: The specified envContrast should be matrix fromat!\n");
					}
					
					envColContrastNames <- colnames(envContrast);
					
					if(!(is.null(envColContrastNames)))
					{
						if(length(envLevels) != length(envColContrastNames))
						{
							stop("\tError: The levels of specified envContrast is not equal to the levels of env factor on model!\n");
						}
						
						if(!all(envLevels %in% envColContrastNames))
						{
							stop("\tError: The column names of specified envContrast is not same as the names of env factor on model!\n");
						}
					} else
					{
						stop("\tError: The column names of specified envContrast could not be null and should be same as environmental level names!\n");
					}
					# checking the env coefficients to zero!
					if(!all(rowSums(envContrast) == 0))
					{
						stop("\tError: Sum of each row should be equal to zero on envContrast!\n");
					}
				}
			} else
			{
				if(!is.null(envContrast))
				{
					warning("\tWarning: The specified envContrast could not be used when the env factor in the model is random!\n");
				}
			}	
		}
		
		for(i in 1:respvar.number)
		{
			result[[i]] <- list();
			result[[i]]$respvar <- msaOutput$output[[i]]$respvar;
			
			is.envFixed <- msaOutput$output[[i]]$envfixed;
			result[[i]]$envfixed <- is.envFixed;
			model <- msaOutput$output[[i]]$model;
			result[[i]]$model <- model;
			
			genoFactorName <- names(msaOutput$output[[i]]$model@frame)[2];
			envFactorName <- names(msaOutput$output[[i]]$model@frame)[3];
			trmtLevels <- levels(msaOutput$output[[i]]$model@frame[,2]);
			tempGenoContrast <- genoContrast;
			trmtContrast <- as.data.frame(t(tempGenoContrast));
			trmtContrast <- list(genoFactorName = trmtContrast);
			names(trmtContrast) <- genoFactorName;
			if(is.envFixed)
			{
				if(is.null(envContrast))
				{	
					result[[i]]$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
					result[[i]]$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast, across = envFactorName);
				} else
				{
					result[[i]]$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
					tempEnvContrast = envContrast;
					tempEnvContrast = as.data.frame(t(tempEnvContrast));
					trmtContrast$tempEnvContrast <- tempEnvContrast;
					names(trmtContrast) <- c(genoFactorName, envFactorName);
					result[[i]]$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast);
				}
			} else
			{
				result[[i]]$contrastOnGeno <- testInteractions(model, custom=trmtContrast);
				result[[i]]$contrastAcrossEnv <- NA;
			}
		}
		
	} else if(contrastOpt == "Default")
	{
		if(msaOutput$populationType != "PyramidedLine")
		{
			stop("\tError: The Default value of contrastOpt argument is only implemented on PyramidedLine multi-site analysis outcomes!\n");
		}else
		{
			if(!is.null(genoContrast))
			{
				warning("\tWarning: The genoContrast would be omitted when the contrastOpt is Default!\n");
			}
			
			#-- checking the msaOutput 
			biGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T));
			triGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC","Cc","cc")), sep="", lex.order = T));
			quadraGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC","Cc","cc"), D = c("DD","Dd","dd")), sep="", lex.order = T));
		
			respvar.number <- length(msaOutput$output);
			for(i in 1 : respvar.number)
			{
				is.envFixed <- msaOutput$output[[i]]$envfixed;
				model <- msaOutput$output[[i]]$model;
				trmtLevels <- levels(model@frame[,2]);
				envLevels <- levels(model@frame[,3]);
				
				trmtLevels.length = length(trmtLevels);
				if(trmtLevels.length == 9)
				{
					if(!(all(biGeneLevels %in% trmtLevels)))
					{
						stop(paste("\tError: There are not default coded genotye levels of bigenes design on ", 
										" of ", msaOutput$output[[i]]$respvar, 
										". Cannot proceed with default contrast analysis!\n", sep = ""));
					}
				} else if(trmtLevels.length == 27)
				{
					if(!(all(triGeneLevels %in% trmtLevels)))
					{
						stop(paste("\tError: There are not default coded genotye levels of trigenes design on ", 
										" of ", msaOutput$output[[i]]$respvar, 
										". Cannot proceed with default contrast analysis!\n", sep = ""));
					}
				} else if(trmtLevels.length == 81)
				{
					if(!(all(quadraGeneLevels %in% trmtLevels)))
					{
						stop(paste("\tError: There are not default coded genotye levels of quadragenes design on ", 
										" of ", msaOutput$output[[i]]$respvar, 
										". Cannot proceed with default contrast analysis!\n", sep = ""));
					}
				} else
				{
					stop(paste("\tError: There are not default genotype levels of bi-, tri- or quadra- genes design on ", 
									msaOutput$output[[i]]$respvar, 
									". Cannot proceed with default contrast analysis!\n", sep = ""));
				}
				
				if(is.envFixed)
				{
					if(!is.null(envContrast))
					{	
						if(!is.matrix(envContrast))
						{
							stop("\tError: The specified envContrast should be matrix fromat!\n");
						}
						
						envColContrastNames <- colnames(envContrast);
						
						if(!(is.null(envColContrastNames)))
						{
							if(length(envLevels) != length(envColContrastNames))
							{
								stop("\tError: The levels of specified envContrast is not equal to the levels of env factor on model!\n");
							}
							
							if(!all(envLevels %in% envColContrastNames))
							{
								stop("\tError: The column names of specified envContrast is not same as the names of env factor on model!\n");
							}
						} else
						{
							stop("\tError: The column names of specified envContrast could not be null and should be same as environmental level names!\n");
						}
						# checking the env coefficients to zero!
						if(!all(rowSums(envContrast) == 0))
						{
							stop("\tError: Sum of each row should be equal to zero on envContrast!\n");
						}
					}
				} else
				{
					if(!is.null(envContrast))
					{
						warning("\tWarning: The specified envContrast could not be used when the env factor in the model is random!\n");
					}
				}	
			}
			
			# --- compute the contrast --- #
			for(i in 1 : respvar.number)
			{
				result[[i]] <- list();
				result[[i]]$respvar <- msaOutput$output[[i]]$respvar;		
				is.envFixed <- msaOutput$output[[i]]$envfixed;
				result[[i]]$envfixed <- is.envFixed;
				model <- msaOutput$output[[i]]$model;
				result[[i]]$model <- model;
				genoFactorName <- names(model@frame)[2];
				envFactorName <- names(model@frame)[3];
				trmtLevels <- levels(model@frame[,2]);
				envLevels <- levels(model@frame[,3]);
				trmtLevels.length <- length(trmtLevels);
				genNum <- 2;
				if(trmtLevels.length == 27)
				{
					genNum = 3;
				} else if(trmtLevels.length == 81)
				{
					genNum = 4;
				} else 
				{
					genNum = 2;
				}
				trmtContrast <- getDefaultGenesContrast(genNum);	
							
				tempGenoContrast <- trmtContrast;
				trmtContrast <- as.data.frame(t(tempGenoContrast));
				trmtContrast <- list(genoFactorName = trmtContrast);
				names(trmtContrast) <- genoFactorName;
				if(is.envFixed)
				{
					if(is.null(envContrast))
					{	
						result[[i]]$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
						result[[i]]$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast, across = envFactorName);
					} else
					{
						result[[i]]$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
						tempEnvContrast = envContrast;
						tempEnvContrast = as.data.frame(t(tempEnvContrast));
						trmtContrast$tempEnvContrast <- tempEnvContrast;
						names(trmtContrast) <- c(genoFactorName, envFactorName);
						result[[i]]$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast);
					}
				} else
				{
					result[[i]]$contrastOnGeno <- testInteractions(model, custom=trmtContrast);
					result[[i]]$contrastAcrossEnv <- NA;
				}
			}
		}
	}
	
	detach("package:phia");
	detach("package:lme4")
	return(list(
					output=result,
					contrastType = contrastOpt
			));
}

