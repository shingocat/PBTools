###############################################################################
# TODO: Perform contrast analysis of single site analysis on Introgression line
# 
# ssaOutput - an output of single site analysis on SSSL, Pyramided Line and Introgression Line
# contrastOpt - contrast option; c("Control", "Custom", "Default") and Default value only work for pyramided line right now.
# genoContrast - NULL or a list of contrast matrix where it is corresponding to site number; required if contrastOpt == "Custom"
# controlLevel - NULL or control level
# acrossEnv - if it is ture, that one level of recurrent parent or length of one on genoContrast list is used to each site.
# 
#
# Author: mqin
# Date: Sep 10, 2014
# FileName: contrastAnalysis.pyramidedline.R
###############################################################################

contrastAnalysisOnIL.ssa <- function(
		ssaOutput, 
		contrastOpt = c("RecurrentParent", "Custom", "Default"), 
		recurrentParent = NULL, 
		genoContrast = NULL, 
		acrossEnv = TRUE,
		alpha = 0.05
	) UseMethod("contrastAnalysisOnIL.ssa");
contrastAnalysisOnIL.ssa.default <- function(
		ssaOutput, 
		contrastOpt = c("RecurrentParent", "Custom", "Default"), 
		recurrentParent = NULL, 
		genoContrast = NULL, 
		acrossEnv = TRUE,
		alpha = 0.05
) 
{	
	if(!is.list(ssaOutput))
	{
		stop("Error: The ssaOutput argument should be a list!");
	} else
	{
		if(is.null(ssaOutput$analysisType))
		{
			stop("Error: The ssaOutput argument should be a outcomes of sigle site Analysis!");
		} else if(!(is.null(ssaOutput$analysisType)) && ssaOutput$analysisType != "SingleSiteAnalysis")
		{
			stop("Error: The ssaOutput argument should be a outcomes of sigle site Analysis!");		
		} else
		{
			contrastOpt <- match.arg(contrastOpt);
			
			if(contrastOpt == "RecurrentParent")
			{
				if(is.null(recurrentParent))
				{
					stop("Error: The recurrentParent argument should not be null when contrastOpt is RecurrentParent.")
				} 
				
				# -- checking the recurrent levels whether exist on each site of each response variable
				# -- for storing the TRUE or FALSE of each conditions, but now it changed design not to stop execution!
				checkOutcomes <- c(); 
				# -- for storing the error message of each conditions
				errorMessages <- c();
				respvar.number <- length(ssaOutput$output);
				for (i in (1 : respvar.number))
				{	
					site.number <- length(ssaOutput$output[[i]]$site);
					if(acrossEnv)
					{
						if(length(recurrentParent) != 1)
						{
							stop("Error: The recurrentParent argument should not be more than one level when acrossEnv argument is TRUE!");
						}
					} else
					{
						if(length(recurrentParent) != site.number)
						{
							stop("Error: The length of recurrentParent argument is not same as the tested environment numbers when acrossEnv argument is FALSE!");
						}
					}
					for( j in (1 : site.number))
					{
						if(ssaOutput$output[[i]]$site[[j]]$responseRate < 0.80) # -- this coditions showed warning messages only not stop any more.
						{
							checkOutcomes <- c(checkOutcomes, F);
							errorMessages <- c(errorMessages, paste("\tWarning, there are less than 0.80 observations on ", 
											ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
											" .\n" , sep = "" ));
							warning(paste("\tWarning, there are less than 0.80 observations on ", 
											ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
											". Cannot proceed with contrast analysis!\n" , sep = "" ));
							next;
						}
						if(ssaOutput$output[[i]]$site[[j]]$lmerRun == "ERROR") # -- this conditions showed warning messages only not stop any more
						{
							checkOutcomes <- c(checkOutcomes, F);
							errorMessages <- c(errorMessages, paste("\tError, there are some errors occured on ", 
											ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
											" .\n" , sep = "" ));
							warning(paste("\tWarning, there are some errors of lmer occured on ", 
											ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
											". Cannot proceed with contrast analysis!\n" , sep = "" ));
							next;
						}
						
						model <- ssaOutput$output[[i]]$site[[j]]$model;
						trmtLevels <- levels(model@frame[,2]);
						if(!(recurrentParent %in% trmtLevels)) # -- this conditions is should be stop.
						{
							checkOutcomes <- c(checkOutcomes,F);
							errorMessages <- c(errorMessages, paste("\tError, there are no specified recurrent parent tested on ", 
											ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
											" .\n" , sep = "" ));
							stop(paste("\tError, there are no specified recurrent parent tested on ", 
											ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
											" .\n" , sep = "" ));
							#next;
						}
					}
				}
#				if(length(checkOutcomes) > 0)
#				{
#					warning(errorMessages);
#				}
				
				# get the number of response variable of single site analysis
				respvar.number <- length(ssaOutput$output);
				result <- list();
				library(multcomp);
				for(i in (1 : respvar.number))
				{
					# define each respone variable contrast outcome structure.
					result[[i]] <- list();
					result[[i]]$respvar <- ssaOutput$output[[i]]$respvar;
					result[[i]]$site <- list();
					# Get the number of site of each response variable
					site.number <- length(ssaOutput$output[[i]]$site);
					for( j in (1 : site.number))
					{
						result[[i]]$site[[j]] <- list();
						result[[i]]$site[[j]]$env <- ssaOutput$output[[i]]$site[[j]]$env;
						
						if(ssaOutput$output[[i]]$site[[j]]$responseRate < 0.80) # -- this coditions showed warning messages only not stop any more.
						{
							result[[i]]$site[[j]]$model <- NA;
							result[[i]]$site[[j]]$messages <- paste("\tWarning, there are less than 0.80 observations on ", 
									ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
									". Cannot proceed with contrast analysis!\n" , sep = "" );
							result[[i]]$site[[j]]$error <- TRUE;
							result[[i]]$site[[j]]$outcome <- NA;
#							warning(paste("\tWarning, there are less than 0.80 observations on ", 
#											ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
#											". Cannot proceed with contrast analysis!\n" , sep = "" ));
							next;
						}
						if(ssaOutput$output[[i]]$site[[j]]$lmerRun == "ERROR") # -- this conditions showed warning messages only not stop any more
						{
							result[[i]]$site[[j]]$model <- NA;
							result[[i]]$site[[j]]$messages <- paste("\tWarning, there are some errors of lmer occured on ", 
									ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
									". Cannot proceed with contrast analysis!\n" , sep = "" );
							result[[i]]$site[[j]]$error <- TRUE;
							result[[i]]$site[[j]]$outcome <- NA;
#							warning(paste("\tWarning, there are some errors of lmer occured on ", 
#											ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
#											". Cannot proceed with contrast analysis!\n" , sep = "" ));
							next;
						}
						
						
						model <- ssaOutput$output[[i]]$site[[j]]$model;
						trmtLevels <- levels(model@frame[,2]);
						
						result[[i]]$site[[j]]$messages <- paste("There are nothing wrong on this ", 
								ssaOutput$output[[i]]$site[[j]]$env, " on ", ssaOutput$output[[i]]$respvar, ".\n", sep = "");
						result[[i]]$site[[j]]$error <- FALSE;
						result[[i]]$site[[j]]$model <- model;
						temp <- compareControlContrast(model, control = recurrentParent, alpha);
						signif <- temp[temp[,"Lower"] <= 0 & temp[,"Upper"] >= 0,];
						
						if (nrow(signif) != 0) {
							rownames(signif) <- 1:nrow(signif);
							#cat("SIGNIFICANT LINEAR CONTRAST AT ALPHA = ", alpha, ":\n", sep = "");
							#print(signif)
						} else { 
							#cat("At alpha = ", alpha, ", no significant contrast.\n", sep = "");
						}
						result[[i]]$site[[j]]$outcome <- signif;
					}
				}
				detach("package:multcomp");
				
			} else if(contrastOpt == "Custom")
			{
				if(is.null(genoContrast))
				{
					stop("Error: The genoContrast should not be null when contrastOpt is Custom!");
				}
				if(!is.list(genoContrast))
				{
					stop("Error: The genoContrast should be of data type of list.");
				}
				genoContrastLength <- length(genoContrast);
				respvar.number <- length(ssaOutput$output);
				# checking wheather the genoContrast length is the same  as the number of environment.
				# -- for storing the TRUE or FALSE of each conditions
				checkOutcomes <- c(); 
				# -- for storing the error message of each conditions
				errorMessages <- c();
				for(i in 1 : respvar.number)
				{
					site.number <- length(ssaOutput$output[[i]]$site);
					if(acrossEnv)
					{
						if(genoContrastLength != 1)
						{
							stop("Error: The length of genoContrast argument should be one when the acrossEnv argument is TRUE.\n");
						}
						
						contrastLevels <- colnames(genoContrast[[1]]);
						for(j in 1 : site.number)
						{
							if(ssaOutput$output[[i]]$site[[j]]$responseRate < 0.80) # -- This conditions showed warning messages not stop excution
							{
								checkOutcomes <- c(checkOutcomes, F);
								errorMessages <- c(errorMessages, paste("\tError, there are less than 0.80 observations on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												" .\n" , sep = "" ));
								warning(paste("\tWarning, there are less than 0.80 observations on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												". Cannot proceed with contrast analysis!\n" , sep = "" ));
								next;
							}
							if(ssaOutput$output[[i]]$site[[j]]$lmerRun == "ERROR") # -- This conditions showed warning messages not stop excution
							{
								checkOutcomes <- c(checkOutcomes, F);
								errorMessages <- c(errorMessages, paste("\tError, there are some errors occured on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												" .\n" , sep = "" ));
								warning(paste("\tWarning, there are some errors of lmer model occured on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												". Cannot proceed with contrast analysis!\n" , sep = "" ))
								next;
							}
							
							model <- ssaOutput$output[[i]]$site[[j]]$model;
							trmtLevels <- levels(model@frame[,2]);
							if(length(contrastLevels) != length(trmtLevels)) 
							{
								checkOutcomes <- c(checkOutcomes,F);
								errorMessages <- c(errorMessages, paste("\tError, The levels of specified genoContrast is not same as genotype levels on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												".\n" , sep = "" ));
								stop(paste("\tError, The levels of specified genoContrast is not same as genotype levels on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												".\n" , sep = "" ));
								#next;
							}
							if(!all( contrastLevels %in% trmtLevels)) # -- This conditions stop excution
							{
								checkOutcomes <- c(checkOutcomes,F);
								errorMessages <- c(errorMessages, paste("\tError, there are no specified genoContrast levels tested on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												" .\n" , sep = "" ));
								stop(paste("\tError, there are no specified genoContrast levels tested on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												" .\n" , sep = "" ));
								#next;
							}
						}
						
					} else
					{
						if(genoContrastLength != site.number)
						{
							stop(paste("Error: The length of genoContrast argument should be as the same number of sites on ",
											ssaOutput$output[[i]]$respvar, ". \n", sep =""));
						}
						
						for(j in 1 : site.number)
						{
							if(ssaOutput$output[[i]]$site[[j]]$responseRate < 0.80) # -- This conditions showed warning messages not stop excution
							{
								checkOutcomes <- c(checkOutcomes, F);
								errorMessages <- c(errorMessages, paste("\tWarning, there are less than 0.80 observations on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												". Cannot proceed with contrast analysis!\n" , sep = "" ));
								warning(paste("\tError, there are less than 0.80 observations on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												". Cannot proceed with contrast analysis!\n" , sep = "" ));
								next;
							}
							if(ssaOutput$output[[i]]$site[[j]]$lmerRun == "ERROR") # -- This conditions showed warning messages not stop excution
							{
								checkOutcomes <- c(checkOutcomes, F);
								errorMessages <- c(errorMessages, paste("\tWarning, there are some errors occured on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												" .\n" , sep = "" ));
								warning(paste("\tWarning, there are some errors of lmer model occured on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												". Cannot proceed with contrast analyasis!\n" , sep = "" ));
								next;
							}
							
							model <- ssaOutput$output[[i]]$site[[j]]$model;
							trmtLevels <- levels(model@frame[,2]);
							contrastLevels <- colnames(genoContrast[[j]]);
							if(length(contrastLevels) != length(trmtLevels))
							{
								checkOutcomes <- c(checkOutcomes,F);
								errorMessages <- c(errorMessages, paste("\tError, The levels of specified genoContrast is not same as genotype levels on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												".\n" , sep = "" ));
								stop(paste("\tError, The levels of specified genoContrast is not same as genotype levels on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												".\n" , sep = "" ));
							#next;
							}
							if(!all(contrastLevels %in% trmtLevels))
							{
								checkOutcomes <- c(checkOutcomes,F);
								errorMessages <- c(errorMessages, paste("\tError, there are no specified genoContrast levels tested on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												" .\n" , sep = "" ));
								stop(paste("\tError, there are no specified genoContrast levels tested on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												" .\n" , sep = "" ));
							#next;
							}
						}
					}
					
				}
				
#				if(length(checkOutcomes) > 0)
#				{
#					stop(errorMessages);
#				}
				
				result <- list();
				library(phia);
				library(lme4);
				respvar.number <- length(ssaOutput$output);
				for( i in 1 : respvar.number)
				{
					result[[i]] <- list();
					result[[i]]$respvar <- ssaOutput$output[[i]]$respvar;
					result[[i]]$site <- list();
					site.number <- length(ssaOutput$output[[i]]$site);
					for(j in 1 : site.number)
					{	
						result[[i]]$site[[j]] <- list();
						result[[i]]$site[[j]]$env <- ssaOutput$output[[i]]$site[[j]]$env;
						if(acrossEnv)
						{	
							if(ssaOutput$output[[i]]$site[[j]]$responseRate < 0.80) # -- This conditions showed warning messages not stop excution
							{
								result[[i]]$site[[j]]$model <- NA;
								result[[i]]$site[[j]]$messages <- paste("\tWarning, there are less than 0.80 observations on ", 
										ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
										". Cannot proceed with contrast analysis!\n" , sep = "" );
								result[[i]]$site[[j]]$error <- TRUE;
								result[[i]]$site[[j]]$outcome <- NA;
#								warning(paste("\tWarning, there are less than 0.80 observations on ", 
#												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
#												". Cannot proceed with contrast analysis!\n" , sep = "" ));
								next;
							}
							if(ssaOutput$output[[i]]$site[[j]]$lmerRun == "ERROR") # -- This conditions showed warning messages not stop excution
							{
								result[[i]]$site[[j]]$model <- NA;
								result[[i]]$site[[j]]$messages <- paste("\tWarning, there are some errors of lmer model occured on ", 
										ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
										". Cannot proceed with contrast analysis!\n" , sep = "" );
								result[[i]]$site[[j]]$error <- TRUE;
								result[[i]]$site[[j]]$outcome <- NA;
#								warning(paste("\tWarning, there are some errors of lmer model occured on ", 
#												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
#												". Cannot proceed with contrast analysis!\n" , sep = "" ));
								next;
							}
							result[[i]]$site[[j]]$messages <- paste("There are nothing wrong on this ", 
									result$output[[i]]$site[[j]]$env, " on ", result$output[[i]]$respvar, ".\n", sep = "");
							result[[i]]$site[[j]]$error <- FALSE;
							model <- ssaOutput$output[[i]]$site[[j]]$model;
							result[[i]]$site[[j]]$model <- model;
							genoContrastTemp <- t(genoContrast[[1]]);
							genoContrastTemp <- as.data.frame(genoContrastTemp);
							genoFactorName <- names(model@frame)[2];
							genoContrastTemp <- list(genoContrastTemp);
							names(genoContrastTemp) <- genoFactorName;
							temp <- testInteractions(model, custom=genoContrastTemp);
							result[[i]]$site[[j]]$outcome <- temp;
						} else
						{	
							if(ssaOutput$output[[i]]$site[[j]]$responseRate < 0.80) # -- This conditions showed warning messages not stop excution
							{
								result[[i]]$site[[j]]$model <- NA;
								result[[i]]$site[[j]]$messages <- paste("\tWarning, there are less than 0.80 observations on ", 
										ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
										". Cannot proceed with contrast analysis!\n" , sep = "" );
								result[[i]]$site[[j]]$error <- TRUE;
								result[[i]]$site[[j]]$outcome <- NA;
#								warning(paste("\tWarning, there are less than 0.80 observations on ", 
#												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
#												". Cannot proceed with contrast analysis!\n" , sep = "" ));
								next;
							}
							if(ssaOutput$output[[i]]$site[[j]]$lmerRun == "ERROR") # -- This conditions showed warning messages not stop excution
							{
								result[[i]]$site[[j]]$model <- NA;
								result[[i]]$site[[j]]$messages <- paste("\tWarning, there are some errors of lmer model occured on ", 
										ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
										". Cannot proceed with contrast analysis!\n" , sep = "" );
								result[[i]]$site[[j]]$error <- TRUE;
								result[[i]]$site[[j]]$outcome <- NA;
#								warning(paste("\tWarning, there are some errors of lmer model occured on ", 
#												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
#												". Cannot proceed with contrast analysis!\n" , sep = "" ));
								next;
							}
							result[[i]]$site[[j]]$messages <- paste("There are nothing wrong on this ", 
									ssaOutput$output[[i]]$site[[j]]$env, " on ", ssaOutput$output[[i]]$respvar, ".\n", sep = "");
							result[[i]]$site[[j]]$error <- FALSE;
							model <- ssaOutput$output[[i]]$site[[j]]$model;
							result[[i]]$site[[j]]$model <- model;
							genoContrastTemp <- t(genoContrast[[j]]);
							genoContrastTemp <- as.data.frame(genoContrastTemp);
							genoFactorName <- names(model@frame)[2];
							genoContrastTemp <- list(genoContrastTemp);
							names(genoContrastTemp) <- genoFactorName;
							temp <- testInteractions(model, custom=genoContrastTemp);
							result[[i]]$site[[j]]$outcome <- temp;
						}
					}
				}
				detach("package:phia");
				detach("package:lme4");
			} else if (contrastOpt == "Default")
			{
				if(ssaOutput$populationType != "PyramidedLine")
				{
					stop("\tError: The Default value of contrastOpt argument is only implemented on PyramidedLine single site analysis outcomes!\n ");
				} else
				{
					#-- checking the ssaOutput 
					checkingOutcome <- c();
					checkingMessage <- c();
					biGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T));
					triGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC","Cc","cc")), sep="", lex.order = T));
					quadraGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC","Cc","cc"), D = c("DD","Dd","dd")), sep="", lex.order = T));
					
					respvar.number <- length(ssaOutput$output);
					for(i in 1 : respvar.number)
					{
						site.number <- length(ssaOutput$output[[i]]$site);
						for( j in 1:site.number)
						{
							if(ssaOutput$output[[i]]$site[[j]]$responseRate < 0.80) # -- The condition showed warning messages not stop execution.
							{
#								checkingOutcome <- c(checkingOutcome, FALSE);
#								checkingMessage <- c(checkingMessage, paste("\tError: There are too many missing observations on ", 
#												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, ". ",
#												"Cannot proceed with contrast analysis!\n", sep = ""));
								warning(paste("\tError: There are too many missing observations on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, ". ",
												"Cannot proceed with contrast analysis!\n", sep = ""))
								next;
							}
							if(ssaOutput$output[[i]]$site[[j]]$lmerRun == "ERROR") # -- The condition showed warning messages not stop execution.
							{
#								checkingOutcome <- c(checkingOutcome, FALSE);
#								checkingMessage <- c(checkingMessage, paste("\tError: There are some error of lmer model on ", 
#												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, ". ",
#												"Cannot proceed with contrast analysis!\n", sep = ""));
								warning(paste("\tError: There are some error of lmer model on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, ". ",
												"Cannot proceed with contrast analysis!\n", sep = ""));
								next;
							}
							
							model <- ssaOutput$output[[i]]$site[[j]]$model;
							trmtLevels <- levels(model@frame[,2]);
							
							trmtLevels.length = length(trmtLevels);
							if(trmtLevels.length == 9)
							{
								if(!all(biGeneLevels %in% trmtLevels))
								{
									checkingOutcome <- c(checkingOutcome, FALSE);
									checkingMessage <- c(checkingMessage, paste("\tError: There are not default coded genotye levels of bigenes design on ", 
													ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
													". Cannot proceed with default contrast analysis!\n", sep = ""));
									stop(paste("\tError: There are not default coded genotye levels of bigenes design on ", 
													ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
													". Cannot proceed with default contrast analysis!\n", sep = ""))
									#next;
								}
							} else if(trmtLevels.length == 27)
							{
								if(!all(triGeneLevels %in% trmtLevels))
								{
									checkingOutcome <- c(checkingOutcome, FALSE);
									checkingMessage <- c(checkingMessage, paste("\tError: There are not default coded genotye levels of trigenes design on ", 
													ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
													". Cannot proceed with default contrast analysis! \n", sep = ""));
									stop(paste("\tError: There are not default coded genotye levels of trigenes design on ", 
													ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
													". Cannot proceed with default contrast analysis!\n", sep = ""));
									#next;
								}
							} else if(trmtLevels.length == 81)
							{
								if(!all(quadraGeneLevels %in% trmtLevels))
								{
									checkingOutcome <- c(checkingOutcome, FALSE);
									checkingMessage <- c(checkingMessage, paste("\tError: There are not default coded genotye levels of quadragenes design on ", 
													ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
													". Cannot proceed with default contrast analysis!\n", sep = ""));
									stop(paste("\tError: There are not default coded genotye levels of quadragenes design on ", 
													ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
													". Cannot proceed with default contrast analysis!\n", sep = ""));
									#next;
								}
							} else
							{
								checkingOutcome <- c(checkingOutcome, FALSE);
								checkingMessage <- c(checkingMessage, paste("\tError: There are not default genotype levels of bi-, tri- or quadra- genes design on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												". Cannot proceed with default contrast analysis!\n", sep = ""));
								stop(paste("\tError: There are not default genotype levels of bi-, tri- or quadra- genes design on ", 
												ssaOutput$output[[i]]$site[[j]]$env, " of ", ssaOutput$output[[i]]$respvar, 
												". Cannot proceed with default contrast analysis!\n", sep = ""));
								#next;
							}
						}
					}
#					if(length(checkingOutcome) > 0)
#					{
#						stop(checkingMessage);
#					}
					
					result <- list();
					library(phia);
					library(lme4);
					for( i in 1:respvar.number)
					{
						result[[i]] <- list();
						result[[i]]$respvar <- ssaOutput$output[[i]]$respvar;
						result[[i]]$site <- list();
						site.number <- length(ssaOutput$output[[i]]$site);
						for( j in 1 : site.number)
						{
							result[[i]]$site[[j]] <- list();
							result[[i]]$site[[j]]$env <- ssaOutput$output[[i]]$site[[j]]$env;
							result[[i]]$site[[j]]$messages <- paste("There are nothing wrong on this ", 
									ssaOutput$output[[i]]$site[[j]]$env, " on ", ssaOutput$output[[i]]$respvar, ".\n", sep = "");
							result[[i]]$site[[j]]$error <- FALSE;
							model <- ssaOutput$output[[i]]$site[[j]]$model;
							result[[i]]$site[[j]]$model <- model;
							trmtLevels <- levels(model@frame[,2]);
							trmtLevelsLength <- length(trmtLevels);
							genoContrast <- NA;
							if(trmtLevelsLength == 9)
							{
								genoContrast <- getDefaultGenesContrast(2);
							} else if (trmtLevelsLength == 27)
							{
								genoContrast <- getDefaultGenesContrast(3);
							} else if (trmtLevelsLength == 81)
							{
								genoContrast <- getDefaultGenesContrast(4);
							}
							genoContrastTemp <- t(genoContrast);
							genoContrastTemp <- as.data.frame(genoContrastTemp);
							genoFactorName <- names(model@frame)[2];
							genoContrastTemp <- list(genoContrastTemp);
							names(genoContrastTemp) <- genoFactorName;
							temp <- testInteractions(model, custom=genoContrastTemp);
							result[[i]]$site[[j]]$outcome <- temp;
							
						}
					}
					detach("package:phia");
					detach("package:lme4");
				} # end statement of if(ssaOutput$populationType != "PyramidedLine")
			} else 
			{
				stop("Error: The value of contrastOpt argument should be one of \"Control, Custom, BiGenes, TriGenes, QuadraGenes\" value.");
			}
		}
	}
	
	return(list(output = result,
					contrastType = contrastOpt
			));
}