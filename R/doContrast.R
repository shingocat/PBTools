###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Nov 25, 2014
# FileName: doContrast.R
###############################################################################


doContrast <- function(
		data,
		contrastOpt = c("RecurrentParent", "Custom", "Default"),
		recurrentParent = NULL,
		genoContrast = NULL,
		acrossEnv = TRUE,
		envContrast = NULL,
		alpha = 0.05
)
{	
	UseMethod("doContrast");
}

doContrast.SingleEnvAnalysis <- function(
		data,
		contrastOpt = c("RecurrentParent", "Custom", "Default"),
		recurrentParent = NULL,
		genoContrast = NULL,
		acrossEnv = TRUE,
		alpha = 0.05
)
{
	if(!inherits(data, "SingleEnvAnalysis"))
		stop("\tError: The data should be of class SingleSiteAnalysis!\n");
	contrastOpt <- match.arg(contrastOpt);
	if(contrastOpt == "RecurrentParent")
	{
		if(missing(recurrentParent))
		{
			stop("\tError: The recurrentParent argument should not be null when contrastOpt is RecurrentParent.")
		}
		# -- checking the recurrent levels whether exist on each site of each response variable
		respvar.number <- length(data$traits);
		for (i in (1 : respvar.number))
		{	
			site.number <- length(data$traits[[i]]$analysis$sea$envs);
			trait.name <- data$traits[[i]]$name;
			
			if(acrossEnv)
			{
				if(length(recurrentParent) != 1)
				{
					stop("\tError: The recurrentParent argument should not be more than one level when acrossEnv argument is TRUE!\n");
				}
			} else
			{
				if(length(recurrentParent) != site.number)
				{
					stop("\tError: The length of recurrentParent argument is not same as the tested environment numbers when acrossEnv argument is FALSE!\n");
				}
			}
			for( j in (1 : site.number))
			{
				env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
				# -- this coditions showed warning messages only not stop any more.---#
				if(data$traits[[i]]$envs[[j]]$restricted$isTRUE) 
				{
					warning(paste("\tWarning, the missing rate is larger than specified missing rate,", data$traits[[i]]$envs$restricted$missingRate,
									" of ", trait.name, " on ", env.name, 
									".\n\tCannot proceed with contrast analysis!\n" , sep = "" ));
					next;
				}
				# -- this conditions showed warning messages only not stop any more ---#
				if(!inherits(data$traits[[i]]$analysis$sea$envs[[j]]$model, "lmerMod")) 
				{
					warning(paste("\tWarning, there are some errors of lmer occured on ", 
									env.name, " of ", trait.name, 
									". Cannot proceed with contrast analysis!\n" , sep = "" ));
					next;
				}
				
				model <- data$traits[[i]]$analysis$sea$envs[[j]]$model;
				trmtLevels <- levels(model@frame[,2]);
				# -- checking whether all envs have the same recurrent parent. this conditions is should be stop.---#
				if(!(recurrentParent %in% trmtLevels))
				{
					stop(paste("\tError, there are no specified recurrent parent tested on ", 
									env.name, " of ", trait.name, 
									" .\n" , sep = "" ));
					#next;
				}
			}
		}
		
		#--- get the number of response variable of single site analysis---#
		library(multcomp)
		for(i in 1:respvar.number)
		{
			trait.name <- data$traits[[i]]$name;
			#--- define each respone variable contrast outcome structure under analysis$sea$envs[[j]]$contrast---#
			site.number <- length(data$traits[[i]]$analysis$sea$envs)
			for( j in 1:site.number){
				env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
				data$traits[[i]]$analysis$sea$envs[[j]]$contrast <- list();
				#---This condition showed warning meseages instead of stop the program---#
				if(data$traits[[i]]$envs[[j]]$restricted$isTRUE)
				{
					#data$traits[[i]]$analysis$sea$envs[[j]]$contrast$model <- NA;
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messange <-  paste("\tWarning, the missing rate is larger than specified missing rate on ", 
							env.name, " of ", trait.name, 
							". Cannot proceed with contrast analysis!\n" , sep = "" );
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$error <- TRUE;
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome <- NA;
					next;
				}
				#---This condition showed warning messages instead of stop the program---#
				if(!inherits(data$traits[[i]]$analysis$sea$envs[[j]]$model, "lmerMod"))
				{
					#data$traits[[i]]$analysis$sea$envs[[j]]$contrast$model <- NA;
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messages <- paste("\tWarning, there are some errors of lmer occured on ", 
							env.name, " of ", trait.name, 
							". Cannot proceed with contrast analysis!\n" , sep = "" );
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$error <- TRUE;
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome <- NA;
					next;
				}
				
				model <- data$traits[[i]]$analysis$sea$envs[[j]]$model;
				trmtLevels <- levels(model@frame[,2]);
				data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messages <- paste("There are nothing wrong on this ", 
						env.name, " on ", trait.name, ".\n", sep = "");
				data$traits[[i]]$analysis$sea$envs[[j]]$contrast$error <- FALSE;
				temp <- compareControlContrast(model, control = recurrentParent, alpha);
# 				signif <- temp[temp[,"Lower"] <= 0 & temp[ ,"Upper"] >= 0];
# 				if(nrow(signif != 0))
# 				{
# 					rownames(signif) <- 1:nrow(signif);
# 				}
        data$traits[[i]]$analysis$sea$envs[[j]]$contrast$type <- contrastOpt;
        data$traits[[i]]$analysis$sea$envs[[j]]$contrast$recurrent <- recurrentParent;
        data$traits[[i]]$analysis$sea$envs[[j]]$contrast$alpha <- alpha;
				data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome <- temp;
			}
		}
		detach("package:multcomp");
	} else if (contrastOpt == "Custom")
	{
		if(missing(genoContrast))
		{
			stop("\tError: The genoContrast argument should not be null when contrastOpt is Custom!\n");
		}
		if(!is.list(genoContrast))
		{
			stop("\tError: The genoContrast argument should be of data type of list.\n");
		}
		genoContrastLength <- length(genoContrast);
		respvar.number <- length(data$traits);
		#--- checking whether the genoContrast length is the same as fthe number of environment.---#
		for(i in 1:respvar.number)
		{
			site.number <- length(data$traits[[i]]$analysis$sea$envs);
			trait.name <- data$traits[[i]]$name;
			if(acrossEnv)
			{
				if(genoContrastLength != 1)
				{
					stop("\tError: The length of genoContrast argument should be one when the acrossEnv argument is TRUE.\n");
				}
				contrastLevels <- colnames(genoContrast[[1]]);
				for(j in 1:site.number)
				{
					env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
					#--- This conditions showed warning messages not stop excution---#
					if(data$traits[[i]]$analysis$sea$envs[[j]]$restricted$isTRUE)
					{
						checkOutcomes <- c(checkOutcomes, F);
						errorMessages <- c(errorMessages, paste("\tError, there are missing rate larger than specified one on ", 
										env.name, " of ", trait.name, 
										" .\n" , sep = "" ));
						warning(paste("\tWarning, there are missing rate larger than specified one on ", 
										env.name, " of ", trait.name, 
										". Cannot proceed with contrast analysis!\n" , sep = "" ));
						next;
					}
					# -- This conditions showed warning messages not stop excution
					if(data$traits[[i]]$analysis$sea$envs[[j]]$lmerRun == "ERROR")
					{
						checkOutcomes <- c(checkOutcomes, F);
						errorMessages <- c(errorMessages, paste("\tError, there are some errors occured on ", 
										env.name, " of ", trait.name, 
										" .\n" , sep = "" ));
						warning(paste("\tWarning, there are some errors of lmer model occured on ", 
										env.name, " of ", trait.name, 
										". Cannot proceed with contrast analysis!\n" , sep = "" ))
						next;
					}
					
					model <- data$traits[[i]]$analysis$sea$envs[[j]]$model;
					trmtLevels <- levels(model@frame[ , 2]);
					#--- This conditions stop excution---#
					if(length(contrastLevels) != length(trmtLevels))
					{
						checkOutcomes <- c(checkOutcomes,F);
						errorMessages <- c(errorMessages, paste("\tError, The levels of specified genoContrast is not same as genotype levels on ", 
										env.name, " of ", trait.name, 
										".\n" , sep = "" ));
						stop(paste("\tError, The levels of specified genoContrast is not same as genotype levels on ", 
										env.name, " of ", trait.name, 
										".\n" , sep = "" ));
					}
					#--- This conditions stop excution---#
					if(!all( contrastLevels %in% trmtLevels)) 
					{
						checkOutcomes <- c(checkOutcomes,F);
						errorMessages <- c(errorMessages, paste("\tError, there are no specified genoContrast levels tested on ", 
										env.name, " of ", trait.name, 
										" .\n" , sep = "" ));
						stop(paste("\tError, there are no specified genoContrast levels tested on ", 
										env.name, " of ", trait.name, 
										" .\n" , sep = "" ));
						#next;
					}
					
				}
			} else 
			{
				if(genoContrastLength != site.number)
				{
					stop(paste("\tError: The length of genoContrast argument should be as the same number of envs on ",
									trait.name, ". \n", sep =""));
				}
				
				for( j in 1:site.number)
				{
					env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
					#--- This conditions showed warning messages not stop excution---#
					if(data$traits[[i]]$envs[[j]]$restricted$isTRUE)
					{
						checkOutcomes <- c(checkOutcomes, F);
						errorMessages <- c(errorMessages, paste("\tWarning, there are missing value larger than specified one on ", 
										env.name, " of ", trait.name, 
										". Cannot proceed with contrast analysis!\n" , sep = "" ));
						warning(paste("\tWarning, there are missing value larger than specified one on ", 
										env.name, " of ", trait.name, 
										". Cannot proceed with contrast analysis!\n" , sep = "" ));
						next;
					}
					#--- This conditions showed warning messages not stop excution---#
					if(data$triats[[i]]$analysis$sea$envs[[j]]$lmerRun == "ERROR") 
					{
						checkOutcomes <- c(checkOutcomes, F);
						errorMessages <- c(errorMessages, paste("\tWarning, there are some errors occured on ", 
										env.name, " of ", trait.name, 
										" .\n" , sep = "" ));
						warning(paste("\tWarning, there are some errors of lmer model occured on ", 
										env.name, " of ", trait.name, 
										". Cannot proceed with contrast analyasis!\n" , sep = "" ));
						next;
					}
					
					model <- data$traits[[i]]$analysis$sea$envs[[j]]$model;
					trmtLevels <- levels(mdoel@frame[,2]);
					contrastLevels <- colnames(genoContrast[[j]]);
					if(length(contrastLevels) != length(trmtLevels))
					{
						checkOutcomes <- c(checkOutcomes,F);
						errorMessages <- c(errorMessages, paste("\tError, The levels of specified genoContrast is not same as genotype levels on ", 
										env.name, " of ", trait.name, 
										".\n" , sep = "" ));
						stop(paste("\tError, The levels of specified genoContrast is not same as genotype levels on ", 
										env.name, " of ", trait.name, 
										".\n" , sep = "" ));
					}
					
					if(!all(contrastLevels %in% trmtLevels))
					{
						checkOutcomes <- c(checkOutcomes,F);
						errorMessages <- c(errorMessages, paste("\tError, there are no specified genoContrast levels tested on ", 
										env.name, " of ", trait.name, 
										" .\n" , sep = "" ));
						stop(paste("\tError, there are no specified genoContrast levels tested on ", 
										env.name, " of ", trait.name, 
										" .\n" , sep = "" ));
					}
				}
			}
		} #---end of statement for( i in 1:respvar.number)---#
		
		#---compute contrast comparing---#
		library(phia);
		library(lme4);
		for(i in 1:respvar.number)
		{
			trait.name <- data$traits[[i]]$name;
			site.number <- length(data$traits[[i]]$analysis$sea$envs);
			for(j in 1:site.number)
			{
				data$traits[[i]]$analysis$sea$envs[[j]]$contrast <- list();
				env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
				if(acrossEnv)
				{
					# -- This conditions showed warning messages not stop excution---#
					if(data$traits[[i]]$envs[[j]]$restricted$isTRUE)
					{
						data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messages <- paste("\tWarning, there are missing rate larger than specified one on ", 
								env.name, " of ", trait.name, 
								". Cannot proceed with contrast analysis!\n" , sep = "" );
						data$traits[[i]]$analysis$sea$envs[[j]]$contrast$error <- TRUE;
						data$traits[[i]]$analysis$sea$envs[[j]]$outcome <- NA;
						next;
					}
					
					# -- This conditions showed warning messages not stop excution---#
					if(data$traits[[i]]$analysis$sea$envs[[j]]$lmerRun == "ERROR")
					{
						data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messages <- paste("\tWarning, there are some errors of lmer model occured on ", 
								env.name, " of ", trait.name, 
								". Cannot proceed with contrast analysis!\n" , sep = "" );
						data$traits[[i]]$analysis$sea$site[[j]]$contrast$error <- TRUE;
						data$traits[[i]]$analysis$sea$site[[j]]$contrast$outcome <- NA;
						next;
					}
					
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messages <- paste("There are nothing wrong on this ", 
							env.name, " on ", trait.name, ".\n", sep = "");
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$error <- FALSE;
					model <- data$traits[[i]]$analysis$sea$envs[[j]]$model;
					genoContrastTemp <- t(genoContrast[[1]]);
					genoContrastTemp <- as.data.frame(genoContrastTemp);
					genoFactorName <- names(model@frame)[2];
					genoContrastTemp <- list(genoContrastTemp);
					names(genoContrastTemp) <- genoFactorName;
					temp <- testInteractions(model, custom=genoContrastTemp);
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome <- temp;
					
				} else
				{
					# -- This conditions showed warning messages not stop excution---#
					if(data$traits[[i]]$envs[[j]]$restricted$isTRUE)
					{
						data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messages <- paste("\tWarning, there are missing rate larger than specified one on ", 
								env.name, " of ", trait.name, 
								". Cannot proceed with contrast analysis!\n" , sep = "" );
						data$traits[[i]]$analysis$sea$envs[[j]]$contrast$error <- TRUE;
						data$traits[[i]]$analysis$sea$envs[[j]]$outcome <- NA;
						next;
					}
					
					# -- This conditions showed warning messages not stop excution---#
					if(data$traits[[i]]$analysis$sea$envs[[j]]$lmerRun == "ERROR")
					{
						data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messages <- paste("\tWarning, there are some errors of lmer model occured on ", 
								env.name, " of ", trait.name, 
								". Cannot proceed with contrast analysis!\n" , sep = "" );
						data$traits[[i]]$analysis$sea$site[[j]]$contrast$error <- TRUE;
						data$traits[[i]]$analysis$sea$site[[j]]$contrast$outcome <- NA;
						next;
					}
					
					data$traits[[i]]$analysis$sea$site[[j]]$contrast$messages <- paste("There are nothing wrong on this ", 
							env.name, " on ", trait.name, ".\n", sep = "");
					data$traits[[i]]$analysis$sea$site[[j]]$contrast$error <- FALSE;
					model <- data$traits[[i]]$analysis$sea$envs[[j]]$model;
					genoContrastTemp <- t(genoContrast[[j]]);
					genoContrastTemp <- as.data.frame(genoContrastTemp);
					genoFactorName <- names(model@frame)[2];
					genoContrastTemp <- list(genoContrastTemp);
					names(genoContrastTemp) <- genoFactorName;
					temp <- testInteractions(model, custom=genoContrastTemp);
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$type <- contrastOpt;
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome <- temp;
				}#--- end statment of if(acrossEnv)---#
				
			}#--- end statment of for(j in 1:site.number---#
		} #--- end statment of for(i in 1:respvar.number---#
		detach("package:phia");
		detach("package:lme4");
	} else if (contrastOpt == "Default")
	{
		if(data$pop.type != "PL")
		{
			stop("\tError: The Default value of contrastOpt argument is only implemented on PyramidedLine single site analysis outcomes!\n ");
		} else
		{
			#---checking the data---#
			biGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T));
			triGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC","Cc","cc")), sep="", lex.order = T));
			quadraGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC","Cc","cc"), D = c("DD","Dd","dd")), sep="", lex.order = T));
			
			respvar.number <- length(data$traits);
			for(i in 1:respvar.number)
			{	
				trait.name <- data$traits[[i]]$name;
				site.number <- length(data$traits[[i]]$analysis$sea$envs);
				for(j in 1:site.number)
				{
					env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
					#--- This conditions show warning messages not stop execution---#
					if(data$traits[[i]]$envs[[j]]$restricted$isTRUE)
					{
						warning(paste("\tError: There are too many missing observations on ", 
										env.name, " of ", triat.name, ". ",
										"Cannot proceed with contrast analysis!\n", sep = ""))
						next;
					}
					# -- The condition showed warning messages not stop execution.---#
					if(data$traits[[i]]$analysis$sea$envs[[j]]$lmerRun == "ERROR") 
					{
						warning(paste("\tError: There are some error of lmer model on ", 
										env.name, " of ", trait.name, ". ",
										"Cannot proceed with contrast analysis!\n", sep = ""));
						next;
					}
					
					model <- data$traits[[i]]$analysis$sea$envs[[j]]$model;
					trmtLevels <- levels(model$frame[,2]);
					trmtLevels.length <- length(trmtLevels);
					if(data$gene.numer == 2)
					{
						if(!all(biGeneLevels %in% trmtLevels))
						{
							stop(paste("\tError: There are not default coded genotye levels of bigenes design on ", 
											env.name, " of ", trait.name, 
											". Cannot proceed with default contrast analysis!\n", sep = ""))
							#next;
						}
					} else if(data$gene.number == 3)
					{
						if(!all(triGeneLevels %in% trmtLevels))
						{
							stop(paste("\tError: There are not default coded genotye levels of triGene design on ", 
											env.name, " of ", trait.name, 
											". Cannot proceed with default contrast analysis!\n", sep = ""))
						}
					} else if(data$gene.number == 4)
					{
						if(!all(quadraGeneLevels %in% trmtLevels))
						{
							stop(paste("\tError: There are not default coded genotye levels of quadraGeneLevels design on ", 
											env.name, " of ", trait.name, 
											". Cannot proceed with default contrast analysis!\n", sep = ""))
						}
					} else
					{
							stop(paste("\tError: There are not default coded genotye levels of bigenes, triGene or quadraGene design on ", 
											env.name, " of ", trait.name, 
											". Cannot proceed with default contrast analysis!\n", sep = ""))
					}
					
				}#---end of statment of for(j in 1:site.number)---#
			}#---end of statment of for(i in 1:respvar.number)
			library(phia);
			library(lme4);
			for(i in 1:respvar.number)
			{
				trait.name <- data$traits[[i]]$name;
				site.number <- length(data$traits[[i]]$analysis$envs);
				data$traits[[i]]$analysis$sea$envs[[j]]$contrast <- list();
				for(j in 1:site.number)
				{
					env.name <- data$traits[[i]]$analysis$sea$envs[[j]]$name;
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$messages <- paste("There are nothing wrong on this ", 
							env.name, " on ", trait.name, ".\n", sep = "");
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$error <- FALSE;
					model <- data$traits[[i]]$analysis$sea$envs[[j]]$model;
					trmtLevels <- levels(model@frame[,2]);
					trmtLevelsLength <- length(trmtLevels);
					genoContrast <- NA;
          gene.number <- 2;
					if(trmtLevelsLength == 9)
					{
						genoContrast <- getDefaultGenesContrast(2);
					} else if (trmtLevelsLength == 27)
					{
            gene.number <- 3;
						genoContrast <- getDefaultGenesContrast(3);
					} else if (trmtLevelsLength == 81)
					{
            gene.number <- 4;
						genoContrast <- getDefaultGenesContrast(4);
					}
					genoContrastTemp <- t(genoContrast);
					genoContrastTemp <- as.data.frame(genoContrastTemp);
					genoFactorName <- names(model@frame)[2];
					genoContrastTemp <- list(genoContrastTemp);
					names(genoContrastTemp) <- genoFactorName;
					temp <- testInteractions(model, custom=genoContrastTemp);
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$type <- contrastOpt;
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$gene.number <- gene.number;
					data$traits[[i]]$analysis$sea$envs[[j]]$contrast$outcome <- temp;
				}
			}
			detach("package:phia");
			detach("package:lme4");
		} #--- end statment of if(data$pop.type != "PL")---#
	}#---end of statment of if(contrastOpt == "recurrentParent")---#
	
	return(data);
}


doContrast.MultiEnvAnalysis <- function(
		data,
		contrastOpt = c("RecurrentParent", "Custom", "Default"), 
		recurrentParent = NULL, 
		genoContrast = NULL, 
		envContrast = NULL
)
{
	if(missing(data))
		stop("\tError: The argument of data could not be null!\n");
	if(!inherits(data, MultiEnvAnalysis))
		stop("\tError: The argument of data should be of class MultiSiteAnalysis!\n");
	contrastOpt <- match.arg(contrastOpt);
	
	library(phia);
	library(lme4);
	
	if(contrastOpt == "RecurrentParent")
	{
		if(missing(recurrentParent))
			stop("\tError: The argument of recurrentParent should be null when contrastOpt is selected RecurrentParent!\n");
		if(!is.vector(recurrentParent, mode = "character"))
			stop("\tError: The argument of recurrentParent should be vector when contrastOpt is selected RecurrentParent!\n");
		if(length(recurrentParent) != 1)
			stop("\tError: The argument of recurrentParent should be of length one when contrastOpt is selected RecurrentParent!\n");
		respvar.number <- length(data$traits);
		#--- checking conditions---#
		for(i in 1:respvar.number)
		{
			is.envFixed <- data$traits[[i]]$msa$envFixed;
			model <- data$traits[[i]]$msa$model;
			trmtLevels <- levels(model@frame[ , 2]);
			envLevels <- levels(model@frame[ , 3]);
			#--- checking whether the specified recurrent parent exist!---#
			if(!(recurrentParent %in% trmtLevels))
				stop("\tError: The specified recurrent parent does not include in the model!\n");
			if(is.envFixed)
			{
				if(!missing(envContrast))
				{
					if(!is.matrix(envContrast))
					{
						stop("\tError: The specified envContrast should be matrix format!\n");
					}
					
					envColContrastNames <- colnames(envContrast);
					if(length(envLevels) != length(envColContrastNames))
						stop("\tError: The levels of specified envContrast is not equal to the levels of env factor on model!\n");
					if(!(is.null(envColContrastNames)))
					{
						if(!all(envLevels %in% encColContrastNames))
							stop("\tError: The column names of specified envContrast is not same as the names of env factor on model!\n");
					} else 
					{
						stop("\tError: The column names of specified envContrast could not be null and should be same as environment level names!\n");
					}
					
					#---checking the env coefficients to zero@---#
					if(!all(rowSums(envContrast) == 0))
						stop("\tError: sum of each row should be equal to zero on envContrast!\n");
				}
			} else 
			{
				if(!missing(envContrast))
					warning("\tWarning: The specified envContrast could not be used when the env factor in the model is random!\n ");
			}
		} #--- end of statement for(i in 1:respvar.number) #checking conditions----#
		
		#---compute contrast comparing on recurrent parent---#
		for(i in 1:respvar.number)
		{
			is.envFixed <- data$traits[[i]]$analysis$msa$envFixed;
			model <- data$traits[[i]]$analysis$msa$model;
			data$traits[[i]]$analysis$msa$contrast <- list();	
			
			genoFactorName <- names(model@frame)[2];
			envFactorName <- names(model$frame)[3];
			trmtLevels <- levels(model@frame[ ,2]);
			trmtLevels.without.rp <- trmtLevels[(trmtLevels != recurrentParent)];
			trmtContrast <- diag(1, nrow = length(trmtLevels.without.rp), ncol=length(trmtLevels.without.rp));
			trmtContrast <- cbind(trmtContrast, -1);
			colnames(trmtContrast) <- c(trmtLevels.without.rp, recurrentParent);
			rownames(trmtContrast) <- interaction(trmtLevels.without.rp, as.factor(paste(" - ", recurrentParent, sep = "")), sep = "");
			trmtContrast <- as.data.frame(t(trmtContrast));
			trmtContrast <- list(genoFactorName = trmtContrast);
			names(trmtContrast) <- genoFactorName;
			if(is.envFixed)
			{
				if(is.null(envContrast))
				{
					data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom = trmtContrast);
					data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- testInteraction(model, custom = trmtContrast, across = envFactorName);
				} else 
				{
					data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom = trmtContrast);
					tempEnvContrast <- envContrast;
					tempEnvContrast <- as.data.frame(t(tempEnvContrast));
					trmtContrast$tempEnvContrast <- tempEnvContrast;
					names(trmtContrast) <- c(genoFactorName, envFactorName);
					data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- testInteraction(model, custom = trmtContrast);
				}
			} else
			{
				data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom = trmtContrast);
				data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- NA;
			}
			
		}
		
	} else if(contrastOpt == "Custom")
	{
		if(missing(genoContrast))
			stop("\tError: The argument of genoContrast could not be null when contrastOpt is Custom!\n");
		if(!is.matrix(genoContrast))
			stop("\tError: The genoContrast should be matrix format!\n");
		respvar.number <- length(data$traits);
		#--- checking conditions---#
		for(i in 1:respvar.number)
		{
			is.envFixed <- data$traits[[i]]$analysis$msa$envFixed;
			model <- data$traits[[i]]$analysis$msa$model;
			trmtLevels <- levels(model$frame[,2]);
			envLevels <- levels(model$frame[,3]);
			genoColContrastName <- colnames(genoContrast);
			
			#---checking whether the specified genoContrast levels exist!---#
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
		
		#---compute contrast comparing---#
		for(i in 1:respvar.number)
		{
			is.envFixed <- data$traits[[i]]$analysis$msa$envFixed;
			model <- data$traits[[i]]$analysis$msa$model;
			data$traits[[i]]$analysis$msa$contrast <- list();
			
			genoFactorName <- names(model@frame)[2];
			envFactorName <- names(model@frame)[3];
			trmtLevels <- levels(model@frame[,2]);
			tempGenoContrast <- genoContrast;
			trmtContrast <- as.data.frame(t(tempGenoContrast));
			trmtContrast <- list(genoFactorName = trmtContrast);
			names(trmtContrast) <- genoFactorName;
			if(is.envFixed)
			{
				if(is.null(envContrast))
				{	
					data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
					data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast, across = envFactorName);
				} else
				{
					data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
					tempEnvContrast = envContrast;
					tempEnvContrast = as.data.frame(t(tempEnvContrast));
					trmtContrast$tempEnvContrast <- tempEnvContrast;
					names(trmtContrast) <- c(genoFactorName, envFactorName);
					data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast);
				}
			} else
			{
				data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom=trmtContrast);
				data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- NA;
			}
		}
	} else if(contrastOpt == "Default")
	{
		if(data$pop.type != "PL")
		{
			stop("\tError: The Default value of contrastOpt argument is only implemented on PyramidedLine multi-site analysis outcomes!\n");
		} else 
		{
			if(!missing(genoContrast))
				warning("\tWarning: The genoContrast would be omitted when the contrastOpt is Default!\n");
			
			#---checking the data---#
			biGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T));
			triGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC","Cc","cc")), sep="", lex.order = T));
			quadraGeneLevels <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC","Cc","cc"), D = c("DD","Dd","dd")), sep="", lex.order = T));
			
			respvar.number <- length(data$traits);
			for( i in 1:respvar.number)
			{
				is.envFixed <- data$traits[[i]]$analysis$msa$envFixed;
				model <- data$traits[[i]]$analysis$msa$model;
				trmtLevels <- levels(model@frame[ , 2]);
				envLevels <- levels(model@frame[ , 3]);
				trait.name <- data$traits[[i]]$name;
				
				trmtLevels.length = length(trmtLevels);
				if(trmtLevels.length == 9)
				{
					if(!(all(biGeneLevels %in% trmtLevels)))
					{
						stop(paste("\tError: There are not default coded genotye levels of bigenes design on ", 
										" of ", trait.name, 
										". Cannot proceed with default contrast analysis!\n", sep = ""));
					}
				} else if(trmtLevels.length == 27)
				{
					if(!(all(triGeneLevels %in% trmtLevels)))
					{
						stop(paste("\tError: There are not default coded genotye levels of trigenes design on ", 
										" of ", trait.name, 
										". Cannot proceed with default contrast analysis!\n", sep = ""));
					}
				} else if(trmtLevels.length == 81)
				{
					if(!(all(quadraGeneLevels %in% trmtLevels)))
					{
						stop(paste("\tError: There are not default coded genotye levels of quadragenes design on ", 
										" of ", trait.name, 
										". Cannot proceed with default contrast analysis!\n", sep = ""));
					}
				} else
				{
					stop(paste("\tError: There are not default genotype levels of bi-, tri- or quadra- genes design on ", 
									trait.name, 
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
			} #--- end of statement for( i in 1:respvar.number) ---#
			
			#--- compute the contrast ---#
			
			for(i in 1 : respvar.number)
			{
				is.envFixed <- data$traits[[i]]$analsyis$msa$envFixed;
				model <- data$traits[[i]]$analysis$msa$model;
				data$traits[[i]]$analysis$msa$contrast <- list();
				
				genoFactorName <- names(model@frame)[2];
				envFactorName <- names(model@frame)[3];
				trmtLevels <- levels(model@frame[,2]);
				envLevels <- levels(model@frame[,3]);
				trmtLevels.length <- length(trmtLevels);
				genNum <- data$gene.num;
				
				trmtContrast <- getDefaultGenesContrast(genNum);	
				
				tempGenoContrast <- trmtContrast;
				trmtContrast <- as.data.frame(t(tempGenoContrast));
				trmtContrast <- list(genoFactorName = trmtContrast);
				names(trmtContrast) <- genoFactorName;
				if(is.envFixed)
				{
					if(is.null(envContrast))
					{	
						data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
						data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast, across = envFactorName);
					} else
					{
						data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom= trmtContrast);
						tempEnvContrast = envContrast;
						tempEnvContrast = as.data.frame(t(tempEnvContrast));
						trmtContrast$tempEnvContrast <- tempEnvContrast;
						names(trmtContrast) <- c(genoFactorName, envFactorName);
						data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- testInteractions(model, custom=trmtContrast);
					}
				} else
				{
					data$traits[[i]]$analysis$msa$contrast$contrastOnGeno <- testInteractions(model, custom=trmtContrast);
					data$traits[[i]]$analysis$msa$contrast$contrastAcrossEnv <- NA;
				}
			}
		}#--- end of statement if(data$pop.type != "PL);
	} #--- end of statement if(contrastOpt == "RecurrentParent")
	
	detach("package:phia");
	detach("package:lme4");
	
	return(data);
}
