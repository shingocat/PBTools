#-------------------------------------------------
# This function performs Stability Analysis
# Author: Alaine A. Gulles
#-------------------------------------------------

stability.analysis <- function(data, respvar, geno, env, method = c("regression", "shukla")) UseMethod("stability.analysis")

stability.analysis.default <- function(data, respvar, geno, env, method = c("regression", "shukla")) {

	method <- match.arg(method)
	if (!is.na(match("RCropStatEnv", search()))) {
		prev.option <- options()$show.error.messages
		options(show.error.messages = FALSE)
	}

	if (is.na(match(respvar, names(data))) ||  
	    is.na(match(geno, names(data))) || 
	    is.na(match(env, names(data)))) { stop("At least one variable name does not match a column in the data frame.") }

	data[,match(geno, names(data))] <- factor(data[,match(geno, names(data))])
	data[,match(env, names(data))] <- factor(data[,match(env, names(data))])
	
	if (method == "regression") {
		result <- list()
		for (i in (1:length(respvar))) {
			result[[i]] <- list()
			site.means <- data.frame(levels(data[, match(env, names(data))]), as.data.frame(tapply(data[, match(respvar[i], names(data))], data[, match(env, names(data))], mean, na.rm = TRUE)))
			colnames(site.means) <- c(env, "site.index")
			rownames(site.means) <- NULL
			temp.data <- subset(data, select = c(geno, env, respvar[i]))
			temp.data <- subset(temp.data, subset = (is.na(temp.data[,match(respvar[i], names(temp.data))]) == F))
			trt.nlevels <- nlevels(temp.data[,match(geno, names(temp.data))]) 
			ge.means.wide <- reshape(temp.data, v.names = respvar[i], idvar = env, timevar = geno, direction = "wide")
			colnames(ge.means.wide) <- gsub(paste(respvar[i], ".", sep = ""), "", colnames(ge.means.wide))
			data.all <- merge(site.means, ge.means.wide, by = env)
			slope <- as.matrix(c(1:(6*trt.nlevels)), nrow = trt.nlevels, ncol = 6)
			dim(slope) <- c(trt.nlevels, 6)

			if (var(data.all$site.index) != 0) {
				for (j in (1:trt.nlevels)) {
					if (length(na.omit(data.all[[j+2]])) < 3) { slope[j,] <- rep(NA, 6) } else {
						model <- lm(data.all[[j+2]] ~ site.index, data.all)
						temp <- summary(model)
						slope[j,] <- c(temp$coef[2,], anova(model)[1, "Mean Sq"], anova(model)[2, "Mean Sq"])
					}
				}
				rownames(slope) <- levels(temp.data[, match(geno,names(temp.data))])
				colnames(slope) <- c("Slope", "SE", "t-value", "Prob", "MSReg", "MSDev")
				slope <- data.frame(slope)
				slope <- subset(slope, subset = (is.na(SE)==F))
			} else { slope <- NULL }
			result[[i]]$respvar <- respvar[i]
			result[[i]]$stability <- slope 
		}
		return(list(result, method = "Stability Analysis using Finlay-Wilkinson Model"))
	}

	if (method == "shukla") {
		if (!is.na(match("package:lme4", search()))) { is.pkgload <- TRUE; detach("package:lme4")} else is.pkgload <- FALSE
		library(nlme)
		result <- list()
		for (i in (1:length(respvar))) {
			result[[i]] <- list()
			result[[i]]$respvar <- respvar[i]
			temp.data <- subset(data, select = c(geno, env, respvar[i]))
			temp.data <- subset(temp.data, subset = (is.na(temp.data[,match(respvar[i], names(temp.data))]) == F))
			#myformula1 <- paste(respvar[i], " ~ 1 + ", geno, " + ", env, sep = "")
			#myformula2 <- paste("varIdent(form = ~ 1 |", geno,")", sep = "")
			command <- paste("gls(formula(", respvar[i]," ~ 1 + ",geno," + ", env,"), weights = varIdent(form = ~ 1 |", geno,"), data = temp.data)", sep = "")
			if (class(try(parse(text = command), TRUE)) == "try-error") {
				result[[i]]$error <- geterrmessage()
			} else {
				model1 <- eval(parse(text = command))
				if (class(try(parse(text = "intervals(model1)"), TRUE)) == "try-error"){ result[[i]]$error <- geterrmessage()  
			}else {
					int <- intervals(model1)
					par <- as.data.frame(int$varStruct)
					par <- rbind(c(1,1,1), par)
					rownames(par) <- levels(temp.data[,match(geno, names(temp.data))])
					sigma <- as.data.frame(int$sigma)
					par$lower <- par$lower*sigma[rownames(sigma) == "lower",]
					par$est.   <- par$est.*sigma[rownames(sigma) == "est.",]
					par$upper <- par$upper*sigma[rownames(sigma) == "upper",]
					result[[i]]$stability <- par
				}

			}
		}
		detach("package:nlme")
		#if (is.pkgload) library(lme4)
		return(list(result, method = "Stability Analysis using Shukla's model"))
	}
}