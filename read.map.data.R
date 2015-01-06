###############################################################################
# TODO: reading map data into R. The map file should be of mandatory format.
#		first column is chromosome[Integer value]
#		second column is position[numeric value]
#		third column is marker[character value]
#
# [Arguments]
#	source - the source of map data;
#	marker - the indicator of marker variable in the map data;
#	chr - the indicator of chromosome variable in the map data;
#	pos	- the indicator of position varaible in the map data;
# 	units - the units of map data, only accept cM(centiMorgan) or bp(base pair)
#
# Author: mqin
# Date: Dec 26, 2014
# FileName: read.map.data.R
###############################################################################

read.map.data <- function(
		source,
		marker,
		chr,
		pos,
		units = c("cM", "bp"),
		sep,
		header = TRUE,
		...
)
{
	if(missing(source))
		stop("\tError: The argument of source could not be null!\n");
	UseMethod("read.map.data");
}

#--- s3 method for reading from file---#
read.map.data.character <- function(
		source,
		marker,
		chr,
		pos,
		units = c("cM", "bp"),
		sep,
		header = TRUE,
		...
)
{
	isMarkAssign <- TRUE;
	isChrAssign <- TRUE;
	isPosAssign <- TRUE;
	
	#	refacting 
	if(missing(marker))
	{
#		stop("\tError: The argument of marker could not be null!\n");
		marker <- 3;
		isMarkAssign <- FALSE;
	}
	if(missing(chr))
	{
#		stop("\tError: The argument of chr could not be null!\n");
		chr <- 1;
		isChrAssign <- FALSE;
	}
	if(missing(pos))
	{
#		stop("\tError: The argument of pos could not be null!\n");
		pos <- 2;
		isPosAssign <- FALSE;
	}
	
	if(!file.exists(source))
		stop("\tError: The specified file does not exist!\n");
	if(missing(sep))
	{
		source <- try(read.csv(file=source, header=header, check.names = FALSE), silent=TRUE);
	} else
	{ 
		source <- try(read.table(file = source, header = header, check.names = FALSE, sep = sep), silent=TRUE);
	}
	if(inherits(source, "try-error"))
		stop("\tError: There are some promblem occurred when read the file!\n");
	units <- match.arg(units);
	data <- inner.read.map.data(source, marker, chr, pos, units,isMarkAssign, isChrAssign, isPosAssign);
	return(data);
}

#--- s3 method for reading from data.frame---#
read.map.data.data.frame <- function(
		source,
		marker,
		chr,
		pos,
		units = c("cM", "bp"),
		...
)
{
	isMarkAssign <- TRUE;
	isChrAssign <- TRUE;
	isPosAssign <- TRUE;
	
	if(missing(marker))
	{
#		stop("\tError: The argument of marker could not be null!\n");
		marker <- 3;
		isMarkAssign <- FALSE;
	}
	if(missing(chr))
	{
#		stop("\tError: The argument of chr could not be null!\n");
		chr <- 1;
		isChrAssign <- FALSE;
	}
	if(missing(pos))
	{
#		stop("\tError: The argument of pos could not be null!\n");
		pos <- 2;
		isPotAssign <- FALSE;
	}
	
	units <- match.arg(units);
	data <- inner.read.map.data(source, marker, chr, pos, units,isMarkAssign, isChrAssign, isPosAssign);
	return(data);
}

inner.read.map.data <- function
(
	source,
	marker,
	chr,
	pos,
	units,
	isMarkAssign,
	isChrAssign,
	isPosAssign
)
{
	#--- checking all the specified indicator in the source---#
	col.names <- colnames(source);
	if(isMarkAssign & !(marker %in% col.names))
		stop("\tError: The sepcified marker varaible does not in the source!\n");
	if(isChrAssign & !(chr %in% col.names))
		stop("\tError: The sepcified chr variable does not in the source!\n");
	if(isPosAssign & !(pos %in% col.names))
		stop("\tError: The sepcified pos variable does not in the source!\n");
	source <- apply(source,2,trimStrings);
	source <- as.data.frame(source);
	source[[marker]] <- as.factor(source[,marker]);
	source[[chr]] <- as.factor(source[,chr]);
	source[[pos]] <- as.numeric(as.character(source[,pos]));
	marker.number <- nlevels(as.factor(source[,marker]));
	chr.number <- nlevels(as.factor(source[,chr]));
	chr.names <- levels(source[,chr]);
	marker.number.by.chr <- aggregate(source[,marker], by = list(source[,chr]), FUN = function(x){return(length(x))});
	colnames(marker.number.by.chr) <- c(chr,marker);
	max.length.by.chr <- aggregate(source[,pos], by = list(source[,chr]), FUN = function(x){x <- as.numeric(as.character(x));return(max(x))});
	colnames(max.length.by.chr) <- c(chr,"Length");
	summary.out <- merge(max.length.by.chr, marker.number.by.chr, by = chr, all = TRUE);
	colnames(summary.out) <- c("Chr.", "Length", "Markers");
	outcomes <- list();
	outcomes$marker.label <- marker;
	outcomes$chr.label <- chr;
	outcomes$pos.label <- pos;
	outcomes$units <- units;
	outcomes$marker.number <- marker.number;
	outcomes$chr.number <- chr.number;
	outcomes$summary.out <- summary.out;
	outcomes$data <- source;
	outcomes$chrs <- list();
	for(i in 1:chr.number)
	{
		outcomes$chrs[[i]] <- list();
		outcomes$chrs[[i]]$name <- chr.names[i];
		outcomes$chrs[[i]]$map <- subset(source[,c(marker,pos)], subset = (source[,chr] == chr.names[i]));
	}
	class(outcomes) <- "MapData";
	return(outcomes);
}

print.MapData <- function(data)
{
	cat("Map Data:\n");
	cat("units:", data$units, ".\n", sep = "");
	cat("There are ", data$marker.number, " Markers.\n", sep = "");
	cat("There are ", data$chr.number, " chromosomes.\n", sep = "");
	cat(rep("-",times=50), sep = "");
	cat("\n");
	cat("Map Summary:\n");
	print(data$summary.out, row.names=FALSE);
	cat("\n");
}