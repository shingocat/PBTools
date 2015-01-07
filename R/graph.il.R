###############################################################################
# TODO: Plot Introgression Lines 
#
# [Arguments]
#	geno, 
#	map, 
# 	f.ind - [UNUSED] Default 2, Position of first individual used for plotting.
# 	l.ind -[UNUSED]  Default 2, Position of last individual used for plotting.
#	f.lr - Default 2, Width between color bars and labels (individual names).
#	f.tb - Default 2, Height between color bars and top and bottom labels (chromosome names and marker names).
#   d.h - Default 0.01, Horizontal distance between chromosomes.
#   d.v	- Default 0.01, Vertival distance between chromosomes.
#   p.t -  Default F, logical. If TRUE marker tickmarks are plotted.
#	p.s - Default F, logical. Not Sure What it is used for.
# 	d.t -  Default 50, Length of marker tickmarks.
#   d.map - Default 100, Distance between marker labels.
# 	c.nme - Default "", Character string specifying the names of chromosomes.
#   i.nme - Default "", Character string specifying the names of individuals.
#   color -  Character string specifying the colors to be used for plotting. Number of colors must be in accordance with z.min and z.max.
#   z.min -  Default 0, Minimum number of colors used for plotting, always set as 0.
#   z.max -  Default 12, Maximum number of colors used for plotting.
#   plt.pdf - Default F, logical. If TRUE a pdf file is generated.
#	plt.png - Default F, logical. If TRUE a png file is generated.
#   plt.name - Default il.plot, Character string indicating the name of the file to be generated.
#   plt.width - Default 10, Width of generated pdf file.
#   plt.height - Default 10, Height of generated pdf file.
#   plt.ptsize - Default 10, Point size of generated pdf file.

# Author: SelectionTools R package and modify by mao qin
# Date: Jan 4, 2015
# FileName: graph.il.R
###############################################################################


graph.il <- function
		(
		geno,
		map, 
		f.lr = 1, 
		f.tb = 1, 
		d.h = 0.02, 
		d.v = 0.01, 
		p.t = F, 
		p.s = F,
		d.t = 50, 
		d.map = 100, 
		c.nme = "", 
		i.nme = "", 
		color =c("yellow", "blue", "red","green", "wheat", "skyblue", "tomato", "palegreen", "yellow4", "darkblue", "darkblue", "darkgreen"), 
		z.min = 0, 
		z.max = 12, 
		cex.chrom = 1, 
		cex.ind = 1, 
		plt.png = F, 
		plt.pdf = F, 
		plt.name = "il.plot", 
		plt.width = 10, 
		plt.height = 10, 
		plt.ptsize = 10
)
{
	if(!inherits(geno, "GenotypicData"))
		stop("\tError: The geno argument should be of class GenotypicData!\n");
	if(!inherits(map, "MapData"))
		stop("\tError: The map argument should be of class MapData!\n");
	if(missing(f.lr))
		f.lr <- 2;
	if(missing(f.tb))
		f.tb <- 2;
	if(missing(d.h))
		d.h <- 0.01;
	if(missing(d.v))
		d.v <- 0.01;
	if(missing(p.t))
		p.t <- F;
	if(missing(p.s))
		p.s <- F;
	if(missing(d.t))
		d.t <- 50;
	if(missing(d.map))
		d.map <- 100;
	if(missing(c.nme))
		c.nme <- "";
	if(missing(i.nme))
		i.nme <- "";
	if(missing(color))
		color <- c("yellow", "blue", "red","green", "wheat", "skyblue", "tomato", "palegreen", "yellow4", "darkblue", "darkblue", "darkgreen"); 
	if(missing(z.min))
		z.min <- 0;
	if(missing(z.max))
		z.max <- 12;
	if(missing(cex.chrom))
		cex.chrom <- 1;
	if(missing(cex.ind))
		cex.ind <- 1;
	if(missing(plt.png))
		plt.png <- F;
	if(missing(plt.pdf))
		plt.pdf <- F;
	if(missing(plt.name))
		plt.name <- "il.plot";
	if(missing(plt.width))
		plt.width <- 10;
	if(missing(plt.height))
		plt.height	<- 10;
	if(missing(plt.ptsize))
		plt.ptsize <- 10;
	
	#reformating data into original data structure---#
	#--- the genotypic data---#
	g <- geno$processed$data;
	geno.names <- g[ , 1];
	marker.names <- colnames(g)[-1];
	g <- t(g[ , -1]);
	g <- as.data.frame(g);
	g <- apply(g, 2, function(x){names(x) = NULL; as.factor(x)});
	g <- as.data.frame(g);
	names(g) <- geno.names;
	rownames(g) <- marker.names;
	
	
	#--- the map data---#
	m <- map$data;
	marker.label <- map$marker.label;
	chr.label <- map$chr.label;
	pos.label <- map$pos.label;
	chr.number <- nlevels(as.factor(m[,chr.label]));
	chr.names <- as.character(levels(as.factor(m[,chr.label])));
	m <- m[,c(chr.label, pos.label, marker.label)]
	m[,1] <- as.integer(as.character(m[,1]));
	m[,2] <- as.numeric(as.character(m[,2]));
	
	close.screen(all = TRUE)
	if (TRUE == plt.pdf) {
		pdf(file = paste(plt.name,".pdf",sep=""), width = plt.width, height = plt.height, 
				pointsize = plt.ptsize);
	}else if (TRUE == plt.png) {
		png(file = paste(plt.name,".png",sep=""), width = plt.width, height = plt.height, 
				pointsize = plt.ptsize, res = 300, units = "in");
	}
	
	par(mar = c(0, 0, 0, 0));
	#--- number of individuals---#
	n.ind <- ncol(g)
	clr <- color
	n.chrom <- max(m[, 1])
	#--- specifying the chromosome name---#
	if ((1 == length(c.nme)) && ("" == c.nme[1])) {
		c.nme <- paste("Chr. ", 1:n.chrom, sep = "")
	}
	#--- specifying the names of individuals---#
	if ((1 == length(i.nme)) && ("" == i.nme[1])) {
		i.nme <- names(g)
	}
	#--- a loop for storing each chromosome ---#
	{
		#--- no. of marker of each chromosome---#
		m.chrom <- 1:n.chrom
		#--- lenght of each chromosome---#
		l.chrom <- 1:n.chrom
		#--- break of each chromosome---#
		b.chrom <- vector("list", n.chrom)
		#--- position of each chromosome---#
		p.chrom <- vector("list", n.chrom)
		for (c in 1:n.chrom) {
			#---index of each chromosome in map data ---#
			idx <- m[, 1] == c;
			m.chrom[c] <- sum(idx)
			#--- each marker position of this chromosome--#
			m.pos <- m[idx, ][, 2]
			l.chrom[c] <- max(m.pos) - min(m.pos)
			#--- marker interval of each chromosome---#
			m.brk <- 0:length(m.pos)
			
			for (j in 2:length(m.pos)) {
				m.brk[j] <- m.pos[j - 1] + (m.pos[j] - m.pos[j - 
											1])/2
			}
			m.brk[1 + length(m.pos)] <- m.pos[length(m.pos)]
			b.chrom[[c]] <- m.brk
			p.chrom[[c]] <- m.pos
		}
	}
	#--- number of drawing column, adding one more column for individuals names--#
	ndc <- 1 + n.chrom
	
	f.disp <- c(f.lr * mean(l.chrom), l.chrom)
	f.disp <- f.disp/sum(f.disp)
	s.disp <- f.disp
	for (kk in 2:ndc) s.disp[kk] <- s.disp[kk - 1] + s.disp[kk]
	s.disp <- c(0, s.disp)
	dcol <- matrix(0, ncol = 2, nrow = ndc)
	for (kk in 1:ndc) {
		dcol[kk, 1] = s.disp[kk]
		dcol[kk, 2] = s.disp[kk + 1]
	}
	ndr <- 2 + n.ind
	f.disp <- c(f.tb, rep(1, n.ind), f.tb)
	f.disp <- f.disp/sum(f.disp)
	s.disp <- f.disp
	for (kk in 2:ndr) s.disp[kk] <- s.disp[kk - 1] + s.disp[kk]
	s.disp <- c(0, s.disp)
	drow <- matrix(0, ncol = 2, nrow = ndr)
	for (kk in 1:ndr) {
		drow[kk, 1] = s.disp[kk]
		drow[kk, 2] = s.disp[kk + 1]
	}
	figs <- matrix(ncol = 4, nrow = ndr * ndc)
	for (i in 1:ndr) for (c in 1:ndc) {
			idx <- (i - 1) * ndc + c
			figs[idx, ] <- c(dcol[c, ], drow[i, ])
		}
	sve <- 1 - figs[, 3]
	figs[, 3] <- 1 - figs[, 4]
	figs[, 4] <- sve
	figs[, 1] <- figs[, 1] + d.h
	figs[, 2] <- figs[, 2] - d.h
	figs[, 3] <- figs[, 3] + d.v
	figs[, 4] <- figs[, 4] - d.v
	figs <- round(figs, digits = 5)
	neu <- split.screen(figs)
	#--- plot chromosome info---#
	for (c in 1:n.chrom) {
		scr <- (c + 1)
		screen(scr)
		plot(0, 0, t = "n", axes = F, xlab = "", ylab = "")
		text(0, 0, c.nme[c], cex = cex.chrom)
		scr <- (ndr - 1) * ndc + (c + 1)
		screen(scr)
		mi <- min(b.chrom[[c]])
		ma <- max(b.chrom[[c]])
		if ((p.t == T) || (p.s == T)) {
			plot(0, 0, axes = F, t = "n", xlab = "", ylab = "", 
					xlim = c(mi, ma), ylim = c(-d.t, d.t))
			lines(c(mi, ma), c(0, 0))
		}
		if (p.t == T) {
			for (mm in 1:length(p.chrom[[c]])) lines(c(p.chrom[[c]][mm], 
								p.chrom[[c]][mm]), c(0, d.t/2))
		}
		if (p.s == T) {
			t <- seq(0, ma, d.map)
			text(t, rep(-d.t/2, length(t)), t)
		}
	}
	#--- plot individuals name---#
	for (i in 1:n.ind) {
		scr <- (i) * ndc + (1)
		screen(scr)
		plot(0, 0, t = "n", axes = F, xlab = "", ylab = "")
		text(0, 0, i.nme[i], cex = cex.ind)
	}
	#--- plot chromosome bar info ---#
	for (i in 1:n.ind) for (c in 1:n.chrom) {
			scr <- (i) * ndc + (c + 1)
			screen(scr)
			p.m <- matrix(0, nrow = m.chrom[c], ncol = 2, byrow = F)
			idx <- m[, 1] == c
			m.nme <- as.character(m[idx, ][, 3])
			idx2 <- row.names(g) %in% m.nme
			p.g <- as.character(g[idx2, i])
			splt <- strsplit(p.g, "/")
			for (p in 1:m.chrom[c]) {
				p.m[p, 1] <- as.numeric(splt[[p]][1])
				p.m[p, 2] <- as.numeric(splt[[p]][2])
			}
			image(x = b.chrom[[c]], y = c(0, 1, 2), z = p.m, axes = 0, 
					zlim = c(z.min, z.max), col = clr)
			box()
		}
	close.screen(all = TRUE)
	if ((TRUE == plt.pdf) | (TRUE == plt.png)) 
		dev.off()
}

