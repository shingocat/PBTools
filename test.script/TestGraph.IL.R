###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Jan 4, 2015
# FileName: TestGraph.IL.R
###############################################################################

geno <- read.geno.data("E://Demo//SelectionTools-examples//il.demo.geno", sep="\t", dp.code = 3, ht.code = 2, rp.code = 1,na.code=".", BCn = 2, Fn = 3)
map <- read.map.data("E://Demo//SelectionTools-examples//il.demo.map",sep=",", marker= "marker", chr="chr", pos="pos")

graph.il(geno, map,
		f.lr = 0.5, 
		f.tb = 1, 
		d.h=0.001, 
		d.v =0.001, 
		p.t=F, 
		p.s=F, 
		d.t=50, 
		d.map=100, 
		c.nme="", 
		i.nme="", 
		color=c("yellow", "blue", "red","green", "wheat", "skyblue", "tomato", "palegreen", "yellow4", "darkblue", "darkblue", "darkgreen"), 
		z.min=0, 
		z.max=12, 
		cex.chrom=1, 
		cex.ind = 1, 
		plt.png=T, 
		plt.pdf=F, 
		plt.name = "GGT", 
		plt.width=10, 
		plt.height=10, 
		plt.ptsize =10)

p <- read.table("E:\\Demo\\SelectionTools-examples\\OUTPUT-FILES\\tmp-76810.mpo")
p <- t(p)
p <- as.data.frame(p)
geno.c <- read.geno.data(p, dp.code = "2/2", ht.code = "1/2", rp.code = "1/1",na.code=".", BCn = 2, Fn = 3)

m <- read.table("E:\\Demo\\SelectionTools-examples\\OUTPUT-FILES\\tmp-76810.mmp")
map.c <- read.map.data(m, marker="V3", chr="V1", pos="V2")
graph.il(geno.c, map.c,
		f.lr = 0.5, 
		f.tb = 1, 
		d.h=0.001, 
		d.v =0.001, 
		p.t=F, 
		p.s=F, 
		d.t=50, 
		d.map=100, 
		c.nme="", 
		i.nme="", 
		color=c("yellow", "blue", "red","green", "wheat", "skyblue", "tomato", "palegreen", "yellow4", "darkblue", "darkblue", "darkgreen"), 
		z.min=0, 
		z.max=12, 
		cex.chrom=1, 
		cex.ind = 1, 
		plt.png=F, 
		plt.pdf=F, 
		plt.name = "GGT", 
		plt.width=10, 
		plt.height=10, 
		plt.ptsize =10)