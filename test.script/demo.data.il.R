###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Sep 25, 2014
# FileName: demo.data.il.R
###############################################################################


marker.names <- paste("M", 1:100, sep = "");
marker.names.200 <- paste("M", 1:200, sep = "");

line.names <- paste("L", 1:100, sep="")

bc.genotype <- c(2,1); # 1 for homozygous on recurrent parent, 2 for heterozygous on recurrent parent.
f.genotype <- c(3,2,1); # 3 for donor parent;
base.Aa.freq <- 1/2. # on bc generation 1;
for(i in 1:10) # for 10 bc generation
{	
	if(i == 1) next;
	base.Aa.freq <- c(base.Aa.freq, (1/2)^i);
}
base.aa.freq <- 1 - base.Aa.freq;
bc.generation.freq <- cbind(Aa = base.Aa.freq, aa=base.aa.freq);


# for bc1
geno.data.bc1 <- matrix(NA, nrow = 100, ncol = 100);
for(i in 1 : length(marker.names))
{
	geno.data.bc1[,i] <- sample(bc.genotype, size = 100, replace = TRUE, prob=bc.generation.freq[1,]);
}
rownames(geno.data.bc1) <- c(line.names);
colnames(geno.data.bc1) <- marker.names;
geno.data.bc1 <- apply(geno.data.bc1, 2, factor);
summary(geno.data.bc1)

# for bc2
geno.data.bc2 <- matrix(NA, nrow = 100, ncol = 100);
for(i in 1 : length(marker.names))
{
	geno.data.bc2[,i] <- sample(bc.genotype, size = 100, replace = TRUE, prob=bc.generation.freq[2,]);
}
rownames(geno.data.bc2) <- c(line.names);
colnames(geno.data.bc2) <- marker.names;
geno.data.bc2 <- apply(geno.data.bc2, 2, factor);


# for bc5 
geno.data.bc5 <- matrix(NA, nrow = 100, ncol = 100);
for(i in 1 : length(marker.names))
{
	geno.data.bc5[,i] <- sample(bc.genotype, size = 100, replace = TRUE, prob=bc.generation.freq[5,]);
}
rownames(geno.data.bc5) <- c(line.names);
colnames(geno.data.bc5) <- marker.names;
geno.data.bc5 <- apply(geno.data.bc5, 2, factor);
summary(geno.data.bc5)

# for bc5f2
geno.data.bc5f2 <- matrix(NA, nrow= 100, ncol = 200);
for(i in 1:length(marker.names.200))
{
	geno.data.bc5f2[,i] <- sample(f.genotype, size= 100, replace = TRUE, prob=getExpectGenotypicFreq(5,2,2));
}
rownames(geno.data.bc5f2) <- c(line.names);
colnames(geno.data.bc5f2) <- c(marker.names.200);
geno.data.bc5f2 <- apply(geno.data.bc5f2, 2, factor);
summary(geno.data.bc5f2);

# for only 
geno.data.30.il <- matrix(NA, nrow = 30, ncol = 60);
for(i in 1:60)
{
	geno.data.30.il[,i] <- sample(f.genotype[c(1,3)], size= 30, replace = TRUE, prob=getExpectGenotypicFreq(5,2,2)[c(1,3)]);
}
rownames(geno.data.30.il) <- c(paste("L", 1:29, sep = ""),"Recurrent");
colnames(geno.data.30.il) <- paste("M", 1:60, sep = "");
geno.data.30.il[30,] <- 1; ## all marker to 1;
geno.data.30.il[c(1:29),1] <- c(3, rep(NA, 28)); # only one line genotyped of donor genotype;
geno.data.30.il[c(1:29),2] <- c(3);  # all the marker genotype are donor genotype;
geno.data.30.il[c(1:29), 3] <- c(rep(3,15),rep(1,14)); # half of marker genotype are donor genotype
geno.data.30.il[c(1:29),11] <- c(3, rep(NA, 28)); # only one line genotyped of donor genotype;
geno.data.30.il[c(1:29),21] <- c(3);  # all the marker genotype are donor genotype;
geno.data.30.il[c(1:29), 27] <- c(rep(3,15),rep(1,14)); # half of marker genotype are donor genotype
geno.data.30.il <- apply(geno.data.30.il, 2, factor);

# for bc3f1 using chisq
geno.data.40.il <- matrix(NA, nrow = 40, ncol = 90);
for(i in 1:90)
{
	geno.data.40.il[,i] <- sample(f.genotype[c(1,3)], size= 40, replace = TRUE, prob=getExpectGenotypicFreq(3,1,2)[c(1,3)]);
}
colnames(geno.data.40.il) <- paste("M", 1:90, sep = "");
rownames(geno.data.40.il) <- c(paste("L", 1:40, sep = ""));
geno.data.40.il[c(1:40),1] <- c(3,rep(NA,39)); # only one donor genotype
geno.data.40.il[c(1:40),2] <- NA; # all genotype missing
geno.data.40.il[c(1:40),3] <- c(1); # all recurrent genotype
geno.data.40.il[c(1:40),4] <- c(3); # all donor genotype
geno.data.40.il[c(1:40),5] <- c(3,3,1,1,rep(NA,36)); # two donor and recurrent
geno.data.40.il[c(1:40),6] <- c(3,3,3,1,1,1, rep(NA, 34)); # three donor and recurrent
geno.data.40.il[c(1:40),7]<- c(3,3,3,3,1,1,rep(NA,34)); # four donor and two recurrent
geno.data.40.il[c(1:40),8]<- c(3,3,3,rep(1,10),rep(NA,27)); # four donor and two recurrent
geno.data.40.il[c(1:40),9]<- c(3,3,1,1,1,rep(NA,35)); # four donor and two recurrent
geno.data.40.il[c(1:40),10]<- c(3,3,3,rep(1,15),rep(NA,22)); # four donor and two recurrent

#########################################################
### single marker analysis
#########################################################
Var   = c(58.5,59.5,77.8,80.9,84.0,83.6,70.1,68.3,69.8,69.8,56.0,54.5, 50.7,49.3,63.8,65.8,56.6,57.5,77.8,79.2,69.9,69.2,62.1,64.5)

Mgeno = c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)
Geno  = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12)
dat<-cbind(Mgeno,Geno,Var)
dat <- as.data.frame(dat);
str(dat)
dat$Mgeno <- factor(dat$Mgeno)
dat$Geno <- factor(dat$Geno)

dat$rep<-as.factor(rep(c(1,2),12))
dat$b<-as.factor(rep(c(1,1,2,2,3,3,4,4),3))

lines <- rep(1:10, each= 4);
block <- 1:4;
M1 <- rep(c(1,3), each = 4);
M2 <- rep(c(3,1), each = 4);
M3 <- rep(c(1,3), each = 20);
var <- rnorm(40, mean = 100, sd = 5)
demo.data.sm <- cbind(lines, M1, M2, M3, block,var)
demo.data.sm <- as.data.frame(demo.data.sm);
demo.data.sm$lines <- factor(demo.data.sm$lines);
demo.data.sm$block <- factor(demo.data.sm$block);
demo.data.sm$M1 <- factor(demo.data.sm$M1);
demo.data.sm$M2 <- factor(demo.data.sm$M2);
demo.data.sm$M3 <- factor(demo.data.sm$M3);


a <- read.table(header = TRUE, text = "
				a b
				1 2
				3 4
				0 NA
				")
b <- read.table(header = TRUE, text = "
				a b
				1 2
				1 0
				0 NA
				")












