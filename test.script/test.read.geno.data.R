###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Oct 8, 2014
# FileName: test.read.geno.data.R
###############################################################################


data <- try(read.csv(file="e://sampledata//geno.data.30.il.csv",header = T),silent=TRUE);
data
data.geno <- data[,1]
data.without.geno <- data[,-1]
data.change <- data;
data.change[1,2] <-4
marker.code <- c(1:3,NA)
all(apply(data[,-1], 2, function(x) all(x %in% marker.code)))
all(apply(data.change[,-1], 2, function(x) all(x %in% marker.code)))
out <- apply(data.change[,-1], 2, function(x) all(x %in% marker.code))
names(out[!out])

debug(read.geno.data)
geno.data.testing <- try(read.geno.data("e://data//sampledata//geno.data.30.il.csv", dp.code = 3, rp.code = 1, ht.code = 2, na.code = NA))
geno.data.testing
undebug(read.geno.data)


library(car)

x <- rep(1:3, 3)
x <- c(x, NA, 1:3)

y <- recode(x, "3=2;2=1;1=0;")
dp.code = 3;
rp.code = 1;
ht.code = 2;
apply(data.without.geno, 2, recode, paste(dp.code , "=", "2;", rp.code,"=", "0;", ht.code,"=", "1;"))

apply(data.without.geno, 2, function(x)
			{
				dp.count <- length(which(x==3));
				rp.count <- length(which(x==1));
				ht.count <- length(which(x==2));
				na.count <- length(which(is.na(x)));
				total <- length(x);
				missing.rate <- na.count / total;
				result <- c(dp.count, rp.count, ht.count, na.count, total, missing.rate);
				names(result) <- c("dp.count", "rp.count", "ht.count", "na.count", "total", "missing.rate");
				return(result);
			})



debug(restrict.geno.data)
geno.data.testing.restrict <- restrict.geno.data(geno.data.testing, 
		missing.rate = 0.80,
		cor.rate = 0.80,
		mono.reduced = TRUE,
		donor.minicount = 2
		)
geno.data.testing.restrict
undebug(restrict.geno.data)


marker.geno <- matrix(sample(c(0,2,NA), replace = T, size = 200, prob= c(0.8, 0.19, 0.01)),, nrow = 40)
colnames(marker.geno) <- paste("M", 1:5, sep = "")
marker.geno.cor <- cor(marker.geno, use="pairwise.complete.obs")
cor(marker.geno[,"M1"], marker.geno[,"M5"], use="complete.obs")
diag(marker.geno.cor) = 0;
iden = apply(marker.geno.cor, 1, function(x) which(x > 0.10));
lens <- unlist(lapply(iden, length));
idx_max = which.max(lens);
dis = NULL;
while(length(iden))
{
	if(lens[idx_max]>1){
#       dis1=snps[idx_max]
		dis=c(dis,idx_max);
		for(i in iden[[idx_max]])
			iden[[i]]=setdiff(iden[[i]],idx_max);
		iden[[idx_max]]=setdiff(iden[[idx_max]],iden[[idx_max]]);
		lens=unlist(lapply(iden,length))
		idx_max=which.max(lens);
	}
	else{
		rest=which(lens==1)
		pairs=data.frame(marker_1=rest,marker_2=0,stringsAsFactors=FALSE);
		for(i in 1:length(rest)){
			marker_2=iden[[rest[i]]];
			pairs[i,2]=ifelse(marker_2 %in% pairs[1:(i-1),1],NA,marker_2);
		}
		pairs=na.omit(pairs)
		for(i in 1:nrow(pairs)){
			inf1 = round()
#			inf1=round(report[report$SNPname==snps[pairs[i,1]],4:5],digits=3)
#			inf2=round(report[report$SNPname==snps[pairs[i,2]],4:5],digits=3)
			dis1=ifelse(inf1[,2]>inf2[,2],pairs[i,1],ifelse(inf1[,1]<inf2[,1],pairs[i,1],pairs[i,2]))
			dis=c(dis,dis1)
		}
		iden=NULL;
	}
}


