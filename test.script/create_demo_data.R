###############################################################################
# TODO: Add comment
# 
# Author: mqin
# Date: Sep 15, 2014
# FileName: create_demo_data.R
###############################################################################


sssl.geno.with.recurrent <- c("HJX", levels(interaction("S", 1:20, sep="")));
sssl.geno.without.recurrent <- c(levels(interaction("S", 1:20, sep="")));
block <- c("Block1","Block2","Block3","Block4");
# For yield 
x <- c(80:95);
yield1 <- sample(x, size = 84, replace = T) + rnorm(84, mean = 0, sd = 2);
yield2 <- sample(x, size = 80, replace = T) + rnorm(80, mean = 0, sd = 2);
yield3 <- sample(x, size = 84, replace = T) + rnorm(84, mean = 4, sd = 2);
yield4 <- sample(x, size = 84, replace = T) + rnorm(84, mean = 6, sd = 2);
yield5 <- sample(x, size = 84, replace = T) + rnorm(84, mean = 6, sd = 2);
yield6 <- sample(x, size = 84, replace = T) + rnorm(84, mean = 6, sd = 2);
yield7 <- sample(x, size = 84, replace = T) + rnorm(84, mean = 6, sd = 2);
for( i in 1:length(yield4))
{
	if(i %% 2 == 0)
	{
		yield4[i] <- NA;
	}
}
yield6[1:4] <- NA;
yield7[1:2] <- NA;
yield7[20:21] <- NA;
yield7[30:31] <- NA;
y <- c(7:15);
# For Tiller Num
tn1 <- sample(y, size = 84, replace = T) + rnorm(84, mean = 0, sd = 2);
tn2 <- sample(y, size = 80, replace = T) + rnorm(80, mean = 0, sd = 2);
tn3 <- sample(y, size = 84, replace = T) + rnorm(84, mean = 4, sd = 2);
tn4 <- sample(y, size = 84, replace = T) + rnorm(84, mean = 6, sd = 2);
tn5 <- sample(y, size = 84, replace = T) + rnorm(84, mean = 6, sd = 2);
tn6 <- sample(y, size = 84, replace = T) + rnorm(84, mean = 6, sd = 2);
tn7 <- sample(y, size = 84, replace = T) + rnorm(84, mean = 6, sd = 2);

for( i in 1:length(tn4))
{
	if(i %% 2 == 0)
	{
		tn4[i] <- NA;
	}
}
tn6[1:4] <- NA;
tn7[1:2] <- NA;
tn7[20:21] <- NA;
tn7[30:31] <- NA;
sssl.data1 <- data.frame(Genotype=rep(sssl.geno.with.recurrent, each=4), Block = rep(block), Yield = yield1, TN = tn1 );
sssl.data1$Env <- "E1";
sssl.data2 <- data.frame(Genotype=rep(sssl.geno.without.recurrent, each=4), Block = rep(block), Yield = yield2, TN = tn2 );
sssl.data2$Env <- "E2";
sssl.data3 <- data.frame(Genotype=rep(sssl.geno.with.recurrent, each=4), Block = rep(block), Yield = yield3, TN = tn3 );
sssl.data3$Env <- "E3";
sssl.data4 <- data.frame(Genotype=rep(sssl.geno.with.recurrent, each=4), Block = rep(block), Yield = yield4, TN = tn4 );
sssl.data4$Env <- "E4";
sssl.data5 <- data.frame(Genotype=rep(sssl.geno.with.recurrent, each=4), Block = rep(block), Yield = yield5, TN = tn5 );
sssl.data5$Env <- "E5";
sssl.data6 <-  data.frame(Genotype=rep(sssl.geno.with.recurrent, each=4), Block = rep(block), Yield = yield6, TN = tn6 );
sssl.data6$Env <- "E6";
sssl.data7 <-  data.frame(Genotype=rep(sssl.geno.with.recurrent, each=4), Block = rep(block), Yield = yield7, TN = tn7 );
sssl.data7$Env <- "E7";

sssl.data.without.rp <- rbind(sssl.data1, sssl.data2, sssl.data3);
sssl.data.with.na <- rbind(sssl.data1, sssl.data4, sssl.data5);
sssl.data.balanced.3env <- rbind(sssl.data1, sssl.data3, sssl.data5)
sssl.data.missing.2 <- rbind(sssl.data1, sssl.data3, sssl.data7)
sssl.data.with.name.without.values <- rbind(sssl.data1, sssl.data6, sssl.data7)

pl.9geno <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T));
py.6geno <- levels(interaction(list(A = c("AA","aa"),B = c("BB","Bb","bb")), sep="", lex.order = T));
py.27geno <- levels(interaction(list(A = c("AA","Aa","aa"),B = c("BB","Bb","bb"), C = c("CC", "Cc","cc")), sep="", lex.order = T));
py.yield1 <- sample(x, size = 36, replace = T) + rnorm(36, mean = 0, sd = 2);
py.yield2 <- sample(x, size = 36, replace = T) + rnorm(36, mean = 3, sd = 2);
py.yield3 <- sample(x, size = 36, replace = T) + rnorm(36, mean = 6, sd = 2);
py.yield4 <- sample(x, size = 24, replace = T) + rnorm(24, mean = 2, sd = 2);
py.yield5 <- sample(x, size = 108, replace = T) + rnorm(108, mean = 3, sd = 2);

py.tn1 <- sample(y, size = 36, replace = T) + rnorm(36, mean = 0, sd = 2);
py.tn2 <- sample(y, size = 36, replace = T) + rnorm(36, mean = 2, sd = 2);
py.tn3 <- sample(y, size = 36, replace = T) + rnorm(36, mean = 4, sd = 2);
py.tn4 <- sample(y, size = 24, replace = T) + rnorm(24, mean = 2, sd = 2);
py.tn5 <- sample(y, size = 108, replace = T) + rnorm(108, mean = 2, sd = 2);

py.data1 <- data.frame(Genotype = rep(pl.9geno, each = 4), Block = rep(block), Yield = py.yield1, TN = py.tn1);
py.data1$Env <- "E1";
py.data2 <- data.frame(Genotype = rep(pl.9geno, each = 4), Block = rep(block), Yield = py.yield2, TN = py.tn2);
py.data2$Env <- "E2";
py.data3 <- data.frame(Genotype = rep(pl.9geno, each = 4), Block = rep(block), Yield = py.yield3, TN = py.tn3);
py.data3$Env <- "E3";
py.data4 <- data.frame(Genotype = rep(py.6geno, each = 4), Block = rep(block), Yield = py.yield4, TN = py.tn4);
py.data4$Env <- "E4";
py.data5 <- data.frame(Genotype = rep(py.27geno, each = 4), Block = rep(block), Yield = py.yield5, TN = py.tn5);
py.data5$Env <- "E5";
py.data.balanced <- rbind(py.data1, py.data2, py.data3); # balanced data;
py.data.unbalanced <- rbind(py.data1, py.data2, py.data3, py.data4); # unbalanced data;
py.data.balanced.3.genes <- rbind(py.data1, py.data2, py.data5); # former is two genes, later is three genes;










