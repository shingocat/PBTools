debug(read.pheno.data)
x <- read.pheno.data(phenodata="e://data//sampledata//bigenes.csv", type = "RAW", pop.type="PL", gene.num = 2, exptl.design="RCB", resp.var=c("Heading", "PH"), geno="Genotype", block = "Block", env="Env", na.code = NA)
read.pheno.data(phenodata="e://data//Small01.csv", type = "RAW", pop.type = "IL", exptl.design="AugLS", resp.var=c("DFL", "PH"), geno="Gname", row = "RowBlock", column ="ColBlock", env="Season", na.code = NA)
y <- read.pheno.data("e://data//xiaolong//TQ_IL_Pheno.csv", type ="MEAN", pop.type = "IL", geno = "L-Code", resp.var=c("PH", "HD", "PN", "PL"), na.code = ".")
TQ_IL <- read.csv("e://data//xiaolong//TQ_IL_Pheno.csv", header = TRUE, check.names = FALSE)
read.pheno.data(TQ_IL, type ="MEAN", pop.type ="IL", geno = "L-Code", resp.var=c("PH", "HD", "PN", "PL"), na.code = ".")

undebug(read.pheno.data)
debug(restrict.pheno.data)
y <- restrict.pheno.data(x);
undebug(restrict.pheno.data)
debug(doSSA)
z <- doSSA(y)
undebug(doSSA)





























