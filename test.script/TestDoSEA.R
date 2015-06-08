x <- read.pheno.data("e://data//sampledata//small01.csv", 
		type="RAW", pop.type = "IL", na.code = NA, 
		geno = "Gname", exptl.design = "RCB", 
		resp.var = c("DFL", "PH"), block = "Rep", 
		env = "Season")
y <- doSEA(x)


x <- read.pheno.data("d://data//test//bigenesT1.csv", 
		type="RAW", gene.num=2, pop.type = "PL", na.code = NA, 
		geno = "Genotype", exptl.design = "RCB", 
		resp.var = c("Heading"), block = "Block"
		
		)
x.sea <- doSEA(x)
