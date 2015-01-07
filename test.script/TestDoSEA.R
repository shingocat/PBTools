x <- read.pheno.data("e://data//sampledata//small01.csv", 
		type="RAW", pop.type = "IL", na.code = NA, 
		geno = "Gname", exptl.design = "RCB", 
		resp.var = c("DFL", "PH"), block = "Rep", 
		env = "Season")
y <- doSEA(x)
