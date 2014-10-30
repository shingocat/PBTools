## This function is used to get default bi-genes, tri-genes and quadra-genes 
## contrast!
## It is by hard coding in R;
##
getDefaultGenesContrast <- function(genNumber = 2) UseMethod("getDefaultGenesContrast");
getDefaultGenesContrast.default <- function(geneNumber = 2)
{
	if(length(geneNumber) > 1)
	{
		stop("Error: The argument of geneNumber should be of length one.");
	}
	if(!is.numeric(geneNumber))
	{
		stop("Error: The argument of geneNumber should be numeric!");
	}
	
	if(geneNumber < 2 || geneNumber > 4)
	{
		stop ("Error: The value of genNumber argument should be between 2 and 4!");
	}
	default.contrast <- matrix(NA);
	if(geneNumber == 2)
	{
		default.contrast <- matrix(c(
						0.25,  0,  0.25,	0,	0,	0,	-0.25,	0,	-0.25,
						0.25,	0,	-0.25,	0,	0,	0,	0.25,	0,	-0.25,
						-0.25,	0,	-0.25,	0.5,	0,	0.5,	-0.25,	0,	-0.25,
						-0.25,	0.5,	-0.25,	0,	0,	0,	-0.25,	0.5,	-0.25,
						0.25,	0,	-0.25,	0,	0,	0,	-0.25,	0,	0.25,
						-0.25,	0.5,	-0.25,	0,	0,	0,	0.25,	-0.5,	0.25,
						-0.25,	0,	0.25,	0.5,	0,	-0.5,	-0.25,	0,	0.25,
						0.25,	-0.5,	0.25,	-0.5,	1,	-0.5,	0.25,	-0.5,	0.25
				),
				nrow=8,ncol=9,byrow=T, 
				dimnames=list(c("a1","a2","d1","d2","a1a2","a1d2","d1a2","d1d2"),
						c("AABB", "AABb",	"AAbb",	"AaBB",	"AaBb",	"Aabb",	"aaBB",	"aaBb",	"aabb"
						)));
	} else if(geneNumber == 3)
	{
		
	} else if(geneNumber == 4)
	{
		
	}
	return (default.contrast);
}