## This function is used to get default bi-genes, tri-genes and quadra-genes 
## contrast!
## It is by hard coding in R;
## modify on Jan 9,2015
## implement not hard coding of contrast matrix.
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
    #--- comment this on Jan 9, 2015, change to use new function to get contrast matrix---#
# 		default.contrast <- matrix(c(
# 						0.25,  0,  0.25,	0,	0,	0,	-0.25,	0,	-0.25,
# 						0.25,	0,	-0.25,	0,	0,	0,	0.25,	0,	-0.25,
# 						-0.25,	0,	-0.25,	0.5,	0,	0.5,	-0.25,	0,	-0.25,
# 						-0.25,	0.5,	-0.25,	0,	0,	0,	-0.25,	0.5,	-0.25,
# 						0.25,	0,	-0.25,	0,	0,	0,	-0.25,	0,	0.25,
# 						-0.25,	0.5,	-0.25,	0,	0,	0,	0.25,	-0.5,	0.25,
# 						-0.25,	0,	0.25,	0.5,	0,	-0.5,	-0.25,	0,	0.25,
# 						0.25,	-0.5,	0.25,	-0.5,	1,	-0.5,	0.25,	-0.5,	0.25
# 				),
# 				nrow=8,ncol=9,byrow=T, 
# 				dimnames=list(c("a1","a2","d1","d2","a1a2","a1d2","d1a2","d1d2"),
# 						c("AABB", "AABb",	"AAbb",	"AaBB",	"AaBb",	"Aabb",	"aaBB",	"aaBb",	"aabb"
# 						)));
    default.contrast <- getDefGeneOrthogonalContrastMatrix(2)$orthogonal.contrast.matrix;
	} else if(geneNumber == 3)
	{
		default.contrast <- getDefGeneOrthogonalContrastMatrix(3)$orthogonal.contrast.matrix;
	} else if(geneNumber == 4)
	{
		default.contrast <- getDefGeneOrthogonalContrastMatrix(4)$orthogonal.contrast.matrix;
	}
	default.contrast <- as.matrix(default.contrast);
	return (default.contrast);
}