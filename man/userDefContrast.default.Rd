\name{userDefContrast.default}
\alias{userDefContrast.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
User Define Contrast Default Function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
userDefContrast.default(model, contrast, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{contrast}{
%%     ~~Describe \code{contrast} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (model, contrast, alpha = 0.05) 
%{
%    trmtLevels <- levels(model@frame[, 2])
%    contrMatName <- colnames(contrast)
%    numLevels <- length(trmtLevels)
%    numContrast <- nrow(contrast)
%    if (numLevels != ncol(contrast)) {
%        stop("Error: The number of columns in the contrast matrix is not equal to the levels of the treatment.")
%    }
%    if (any(contrMatName \%in\% trmtLevels == FALSE)) {
%        stop("Error: At least one of the column names in the contrast matrix does not match the levels in the data set.")
%    }
%    if (numContrast >= numLevels) {
%        stop("Error: Too many rows in the contrast matrix.")
%    }
%    if (!isContrast(contrast)) {
%        stop("Error: At least one of the rows in the matrix is not a contrast.")
%    }
%    newContrast <- contrast[, match(trmtLevels, contrMatName)]
%    result <- contrastCompute(model, newContrast, alpha)
%    return(result)
%  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
