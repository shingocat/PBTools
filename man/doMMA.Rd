\name{doMMA}
\alias{doMMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Doing Multiple Marker Analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
doMMA(phenotypicData, genotypicData, resp.var, method = c("LASSO", "RIDGE_REGRESSION", "ELASTIC_NET"), siglevel = 0.05, bootstrap = 10, pval.method = c("median", "fdr", "holm", "QA"), family = c("gaussian", "binomial", "possion"), include.ht = TRUE, nfolds = 3, step = 0.1, max.try = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phenotypicData}{
%%     ~~Describe \code{phenotypicData} here~~
}
  \item{genotypicData}{
%%     ~~Describe \code{genotypicData} here~~
}
  \item{resp.var}{
%%     ~~Describe \code{resp.var} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{siglevel}{
%%     ~~Describe \code{siglevel} here~~
}
  \item{bootstrap}{
%%     ~~Describe \code{bootstrap} here~~
}
  \item{pval.method}{
%%     ~~Describe \code{pval.method} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{include.ht}{
%%     ~~Describe \code{include.ht} here~~
}
  \item{nfolds}{
%%     ~~Describe \code{nfolds} here~~
}
  \item{step}{
%%     ~~Describe \code{step} here~~
}
  \item{max.try}{
%%     ~~Describe \code{max.try} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
