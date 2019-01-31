\name{qDplot}
\alias{qDplot}
\title{Plot the genetic diversity profiles (q=0,1,2)}
\description{This function calculates the diversity profiles and generates a line plot. You can get either one of the diversity index plot, richness (q=0), exponential Shannon entropy (q=1), and heterozygosity related measure of diversity (q=2), or all of them.}
\usage{qDplot(x, q, ncode)}
\arguments{
  \item{x}{ The input genetic files, genepop format}
  \item{q}{The order of Hill numbers,  q =0,1,2, or "all".If q= "all", this will return plot including three diversities.}
  \item{ncode}{ The coding type of your data.}
}
\details{Result returns a diversity plot.}
\references{
Marcon, E., & Herault, B. (2015). entropart: An R package to measure and partition diversity. Journal of Statistical Software, 67(8).

Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{qinxinghu@gmail.com}
\seealso{qD}
\examples{
# example genepop file
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
# plot allelic richness
qDplot(infile,q=0,ncode=3)
#plot all three diversity profiles
qDplot(infile,q="all",ncode=3)}
\keyword{Plot}
\keyword{Diversity profiles}
