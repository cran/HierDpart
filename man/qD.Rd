\name{qD}
\alias{qD}
\title{Function for calculating genetic diversity profile (q=0,1,2)}
\description{This function calculates genetic diversity profile (q=0, 1, 2).}
\usage{qD(x, q, ncode)}
\arguments{
  \item{x}{Input data, genepop format.}
  \item{q}{Hill numbers, q=0,1,2.}
  \item{ncode}{The coding type of the genetic file.}
}
\details{This function returns different genetic diversities, q is the order of Hill number, specifically, richness (q=0), exponential Shannon entropy (q=1), and heterozygosity related measure (q=2, ^2D=1/(1-He)).}
\value{Result returns to the diversity profile (q) as you defined.}
\references{
Marcon, E., & Herault, B. (2015). entropart: An R package to measure and partition diversity. Journal of Statistical Software, 67(8).

Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
}
\author{qinxinghu@gmail.com}
\examples{
# example genepop file
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
qD(infile,q=0,ncode=3)
}
\keyword{Diversity profiles}

