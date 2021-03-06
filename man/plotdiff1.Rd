\name{plotdiff1}
\alias{plotdiff1}
\title{Plot genetic differentiation (Delta D) across loci}
\description{Plot genetic differentiation (Delta D, Gaggiotti OE, et al, 2018) across loci.}
\usage{plotdiff1(x, ncode)}
\arguments{
  \item{x}{Genetic data, should be genepop format}
  \item{ncode}{The coding type of your genetic data.}
  }
\details{Function generates two genetic differentiation (Delta D) plots, a box plot and a scatter plot}
\references{
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.
}
\author{qinxinghu@gmail.com}
\seealso{DeltaDcorplot}
\examples{
##plot genetic differentiation across loci
# example genepop file
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
plotdiff1(infile,ncode=3)}
\keyword{Plot}
\keyword{Plot genetic differentiation}
