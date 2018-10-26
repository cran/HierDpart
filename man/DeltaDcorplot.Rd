\name{DeltaDcorplot}
\alias{DeltaDcorplot}
\title{Plot correlation matrix for pairwise Delta D}
\description{This function plots the pairwise genetic differentiation (Delta D).}
\usage{DeltaDcorplot(x, ncode)}
\arguments{
  \item{x}{The genetic data file, specifically genepop format.}
  \item{ncode}{The coding type of the genetic data.}
}
\details{ This function returns the pairwise genetic differentiation (Delta D) plot.}
\value{The result will return a genetic differentiation (Delta D) correlation plot.}
\references{
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{qinxinghu@gmail.com}
\seealso{HierD; CorDdplot}
\examples{
# example genepop file
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
DeltaDcorplot(infile,ncode=3)
}
\keyword{Correlation plot}
\keyword{Plot genetic differentiation}
