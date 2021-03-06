\name{COR_DeltaDd}
\alias{COR_DeltaDd}
\title{Function for calculating pairwise Delta D and correlation between genetic differentiation (Delta D) and geographic distance}
\description{This function calculates new pairwise genetic differentiation, Delta D (Gaggiotti OE, et al, 2018), as well as the correlation between genetic differentiation (Delta D) and geographic distance.}
\usage{COR_DeltaDd(f, d, ncode, nrepet)}
\arguments{
  \item{f}{A genetic data files, the format is genepop format.}
  \item{d}{Pairwise geographic distance for calculating correlations. If no actual geographic distance, d=NULL. If TRUE, the geographic distance should be a matrix that corresponds to the pairwise genetic differentiation matrix.}
  \item{ncode}{The coding type of your genetic data.}
  \item{nrepet}{Number of permutations.}
}
\details{This function returns to lists of pairwise Delta D matrix, and a matrix of theorical geographical distance if d is NULL, as well as the correlation coefficient.}
\value{
  \item{PairwiseDeltaD}{The pairwise genetic differentiation matrix (Delta D).}
  \item{Dgeo}{The pairwise geographical distance matrix.}
  \item{CorDeltaDd}{Pearson correlation coefficient between genetic differentiation (Delta D) and geographic distance.}
}
\references{
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{qinxinghu@gmail.com}
\seealso{COR_Fst}
\examples{
# example genepop file
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
d1=COR_DeltaDd(infile,d=NULL,ncode=3,nrepet=999)
print(d1)}
\keyword{Correlation between genetic differentiation and geographic distance}

