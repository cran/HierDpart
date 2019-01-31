\name{COR_Fstd}
\alias{COR_Fstd}
\title{Function for calculating pairwise Fst and correlation coefficient between Fst (Weir and Cockerham, 1984) and geographic distance}
\description{This function calculates pairwise genetic differentiation (Fst), as well as the correlation coefficient between genetic differentiation (Fst) and geographic distance.}
\usage{COR_Fstd(x, d, ncode,nrepet)}
\arguments{
  \item{x}{A genetic data files, here the format is genepop format.}
  \item{d}{Pairwise geographic distance for calculating correlations. If no actual geographic distance, d=NULL. The geographic distance should be a matrix that corresponds to the pairwise genetic differentiation matrix.}
  \item{ncode}{The code type of the genotype.}
  \item{nrepet}{Number of permutations.}
}
\details{ This function returns lists of pairwise Fst matrix, and a matrix of theorical geographical distance if d is NULL, as well as the correlation coefficient.}
\value{
  \item{pwFst}{The pairwise genetic differentiation (Fst) matrix.}
  \item{COR_Fstd}{Pearson correlation coefficient between genetic differentiation (Fst) and geographic distance.}
}
\references{
Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5(1), 184-186.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24 (11), 1403-1405.

Weir, B.S. (1996) Genetic Data Analysis II. Sinauer Associates.

Weir B.S. and Cockerham C.C. (1984) Estimating F-Statistics for the Analysis of Population Structure. Evolution 38:1358

Evolution and the Genetics of Populations Volume 2: the Theory of Gene Frequencies, pg 294-295, S. Wright, Univ. of Chicago Press, Chicago, 1969.}
\author{ qinxinghu@gmail.com}
\seealso{COR_detaDd}
\examples{
# example genepop file, the example file is a complicated hierarchical island model
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
#d2=COR_Fstd(infile,d=NULL,ncode=3,nrepet=999)
#print(d2)
}
\keyword{Correlation between genetic differentiation and geographic distance}

