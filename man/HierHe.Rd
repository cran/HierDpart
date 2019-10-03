\name{HierHe}
\alias{HierHe}
\title{Function for calculating hierarchical heterozygosity}
\description{This function calculates and decomposes hierarchical heterozygosity into different levels. We can obtain allelic diversity (He) in any hierarchy.}
\usage{HierHe(x, nreg, r, ncode)}
\arguments{
  \item{x}{Genetic data, the format should be genepop file with proper extention.}
  \item{nreg}{Number of aggregate in the hierarchy.}
  \item{r}{The number of subaggregate in the hierarchy.}
  \item{ncode}{Coding type of genetic data, whether 3, or 2.}
}
\details{In line with other genetic metrics, this function obtains hierarchical expected heterozygosity (He) for metapopulations. The result gives us hierarchical allelic diversity at different levels that you set in the input argument.}
\value{
  \item{HierHe_perloc}{List of hierarchical allelic diversity (heterozygosity) at different levels per locus.}
  \item{HierHr}{Hierarchical allelic diversity (heterozygosity) at aggregate/regional levels.}
    \item{Hrperloc}{Hierarchical allelic diversity (heterozygosity) per locus at aggregate/regional levels.}
   \item{HieHr_overall}{Overall (mean) hierarchical allelic diversity (heterozygosity) at different levels.}
}
\references{
Goudet, J. (2005). Hierfstat, a package for R to compute and test hierarchical F-statistics. Molecular Ecology Notes, 5(1), 184-186.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), 1403-1405.

Nei, M. (1977). F-statistics and analysis of gene diversity in subdivided populations. Annals of human genetics, 41(2), 225-233.
}
\author{qinxinghu@gmail.com}
\seealso{HierAr,HierD,HierFst}
\examples{
# example genepop file
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
HiHe=HierHe(infile,nreg=4,r=c(7,4,2,3),ncode=3)
print(HiHe)
}
\keyword{Hierarchical genetic decomposition}

