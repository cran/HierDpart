\name{HierD}
\alias{HierD}
\title{Function for calculating and decomposing genetic diversity (D, q=1) and differentiation (Delta D)}
\description{ This function calculates and decomposes genetic diversity (Hill number q=1) and differentiation (Delta D) (Oscar, E et al, 2018) into alpha, beta, gamma diversity, as well as Delta D2 and Delta D1 differentiation.}
\usage{HierD(x, nreg, r, ncode)}
\arguments{
  \item{x}{The genotypic data, genepop format}
  \item{nreg}{Number of aggregate in the hierarchical structure.}
  \item{r}{The number of subaggregate in the hierarchy.}
  \item{ncode}{Coding type of genetic data.}
}
\details{ This function calculates the hierarchical genetic differentiation (Delta D) for metapopulations structured in any hierarchy. In this function, you can input genepop file (x) and define structure (nreg, r) in the arguments, and then returns a list of your hierarchical genetic diversity and differentiation (Delta D ). The result contains D_gamma, D_alpha.2, D_alpha.1, D_beta.2, D_beta.1, Differentiation.2, Differentiation.1, which corresponds to the diversity at different hierarchical levels.}
\value{
  \item{Dst}{Result returns to a list of hierarchical diversities (q=1) and genetic differentiations, including D_gamma (total diversity), D_alpha.2 (average total diversity between aggregates (eg. regions) within the ecosystem), D_alpha.1 (average diversity between subaggregates (eg.subpopulations)), D_beta.2 (number of aggregate equivalent), D_beta.1 (number of subaggregate equivalent), Differentiation.2 (genetic differentiation between aggragates), Differentiation.1 (genetic differentiation between subaggregates).}
}
\references{
Gaggiotti, O. E., Chao, A., Peres-Neto, P., Chiu, C. H., Edwards, C., Fortin, M. J., ... & Selkoe, K. A. (2018). Diversity from genes to ecosystems: A unifying framework to study variation across biological metrics and scales. Evolutionary Applications.}
\author{ qinxinghu@gmail.com}
\seealso{HierAr; HierFst; HierHe}
\examples{
# example genepop file
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
hierD=HierD(infile,nreg=4,r=c(7,4,2,3),ncode=3)
print(hierD)}
\keyword{Hierarchical genetic decomposition}
