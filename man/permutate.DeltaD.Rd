\name{permutate.DeltaD}
\alias{permutate.DeltaD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Permutation test for allelic differentiation (Delta D)
%%  ~~function to do ... ~~
}
\description{ This function performs permutations for pairwise comparison of allelic differentiation.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
permutate.DeltaD(file, ncode, permutations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ Genepop file
%%     ~~Describe \code{file} here~~
}
  \item{ncode}{Genotype coding
%%     ~~Describe \code{ncode} here~~
}
  \item{permutations}{ Number of permutations
%%     ~~Describe \code{permutations} here~~
}
}
\details{ This function runs slowly because the permutations work on the random sampling of the individuals from raw data and calculating the statistics every iteration. The permutation test here shuffles the individuals  and samples the individuals to calculate the Delta D every permutation. Significance of the allelic differentiation is obtained from Monte-Carlo test.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{pwDeltaD.obs }{ The observed value of allelic differentiation Delta D}
\item{pwalltests}{The test results of pairwise Delta D}
%% ...
}
\references{
Thioulouse, J., & Dray, S. (2007). Interactive multivariate data analysis in R with the ade4 and ade4TkGUI packages. Journal of Statistical Software, 22(5), 1-14.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
permutate.DeltaD=function(infile,ncode=3,permutations=999)
}

