\name{Hierpertest}
\alias{Hierpertest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test the significance among aggregates (regions, ecosystems...)
%%  ~~function to do ... ~~
}
\description{ Permutation test of the significance among aggregates,e.g., test the significance of the allelic differentiation within ecosystem among different regions
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Hierpertest(f, ncode, level, permutations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ Input file, genepop format
%%     ~~Describe \code{f} here~~
}
  \item{ncode}{ Genotype coding type, whether 2 digital or 3 digital
%%     ~~Describe \code{ncode} here~~
}
  \item{level}{ The hierarchical level you want to test
%%     ~~Describe \code{level} here~~
}
  \item{permutations}{ permutations
%%     ~~Describe \code{permutations} here~~
}
}

\value{
\item{pertest}{Permutation test result for the level you specified}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Anderson, M. J. (2014). Permutational multivariate analysis of variance (PERMANOVA). Wiley StatsRef: Statistics Reference Online, 1-15.
%% ~put references to the literature/web site here ~
}
\author{qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
# example genepop file
f <- system.file('extdata',package='HierDpart')
infile <- file.path(f, "Island.gen")
region1=paste("region",rep(1,time=6))
region2=paste("region",rep(2,time=6))
region3=paste("region",rep(3,time=2))
region4=paste("region",rep(4,time=2))

level1=data.frame(matrix(data=0,nrow=16,ncol=1))
level1[1:6,1]=region1
level1[7:12,1]=region2
level1[13:14,1]=region3
level1[15:16,1]=region4
colnames(level1)=c("region")

hiertest=Hierpertest(infile,ncode=3,level1$region,permutations=999)
hiertest
}

