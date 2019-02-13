\name{pwag.permutest}
\alias{pwag.permutest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Permutation test for allelic differentiation (Delta D) among subaggregates within aggregate
%%  ~~function to do ... ~~
}
\description{ This test examines the allelic differentiation (Delta D) among subaggregates within aggregate. The permutation works on the pairwise allelic differentiation between subaggregates within different aggregates.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pwag.permutest(f, ncode, group, permutations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ Genepop file
%%     ~~Describe \code{f} here~~
}
  \item{ncode}{The number coding genotypes
%%     ~~Describe \code{ncode} here~~
}
  \item{group}{ The group or level of the aggregates
%%     ~~Describe \code{group} here~~
}
  \item{permutations}{ Number of permutations
%%     ~~Describe \code{permutations} here~~
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{pwpermutest.detail}{Permutation test detail of the pairwise difference between aggregates}
\item{pwtest.aggregate}{Pairwise permutation P-value between aggregates}
%% ...
}
\references{
Anderson, M. J. (2014). Permutational multivariate analysis of variance (PERMANOVA). Wiley StatsRef: Statistics Reference Online, 1-15.
%% ~put references to the literature/web site here ~
}
\author{ qinxinghu@gmail.com
%%  ~~who you are~~
}

\examples{
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
pwag.permutest(infile, ncode=3, group=level1$region, permutations=999)
}
