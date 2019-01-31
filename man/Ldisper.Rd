\name{Ldisper}
\alias{Ldisper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test significance of dispersion among subaggregate within aggregate
%%  ~~function to do ... ~~
}
\description{ Test significance of dispersion among subaggregate within aggregate.Here we used dispersion for testing the variance within the aggregate
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Ldisper(f, ncode, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ genepop file
%%     ~~Describe \code{f} here~~
}
  \item{ncode}{ genetype coding type
%%     ~~Describe \code{ncode} here~~
}
  \item{group}{ the level of aggregate you want to test
%%     ~~Describe \code{group} here~~
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{dispersion}{Multivariate homogeneity of group dispersions}
\item{pertestdis}{Permutation test of significance among groups(aggregates)}
%% ...
}
\references{
Anderson, M.J. (2006) Distance-based tests for homogeneity of multivariate dispersions. Biometrics 62, 245-253.

Anderson, M.J., Ellingsen, K.E., McArdle, B.H. (2006) Multivariate dispersion as a measure of beta diversity. Ecology Letters 9, 683-693.

ONeill, M.E. (2000) A Weighted Least Squares Approach to Levene s Test of Homogeneity of Variance. Australian and New Zealand Journal of Statistics 42, 8-100.

}
\author{qinxinghu@gmail.com
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

ldis=Ldisper(infile,ncode=3,group=level1$region)
ldis$dispersion
ldis$pertestdis
}


