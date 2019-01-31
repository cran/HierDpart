\name{plot.aggregates}
\alias{plot.aggregates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot aggregates with dispersions (variances)
%%  ~~function to do ... ~~
}
\description{ Plot the aggregates by groups
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.aggregates(disper, hull = FALSE, ellipse = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{disper}{ The dispersion data
%%     ~~Describe \code{disper} here~~
}
  \item{hull}{ whether plot hull shape
%%     ~~Describe \code{hull} here~~
}
  \item{ellipse}{ whether plot ellipse by group
%%     ~~Describe \code{ellipse} here~~
}
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
plot(ldis$dispersion,hull=FALSE, ellipse=TRUE)

}

