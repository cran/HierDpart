\name{plot_pcoa_aggregate}
\alias{plot_pcoa_aggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Draw the aggregation plot for allelic differentiation
%%  ~~function to do ... ~~
}
\description{ Ordination plot of the allelic differentiation. Draw the aggregation groups according to the hierarchical level.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_pcoa_aggregate(f, ncode, level, label = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{ Input files, genepop format
%%     ~~Describe \code{f} here~~
}
  \item{ncode}{ Genotype coding type
%%     ~~Describe \code{ncode} here~~
}
  \item{level}{Hierarchical level (group)
%%     ~~Describe \code{level} here~~
}
  \item{label}{ Whether lable the aggregates or not, TRUE or FALSE
%%     ~~Describe \code{label} here~~
}
}
\details{ PcoA ordination plot of the allelic differentiation between aggregates. If there is distinguishable difference between aggregates, the plot identifies the different aggregations
%%  ~~ If necessary, more details than the description above ~~
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
plot_pcoa_aggregate(infile, ncode=3, level=level1$region, label = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
